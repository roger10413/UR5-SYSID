# -*- coding: utf-8 -*-
"""
J4、J5 速度調變實驗（robot/ur5_j4_velocity_modulation_ident.py、ur5_j5_velocity_modulation_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -p "test_*.py"

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在開始姿態不動）。兩支程式由 J3 速度調變改寫，以 AXES 設定區分。

測試內容（J4、J5 各一份）
--------
1. 語法正確；姿態表：J4 預設 UP、J5 預設 DIAG，中心 J4 −90°、J5 0°；條件同 J3 速度調變（J2 的 1.5 倍）、調變 3 s。
2. 規劃：三個階段為 1、6、42 段；調變段中心 = 姿態中心 + 錯開量；兩個姿態、兩個方向、三個階段的規劃範圍都在
   同一軸定速 FULL 的範圍（含餘裕）內；第一次 movej 不超過上限。
3. 安全檢查：通過；視窗未設定、太窄、最低速度過低、非整數週期、speedj 加速度上限不足、J4 法蘭線速度超限時 FAIL。
4. 姿態：各自六軸 ±0.5°；互換、J3 程式的 UP、J0 不在 180°、NaN、未知姿態一律拒絕。
5. URScript：縮排規則、只含 ASCII、speedj 只有該軸非零、movej 只改該軸、最後回到開始角度、角度保護方向正確。
6. main()：各種錯誤都在送出前中止；QUICK 端對端自動命名（j4velmod_data_up_p130_pos_quick_…）；UP 提醒把 J2 轉回；
   工具加速度寫進 CSV。
7. 分析（合成資料：J = 0.15、J0 實測的漣波、白雜訊 + 慢變雜訊、伺服延遲 3 ms；重力保持力矩：J4 UP、J5 為 0，
   J4 DIAG 為 2.34·sin(q − (−89.7°)) N·m）：各條件找回 J；J4 DIAG「無重力項」時 1 Hz 的 J 偏低（假慣量約 −0.06），
   J4 UP 與 J5 幾乎不變；target_moment 回歸找回 J_model；QUICK 只列各段、不當機。
"""

import os
import sys
import csv
import math
import time
import types
import socket
import builtins
import tempfile
import unittest
import py_compile
import importlib.util
import contextlib
import io

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROBOT = os.path.join(ROOT, "robot") if os.path.isdir(os.path.join(ROOT, "robot")) else ROOT

RAD = lambda q: [math.radians(x) for x in q]
UP_Q = RAD((180.0, -90.0, 0.0, -90.0, -90.0, 0.0))
HOME_Q = RAD((180.0, -90.0, 90.0, -90.0, -90.0, 0.0))
J3UP_Q = RAD((180.0, -90.0, 0.0, -130.0, -90.0, 0.0))
POSE_Q = {"UP": UP_Q, "DIAG": HOME_Q}
PAYLOAD = 1.3

AXES = {
    4: dict(path=os.path.join(ROBOT, "ur5_j4_velocity_modulation_ident.py"), default="UP", center=-90.0,
            win=(-112.0, -68.0), cv_full=(-110.71, -69.29), radius=0.0823,
            grav={"UP": (0.0, 0.0), "DIAG": (2.34, math.radians(-89.7))}),
    5: dict(path=os.path.join(ROBOT, "ur5_j5_velocity_modulation_ident.py"), default="DIAG", center=0.0,
            win=(-22.0, 22.0), cv_full=(-20.71, 20.71), radius=0.0,
            grav={"UP": (0.0, 0.0), "DIAG": (0.0, 0.0)}),
}


class FakeRTDE:
    q = list(UP_Q)
    sf = 1.0
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getTargetQdd(self): return [0.0] * 6
    def getTargetMoment(self): return [0.0] * 6
    def getJointTemperatures(self): return [30.0] * 6
    def getTargetSpeedFraction(self): return FakeRTDE.sf
    def disconnect(self): pass


def load(n, name, pose=None):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    spec = importlib.util.spec_from_file_location(name, AXES[n]["path"])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PAYLOAD_KG = PAYLOAD
    if pose is not None:
        mod.POSE = pose
    setattr(mod, f"J{n}_SAFE_MIN_DEG", AXES[n]["win"][0])
    setattr(mod, f"J{n}_SAFE_MAX_DEG", AXES[n]["win"][1])
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(default, payload=str(PAYLOAD)):
    return lambda prompt="", *a, **k: payload if "負載質量" in str(prompt) else default


def plan(m, stage, d, pose):
    conds, n = m.stage_config(stage)
    return m.plan_segments(conds, n, d, m.V_CRUISE, m.T_MOD, m.SETTLE_TIME, m.RAMP_ACCEL, m.STOPJ_DECEL,
                           m.OVERRUN_DEG, m.POSES[pose]["center_deg"], m.GEAR_RATIO, m.DT)


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


class Base(unittest.TestCase):
    N = None

    @property
    def cfg(self):
        return AXES[self.N]

    def setUp(self):
        self._sock, self._input = socket.socket, builtins.input
        socket.socket = no_network
        builtins.input = answers(self.cfg["default"])
        FakeRTDE.q, FakeRTDE.sf = list(POSE_Q[self.cfg["default"]]), 1.0

    def tearDown(self):
        socket.socket, builtins.input = self._sock, self._input
        FakeRTDE.q, FakeRTDE.sf = list(UP_Q), 1.0
        sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE

    def m(self, tag, pose=None):
        return load(self.N, f"jv{self.N}_{tag}", pose)

    def safety(self, m, segs, pose, **kw):
        q = m.POSES[pose]["q_deg"]
        a = dict(v_cruise=m.V_CRUISE, t_mod=m.T_MOD, qd_max=m.QD_MAX, qdd_max=m.QDD_CHECK_LIMIT,
                 v_min_allowed=m.V_MIN_ALLOWED, ramp_accel=m.RAMP_ACCEL, mod_accel=m.MOD_ACCEL_LIMIT,
                 decel=m.STOPJ_DECEL, movej_vel=m.MOVEJ_VEL, movej_accel=m.MOVEJ_ACCEL,
                 safe_min_deg=self.cfg["win"][0], safe_max_deg=self.cfg["win"][1],
                 margin_deg=m.SAFETY_MARGIN_DEG, q_pose_deg=q, end_deg=q[self.N])
        a.update(kw)
        return m.safety_check(segs, **a)


class PlanCases:
    def test_compile(self):
        py_compile.compile(self.cfg["path"], doraise=True)

    def test_settings(self):
        m = self.m("set")
        n = self.N
        self.assertEqual(m.JOINT_INDEX, n)
        self.assertEqual(m.POSE, self.cfg["default"])
        self.assertEqual(set(m.POSES), {"UP", "DIAG"})
        for pose, q in POSE_Q.items():
            self.assertEqual(m.POSES[pose]["center_deg"], self.cfg["center"])
            for a, b in zip(m.POSES[pose]["q_deg"], q):
                self.assertAlmostEqual(math.radians(a), b, places=9, msg=pose)
        self.assertEqual(m.V_CRUISE, 0.12)
        self.assertEqual(m.CONDITIONS, [(1.0, 0.09), (2.0, 0.09), (3.0, 0.09), (2.0, 0.045), (1.0, 0.045), (3.0, 0.06)])
        self.assertEqual(m.T_MOD, 3.0)
        acc = max(2 * math.pi * f * a for f, a in m.CONDITIONS)
        self.assertGreater(m.MOD_ACCEL_LIMIT, acc)
        self.assertLessEqual(m.MOD_ACCEL_LIMIT, m.QDD_CHECK_LIMIT)
        self.assertEqual([s for s, _, _ in m.VECTOR_SIGNALS], ["tool_acc"])

    def test_stages_and_ranges(self):
        m = self.m("plan")
        n = self.N
        for pose in ("UP", "DIAG"):
            c = m.POSES[pose]["center_deg"]
            for stage, n_seg in (("QUICK", 1), ("MID", 6), ("FULL", 42)):
                for d in (+1, -1):
                    segs = plan(m, stage, d, pose)
                    msg = f"J{n} {pose} {stage} {d:+d}"
                    self.assertEqual(len(segs), n_seg, msg)
                    for s in segs:
                        mid = 0.5 * (s["mod_start_deg"] + s["mod_end_deg"])
                        self.assertAlmostEqual(mid, c + d * s["offset_deg"], places=9)
                    ok, lines, (lo, hi) = self.safety(m, segs, pose)
                    self.assertTrue(ok, msg + "\n" + "\n".join(lines))
                    self.assertGreaterEqual(lo, self.cfg["cv_full"][0] - 0.01, msg)
                    self.assertLessEqual(hi, self.cfg["cv_full"][1] + 0.01, msg)
                    first = abs(segs[0]["start_deg"] - m.POSES[pose]["q_deg"][n])
                    self.assertLess(first, 19.0, msg)

    def test_safety_failures(self):
        m = self.m("safety")
        for pose in ("UP", "DIAG"):
            segs = plan(m, "FULL", +1, pose)
            self.assertFalse(self.safety(m, segs, pose, safe_min_deg=None)[0], "視窗未設定")
            self.assertFalse(self.safety(m, segs, pose, safe_min_deg=self.cfg["center"] - 5.0)[0], "視窗太窄")
            self.assertFalse(self.safety(m, segs, pose, v_min_allowed=0.07)[0], "最低速度過低")
            self.assertFalse(self.safety(m, segs, pose, t_mod=3.3)[0], "非整數週期")
            self.assertFalse(self.safety(m, segs, pose, mod_accel=1.5)[0], "speedj 加速度上限小於調變所需")
            ok, lines, _ = self.safety(m, segs, pose)
            k = next(i for i, l in enumerate(lines) if l.startswith("法蘭線速度"))
            self.assertIn(f"法蘭離 J{self.N} 軸 {self.cfg['radius']:.3f} m", lines[k])
            self.assertIn(f"{self.cfg['radius'] * 0.21:.4f} m/s", lines[k])
            if self.N == 4:
                self.assertFalse(self.safety(m, segs, pose, flange_speed_max=0.01)[0], "法蘭線速度超限")
        c = self.cfg["center"]
        self.assertFalse(m.check_start_position(c, c - 18.0, 40.0, None, None)[0])
        self.assertFalse(m.check_start_position(c + 30.0, c - 18.0, 40.0, *self.cfg["win"])[0], "目前角度在視窗外")
        self.assertFalse(m.check_start_position(c, c - 18.0, 10.0, *self.cfg["win"])[0], "第一次 movej 過遠")
        self.assertTrue(m.check_start_position(c, c - 18.0, 40.0, *self.cfg["win"])[0])

    def test_check_pose(self):
        m = self.m("pose")
        n = self.N
        self.assertTrue(m.check_pose(UP_Q, "UP")[0])
        self.assertTrue(m.check_pose(HOME_Q, "DIAG")[0])
        self.assertTrue(m.check_pose([UP_Q[0] + 2 * math.pi] + UP_Q[1:], "UP")[0], "±360° 等價")
        self.assertFalse(m.check_pose(HOME_Q, "UP")[0])
        self.assertFalse(m.check_pose(UP_Q, "DIAG")[0])
        self.assertFalse(m.check_pose(J3UP_Q, "UP")[0], "J3 程式的 UP（J3 = −130°）")
        self.assertFalse(m.check_pose([math.radians(170.0)] + UP_Q[1:], "UP")[0], "J0 不在 180°")
        q = list(HOME_Q); q[n] += math.radians(0.6)
        self.assertFalse(m.check_pose(q, "DIAG")[0], f"J{n} 不在開始角度")
        self.assertFalse(m.check_pose([float("nan")] + UP_Q[1:], "UP")[0])
        self.assertFalse(m.check_pose(None, "UP")[0])
        self.assertFalse(m.check_pose(UP_Q, "EXT45")[0], "未知姿態")


class URScriptCases:
    def test_structure(self):
        m = self.m("script")
        n = self.N
        for pose in ("UP", "DIAG"):
            base, end = POSE_Q[pose], m.POSES[pose]["q_deg"][n]
            for d in (+1, -1):
                segs = plan(m, "MID", d, pose)
                text = m.build_urscript(segs, d, n, base, m.V_CRUISE, m.RAMP_ACCEL, m.MOD_ACCEL_LIMIT,
                                        m.STOPJ_DECEL, m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=end)
                self.assertTrue(urscript_indent_ok(text))
                text.encode("ascii")
                self.assertTrue(text.startswith(f"def j{n}_velocity_modulation_program():"))
                for line in text.splitlines():
                    if "speedj(" in line:
                        v = line.split("[")[1].split("]")[0].split(",")
                        self.assertEqual([x.strip() for k, x in enumerate(v) if k != n], ["0.0"] * 5, line)
                    if "while" in line:
                        self.assertIn(f"q_now[{n}] {'<' if d > 0 else '>'} ", line)
                movej = [l for l in text.splitlines() if "movej(" in l]
                self.assertEqual(len(movej), len(segs) + 1)
                for line in movej:
                    vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    for j in range(6):
                        if j != n:
                            self.assertAlmostEqual(vals[j], base[j], places=5, msg=f"只改 J{n}")
                last = [float(x) for x in movej[-1].split("[")[1].split("]")[0].split(",")]
                self.assertAlmostEqual(math.degrees(last[n]), end, places=4)


class MainCases:
    def _run(self, m, stage="QUICK", direction=+1):
        scripts, files = [], {}
        m.send_urscript = lambda ip, text, port=30002: scripts.append(text)
        m.wait_for_motion_complete = lambda logger, planned_s, **k: ("done", planned_s, planned_s)
        out = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            m.TEST_STAGE, m.DIRECTION = stage, direction
            with contextlib.redirect_stdout(out):
                m.main()
            for root, _, names in os.walk(tmp):
                for nm in names:
                    with open(os.path.join(root, nm), encoding="utf-8") as f:
                        files[nm] = f.read()
        self.last_out = out.getvalue()
        return scripts, files

    def test_aborts_before_sending(self):
        n, default = self.N, self.cfg["default"]
        other = "DIAG" if default == "UP" else "UP"
        cases = [("window", "安全檢查未通過"), ("pose", "姿態檢查未通過"), ("j3up", "姿態檢查未通過"),
                 ("name", "已取消"), ("payload", "負載確認未通過"), ("slider", "速度滑桿檢查未通過"),
                 ("rtde", "RTDE 記錄未能啟動")]
        for k, (case, reason) in enumerate(cases):
            m = self.m(f"abort{k}")
            if case == "window":
                setattr(m, f"J{n}_SAFE_MIN_DEG", None)
            elif case == "pose":
                FakeRTDE.q = list(POSE_Q[other])
            elif case == "j3up":
                FakeRTDE.q = list(J3UP_Q)
            elif case == "name":
                builtins.input = answers("yes")
            elif case == "payload":
                builtins.input = answers(default, payload="2.30")
            elif case == "slider":
                FakeRTDE.sf = 0.5
            elif case == "rtde":
                class Failing:
                    def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
                sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
            scripts, _ = self._run(m)
            self.assertEqual(scripts, [], case)
            self.assertIn(reason, self.last_out, case)
            FakeRTDE.q, FakeRTDE.sf = list(POSE_Q[default]), 1.0
            builtins.input = answers(default)
            sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE

    def test_quick_end_to_end(self):
        n = self.N
        for pose in ("UP", "DIAG"):
            m = self.m(f"e2e_{pose}", pose)
            FakeRTDE.q = list(POSE_Q[pose])
            builtins.input = answers(pose)
            for d, label in ((+1, "pos"), (-1, "neg")):
                scripts, files = self._run(m, direction=d)
                self.assertEqual(len(scripts), 1, f"{pose} {label}")
                self.assertTrue(scripts[0].startswith(f"def j{n}_velocity_modulation_program():"))
                self.assertTrue(any(nm.startswith(f"j{n}velmod_data_{pose.lower()}_p130_{label}_quick_") for nm in files),
                                list(files))
                info = next(v for nm, v in files.items() if nm.startswith(f"j{n}velmod_info_"))
                self.assertIn(f"姿態: {pose}", info)
                self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
                if pose == "UP":
                    self.assertIn("J2 轉回 90°", self.last_out)
                else:
                    self.assertNotIn("[提醒]", self.last_out)

    def test_records_tool_accelerometer(self):
        n, default = self.N, self.cfg["default"]
        m = self.m("toolacc", default)
        FakeRTDE.getActualToolAccelerometer = lambda self: [0.1, -9.8, 0.2]
        try:
            scripts, files = self._run(m)
        finally:
            del FakeRTDE.getActualToolAccelerometer
        self.assertEqual(len(scripts), 1)
        data = next(v for nm, v in files.items() if nm.startswith(f"j{n}velmod_data_"))
        rows = list(csv.DictReader(io.StringIO(data)))
        self.assertTrue(rows)
        self.assertEqual([float(rows[0][f"tool_acc_{a}"]) for a in "xyz"], [0.1, -9.8, 0.2])
        info = next(v for nm, v in files.items() if nm.startswith(f"j{n}velmod_info_"))
        self.assertIn(f"工具加速度（成功讀取列數）: tool_acc {len(rows)}\n", info)


# ------------------------------------------------------------
# 合成資料（同 test_j3_velmod_offline.py，改成任意軸）
# ------------------------------------------------------------

J_TRUE, TC_TRUE, B_TRUE, C_TRUE = 0.15, 2.5, 6.0, -2.0
LSB = 1.192e-5
RIPPLE = {+1: [(2, -0.212, 0.579), (5, -0.164, -1.051), (10, 0.453, 0.077), (12, -0.54, 0.738)],
          -1: [(2, -0.315, 0.027), (5, 0.133, 1.068), (10, -0.385, 0.005), (12, 0.48, -0.734)]}


def movej_profile(dist, vmax, acc, h):
    D = abs(dist)
    if D < 1e-12:
        return np.zeros(0)
    if D <= vmax ** 2 / acc:
        T = 2 * math.sqrt(D / acc); vpk = acc * T / 2
    else:
        T = D / vmax + vmax / acc; vpk = vmax
    t = np.arange(0, T, h)
    v = np.clip(np.minimum.reduce([acc * t, np.full_like(t, vpk), acc * (T - t)]), 0, None)
    return math.copysign(1.0, dist) * v * D / (np.sum(v) * h)


def synth_csv(m, path, n, d, segs, pose, m_grav, q_zero, servo_ms=3.0, noise=0.6, slow=0.3, seed=0, h=0.001,
              model=None):
    """J(n) 力矩 = J·a + sign(v)·(Tc + B|v| + C v² + 漣波) + M·sin(q − q0) ＋ 雜訊；實際速度比目標晚 servo_ms。"""
    rng = np.random.default_rng(seed)
    sub = int(round(m.DT / h))

    def cyc(x):
        return np.r_[x, np.zeros((-len(x)) % sub)]
    n_pause = int(round(m.PAUSE_TIME / m.DT)) * sub
    n_stop = int(math.ceil(m.V_CRUISE / m.STOPJ_DECEL / m.DT)) * sub
    q0 = POSE_Q[pose][n]
    pieces, pos = [], q0
    for s in segs:
        st = math.radians(s["start_deg"])
        pieces.append(cyc(movej_profile(st - pos, m.MOVEJ_VEL, m.MOVEJ_ACCEL, h)))
        tp = np.arange(s["n_pre"] * sub) * h
        tm = np.arange(s["n_mod"] * sub) * h
        ts = np.arange(n_stop) * h
        seg = [np.zeros(n_pause), d * np.minimum(m.RAMP_ACCEL * tp, m.V_CRUISE),
               d * (m.V_CRUISE + s["amp"] * np.sin(2 * math.pi * s["freq"] * tm)),
               d * np.maximum(m.V_CRUISE - m.STOPJ_DECEL * ts, 0.0), np.zeros(n_pause)]
        pieces += seg
        pos = st + np.sum(np.concatenate(seg)) * h
    pieces.append(cyc(movej_profile(q0 - pos, m.MOVEJ_VEL, m.MOVEJ_ACCEL, h)))
    pieces.append(np.zeros(125 * sub))
    vt = np.concatenate(pieces)
    tf = np.arange(len(vt)) * h
    va = np.interp(tf - servo_ms / 1000.0, tf, vt)
    q = q0 + np.concatenate([[0.0], np.cumsum((va[1:] + va[:-1]) / 2 * h)])
    a = np.gradient(va, h)
    sv = np.sign(np.where(np.abs(va) > 1e-6, va, 0.0))
    phi = m.GEAR_RATIO * q
    fric = np.zeros(len(va))
    for sg in (+1, -1):
        mm = sv == sg
        vm = np.abs(va[mm])
        rip = sum(c * np.cos(k * phi[mm]) + s_ * np.sin(k * phi[mm]) for k, c, s_ in RIPPLE[sg])
        fric[mm] = sg * (TC_TRUE + B_TRUE * vm + C_TRUE * vm ** 2 + rip)
    tau = J_TRUE * a + fric + m_grav * np.sin(q - q_zero)
    ii = np.arange(0, len(vt), sub)
    nn = len(ii)
    ar, rho = np.zeros(nn), math.exp(-m.DT / 0.3)
    e = rng.normal(0, slow * math.sqrt(1 - rho ** 2), nn)
    for k in range(1, nn):
        ar[k] = rho * ar[k - 1] + e[k]
    tau_s = tau[ii] + rng.normal(0, noise, nn) + ar
    q_s = np.round((q[ii] + rng.normal(0, 1e-5, nn)) / LSB) * LSB
    head = ["timestamp", f"actual_q_{n}", f"actual_qd_{n}", f"actual_current_{n}", f"target_qd_{n}", f"joint_temp_{n}"]
    extra = np.zeros((nn, 0))
    if model is not None:
        qdd_t = np.gradient(vt, h)[ii]
        qt = q0 + np.concatenate([[0.0], np.cumsum((vt[1:] + vt[:-1]) / 2 * h)])[ii]
        tmom = model[0] * qdd_t + model[1] * vt[ii] + m_grav * np.sin(qt - q_zero)
        head += [f"target_qdd_{n}", f"target_moment_{n}"]
        extra = np.column_stack([qdd_t, tmom])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(head)
        for k in range(nn):
            w.writerow([round(k * m.DT, 6), q_s[k], va[ii[k]], tau_s[k] / m.KT_OUT, vt[ii[k]], 30.0]
                       + list(extra[k]))


_CACHE = {}


class AnalysisCases:
    # J4：兩個姿態、兩個方向；J5：DIAG 兩個方向 + UP 正轉（J5 的重力為 0，兩個姿態的合成資料只差起點）
    RUNS = {4: (("UP", +1), ("UP", -1), ("DIAG", +1), ("DIAG", -1)), 5: (("DIAG", +1), ("DIAG", -1), ("UP", +1))}

    def result(self, pose, d, stage="FULL", **kw):
        n = self.N
        key = (n, pose, d, stage, tuple(sorted(kw.items())))
        if key not in _CACHE:
            m = load(n, f"jv{n}_ana_{len(_CACHE)}", pose)
            segs = plan(m, stage, d, pose)
            mg, q0 = self.cfg["grav"][pose]
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                synth_csv(m, p, n, d, segs, pose, mg, q0, **kw)
                with contextlib.redirect_stdout(io.StringIO()):
                    analyze = getattr(m, f"analyze_j{n}_velocity_modulation")
                    _CACHE[key] = analyze(p, n, m.KT_OUT, d, m.CONDITIONS, m.V_CRUISE)
        return _CACHE[key]

    def test_recovers_J(self):
        for pose, d in self.RUNS[self.N]:
            r = self.result(pose, d)
            msg = f"J{self.N} {pose} 方向 {d:+d}"
            self.assertEqual(len(r["segments"]), 42, msg)
            self.assertLessEqual(r["skipped"], 1, msg)
            self.assertEqual(sorted(r["conditions"]), list(range(6)), msg)
            for ci, c in r["conditions"].items():
                cm = f"{msg} 條件 {ci + 1}"
                self.assertLess(abs(c["J"] - J_TRUE), 3 * c["se_J"] + 0.03, f"{cm}：J = {c['J']:.4f} ± {c['se_J']:.4f}")
                self.assertAlmostEqual(c["track_gain"], 1.0, delta=0.01, msg=cm)
            self.assertAlmostEqual(r["J_mean"], J_TRUE, delta=0.03, msg=msg)

    def test_gravity_coupling_without_model(self):
        """「無重力項」減「完整」（同一份資料，1 Hz 大振幅）：J4 DIAG 像單擺（K_g ≈ +2.34 N·m/rad），偏低約 −0.06；
        J4 UP 與 J5 沒有重力，兩種寫法只差雜訊（< 0.04）。"""
        k = "無重力項（只看影響）"
        for pose, d in self.RUNS[self.N]:
            c = self.result(pose, d)["conditions"][0]
            diff = c["J_variants"][k] - c["J"]
            msg = f"J{self.N} {pose} {d:+d}：{diff:+.4f}"
            if (self.N, pose) == (4, "DIAG"):
                self.assertLess(diff, -0.03, msg)
                self.assertGreater(diff, -0.10, msg)
            else:
                # 沒有重力時，兩種寫法只差重力項吸收掉的慢變雜訊（約 0.02，遠小於 1 Hz 的標準差約 0.07）
                self.assertLess(abs(diff), 0.04, msg)

    def test_target_moment_model(self):
        r = self.result(self.cfg["default"], -1, model=(0.12, 0.5))
        self.assertEqual(sorted(r["model"]), list(range(6)))
        for ci, mr in r["model"].items():
            self.assertAlmostEqual(mr["J_model"], 0.12, delta=1e-3, msg=f"條件 {ci + 1}")
        self.assertAlmostEqual(r["J_model_mean"], 0.12, delta=1e-3)

    def test_quick_only_lists_segments(self):
        r = self.result(self.cfg["default"], +1, stage="QUICK")
        self.assertEqual(len(r["segments"]), 1)
        self.assertEqual(r["conditions"], {})


class TestJ4VelmodPlan(PlanCases, Base):
    N = 4


class TestJ5VelmodPlan(PlanCases, Base):
    N = 5


class TestJ4VelmodScript(URScriptCases, Base):
    N = 4


class TestJ5VelmodScript(URScriptCases, Base):
    N = 5


class TestJ4VelmodMain(MainCases, Base):
    N = 4


class TestJ5VelmodMain(MainCases, Base):
    N = 5


class TestJ4VelmodAnalysis(AnalysisCases, Base):
    N = 4


class TestJ5VelmodAnalysis(AnalysisCases, Base):
    N = 5


if __name__ == "__main__":
    unittest.main(verbosity=2)
