# -*- coding: utf-8 -*-
"""
J2 速度調變實驗（robot/ur5_j2_velocity_modulation_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在開始姿態不動）。

測試內容
--------
1. 語法正確；姿態表：UP（預設，中心 −5.7°）、DIAG（中心 71°）；條件同 J0 速度調變、調變 3 s。
2. 規劃：三個階段為 1、6、42 段；調變段中心 = 姿態中心 + 錯開量；兩個姿態、兩個方向、三個階段的規劃範圍都在
   J2 定速 FULL 的範圍（UP −26.41°～+15.01°，DIAG 50.29°～92.00°，含餘裕）內；第一次 movej 不超過上限。
3. 安全檢查：通過；J2 視窗未設定、太窄、最低速度過低、非整數週期、法蘭線速度超限、speedj 加速度上限不足時 FAIL。
4. 姿態：UP、DIAG 各自六軸 ±0.5°；互換、J0 不在 180°、NaN、未知姿態一律拒絕。
5. URScript：縮排規則、只含 ASCII、speedj 只有 J2 非零、movej 只改 J2、最後回到開始角度、角度保護方向正確。
6. main()：視窗未設定、姿態不符、未輸入姿態名稱、負載不符、速度滑桿不是 100%、RTDE 失敗時都在送出前中止；
   QUICK 端對端自動命名（j2velmod_data_up_p130_pos_quick_…、j2velmod_data_diag_…）；UP 提醒把 J2 轉回 90°。
7. 分析（合成資料：J = 1.3、重力保持力矩 −23.4·sin(q + 5.73°) N·m（官方模型：UP 為重力零點附近、
   DIAG 約 −21.6～−23.4）、J0 實測的漣波、白雜訊 + 慢變雜訊、伺服延遲 3 ms）：
   兩個姿態、兩個方向各條件都找回 J；「重力改為各段斜率」也找回 J；「無重力項」時 UP 條件 1 的 J 偏高 0.4～0.8
   （−K_g/ω² ≈ +0.59），DIAG 偏差小於 UP；target_moment（含重力）回歸找回 J_model；QUICK 只列各段、不當機。
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
JV = os.path.join(ROBOT, "ur5_j2_velocity_modulation_ident.py")

UP_Q = [math.radians(x) for x in (180.0, -90.0, 0.0, -90.0, -90.0, 0.0)]
HOME_Q = [math.radians(x) for x in (180.0, -90.0, 90.0, -90.0, -90.0, 0.0)]
POSE_Q = {"UP": UP_Q, "DIAG": HOME_Q}
PAYLOAD = 1.3
WIN = {"UP": (-28.0, 17.0), "DIAG": (48.0, 93.0)}               # J2 定速測試用的視窗
J2CV_FULL = {"UP": (-26.41, 15.01), "DIAG": (50.29, 92.0)}      # J2 定速 FULL 的規劃範圍（含餘裕）


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


def load(name, pose="UP"):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    spec = importlib.util.spec_from_file_location(name, JV)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PAYLOAD_KG = PAYLOAD
    mod.POSE = pose
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(default="UP", payload=str(PAYLOAD)):
    return lambda prompt="", *a, **k: payload if "負載質量" in str(prompt) else default


def plan(m, stage, d, pose="UP"):
    conds, n = m.stage_config(stage)
    return m.plan_segments(conds, n, d, m.V_CRUISE, m.T_MOD, m.SETTLE_TIME, m.RAMP_ACCEL, m.STOPJ_DECEL,
                           m.OVERRUN_DEG, m.POSES[pose]["center_deg"], m.GEAR_RATIO, m.DT)


def safety(m, segs, pose="UP", **kw):
    q = m.POSES[pose]["q_deg"]
    a = dict(v_cruise=m.V_CRUISE, t_mod=m.T_MOD, qd_max=m.QD_MAX, qdd_max=m.QDD_CHECK_LIMIT,
             v_min_allowed=m.V_MIN_ALLOWED, ramp_accel=m.RAMP_ACCEL, mod_accel=m.MOD_ACCEL_LIMIT,
             decel=m.STOPJ_DECEL, movej_vel=m.MOVEJ_VEL, movej_accel=m.MOVEJ_ACCEL, safe_min_deg=WIN[pose][0],
             safe_max_deg=WIN[pose][1], margin_deg=m.SAFETY_MARGIN_DEG, q_pose_deg=q, end_deg=q[2])
    a.update(kw)
    return m.safety_check(segs, **a)


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


class Base(unittest.TestCase):
    def setUp(self):
        self._sock, self._input = socket.socket, builtins.input
        socket.socket = no_network
        builtins.input = answers()
        FakeRTDE.q, FakeRTDE.sf = list(UP_Q), 1.0

    def tearDown(self):
        socket.socket, builtins.input = self._sock, self._input
        FakeRTDE.q, FakeRTDE.sf = list(UP_Q), 1.0
        sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE


class TestPlan(Base):
    def test_compile(self):
        py_compile.compile(JV, doraise=True)

    def test_settings(self):
        m = load("jv2_set")
        self.assertEqual(m.JOINT_INDEX, 2)
        self.assertEqual(m.POSE, "UP")
        self.assertEqual(set(m.POSES), {"UP", "DIAG"})
        self.assertEqual(m.POSES["UP"]["center_deg"], -5.7)
        self.assertEqual(m.POSES["DIAG"]["center_deg"], 71.0)
        self.assertEqual(m.V_CRUISE, 0.12)
        self.assertEqual(m.CONDITIONS, [(1.0, 0.06), (2.0, 0.06), (3.0, 0.06), (2.0, 0.03), (1.0, 0.03), (3.0, 0.04)])
        self.assertEqual(m.T_MOD, 3.0)
        acc = max(2 * math.pi * f * a for f, a in m.CONDITIONS)
        self.assertAlmostEqual(acc, 1.131, delta=0.001)
        self.assertGreater(m.MOD_ACCEL_LIMIT, acc)

    def test_stages_and_ranges(self):
        m = load("jv2_plan")
        for pose in ("UP", "DIAG"):
            c = m.POSES[pose]["center_deg"]
            for stage, n_seg in (("QUICK", 1), ("MID", 6), ("FULL", 42)):
                for d in (+1, -1):
                    segs = plan(m, stage, d, pose)
                    msg = f"{pose} {stage} {d:+d}"
                    self.assertEqual(len(segs), n_seg, msg)
                    for s in segs:
                        mid = 0.5 * (s["mod_start_deg"] + s["mod_end_deg"])
                        self.assertAlmostEqual(mid, c + d * s["offset_deg"], places=9)
                        self.assertAlmostEqual(abs(s["mod_end_deg"] - s["mod_start_deg"]),
                                               math.degrees(m.V_CRUISE * m.T_MOD), places=9)
                    ok, lines, (lo, hi) = safety(m, segs, pose)
                    self.assertTrue(ok, msg + "\n" + "\n".join(lines))
                    self.assertGreaterEqual(lo, J2CV_FULL[pose][0] - 0.01, msg)
                    self.assertLessEqual(hi, J2CV_FULL[pose][1] + 0.01, msg)
                    first = abs(segs[0]["start_deg"] - m.POSES[pose]["q_deg"][2])
                    self.assertLess(first, m.MAX_INITIAL_MOVE_DEG, msg)

    def test_safety_failures(self):
        m = load("jv2_safety")
        for pose in ("UP", "DIAG"):
            segs = plan(m, "FULL", +1, pose)
            self.assertFalse(safety(m, segs, pose, safe_min_deg=None)[0], "視窗未設定")
            self.assertFalse(safety(m, segs, pose, safe_min_deg=m.POSES[pose]["center_deg"] - 5.0)[0], "視窗太窄")
            self.assertFalse(safety(m, segs, pose, v_min_allowed=0.07)[0], "最低速度過低")
            self.assertFalse(safety(m, segs, pose, t_mod=3.3)[0], "非整數週期")
            self.assertFalse(safety(m, segs, pose, flange_speed_max=0.05)[0], "法蘭線速度超限")
            self.assertFalse(safety(m, segs, pose, mod_accel=1.0)[0], "speedj 加速度上限小於調變所需")
        ok, lines, _ = safety(m, plan(m, "FULL", +1), "UP")
        k = next(i for i, l in enumerate(lines) if l.startswith("法蘭線速度"))
        self.assertIn("0.0889 m/s", lines[k])                    # 0.494 m × 0.18 rad/s
        self.assertFalse(m.check_start_position(0.0, -23.9, 40.0, None, None)[0])
        self.assertFalse(m.check_start_position(30.0, -23.9, 40.0, *WIN["UP"])[0], "目前角度在視窗外")
        self.assertFalse(m.check_start_position(0.0, -23.9, 20.0, *WIN["UP"])[0], "第一次 movej 過遠")
        self.assertTrue(m.check_start_position(0.0, -23.9, 40.0, *WIN["UP"])[0])

    def test_check_pose(self):
        m = load("jv2_pose")
        self.assertTrue(m.check_pose(UP_Q, "UP")[0])
        self.assertTrue(m.check_pose(HOME_Q, "DIAG")[0])
        self.assertTrue(m.check_pose([UP_Q[0] + 2 * math.pi] + UP_Q[1:], "UP")[0], "±360° 等價")
        self.assertFalse(m.check_pose(HOME_Q, "UP")[0], "診斷原點不是 UP")
        self.assertFalse(m.check_pose(UP_Q, "DIAG")[0])
        self.assertFalse(m.check_pose([math.radians(170.0)] + UP_Q[1:], "UP")[0], "J0 不在 180°")
        self.assertFalse(m.check_pose(UP_Q[:2] + [math.radians(0.6)] + UP_Q[3:], "UP")[0], "J2 超過 0.5°")
        self.assertFalse(m.check_pose([float("nan")] + UP_Q[1:], "UP")[0])
        self.assertFalse(m.check_pose(None, "UP")[0])
        self.assertFalse(m.check_pose(UP_Q, "EXT45")[0], "未知姿態")


class TestURScript(Base):
    def test_structure(self):
        m = load("jv2_script")
        for pose in ("UP", "DIAG"):
            base, end = POSE_Q[pose], m.POSES[pose]["q_deg"][2]
            for d in (+1, -1):
                segs = plan(m, "MID", d, pose)
                text = m.build_urscript(segs, d, 2, base, m.V_CRUISE, m.RAMP_ACCEL, m.MOD_ACCEL_LIMIT,
                                        m.STOPJ_DECEL, m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=end)
                self.assertTrue(urscript_indent_ok(text))
                text.encode("ascii")
                self.assertTrue(text.startswith("def j2_velocity_modulation_program():"))
                for line in text.splitlines():
                    if "speedj(" in line:
                        v = line.split("[")[1].split("]")[0].split(",")
                        self.assertEqual([x.strip() for k, x in enumerate(v) if k != 2], ["0.0"] * 5, line)
                    if "while" in line:
                        self.assertIn(f"q_now[2] {'<' if d > 0 else '>'} ", line)
                movej = [l for l in text.splitlines() if "movej(" in l]
                self.assertEqual(len(movej), len(segs) + 1)
                for line in movej:
                    vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    for j in (0, 1, 3, 4, 5):
                        self.assertAlmostEqual(vals[j], base[j], places=5, msg="只改 J2")
                last = [float(x) for x in movej[-1].split("[")[1].split("]")[0].split(",")]
                self.assertAlmostEqual(math.degrees(last[2]), end, places=4)


class TestMain(Base):
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
                for n in names:
                    with open(os.path.join(root, n), encoding="utf-8") as f:
                        files[n] = f.read()
        self.last_out = out.getvalue()
        return scripts, files

    def test_aborts_before_sending(self):
        cases = [("window", "安全檢查未通過"), ("pose", "姿態檢查未通過"), ("name", "已取消"),
                 ("payload", "負載確認未通過"), ("slider", "速度滑桿檢查未通過"), ("rtde", "RTDE 記錄未能啟動")]
        for k, (case, reason) in enumerate(cases):
            m = load(f"jv2_abort{k}")
            m.J2_SAFE_MIN_DEG, m.J2_SAFE_MAX_DEG = WIN["UP"]
            if case == "window":
                m.J2_SAFE_MIN_DEG = None
            elif case == "pose":
                FakeRTDE.q = list(HOME_Q)
            elif case == "name":
                builtins.input = answers("yes")
            elif case == "payload":
                builtins.input = answers(payload="2.30")
            elif case == "slider":
                FakeRTDE.sf = 0.5
            elif case == "rtde":
                class Failing:
                    def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
                sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
            scripts, files = self._run(m)
            self.assertEqual(scripts, [], case)
            self.assertIn(reason, self.last_out, case)
            FakeRTDE.q, FakeRTDE.sf = list(UP_Q), 1.0
            builtins.input = answers()
            sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE

    def test_quick_end_to_end(self):
        for pose in ("UP", "DIAG"):
            m = load(f"jv2_e2e_{pose}", pose)
            m.J2_SAFE_MIN_DEG, m.J2_SAFE_MAX_DEG = WIN[pose]
            FakeRTDE.q = list(POSE_Q[pose])
            builtins.input = answers(pose)
            for d, label in ((+1, "pos"), (-1, "neg")):
                scripts, files = self._run(m, direction=d)
                self.assertEqual(len(scripts), 1, f"{pose} {label}")
                self.assertTrue(any(n.startswith(f"j2velmod_data_{pose.lower()}_p130_{label}_quick_") for n in files),
                                list(files))
                info = next(v for n, v in files.items() if n.startswith("j2velmod_info_"))
                self.assertIn(f"姿態: {pose}", info)
                self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
                if pose == "UP":
                    self.assertIn("J2 轉回 90°", self.last_out)
                else:
                    self.assertNotIn("[提醒]", self.last_out)


# ------------------------------------------------------------
# 合成資料
# ------------------------------------------------------------

J_TRUE, TC_TRUE, B_TRUE, C_TRUE = 1.3, 8.0, 25.0, -8.0
M_GRAV, Q_ZERO = -23.4, math.radians(-5.73)           # 重力保持力矩 = M·sin(q − q0)：官方模型 J2 重力零點 −5.73°、最大 23.4 N·m
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


def synth_csv(m, path, d, segs, pose="UP", servo_ms=3.0, noise=0.6, slow=0.3, seed=0, h=0.001, model=None):
    """J2 力矩 = J·a + sign(v)·(Tc + B|v| + C v² + 漣波) + M·sin(q − q0) ＋ 雜訊；實際速度比目標晚 servo_ms。
    model = (J_model, B_model) 時另寫 target_qdd_2 與 target_moment_2 = J_model·qdd + B_model·qd + 重力。"""
    rng = np.random.default_rng(seed)
    sub = int(round(m.DT / h))

    def cyc(x):
        return np.r_[x, np.zeros((-len(x)) % sub)]
    n_pause = int(round(m.PAUSE_TIME / m.DT)) * sub
    n_stop = int(math.ceil(m.V_CRUISE / m.STOPJ_DECEL / m.DT)) * sub
    q0 = POSE_Q[pose][2]
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
    grav = M_GRAV * np.sin(q - Q_ZERO)
    tau = J_TRUE * a + fric + grav
    ii = np.arange(0, len(vt), sub)
    n = len(ii)
    ar, rho = np.zeros(n), math.exp(-m.DT / 0.3)
    e = rng.normal(0, slow * math.sqrt(1 - rho ** 2), n)
    for k in range(1, n):
        ar[k] = rho * ar[k - 1] + e[k]
    tau_s = tau[ii] + rng.normal(0, noise, n) + ar
    q_s = np.round((q[ii] + rng.normal(0, 1e-5, n)) / LSB) * LSB
    head = ["timestamp", "actual_q_2", "actual_qd_2", "actual_current_2", "target_qd_2", "joint_temp_2"]
    extra = np.zeros((n, 0))
    if model is not None:
        qdd_t = np.gradient(vt, h)[ii]
        qt = q0 + np.concatenate([[0.0], np.cumsum((vt[1:] + vt[:-1]) / 2 * h)])[ii]
        tmom = model[0] * qdd_t + model[1] * vt[ii] + M_GRAV * np.sin(qt - Q_ZERO)
        head += ["target_qdd_2", "target_moment_2"]
        extra = np.column_stack([qdd_t, tmom])
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(head)
        for k in range(n):
            w.writerow([round(k * m.DT, 6), q_s[k], va[ii[k]], tau_s[k] / m.KT_OUT, vt[ii[k]], 30.0]
                       + list(extra[k]))


class TestAnalysis(unittest.TestCase):
    _cache = {}

    @classmethod
    def result(cls, key, d=+1, pose="UP", stage="FULL", **kw):
        if key not in cls._cache:
            m = load("jv2_ana_" + key, pose)
            segs = plan(m, stage, d, pose)
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                synth_csv(m, p, d, segs, pose=pose, **kw)
                with contextlib.redirect_stdout(io.StringIO()):
                    cls._cache[key] = (m, m.analyze_j2_velocity_modulation(p, 2, m.KT_OUT, d, m.CONDITIONS,
                                                                            m.V_CRUISE))
        return cls._cache[key]

    def test_recovers_J_with_gravity(self):
        for pose in ("UP", "DIAG"):
            for d in (+1, -1):
                m, r = self.result(f"{pose}_{'pos' if d > 0 else 'neg'}", d=d, pose=pose)
                msg = f"{pose} 方向 {d:+d}"
                self.assertEqual(len(r["segments"]), 42, msg)
                # DIAG 正轉：各段都在 90° 以下，最後 movej 回 90° 也是正方向，回程可能被算成「無法配對 1 段」（不影響分析）
                self.assertLessEqual(r["skipped"], 1 if (pose, d) == ("DIAG", +1) else 0, msg)
                self.assertEqual(sorted(r["conditions"]), list(range(6)), msg)
                for ci, c in r["conditions"].items():
                    cm = f"{msg} 條件 {ci + 1}"
                    self.assertLess(abs(c["J"] - J_TRUE), 3 * c["se_J"] + 0.03,
                                    f"{cm}：J = {c['J']:.4f} ± {c['se_J']:.4f}")
                    k = "重力改為各段斜率"
                    self.assertLess(abs(c["J_variants"][k] - J_TRUE), 3 * c["se_J_variants"][k] + 0.05, f"{cm} {k}")
                    self.assertAlmostEqual(c["track_gain"], 1.0, delta=0.01, msg=cm)
                self.assertAlmostEqual(r["J_mean"], J_TRUE, delta=0.03, msg=msg)

    def test_gravity_coupling_without_model(self):
        """不放重力項時：UP 條件 1（1 Hz、0.06）的 J 偏高約 +0.59（K_g ≈ −23.4 N·m/rad，−K_g/ω²）；
        DIAG 前臂接近水平，重力對角度的斜率小，偏差比 UP 小。
        （條件 5 振幅只有一半，沒扣掉的重力變化讓這個變體的標準差很大，不拿來比。）"""
        _, up = self.result("UP_pos", d=+1, pose="UP")
        _, dg = self.result("DIAG_pos", d=+1, pose="DIAG")
        k = "無重力項（只看影響）"
        b_up = up["conditions"][0]["J_variants"][k] - J_TRUE
        b_dg = dg["conditions"][0]["J_variants"][k] - J_TRUE
        self.assertGreater(b_up, 0.4, f"UP 條件 1：偏差 {b_up:.3f}")
        self.assertLess(b_up, 0.8, f"UP 條件 1：偏差 {b_up:.3f}")
        self.assertLess(abs(b_dg), abs(b_up), f"DIAG 條件 1：偏差 {b_dg:.3f}")

    def test_target_moment_model(self):
        _, r = self.result("model", d=-1, pose="UP", model=(2.9, 0.5))
        self.assertEqual(sorted(r["model"]), list(range(6)))
        for ci, mr in r["model"].items():
            self.assertAlmostEqual(mr["J_model"], 2.9, delta=1e-3, msg=f"條件 {ci + 1}")
        self.assertAlmostEqual(r["J_model_mean"], 2.9, delta=1e-3)

    def test_quick_only_lists_segments(self):
        _, r = self.result("quick", stage="QUICK")
        self.assertEqual(len(r["segments"]), 1)
        self.assertEqual(r["conditions"], {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
