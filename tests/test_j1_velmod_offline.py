# -*- coding: utf-8 -*-
"""
J1 速度調變實驗（robot/ur5_j1_velocity_modulation_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在 UP 姿態不動）。

測試內容
--------
1. 語法正確。
2. 規劃：三個階段為 1、6、42 段；調變段中心 = −92° + 錯開量；兩個方向、三個階段的規劃範圍都在
   J1 定速 FULL 的範圍（−112.71°～−71.29°，含餘裕）內，也在 J1 定速用的視窗 −114°～−70° 內。
3. 安全檢查：通過；J1 視窗未設定、太窄、最低速度過低、非整數週期、法蘭線速度超限時 FAIL。
4. 姿態：只接受 UP 六軸 ±0.5°；診斷原點、J0 不在 180°、NaN、未知姿態一律拒絕。
5. URScript：縮排規則、只含 ASCII、speedj 只有 J1 非零、movej 只改 J1、最後回到 −90°、角度保護方向正確。
6. main()：視窗未設定、姿態不符、未輸入姿態名稱、負載不符、速度滑桿不是 100% 時都在送出前中止；
   QUICK 端對端自動命名（j1velmod_data_up_p130_pos_quick_…），資料存在 OUTPUT_DIR。
7. 分析（合成資料：J = 6.0、重力 −66.8·sin(q + 92°) N·m（手臂朝上像倒立的擺，正負號同 11/5 實測）、
   J0 實測的漣波、白雜訊 + 慢變雜訊、伺服延遲 3 ms）：
   兩個方向各條件都找回 J；「重力改為各段斜率」也找回 J；「無重力項」時 1 Hz 的 J 明顯偏高（多 1.0 以上）；
   target_moment（含重力）回歸找回 J_model；QUICK 只列各段、不當機。
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
JV = os.path.join(ROBOT, "ur5_j1_velocity_modulation_ident.py")

UP_Q = [math.radians(x) for x in (180.0, -90.0, 0.0, -90.0, -90.0, 0.0)]
HOME_Q = [math.radians(x) for x in (180.0, -90.0, 90.0, -90.0, -90.0, 0.0)]
PAYLOAD = 1.3
WIN = (-114.0, -70.0)                                 # J1 定速用的視窗
J1CV_FULL = (-112.71, -71.29)                         # J1 定速 FULL 的規劃範圍（含餘裕）


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


def load(name):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    spec = importlib.util.spec_from_file_location(name, JV)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PAYLOAD_KG = PAYLOAD
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(default="UP", payload=str(PAYLOAD)):
    return lambda prompt="", *a, **k: payload if "負載質量" in str(prompt) else default


def plan(m, stage, d):
    conds, n = m.stage_config(stage)
    return m.plan_segments(conds, n, d, m.V_CRUISE, m.T_MOD, m.SETTLE_TIME, m.RAMP_ACCEL, m.STOPJ_DECEL,
                           m.OVERRUN_DEG, m.POSES["UP"]["center_deg"], m.GEAR_RATIO, m.DT)


def safety(m, segs, **kw):
    a = dict(v_cruise=m.V_CRUISE, t_mod=m.T_MOD, qd_max=m.QD_MAX, qdd_max=m.QDD_CHECK_LIMIT,
             v_min_allowed=m.V_MIN_ALLOWED, ramp_accel=m.RAMP_ACCEL, mod_accel=m.MOD_ACCEL_LIMIT,
             decel=m.STOPJ_DECEL, movej_vel=m.MOVEJ_VEL, movej_accel=m.MOVEJ_ACCEL, safe_min_deg=WIN[0],
             safe_max_deg=WIN[1], margin_deg=m.SAFETY_MARGIN_DEG, q_pose_deg=m.POSES["UP"]["q_deg"], end_deg=-90.0)
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

    def test_stages_and_ranges(self):
        m = load("jv_plan")
        for stage, n_seg in (("QUICK", 1), ("MID", 6), ("FULL", 42)):
            for d in (+1, -1):
                segs = plan(m, stage, d)
                self.assertEqual(len(segs), n_seg)
                for s in segs:
                    mid = 0.5 * (s["mod_start_deg"] + s["mod_end_deg"])
                    self.assertAlmostEqual(mid, -92.0 + d * s["offset_deg"], places=9)
                    self.assertAlmostEqual(abs(s["mod_end_deg"] - s["mod_start_deg"]),
                                           math.degrees(m.V_CRUISE * m.T_MOD), places=9)
                ok, lines, (lo, hi) = safety(m, segs)
                self.assertTrue(ok, "\n".join(lines))
                self.assertGreaterEqual(lo, J1CV_FULL[0], f"{stage} {d:+d}")
                self.assertLessEqual(hi, J1CV_FULL[1], f"{stage} {d:+d}")
                first = abs(segs[0]["start_deg"] - (-90.0))
                self.assertLess(first, m.MAX_INITIAL_MOVE_DEG)

    def test_safety_failures(self):
        m = load("jv_safety")
        segs = plan(m, "FULL", +1)
        self.assertFalse(safety(m, segs, safe_min_deg=None)[0], "視窗未設定")
        self.assertFalse(safety(m, segs, safe_min_deg=-100.0)[0], "視窗太窄")
        self.assertFalse(safety(m, segs, v_min_allowed=0.05)[0], "最低速度過低")
        self.assertFalse(safety(m, segs, t_mod=3.3)[0], "非整數週期")
        self.assertFalse(safety(m, segs, flange_speed_max=0.05)[0], "法蘭線速度超限")
        self.assertFalse(safety(m, segs, mod_accel=0.5)[0], "speedj 加速度上限小於調變所需")
        self.assertFalse(m.check_start_position(-90.0, -104.4, 20.0, None, None)[0])
        self.assertFalse(m.check_start_position(-60.0, -104.4, 40.0, *WIN)[0], "目前角度在視窗外")
        self.assertFalse(m.check_start_position(-90.0, -104.4, 10.0, *WIN)[0], "第一次 movej 過遠")
        self.assertTrue(m.check_start_position(-90.0, -104.4, 20.0, *WIN)[0])

    def test_check_pose(self):
        m = load("jv_pose")
        self.assertTrue(m.check_pose(UP_Q, "UP")[0])
        self.assertTrue(m.check_pose([UP_Q[0] + 2 * math.pi] + UP_Q[1:], "UP")[0], "±360° 等價")
        self.assertFalse(m.check_pose(HOME_Q, "UP")[0], "診斷原點不是 UP")
        self.assertFalse(m.check_pose([math.radians(170.0)] + UP_Q[1:], "UP")[0], "J0 不在 180°")
        self.assertFalse(m.check_pose(UP_Q[:1] + [UP_Q[1] + math.radians(0.6)] + UP_Q[2:], "UP")[0], "超過 0.5°")
        self.assertFalse(m.check_pose([float("nan")] + UP_Q[1:], "UP")[0])
        self.assertFalse(m.check_pose(None, "UP")[0])
        self.assertFalse(m.check_pose(UP_Q, "DIAG")[0], "未知姿態")


class TestURScript(Base):
    def test_structure(self):
        m = load("jv_script")
        for d in (+1, -1):
            segs = plan(m, "MID", d)
            text = m.build_urscript(segs, d, 1, UP_Q, m.V_CRUISE, m.RAMP_ACCEL, m.MOD_ACCEL_LIMIT, m.STOPJ_DECEL,
                                    m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=-90.0)
            self.assertTrue(urscript_indent_ok(text))
            text.encode("ascii")
            self.assertTrue(text.startswith("def j1_velocity_modulation_program():"))
            for line in text.splitlines():
                if "speedj(" in line:
                    v = line.split("[")[1].split("]")[0].split(",")
                    self.assertEqual([x.strip() for k, x in enumerate(v) if k != 1], ["0.0"] * 5, line)
                if "while" in line:
                    self.assertIn(f"q_now[1] {'<' if d > 0 else '>'} ", line)
            movej = [l for l in text.splitlines() if "movej(" in l]
            self.assertEqual(len(movej), len(segs) + 1)
            for line in movej:
                vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                for j in (0, 2, 3, 4, 5):
                    self.assertAlmostEqual(vals[j], UP_Q[j], places=5, msg="只改 J1")
            self.assertAlmostEqual(math.degrees(float(movej[-1].split("[")[1].split(",")[1])), -90.0, places=4)


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
            m = load(f"jv_abort{k}")
            m.J1_SAFE_MIN_DEG, m.J1_SAFE_MAX_DEG = WIN
            if case == "window":
                m.J1_SAFE_MIN_DEG = None
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
        m = load("jv_e2e")
        m.J1_SAFE_MIN_DEG, m.J1_SAFE_MAX_DEG = WIN
        for d, label in ((+1, "pos"), (-1, "neg")):
            scripts, files = self._run(m, direction=d)
            self.assertEqual(len(scripts), 1)
            self.assertTrue(any(n.startswith(f"j1velmod_data_up_p130_{label}_quick_") for n in files), list(files))
            info = next(v for n, v in files.items() if n.startswith("j1velmod_info_"))
            self.assertIn("姿態: UP", info)
            self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
            self.assertIn("請以教導器把 J2 轉回 90°", self.last_out)


# ------------------------------------------------------------
# 合成資料
# ------------------------------------------------------------

J_TRUE, TC_TRUE, B_TRUE, C_TRUE = 6.0, 11.0, 30.0, -10.0
M_GRAV, Q_ZERO = -66.8, math.radians(-92.0)           # 重力力矩 = M·sin(q − q0)：UP 像倒立的擺，斜率約 −67 N·m/rad
#                                                       （11/5 實測 target_moment −1.165 N·m/° = −66.8 N·m/rad）
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


def synth_csv(m, path, d, segs, servo_ms=3.0, noise=0.6, slow=0.3, seed=0, h=0.001, model=None):
    """J1 力矩 = J·a + sign(v)·(Tc + B|v| + C v² + 漣波) + M·sin(q − q0) ＋ 雜訊；實際速度比目標晚 servo_ms。
    model = (J_model, B_model) 時另寫 target_qdd_1 與 target_moment_1 = J_model·qdd + B_model·qd + 重力。"""
    rng = np.random.default_rng(seed)
    sub = int(round(m.DT / h))

    def cyc(x):
        return np.r_[x, np.zeros((-len(x)) % sub)]
    n_pause = int(round(m.PAUSE_TIME / m.DT)) * sub
    n_stop = int(math.ceil(m.V_CRUISE / m.STOPJ_DECEL / m.DT)) * sub
    q0 = UP_Q[1]
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
    head = ["timestamp", "actual_q_1", "actual_qd_1", "actual_current_1", "target_qd_1", "joint_temp_1"]
    extra = np.zeros((n, 0))
    if model is not None:
        qdd_t = np.gradient(vt, h)[ii]
        qt = q0 + np.concatenate([[0.0], np.cumsum((vt[1:] + vt[:-1]) / 2 * h)])[ii]
        tmom = model[0] * qdd_t + model[1] * vt[ii] + M_GRAV * np.sin(qt - Q_ZERO)
        head += ["target_qdd_1", "target_moment_1"]
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
    def result(cls, key, d=+1, stage="FULL", **kw):
        if key not in cls._cache:
            m = load("jv_ana_" + key)
            segs = plan(m, stage, d)
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                synth_csv(m, p, d, segs, **kw)
                with contextlib.redirect_stdout(io.StringIO()):
                    cls._cache[key] = (m, m.analyze_j1_velocity_modulation(p, 1, m.KT_OUT, d, m.CONDITIONS,
                                                                            m.V_CRUISE))
        return cls._cache[key]

    def test_recovers_J_with_gravity(self):
        for d in (+1, -1):
            m, r = self.result("pos" if d > 0 else "neg", d=d)
            msg = f"方向 {d:+d}"
            self.assertEqual(len(r["segments"]), 42, msg)
            self.assertEqual(r["skipped"], 0, msg)
            self.assertEqual(sorted(r["conditions"]), list(range(6)), msg)
            for ci, c in r["conditions"].items():
                cm = f"{msg} 條件 {ci + 1}"
                self.assertLess(abs(c["J"] - J_TRUE), 3 * c["se_J"] + 0.05, f"{cm}：J = {c['J']:.4f} ± {c['se_J']:.4f}")
                k = "重力改為各段斜率"
                self.assertLess(abs(c["J_variants"][k] - J_TRUE), 3 * c["se_J_variants"][k] + 0.08, f"{cm} {k}")
                self.assertAlmostEqual(c["track_gain"], 1.0, delta=0.01, msg=cm)
            self.assertAlmostEqual(r["J_mean"], J_TRUE, delta=0.05, msg=msg)

    def test_gravity_coupling_without_model(self):
        """不放重力項時，1 Hz 的 J 明顯偏高（重力對角度的斜率 K_g ≈ −66.8 N·m/rad，角度振盪的部分約
        −K_g/ω² ≈ +1.7；段內重力的慢變化也混入；11/5 實測 +2.1～+2.4）：重力項不可省略。"""
        m, r = self.result("pos", d=+1)
        for ci in (0, 4):                                   # 1 Hz 的兩個條件
            bias = r["conditions"][ci]["J_variants"]["無重力項（只看影響）"] - J_TRUE
            self.assertGreater(bias, 1.0, f"條件 {ci + 1}：偏差 {bias:.3f}")

    def test_target_moment_model(self):
        _, r = self.result("model", d=-1, model=(8.1, 0.5))
        self.assertEqual(sorted(r["model"]), list(range(6)))
        for ci, mr in r["model"].items():
            self.assertAlmostEqual(mr["J_model"], 8.1, delta=1e-3, msg=f"條件 {ci + 1}")
        self.assertAlmostEqual(r["J_model_mean"], 8.1, delta=1e-3)

    def test_quick_only_lists_segments(self):
        _, r = self.result("quick", stage="QUICK")
        self.assertEqual(len(r["segments"]), 1)
        self.assertEqual(r["conditions"], {})


if __name__ == "__main__":
    unittest.main(verbosity=2)
