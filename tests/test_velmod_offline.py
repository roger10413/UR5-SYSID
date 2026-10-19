# -*- coding: utf-8 -*-
"""
速度調變實驗（robot/ur5_velocity_modulation_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在 J0=180° 不動）、姿態模組設為 None。

測試內容
--------
1. 語法正確。
2. 規劃：三個階段為 1、4、28 段；FULL 條件輪流、起點錯開 (k − 3)/7 個馬達圈（7 的倍數次以外的漣波諧波
   在 7 次重複中平均為 0）；調變段長 V_CRUISE·T_MOD、中心 = 180° + 錯開量；正反轉以 180° 對稱。
3. 安全檢查：三個階段、兩個方向都通過，且規劃範圍在固定角度窗口定速實驗 FULL 的規劃範圍內；
   J0 視窗未設定、視窗太窄、速度超限、speedj 加速度上限小於調變所需、最低速度過低、非整數週期時 FAIL。
4. 起始位置檢查：視窗未設定、目前角度在視窗外、離第一段起點過遠時拒絕。
5. URScript：縮排規則、只含 ASCII、每段兩個有次數上限且有角度保護的迴圈（正轉 "<"、反轉 ">"）、
   movej 只改 J0 且最後回到 180°；照腳本逐週期算出的速度指令 = 巡航 + A·sin(2πft)，
   最低速度不低於下限、加速度低於 speedj 上限、整數週期使調變不改變總位移。
6. main()：視窗未設定、RTDE 失敗、起始位置過遠時都在送出前中止；QUICK 端對端自動命名。
7. 分析（合成資料：J = 1.5、10/15 實測的漣波與起步後行程曲線、白雜訊 + 慢變雜訊、伺服延遲 3 ms）：
   兩個方向都找回 J、B_local、追隨增益與延遲；減速時摩擦較大時各條件 J 依 J − 2c/(πA) 變化、
   J 對 1/A 回歸找回 J 與 c；電流延遲 2 ms 時 J 依 J − B_local·Δ 偏低、由 B_local 隨頻率的變化找回延遲；
   起點錯開使無漣波項模型的 J 不偏（不錯開則 2 Hz 偏低）；QUICK、MID 資料只列各段、不當機。
"""

import os
import re
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
VM = os.path.join(ROBOT, "ur5_velocity_modulation_ident.py")
ANG = os.path.join(ROBOT, "ur5_const_velocity_angle_ident.py")

HOME_Q = [3.14159265, -1.57079633, 1.57079633, -1.57079633, -1.57079633, 0.0]


class FakeRTDE:
    q = list(HOME_Q)
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getJointTemperatures(self): return [30.0] * 6
    def disconnect(self): pass


def load(name, path=VM):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.home_pose = None
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def plan(m, stage, direction, conditions=None, t_mod=None):
    conds, n = m.stage_config(stage)
    return m.plan_segments(conditions or conds, n, direction, m.V_CRUISE, t_mod or m.T_MOD, m.SETTLE_TIME,
                           m.RAMP_ACCEL, m.STOPJ_DECEL, m.OVERRUN_DEG, m.WINDOW_CENTER_DEG, m.GEAR_RATIO, m.DT)


def safety(m, segs, **kw):
    a = dict(v_cruise=m.V_CRUISE, t_mod=m.T_MOD, qd_max=m.QD_MAX, qdd_max=m.QDD_CHECK_LIMIT,
             v_min_allowed=m.V_MIN_ALLOWED, ramp_accel=m.RAMP_ACCEL, mod_accel=m.MOD_ACCEL_LIMIT,
             decel=m.STOPJ_DECEL, movej_vel=m.MOVEJ_VEL, movej_accel=m.MOVEJ_ACCEL, safe_min_deg=65,
             safe_max_deg=300, margin_deg=m.SAFETY_MARGIN_DEG, end_deg=m.END_J0_DEG)
    a.update(kw)
    return m.safety_check(segs, **a)


def build(m, segs, d):
    return m.build_urscript(segs, d, 0, HOME_Q, m.V_CRUISE, m.RAMP_ACCEL, m.MOD_ACCEL_LIMIT, m.STOPJ_DECEL,
                            m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=m.END_J0_DEG)


class Base(unittest.TestCase):
    def setUp(self):
        self._sock = socket.socket
        self._input = builtins.input
        socket.socket = no_network
        builtins.input = lambda *a, **k: "yes"
        FakeRTDE.q = list(HOME_Q)

    def tearDown(self):
        socket.socket = self._sock
        builtins.input = self._input
        FakeRTDE.q = list(HOME_Q)


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


class TestPlan(Base):
    def test_compile(self):
        py_compile.compile(VM, doraise=True)

    def test_stages_order_offsets(self):
        m = load("vm_plan1")
        self.assertEqual([len(plan(m, st, +1)) for st in ("QUICK", "MID", "FULL")], [1, 4, 28])
        segs = plan(m, "FULL", +1)
        self.assertEqual([s["cond"] for s in segs], [0, 1, 2, 3] * 7, "條件應依重複輪流")
        self.assertEqual([s["repeat"] for s in segs], [k for k in range(7) for _ in range(4)])
        rev = 360.0 / m.GEAR_RATIO
        offs = m.start_offsets_deg(7)
        self.assertAlmostEqual(offs[1] - offs[0], rev / 7, places=12)
        self.assertAlmostEqual(sum(offs), 0.0, places=12)
        self.assertEqual(m.start_offsets_deg(1), [0.0])
        for d in (+1, -1):
            segs = plan(m, "FULL", d)
            for ci in range(4):
                phi = np.array([m.GEAR_RATIO * math.radians(s["mod_start_deg"]) for s in segs if s["cond"] == ci])
                for h in (1, 2, 4, 5, 10, 12, 20, 24):
                    self.assertLess(abs(np.mean(np.exp(1j * h * phi))), 1e-9, f"{h} 次漣波應平均掉")
                self.assertAlmostEqual(abs(np.mean(np.exp(7j * phi))), 1.0, places=9)

    def test_geometry_and_symmetry(self):
        m = load("vm_plan2")
        mod_len = math.degrees(m.V_CRUISE * m.T_MOD)
        pre = math.degrees(m.V_CRUISE ** 2 / (2 * m.RAMP_ACCEL) + m.V_CRUISE * m.SETTLE_TIME)
        for d in (+1, -1):
            for s in plan(m, "FULL", d):
                self.assertAlmostEqual(d * (s["mod_end_deg"] - s["mod_start_deg"]), mod_len, places=9)
                self.assertAlmostEqual((s["mod_start_deg"] + s["mod_end_deg"]) / 2,
                                       m.WINDOW_CENTER_DEG + d * s["offset_deg"], places=9)
                self.assertAlmostEqual(d * (s["mod_start_deg"] - s["start_deg"]), pre, places=9)
                self.assertAlmostEqual(s["n_pre"] * m.DT, m.V_CRUISE / m.RAMP_ACCEL + m.SETTLE_TIME, places=9)
                self.assertAlmostEqual(s["n_mod"] * m.DT, m.T_MOD, places=9)
                n_cyc = s["freq"] * s["n_mod"] * m.DT
                self.assertAlmostEqual(n_cyc, round(n_cyc), places=9, msg="調變應為整數個週期")
                self.assertAlmostEqual(d * (s["guard_deg"] - s["mod_end_deg"]), m.OVERRUN_DEG, places=9)
                self.assertGreater(d * (s["end_deg"] - s["mod_end_deg"]), 0)
                self.assertGreater(d * (s["extreme_deg"] - s["guard_deg"]), 0)
        sp, sn = plan(m, "FULL", +1), plan(m, "FULL", -1)
        c = m.WINDOW_CENTER_DEG
        for a, b in zip(sp, sn):
            for key in ("start_deg", "mod_start_deg", "mod_end_deg", "guard_deg", "end_deg", "extreme_deg"):
                self.assertAlmostEqual(a[key] - c, -(b[key] - c), places=9)

    def test_safety(self):
        m = load("vm_safety")
        a = load("vm_safety_ang", ANG)
        env_lo, env_hi = [], []
        for d in (+1, -1):
            levels, order, n_revs = a.stage_config("FULL")
            ta, tb = a.window_bounds_deg(a.WINDOW_CENTER_DEG, n_revs, a.GEAR_RATIO)
            sa = a.plan_segments(a.speed_sequence(levels, order), d, ta, tb, a.SPEEDJ_ACCEL, a.STOPJ_DECEL,
                                 a.SETTLE_TIME, a.OVERRUN_DEG, a.LOOP_TIME_FACTOR, a.DT)
            _, _, (lo, hi) = a.safety_check(sa, a.QD_MAX, a.QDD_CHECK_LIMIT, a.SPEEDJ_ACCEL, a.STOPJ_DECEL,
                                            a.MOVEJ_VEL, a.MOVEJ_ACCEL, 65, 300, a.SAFETY_MARGIN_DEG)
            env_lo.append(lo); env_hi.append(hi)
        for stage in ("QUICK", "MID", "FULL"):
            for d in (+1, -1):
                ok, lines, (lo, hi) = safety(m, plan(m, stage, d))
                self.assertTrue(ok, f"{stage} {d:+d}：" + "\n".join(lines))
                # 已上機的定速實驗 FULL 規劃範圍（兩個方向合併）之內
                self.assertGreaterEqual(lo, min(env_lo), f"{stage} {d:+d}")
                self.assertLessEqual(hi, max(env_hi), f"{stage} {d:+d}")
                ok, _, _ = safety(m, plan(m, stage, d), safe_min_deg=min(env_lo), safe_max_deg=max(env_hi))
                self.assertTrue(ok)
        segs = plan(m, "FULL", +1)
        self.assertFalse(safety(m, segs, safe_min_deg=None, safe_max_deg=None)[0], "J0 視窗未設定時應 FAIL")
        self.assertFalse(safety(m, segs, safe_min_deg=170)[0], "規劃範圍超出 J0 視窗時應 FAIL")
        self.assertFalse(safety(m, segs, qd_max=0.15)[0], "最高速度超限時應 FAIL")
        self.assertFalse(safety(m, segs, mod_accel=1.0)[0], "speedj 加速度上限小於調變所需時應 FAIL")
        self.assertFalse(safety(m, segs, mod_accel=3.5)[0], "speedj 加速度上限超過檢查上限時應 FAIL")
        big = plan(m, "FULL", +1, conditions=[(1.0, 0.10)])
        self.assertFalse(safety(m, big)[0], "最低速度低於下限時應 FAIL")
        odd = plan(m, "FULL", +1, t_mod=5.1)
        self.assertFalse(safety(m, odd, t_mod=5.1)[0], "非整數週期時應 FAIL")

    def test_check_start_position(self):
        m = load("vm_start")
        first = plan(m, "FULL", +1)[0]["start_deg"]
        self.assertFalse(m.check_start_position(180.0, first, 40.0, None, None)[0])
        self.assertTrue(m.check_start_position(180.0, first, 40.0, 65, 300)[0])
        self.assertFalse(m.check_start_position(100.0, first, 40.0, 65, 300)[0], "離第一段起點過遠應拒絕")
        self.assertFalse(m.check_start_position(180.0, first, 40.0, 190, 300)[0], "目前角度在視窗外應拒絕")

    def test_time_estimate(self):
        m = load("vm_time")
        t = {st: m.estimate_total_time(plan(m, st, +1), 180.0, m.STOPJ_DECEL, end_deg=m.END_J0_DEG)
             for st in ("QUICK", "MID", "FULL")}
        self.assertLess(t["QUICK"], t["MID"])
        self.assertLess(t["MID"], t["FULL"])
        self.assertLess(t["FULL"], 400.0, "FULL 每個方向應在 7 分鐘內")
        self.assertGreater(t["FULL"], 28 * (m.T_MOD + 1.0))


class TestURScript(Base):
    def test_structure(self):
        m = load("vm_script")
        for d, op in ((+1, "<"), (-1, ">")):
            segs = plan(m, "FULL", d)
            s = build(m, segs, d)
            self.assertTrue(urscript_indent_ok(s))
            self.assertTrue(all(ord(ch) < 128 for ch in s), "URScript 應只含 ASCII")
            self.assertEqual(s.count("while (n < "), 2 * len(segs), "每段兩個迴圈都要有次數上限")
            self.assertEqual(s.count(f"q_now[0] {op} "), 2 * len(segs), "每段兩個迴圈都要有角度保護")
            self.assertEqual(s.count("stopj("), len(segs))
            self.assertEqual(s.count(f"a={m.MOD_ACCEL_LIMIT}"), len(segs))
            caps = [int(x) for x in re.findall(r"while \(n < (\d+)\)", s)]
            self.assertEqual(caps, [v for sg in segs for v in (sg["n_pre"], sg["n_mod"])])
            guards = [float(x) for x in re.findall(rf"q_now\[0\] {re.escape(op)} (-?\d+\.\d+)\)", s)]
            for g, sg in zip(guards[::2], segs):
                self.assertAlmostEqual(math.degrees(g), sg["guard_deg"], places=4)
            self.assertEqual(guards[::2], guards[1::2])
            # movej 只改 J0：依序到各段起點，最後回到 180°
            j0 = []
            for line in s.splitlines():
                if "movej(" in line:
                    vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    for j in range(1, 6):
                        self.assertAlmostEqual(vals[j], HOME_Q[j], places=5)
                    j0.append(math.degrees(vals[0]))
            self.assertEqual(len(j0), len(segs) + 1)
            for a, sg in zip(j0, segs):
                self.assertAlmostEqual(a, sg["start_deg"], places=4)
            self.assertAlmostEqual(j0[-1], 180.0, places=4)
            self.assertGreater(s.rfind("movej("), s.rfind("stopj("), "回原點應在最後一段 stopj 之後")

    def test_velocity_commands(self):
        """照腳本的 sin 式子逐週期算速度指令，與設計值比較。"""
        m = load("vm_cmd")
        pat = re.compile(r"qd_t = (-?\d+\.\d+) ([+-]) (\d+\.\d+) \* sin\((\d+\.\d+) \* n\)")
        for d in (+1, -1):
            segs = plan(m, "FULL", d)
            mods = pat.findall(build(m, segs, d))
            self.assertEqual(len(mods), len(segs))
            for (c0, pm, amp, w), sg in zip(mods, segs):
                n = np.arange(sg["n_mod"])
                qd = float(c0) + (1.0 if pm == "+" else -1.0) * float(amp) * np.sin(float(w) * n)
                design = d * (m.V_CRUISE + sg["amp"] * np.sin(2 * math.pi * sg["freq"] * m.DT * n))
                self.assertLess(np.max(np.abs(qd - design)), 1e-6)
                self.assertGreaterEqual(np.min(d * qd), m.V_MIN_ALLOWED, "速度不歸零、不反向")
                self.assertLessEqual(np.max(np.abs(qd)), m.QD_MAX)
                self.assertLess(np.max(np.abs(np.diff(qd))) / m.DT, m.MOD_ACCEL_LIMIT)
                self.assertLess(abs(np.sum(qd - d * m.V_CRUISE) * m.DT), 1e-6, "整數週期：調變不改變總位移")


class TestMain(Base):
    def _run(self, m, stage="QUICK", direction=+1):
        sent = []
        m.send_urscript = lambda *a, **k: sent.append(1)
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            m.TEST_STAGE, m.DIRECTION = stage, direction
            with contextlib.redirect_stdout(io.StringIO()):
                m.main()
        return sent

    def test_window_unset_aborts(self):
        m = load("vm_main1")
        m.J0_SAFE_MIN_DEG = m.J0_SAFE_MAX_DEG = None
        self.assertEqual(self._run(m), [], "視窗未設定時不該送出任何軌跡")

    def test_rtde_failure_aborts(self):
        m = load("vm_main2")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m), [], "RTDE 失敗時不該送出任何軌跡")

    def test_far_start_aborts(self):
        m = load("vm_main3")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300
        FakeRTDE.q = [math.radians(100.0)] + HOME_Q[1:]       # J0 離第一段起點很遠
        self.assertEqual(self._run(m), [], "起始位置過遠時不該送出任何軌跡")

    def test_autoname_end_to_end(self):
        m = load("vm_main4")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300
        scripts = []
        m.send_urscript = lambda ip, text, port=30002: scripts.append(text)
        m.wait_for_motion_complete = lambda logger, planned_s, **k: ("done", planned_s, planned_s)
        names = set()
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            for d in (+1, -1):
                m.TEST_STAGE, m.DIRECTION = "QUICK", d
                with contextlib.redirect_stdout(io.StringIO()):
                    m.main()
            for _, _, files in os.walk(tmp):
                names.update(files)
        csvs = sorted(n for n in names if n.endswith(".csv"))
        self.assertEqual(len(csvs), 2, csvs)
        self.assertTrue(any(n.startswith("velmod_data_pos_quick_") for n in csvs))
        self.assertTrue(any(n.startswith("velmod_data_neg_quick_") for n in csvs))
        self.assertEqual(sum(n.startswith("velmod_info_") for n in names), 2)
        self.assertEqual(len(scripts), 2)
        for s in scripts:
            self.assertTrue(s.startswith("def velocity_modulation_program():"))
            self.assertTrue(urscript_indent_ok(s))


# ------------------------------------------------------------
# 合成資料
# ------------------------------------------------------------

J_TRUE, TC_TRUE, B_TRUE, C_TRUE = 1.5, 9.3, 37.4, -16.0
LSB = 1.192e-5                                        # 編碼器解析度 [rad]
# 10/15 殘差分解的漣波（運動方向力矩，馬達角 φ = 101·q）：(次數, cos 係數, sin 係數)
RIPPLE = {+1: [(2, -0.212, 0.579), (4, 0.246, 0.005), (5, -0.164, -1.051), (10, 0.453, 0.077),
               (12, -0.54, 0.738), (20, -0.028, -0.189), (24, 0.105, 0.159)],
          -1: [(2, -0.315, 0.027), (4, -0.07, 0.195), (5, 0.133, 1.068), (10, -0.385, 0.005),
               (12, 0.48, -0.734), (20, -0.006, 0.098), (24, -0.127, -0.046)]}
# 10/15 的起步後行程曲線 [deg] → [N·m]
G_KNOTS = [0, 1, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 30, 40, 50]
G_VALS = [0.52, -0.29, -0.24, -0.21, -0.21, -0.14, -0.06, 0.0, -0.02, 0.01, 0.03, 0.09, 0.12, 0.23, 0.21]


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


def synth_csv(m, path, d, segs, servo_ms=3.0, delay_ms=0.0, dec_extra=0.0, noise=0.4, slow=0.2, drift=0.3,
              seed=0, h=0.001):
    """依 segs 產生 8 ms 取樣的 CSV。運動方向力矩 = J·a + Tc + B v + C v² + 漣波 + 起步後行程 + 漂移
    （+ 減速時多 dec_extra）＋白雜訊與慢變雜訊；實際速度比目標晚 servo_ms；電流可再晚 delay_ms。"""
    rng = np.random.default_rng(seed)
    sub = int(round(m.DT / h))

    def cyc(x):                                       # 補零到整數個控制週期
        return np.r_[x, np.zeros((-len(x)) % sub)]
    n_pause = int(round(m.PAUSE_TIME / m.DT)) * sub
    n_stop = int(math.ceil(m.V_CRUISE / m.STOPJ_DECEL / m.DT)) * sub
    pieces, pos = [], HOME_Q[0]
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
    pieces.append(cyc(movej_profile(math.radians(m.END_J0_DEG) - pos, m.MOVEJ_VEL, m.MOVEJ_ACCEL, h)))
    pieces.append(np.zeros(125 * sub))
    vt = np.concatenate(pieces)
    tf = np.arange(len(vt)) * h
    va = np.interp(tf - servo_ms / 1000.0, tf, vt)
    q = HOME_Q[0] + np.concatenate([[0.0], np.cumsum((va[1:] + va[:-1]) / 2 * h)])
    a = np.gradient(va, h)
    sv = np.sign(np.where(np.abs(va) > 1e-6, va, 0.0))
    moving = sv != 0
    start = np.flatnonzero(moving & ~np.r_[False, moving[:-1]])
    seg_id = np.cumsum(np.isin(np.arange(len(va)), start))
    travel = np.degrees(np.abs(q - q[start][np.maximum(seg_id - 1, 0)]))
    phi = m.GEAR_RATIO * q
    tau_m = np.zeros(len(va))
    for sg in (+1, -1):
        mm = sv == sg
        rip = sum(c * np.cos(k * phi[mm]) + s_ * np.sin(k * phi[mm]) for k, c, s_ in RIPPLE[sg])
        acc_m = sg * a[mm]
        vm = np.abs(va[mm])
        tau_m[mm] = (J_TRUE * acc_m + TC_TRUE + B_TRUE * vm + C_TRUE * vm ** 2 + rip
                     + np.interp(travel[mm], G_KNOTS, G_VALS) + dec_extra * (acc_m < -1e-6) + drift * tf[mm] / tf[-1])
    tau = sv * tau_m
    tau[~moving] = J_TRUE * a[~moving]
    if delay_ms:
        tau = np.interp(tf - delay_ms / 1000.0, tf, tau)
    ii = np.arange(0, len(vt), sub)
    n = len(ii)
    ar, rho = np.zeros(n), math.exp(-m.DT / 0.3)
    e = rng.normal(0, slow * math.sqrt(1 - rho ** 2), n)
    for k in range(1, n):
        ar[k] = rho * ar[k - 1] + e[k]
    tau_s = tau[ii] + rng.normal(0, noise, n) + ar * moving[ii]
    q_s = np.round((q[ii] + rng.normal(0, 1e-5, n)) / LSB) * LSB
    temp = 30.0 + 2.0 * tf[ii] / tf[-1]
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(["timestamp", "actual_q_0", "actual_qd_0", "actual_current_0", "target_qd_0", "joint_temp_0"])
        for k in range(n):
            w.writerow([round(k * m.DT, 6), q_s[k], va[ii[k]], tau_s[k] / m.KT_OUT, vt[ii[k]], temp[k]])


class TestAnalysis(unittest.TestCase):
    """每種情境只產生、分析一次（各約 3 秒），結果在測試間共用。"""
    _cache = {}

    @classmethod
    def result(cls, key, d=+1, stage="FULL", no_offsets=False, **kw):
        if key not in cls._cache:
            m = load("vm_ana_" + key)
            segs = plan(m, stage, d)
            if no_offsets:                            # 對照：每次同一起點（條件仍輪流）
                conds, n = m.stage_config(stage)
                base = m.plan_segments(conds, 1, d, m.V_CRUISE, m.T_MOD, m.SETTLE_TIME, m.RAMP_ACCEL,
                                       m.STOPJ_DECEL, m.OVERRUN_DEG, m.WINDOW_CENTER_DEG, m.GEAR_RATIO, m.DT)
                segs = [dict(s, index=i) for i, s in enumerate(base * n)]
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                synth_csv(m, p, d, segs, **kw)
                with contextlib.redirect_stdout(io.StringIO()):
                    cls._cache[key] = (m, m.analyze_velocity_modulation(p, 0, m.KT_OUT, d, m.CONDITIONS,
                                                                          m.V_CRUISE))
        return cls._cache[key]

    def test_recovers_J_both_directions(self):
        B_loc = B_TRUE + 2 * C_TRUE * 0.12
        for d in (+1, -1):
            m, r = self.result("pos" if d > 0 else "neg", d=d)
            msg = f"方向 {d:+d}"
            self.assertEqual(len(r["segments"]), 28, msg)
            self.assertEqual(r["skipped"], 0, msg)
            self.assertEqual(sorted(r["conditions"]), [0, 1, 2, 3], msg)
            for sg in r["segments"]:
                self.assertAlmostEqual(sg["n_cycles"], m.CONDITIONS[sg["cond"]][0] * m.T_MOD, delta=0.05, msg=msg)
            for ci, c in r["conditions"].items():
                f, amp = m.CONDITIONS[ci]
                cm = f"{msg} 條件 {ci + 1}"
                self.assertEqual(c["n_seg"], 7, cm)
                self.assertAlmostEqual(c["freq_meas"], f, delta=1e-3 * f, msg=cm)
                self.assertAlmostEqual(c["a_amp"], 2 * math.pi * f * amp, delta=0.01 * 2 * math.pi * f * amp, msg=cm)
                self.assertLess(abs(c["J"] - J_TRUE), 3 * c["se_J"] + 0.02, f"{cm}：J = {c['J']:.4f} ± {c['se_J']:.4f}")
                self.assertLess(abs(c["B_local"] - B_loc), 3 * c["se_B_local"] + 0.1, cm)
                self.assertAlmostEqual(c["track_gain"], 1.0, delta=0.01, msg=cm)
                self.assertAlmostEqual(c["track_lag_ms"], 3.0, delta=0.3, msg=cm)
            self.assertAlmostEqual(r["J_mean"], J_TRUE, delta=0.03, msg=msg)
            self.assertLess(abs(r["k_invA"]), 3 * r["se_k_invA"], msg)
            self.assertLess(abs(r["delay_ms"]), 3 * r["se_delay_ms"], msg)

    def test_dec_friction(self):
        """減速時摩擦多 c = 0.3 N·m：各條件 J ≈ J − 2c/(πA)，J 對 1/A 回歸找回 J 與 c。"""
        c_true = 0.3
        m, r = self.result("dec", dec_extra=c_true)
        for ci, c in r["conditions"].items():
            expect = J_TRUE - 2 * c_true / (math.pi * c["a_amp"])
            self.assertLess(abs(c["J"] - expect), 3 * c["se_J"] + 0.03, f"條件 {ci + 1}：{c['J']:.4f} vs {expect:.4f}")
        self.assertGreater(r["J_chi2"], 10 * r["J_dof"], "各條件的 J 應明顯不一致")
        self.assertLess(abs(r["J_inf"] - J_TRUE), 3 * r["se_J_inf"])
        self.assertLess(r["k_invA"], -3 * r["se_k_invA"])
        self.assertAlmostEqual(-math.pi * r["k_invA"] / 2, c_true, delta=0.1)

    def test_current_delay(self):
        """電流比編碼器晚 2 ms：J 偏低約 B_local·Δ，B_local 對 ω² 的斜率找回延遲。"""
        m, r = self.result("delay", delay_ms=2.0)
        B_loc = B_TRUE + 2 * C_TRUE * 0.12
        self.assertAlmostEqual(r["J_mean"], J_TRUE - B_loc * 0.002, delta=0.03)
        self.assertLess(abs(r["delay_ms"] - 2.0), 3 * r["se_delay_ms"])

    def test_offsets_protect_against_ripple(self):
        """不放漣波項時：起點錯開的 2 Hz 條件 J 不偏；每次同一起點則偏離（漣波基頻約 1.93 Hz）。"""
        k = "無漣波項"
        _, r_off = self.result("pos", d=+1)
        _, r_no = self.result("nooff", no_offsets=True)
        for ci in (1, 3):
            c = r_off["conditions"][ci]
            self.assertLess(abs(c["J_variants"][k] - J_TRUE), 3 * c["se_J_variants"][k] + 0.02, f"條件 {ci + 1}")
        worst = max(abs(r_no["conditions"][ci]["J_variants"][k] - J_TRUE) for ci in (1, 3))
        self.assertGreater(worst, 0.08)

    def test_quick_and_mid_only_list_segments(self):
        for stage, n_seg in (("QUICK", 1), ("MID", 4)):
            m, r = self.result(stage.lower(), stage=stage)
            self.assertEqual([sg["cond"] for sg in r["segments"]], list(range(n_seg)), stage)
            self.assertEqual(r["conditions"], {}, stage)
            for sg in r["segments"]:
                self.assertAlmostEqual(sg["track_gain"], 1.0, delta=0.01)


if __name__ == "__main__":
    unittest.main(verbosity=2)
