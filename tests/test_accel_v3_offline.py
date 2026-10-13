# -*- coding: utf-8 -*-
"""
定加速度 v3（ur5_const_accel_ident_v3.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式、假的 rtde_receive（手臂停在 J0=180°）、
姿態模組設為 None。腳本位置自動判斷：有 robot/ 資料夾就用 robot/，否則用專案根目錄。

測試內容
--------
1. 語法正確。
2. 三階段安全檢查通過；峰值速度、加速度超限時 FAIL；起始錯開剛好涵蓋 1 個馬達圈。
3. URScript：縮排規則、只含 ASCII、每次循環兩個 speedj（加速＋定速、減速）、movej 只改 J0 且依次錯開。
4. main()：J0 視窗未設定、RTDE 失敗時在送出前中止；QUICK 端對端自動命名。
5. 分析：合成資料（已知 J、B、Tc，加上鎖定馬達角度的漣波與雜訊）能找回 J、B、Tc（定速段 4 圈與 0 圈）；
   「起始角錯開」比「每次同一起點」的 J 誤差小；換階段後力矩有延遲時，跳過前 40 ms 可去除偏差；
   減速段摩擦比加速段多一個常數時（10/13 實測現象），含截距回歸仍找回 J，且截距被標示為顯著。
6. 起點固定 180°：目前 J0 在 170° 時第一個 movej 到 180°；離起點超過 40° 時中止、不送出。
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
V3 = os.path.join(ROBOT, "ur5_const_accel_ident_v3.py")

HOME_Q = [3.14159265, -1.57079633, 1.57079633, -1.57079633, -1.57079633, 0.0]


class FakeRTDE:
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(HOME_Q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getJointTemperatures(self): return [30.0] * 6
    def disconnect(self): pass


def load(name):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, V3)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.home_pose = None
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


class Base(unittest.TestCase):
    def setUp(self):
        self._sock = socket.socket
        self._input = builtins.input
        socket.socket = no_network
        builtins.input = lambda *a, **k: "yes"

    def tearDown(self):
        socket.socket = self._sock
        builtins.input = self._input


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


class TestPlan(Base):
    def test_compile(self):
        py_compile.compile(V3, doraise=True)

    def test_safety_and_offsets(self):
        m = load("v3_safe")
        th = m.hold_time(m.V_PEAK, m.HOLD_MOTOR_REVS)
        for stage in ("QUICK", "MID", "FULL"):
            acc, n = m.stage_config(stage)
            ok, lines, _ = m.safety_check(acc, m.V_PEAK, n, th, m.QD_MAX, m.QDD_CHECK_LIMIT,
                                          m.MAX_EXCURSION_DEG, m.EXCURSION_MARGIN_FACTOR,
                                          m.MOVEJ_RETURN_ACCEL, m.MOVEJ_RETURN_VEL)
            self.assertTrue(ok, f"{stage}\n" + "\n".join(lines))
        acc, n = m.stage_config("FULL")
        self.assertFalse(m.safety_check(acc, 0.35, n, th, m.QD_MAX, m.QDD_CHECK_LIMIT, m.MAX_EXCURSION_DEG,
                                        1.0, m.MOVEJ_RETURN_ACCEL, m.MOVEJ_RETURN_VEL)[0], "峰值速度超限應 FAIL")
        self.assertFalse(m.safety_check([3.5], m.V_PEAK, n, th, m.QD_MAX, m.QDD_CHECK_LIMIT, m.MAX_EXCURSION_DEG,
                                        1.0, m.MOVEJ_RETURN_ACCEL, m.MOVEJ_RETURN_VEL)[0], "加速度超限應 FAIL")
        offs = m.start_offsets_rad(20)
        step = 2 * math.pi / 101 / 20
        self.assertEqual(len(offs), 20)
        self.assertAlmostEqual(offs[1] - offs[0], step, places=12)
        self.assertAlmostEqual(offs[-1] + step, 2 * math.pi / 101, places=12)
        self.assertAlmostEqual(th * m.V_PEAK, m.HOLD_MOTOR_REVS * 2 * math.pi / 101, places=12)

    def test_urscript(self):
        m = load("v3_script")
        th = m.hold_time(m.V_PEAK, m.HOLD_MOTOR_REVS)
        for d in (+1, -1):
            acc, n = m.stage_config("FULL")
            s = m.build_full_urscript(acc, m.V_PEAK, n, th, d, 0, HOME_Q, m.MOVEJ_RETURN_ACCEL,
                                      m.MOVEJ_RETURN_VEL, m.PAUSE_TIME)
            self.assertTrue(urscript_indent_ok(s))
            self.assertTrue(all(ord(ch) < 128 for ch in s))
            n_cyc = len(acc) * n
            self.assertEqual(s.count("speedj("), 2 * n_cyc)
            self.assertEqual(s.count("movej("), n_cyc + 1)
            j0 = []
            for line in s.splitlines():
                if "movej(" in line:
                    vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    for j in range(1, 6):
                        self.assertAlmostEqual(vals[j], HOME_Q[j], places=5)
                    j0.append(vals[0])
            first = np.array(j0[:n]) - HOME_Q[0]
            self.assertTrue(np.all(np.diff(d * first) > 0), "起始角應沿運動方向依次錯開")
            self.assertAlmostEqual(j0[-1], HOME_Q[0], places=5)


class TestStart(Base):
    def test_check_initial_move(self):
        m = load("v3_start")
        self.assertFalse(m.check_initial_move(170, 180, 40, None, None)[0], "視窗未設定應拒絕")
        self.assertTrue(m.check_initial_move(170, 180, 40, 65, 300)[0])
        self.assertFalse(m.check_initial_move(130, 180, 40, 65, 300)[0], "距離超過上限應拒絕")
        self.assertFalse(m.check_initial_move(60, 80, 40, 65, 300)[0], "目前角度在視窗外應拒絕")

    def _main_with_q0(self, name, q0_deg):
        m = load(name)
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300
        q = list(HOME_Q); q[0] = math.radians(q0_deg)
        FakeRTDE.getActualQ = lambda self: list(q)
        sent = []
        m.send_urscript = lambda ip, text, *a, **k: sent.append(text)
        m.wait_for_motion_complete = lambda logger, planned_s, **k: ("done", planned_s, planned_s)
        try:
            with tempfile.TemporaryDirectory() as tmp:
                m.OUTPUT_DIR, m.TEST_STAGE, m.DIRECTION = tmp, "QUICK", +1
                with contextlib.redirect_stdout(io.StringIO()):
                    m.main()
        finally:
            FakeRTDE.getActualQ = lambda self: list(HOME_Q)
        return sent

    def test_moves_to_fixed_start(self):
        sent = self._main_with_q0("v3_start2", 170.0)
        self.assertEqual(len(sent), 1)
        movej = [ln for ln in sent[0].splitlines() if "movej(" in ln]
        first = float(movej[0].split("[")[1].split(",")[0]); last = float(movej[-1].split("[")[1].split(",")[0])
        self.assertAlmostEqual(math.degrees(first), 180.0, places=4, msg="第一次循環應從 180° 出發")
        self.assertAlmostEqual(math.degrees(last), 180.0, places=4, msg="結束應回到 180°")

    def test_far_start_aborts(self):
        self.assertEqual(self._main_with_q0("v3_start3", 130.0), [], "離起點 50° 應中止、不送出")


class TestMain(Base):
    def _run(self, m, stage="QUICK", d=+1):
        sent = []
        m.send_urscript = lambda *a, **k: sent.append(1)
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR, m.TEST_STAGE, m.DIRECTION = tmp, stage, d
            with contextlib.redirect_stdout(io.StringIO()):
                m.main()
        return sent

    def test_window_unset_aborts(self):
        m = load("v3_m1")
        m.J0_SAFE_MIN_DEG = m.J0_SAFE_MAX_DEG = None
        self.assertEqual(self._run(m), [])

    def test_rtde_failure_aborts(self):
        m = load("v3_m2")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m), [])

    def test_autoname_end_to_end(self):
        m = load("v3_m3")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300
        m.send_urscript = lambda *a, **k: None
        m.wait_for_motion_complete = lambda logger, planned_s, **k: ("done", planned_s, planned_s)
        names = set()
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            for d in (+1, -1):
                m.TEST_STAGE, m.DIRECTION = "QUICK", d
                with contextlib.redirect_stdout(io.StringIO()):
                    m.main()
                time.sleep(1.1)
            for _, _, fs in os.walk(tmp):
                names.update(fs)
        csvs = [n for n in names if n.endswith(".csv")]
        self.assertEqual(len(csvs), 2, csvs)
        self.assertTrue(any("_pos_quick_" in n for n in csvs))
        self.assertTrue(any("_neg_quick_" in n for n in csvs))


def synth_csv(m, path, d, J0, B0, TC0, use_offsets=True, seed=0, lag_s=0.0, dec_extra=0.0):
    """依 v3 軌跡產生合成資料：理想追隨（qd = target），tau = J·a + B·v + Tc + 漣波(φ) + 雜訊。
    lag_s > 0 時，慣性力矩 J·a 以一階延遲（時間常數 lag_s）跟隨，模擬換階段後力矩建立的暫態。
    dec_extra：減速段摩擦額外增加的量 [N·m]（10/13 實測約 1.0），模擬加減速段摩擦不同。"""
    rng = np.random.default_rng(seed)
    ja = 0.0
    acc, n = m.stage_config("FULL")
    th = m.hold_time(m.V_PEAK, m.HOLD_MOTOR_REVS)
    offs = m.start_offsets_rad(n) if use_offsets else [0.0] * n
    q0, dt, t = HOME_Q[0], m.DT, 0.0
    rows = []

    def emit(q, v, a, extra=0.0):
        nonlocal t, ja
        ja = J0 * a if lag_s <= 0 else ja + (J0 * a - ja) * dt / (lag_s + dt)
        phi = (q * m.GEAR_RATIO) % (2 * math.pi)
        ripple = 1.2 * math.cos(5 * phi + 0.7) + 0.8 * math.sin(12 * phi) + 0.4 * math.cos(2 * phi)
        fr = (B0 * abs(v) + TC0 + extra) * (1 if v > 1e-9 else (-1 if v < -1e-9 else 0))
        tau = ja + fr + (ripple if abs(v) > 1e-9 else 0.0) + rng.normal(0, 0.4)
        rows.append({"timestamp": t, "actual_q_0": q, "actual_qd_0": v + rng.normal(0, 0.002),
                     "actual_current_0": tau / m.KT_OUT, "target_qd_0": v})
        t += dt

    for a in acc:
        for off in offs:
            q = q0 + d * off
            for _ in range(int(0.3 / dt)):
                emit(q, 0.0, 0.0)
            v = 0.0
            while v < m.V_PEAK - 1e-12:                      # 加速
                v = min(v + a * dt, m.V_PEAK); q += d * v * dt; emit(q, d * v, d * a)
            for _ in range(int(round(th / dt))):             # 定速
                q += d * v * dt; emit(q, d * v, 0.0)
            while v > 1e-12:                                 # 同 a 減速
                v = max(v - a * dt, 0.0); q += d * v * dt; emit(q, d * v, -d * a, dec_extra)
            for _ in range(int(0.3 / dt)):
                emit(q, 0.0, 0.0)
            back = q - (q0 + d * off)                        # movej 回程（反方向，不參與分析）
            nb = int(abs(back) / (0.15 * dt)) + 1
            for _ in range(nb):
                q -= back / nb; emit(q, -d * 0.15, 0.0)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


class TestAnalysis(unittest.TestCase):
    def test_recovers_J_B_Tc(self):
        m = load("v3_ana")
        J0, B0, TC0 = 2.0, 32.5, 8.9
        for hold in (4, 0):
            m.HOLD_MOTOR_REVS = hold
            for d in (+1, -1):
                with tempfile.TemporaryDirectory() as tmp:
                    p = os.path.join(tmp, "syn.csv")
                    synth_csv(m, p, d, J0, B0, TC0)
                    with contextlib.redirect_stdout(io.StringIO()):
                        r = m.analyze_const_accel_v3(p, 0, m.KT_OUT, d, m.ACCEL_LEVELS, m.V_PEAK, n_boot=100)
                msg = f"定速段 {hold} 圈，方向 {d:+d}"
                self.assertEqual(len(r["cycles"]), 80, msg)
                self.assertAlmostEqual(r["J_origin"], J0, delta=0.05, msg=msg)
                self.assertAlmostEqual(r["B"], B0, delta=0.5, msg=msg)
                self.assertAlmostEqual(r["Tc"], TC0, delta=0.1, msg=msg)
                self.assertAlmostEqual(r["J_intercept"], 0.0, delta=0.05, msg=msg)
                self.assertFalse(r["intercept_significant"], msg)
                if hold > 0:
                    self.assertAlmostEqual(r["hold_minus_f"], 0.0, delta=0.1, msg=msg)
                else:
                    self.assertTrue(math.isnan(r["hold_minus_f"]), "無定速段時 hold_minus_f 應為 NaN")

    def test_dec_friction_offset(self):
        """減速段摩擦多 1.0 N·m：各檔 J 隨 a 偏離，含截距回歸找回 J、截距 ≈ −0.5 且標示為顯著。"""
        m = load("v3_ana4")
        J0, B0, TC0 = 2.0, 32.5, 8.9
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "syn.csv")
            synth_csv(m, p, +1, J0, B0, TC0, dec_extra=1.0)
            with contextlib.redirect_stdout(io.StringIO()):
                r = m.analyze_const_accel_v3(p, 0, m.KT_OUT, +1, m.ACCEL_LEVELS, m.V_PEAK, n_boot=50)
        self.assertGreater(r["J_levels"][-1] - r["J_levels"][0], 0.5, "各檔 J 應隨 a 增加")
        self.assertAlmostEqual(r["J_slope"], J0, delta=0.1)
        self.assertAlmostEqual(r["J_intercept"], -0.5, delta=0.07)
        self.assertTrue(r["intercept_significant"])

    def test_offsets_reduce_ripple_bias(self):
        m = load("v3_ana2")
        J0, B0, TC0 = 2.0, 32.5, 8.9
        err = {}
        for use in (True, False):
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                synth_csv(m, p, +1, J0, B0, TC0, use_offsets=use)
                with contextlib.redirect_stdout(io.StringIO()):
                    r = m.analyze_const_accel_v3(p, 0, m.KT_OUT, +1, m.ACCEL_LEVELS, m.V_PEAK, n_boot=20)
            err[use] = np.max(np.abs(r["J_levels"] - J0))
        self.assertLess(err[True], err[False], f"錯開 {err[True]:.4f} vs 不錯開 {err[False]:.4f}")

    def test_settle_removes_transient_bias(self):
        """換階段後力矩有一階延遲（15 ms）時，跳過前 40 ms 可去除大部分偏差。"""
        m = load("v3_ana3")
        J0, B0, TC0 = 2.0, 32.5, 8.9
        err = {}
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "syn.csv")
            synth_csv(m, p, +1, J0, B0, TC0, lag_s=0.015)
            for st in (0.0, 0.040):
                with contextlib.redirect_stdout(io.StringIO()):
                    r = m.analyze_const_accel_v3(p, 0, m.KT_OUT, +1, m.ACCEL_LEVELS, m.V_PEAK,
                                                 settle_s=st, n_boot=20)
                err[st] = np.max(np.abs(r["J_levels"] - J0))
        self.assertLess(err[0.040], err[0.0], f"跳過 {err[0.040]:.4f} vs 不跳過 {err[0.0]:.4f}")
        self.assertLess(err[0.040], 0.08, f"跳過 40 ms 後各檔 J 最大誤差 {err[0.040]:.4f}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
