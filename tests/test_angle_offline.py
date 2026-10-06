# -*- coding: utf-8 -*-
"""
固定角度窗口定速腳本（robot/ur5_const_velocity_angle_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同 test_offline.py：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在 J0=180° 不動）、姿態模組設為 None。

測試內容
--------
1. 語法正確。
2. 角度窗口長度 = 10 個馬達圈；正反轉的規劃以窗口中心對稱；FULL 為 12 段（先升後降）。
3. 安全檢查：J0 視窗設定時三個階段都通過；視窗未設定、速度超限時 FAIL。
4. URScript：縮排規則、每段一個有次數上限的迴圈、正轉用 "<"、反轉用 ">"、只含 ASCII。
5. main()：視窗未設定、RTDE 失敗、起始位置離第一段起點過遠時，都在送出軌跡前中止。
6. 端對端：CSV 與摘要檔自動命名且不撞名。
7. 分析：用已知 B、Tc 加上鎖定角度的漣波與雜訊合成資料，analyze_angle_window 能找回答案。
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


def load(name):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, ANG)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.home_pose = None
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def plan(m, stage, direction):
    levels, order, n_revs = m.stage_config(stage)
    ta, tb = m.window_bounds_deg(m.WINDOW_CENTER_DEG, n_revs, m.GEAR_RATIO)
    segs = m.plan_segments(m.speed_sequence(levels, order), direction, ta, tb, m.SPEEDJ_ACCEL,
                           m.STOPJ_DECEL, m.SETTLE_TIME, m.OVERRUN_DEG, m.LOOP_TIME_FACTOR, m.DT)
    return ta, tb, segs


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
        py_compile.compile(ANG, doraise=True)

    def test_window_and_symmetry(self):
        m = load("ang_plan")
        ta, tb, sp = plan(m, "FULL", +1)
        _, _, sn = plan(m, "FULL", -1)
        self.assertAlmostEqual(tb - ta, 10 * 360 / 101, places=9)
        self.assertEqual(len(sp), 12)
        self.assertEqual([s["speed"] for s in sp], m.SPEED_LEVELS + m.SPEED_LEVELS[::-1])
        c = m.WINDOW_CENTER_DEG
        for a, b in zip(sp, sn):
            self.assertAlmostEqual(a["start_deg"] - c, -(b["start_deg"] - c), places=9)
            self.assertAlmostEqual(a["extreme_deg"] - c, -(b["extreme_deg"] - c), places=9)
            # 進入窗口前至少有「加速 + 穩定」的距離
            self.assertLess(a["start_deg"], ta)
            self.assertGreater(a["stop_cmd_deg"], tb)

    def test_safety(self):
        m = load("ang_safety")
        for stage in ("QUICK", "MID", "FULL"):
            for d in (+1, -1):
                _, _, segs = plan(m, stage, d)
                ok, lines, _ = m.safety_check(segs, m.QD_MAX, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL,
                                              m.STOPJ_DECEL, m.MOVEJ_VEL, m.MOVEJ_ACCEL, 65, 300,
                                              m.SAFETY_MARGIN_DEG)
                self.assertTrue(ok, f"{stage} {d:+d}：" + "\n".join(lines))
        _, _, segs = plan(m, "FULL", +1)
        ok, _, _ = m.safety_check(segs, m.QD_MAX, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                  m.MOVEJ_VEL, m.MOVEJ_ACCEL, None, None, m.SAFETY_MARGIN_DEG)
        self.assertFalse(ok, "J0 視窗未設定時應 FAIL")
        ok, _, _ = m.safety_check(segs, 0.15, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                  m.MOVEJ_VEL, m.MOVEJ_ACCEL, 65, 300, m.SAFETY_MARGIN_DEG)
        self.assertFalse(ok, "速度超過上限時應 FAIL")
        ok, _, _ = m.safety_check(segs, m.QD_MAX, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                  m.MOVEJ_VEL, m.MOVEJ_ACCEL, 170, 300, m.SAFETY_MARGIN_DEG)
        self.assertFalse(ok, "規劃範圍超出 J0 視窗時應 FAIL")

    def test_urscript(self):
        m = load("ang_script")
        for d, op in ((+1, "<"), (-1, ">")):
            _, _, segs = plan(m, "FULL", d)
            s = m.build_urscript(segs, d, 0, HOME_Q, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                 m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT)
            self.assertTrue(urscript_indent_ok(s))
            self.assertTrue(all(ord(ch) < 128 for ch in s), "URScript 應只含 ASCII")
            self.assertEqual(s.count("while (n < "), len(segs), "每段都要有次數上限")
            self.assertEqual(s.count(f"q_now[0] {op} "), len(segs))
            self.assertEqual(s.count("movej("), len(segs))
            self.assertEqual(s.count("stopj("), len(segs))
            # movej 只改 J0：其他五軸維持 HOME_Q
            for line in s.splitlines():
                if "movej(" in line:
                    vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    for j in range(1, 6):
                        self.assertAlmostEqual(vals[j], HOME_Q[j], places=5)


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
        m = load("ang_main1")
        m.J0_SAFE_MIN_DEG = m.J0_SAFE_MAX_DEG = None
        self.assertEqual(self._run(m), [], "視窗未設定時不該送出任何軌跡")

    def test_rtde_failure_aborts(self):
        m = load("ang_main2")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m), [], "RTDE 失敗時不該送出任何軌跡")

    def test_far_start_aborts(self):
        m = load("ang_main3")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG = 65, 300
        FakeRTDE.q = [math.radians(100.0)] + HOME_Q[1:]       # J0 離窗口很遠
        self.assertEqual(self._run(m), [], "起始位置過遠時不該送出任何軌跡")

    def test_autoname_end_to_end(self):
        m = load("ang_main4")
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
            for _, _, files in os.walk(tmp):
                names.update(files)
        csvs = sorted(n for n in names if n.endswith(".csv"))
        self.assertEqual(len(csvs), 2, csvs)
        self.assertTrue(any("_pos_quick_" in n for n in csvs))
        self.assertTrue(any("_neg_quick_" in n for n in csvs))


class TestAnalysis(unittest.TestCase):
    def test_recovers_B_Tc_with_locked_ripple(self):
        """已知 B=32.5、Tc=8.9，加上跟馬達角度鎖定的漣波（±1.5 N·m 量級）與雜訊，應找回答案。"""
        m = load("ang_ana")
        B0, TC0 = 32.5, 8.9
        rng = np.random.default_rng(0)
        for d in (+1, -1):
            ta, tb, segs = plan(m, "FULL", d)
            rows, t = [], 0.0
            for s in segs:
                v = s["speed"]
                q_deg = s["start_deg"]
                end = s["stop_cmd_deg"]
                while (q_deg - end) * d < 0:
                    q_rad = math.radians(q_deg)
                    phi = (q_rad * m.GEAR_RATIO) % (2 * math.pi)
                    ripple = 1.2 * math.cos(5 * phi + 0.7) + 0.8 * math.sin(12 * phi) + 0.4 * math.cos(2 * phi)
                    tau = B0 * d * v + TC0 * d + ripple + rng.normal(0, 0.4)
                    row = {"timestamp": t, "actual_q_0": q_rad, "actual_qd_0": d * v + rng.normal(0, 0.002),
                           "actual_current_0": tau / m.KT_OUT, "target_qd_0": d * v}
                    rows.append(row)
                    q_deg += d * math.degrees(v * m.DT)
                    t += m.DT
                t += 5.0
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                with open(p, "w", newline="", encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                    w.writeheader(); w.writerows(rows)
                with contextlib.redirect_stdout(io.StringIO()):
                    r = m.analyze_angle_window(p, 0, m.KT_OUT, d, m.SPEED_LEVELS, ta, tb)
            self.assertEqual(len(r["segments"]), 12)
            for sg in r["segments"]:
                self.assertAlmostEqual(sg["revs"], 10.0, delta=0.05)
            self.assertAlmostEqual(r["B"], B0, delta=0.15, msg=f"方向 {d:+d}")
            self.assertAlmostEqual(r["Tc"], TC0, delta=0.02, msg=f"方向 {d:+d}")


if __name__ == "__main__":
    unittest.main(verbosity=2)
