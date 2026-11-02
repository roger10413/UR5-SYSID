# -*- coding: utf-8 -*-
"""
J1 固定角度窗口定速（robot/ur5_j1_const_velocity_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在診斷原點不動）。

測試內容
--------
1. 語法正確。
2. 規劃：窗口 6 個馬達圈、中心 −109°；三個階段、兩個方向的起點與煞車後最遠點都在 −127.71°～−90° 之間
   （不低於診斷原點），正反轉以窗口中心對稱。
3. 安全檢查：J1 視窗設定時三個階段都通過；視窗未設定、速度超限、規劃範圍超出視窗、
   窗口移到手臂下降側（法蘭過低）、法蘭高度下限不是數值時 FAIL。
4. 法蘭高度：診斷原點 0.432 m、J1 −127.71° 時 0.658 m。
5. 姿態檢查：六軸都在容差內才通過；J0、J1 不符、NaN、讀不到時拒絕；±360° 等價。
6. URScript：縮排規則、只含 ASCII、每段一個有次數上限的迴圈、正轉用 "<"、反轉用 ">"、
   movej 只改 J1、最後回到 J1 = −90°。
7. 負載設定與輸入確認。
8. main()：J1 視窗未設定、負載未設定、未輸入 J1、負載不符、RTDE 失敗、姿態不符、第一次 movej 過遠、
   速度滑桿不是 100% 時都在送出前中止；QUICK 兩個方向端對端各送出一次腳本，
   產生自動命名的 CSV（含 target_moment）與摘要檔。
9. 分析（合成資料）：已知 B、Tc、重力保持力矩（官方模型形狀，窗口內由 0 到 −17 N·m）與換算比例，
   加上跟馬達角度鎖定的漣波與雜訊，analyze_angle_window + combine_directions 找回 B、Tc、
   窗口平均重力保持力矩與換算比例。
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
J1C = os.path.join(ROBOT, "ur5_j1_const_velocity_ident.py")

HOME_Q = [3.14159265, -1.57079633, 1.57079633, -1.57079633, -1.57079633, 0.0]
EXT45_Q = [HOME_Q[0]] + [math.radians(x) for x in (-45.0, 45.0, -90.0, -90.0, 0.0)]
PAYLOAD = 1.3
SAFE = (-131.0, -87.0)                               # 測試用 J1 視窗（含 FULL 規劃範圍 −129.71°～−88.00°）
TM_HOME = [0.0, -23.3, 5.0, 1.0, 0.1, 0.0]           # 假的 target_moment（診斷原點）


class FakeRTDE:
    q = list(HOME_Q)
    sf = 1.0
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getTargetQ(self): return list(FakeRTDE.q)
    def getTargetMoment(self): return list(TM_HOME)
    def getJointTemperatures(self): return [30.0] * 6
    def getTargetSpeedFraction(self): return FakeRTDE.sf
    def disconnect(self): pass


def load(name):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, J1C)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PAYLOAD_KG = PAYLOAD
    mod.J1_SAFE_MIN_DEG, mod.J1_SAFE_MAX_DEG = SAFE
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(confirm="J1", payload=str(PAYLOAD)):
    """模擬操作者輸入：問負載質量時回答 payload，要求輸入 J1 時回答 confirm，其他提示回答空字串。"""
    def ans(prompt="", *a, **k):
        p = str(prompt)
        if "負載質量" in p:
            return payload
        if "輸入 J1" in p:
            return confirm
        return ""
    return ans


def plan(m, stage, direction, center=None):
    levels, order, n_revs = m.stage_config(stage)
    c = m.WINDOW_CENTER_DEG if center is None else center
    ta, tb = m.window_bounds_deg(c, n_revs, m.GEAR_RATIO)
    segs = m.plan_segments(m.speed_sequence(levels, order), direction, ta, tb, m.SPEEDJ_ACCEL,
                           m.STOPJ_DECEL, m.SETTLE_TIME, m.OVERRUN_DEG, m.LOOP_TIME_FACTOR, m.DT)
    return ta, tb, segs


def check(m, segs, safe=SAFE, qd_max=None, **kw):
    return m.safety_check(segs, m.QD_MAX if qd_max is None else qd_max, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL,
                          m.STOPJ_DECEL, m.MOVEJ_VEL, m.MOVEJ_ACCEL, safe[0], safe[1], m.SAFETY_MARGIN_DEG,
                          end_deg=m.END_J1_DEG, **kw)


class Base(unittest.TestCase):
    def setUp(self):
        self._sock = socket.socket
        self._input = builtins.input
        socket.socket = no_network
        builtins.input = answers()
        FakeRTDE.q = list(HOME_Q)
        FakeRTDE.sf = 1.0

    def tearDown(self):
        socket.socket = self._sock
        builtins.input = self._input
        FakeRTDE.q = list(HOME_Q)
        FakeRTDE.sf = 1.0


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


def movej_targets(text):
    return [[float(x) for x in l.split("[")[1].split("]")[0].split(",")] for l in text.splitlines() if "movej(" in l]


class TestPlan(Base):
    def test_compile(self):
        py_compile.compile(J1C, doraise=True)

    def test_window_and_range(self):
        m = load("j1c_plan")
        ta, tb, sp = plan(m, "FULL", +1)
        _, _, sn = plan(m, "FULL", -1)
        self.assertAlmostEqual(tb - ta, 6 * 360 / 101, places=9)
        self.assertAlmostEqual((ta + tb) / 2, -109.0, places=9)
        self.assertEqual(len(sp), 12)
        self.assertEqual([s["speed"] for s in sp], m.SPEED_LEVELS + m.SPEED_LEVELS[::-1])
        c = m.WINDOW_CENTER_DEG
        for a, b in zip(sp, sn):
            self.assertAlmostEqual(a["start_deg"] - c, -(b["start_deg"] - c), places=9)
            self.assertAlmostEqual(a["extreme_deg"] - c, -(b["extreme_deg"] - c), places=9)
            self.assertLess(a["start_deg"], ta)
            self.assertGreater(a["stop_cmd_deg"], tb)
        for stage in ("QUICK", "MID", "FULL"):
            for d in (+1, -1):
                _, _, segs = plan(m, stage, d)
                pts = [s["start_deg"] for s in segs] + [s["extreme_deg"] for s in segs]
                self.assertLessEqual(max(pts), -90.0, f"{stage} {d:+d}：不可低於診斷原點")
                self.assertGreaterEqual(min(pts), -127.72, f"{stage} {d:+d}")
        self.assertAlmostEqual(min(s["start_deg"] for s in sp), -127.71, delta=0.01)
        self.assertAlmostEqual(max(s["start_deg"] for s in sn), -90.29, delta=0.01)

    def test_safety(self):
        m = load("j1c_safety")
        for stage in ("QUICK", "MID", "FULL"):
            for d in (+1, -1):
                _, _, segs = plan(m, stage, d)
                ok, lines, _ = check(m, segs)
                self.assertTrue(ok, f"{stage} {d:+d}：" + "\n".join(lines))
        _, _, segs = plan(m, "FULL", +1)
        ok, _, (lo, hi) = check(m, segs)
        self.assertAlmostEqual(lo, -129.71, delta=0.01)
        self.assertAlmostEqual(hi, -88.0, places=9)
        self.assertFalse(check(m, segs, safe=(None, None))[0], "J1 視窗未設定時應 FAIL")
        self.assertFalse(check(m, segs, safe=(None, -87.0))[0])
        self.assertFalse(check(m, segs, qd_max=0.15)[0], "速度超過上限時應 FAIL")
        self.assertFalse(check(m, segs, safe=(-125.0, -87.0))[0], "規劃範圍超出 J1 視窗時應 FAIL")
        self.assertFalse(check(m, segs, flange_z_min=None)[0], "法蘭高度下限不是數值時應 FAIL")
        self.assertFalse(check(m, segs, flange_z_min=float("nan"))[0])
        # 窗口設到手臂下降側：J1 視窗再寬，法蘭高度檢查也要擋下
        _, _, low = plan(m, "FULL", -1, center=-75.0)
        ok, lines, _ = check(m, low, safe=(-131.0, -30.0))
        self.assertFalse(ok)
        k = next(i for i, l in enumerate(lines) if l.startswith("法蘭"))
        self.assertEqual(lines[k + 1].strip(), "[FAIL]")
        self.assertTrue(all("[FAIL]" not in l for l in lines[:k]), "其他項目應通過，只有法蘭高度 FAIL")

    def test_flange_height(self):
        m = load("j1c_fk")
        self.assertAlmostEqual(m.flange_height(-90.0), 0.432, delta=0.001)
        self.assertAlmostEqual(m.flange_height(-127.71), 0.658, delta=0.001)
        self.assertGreater(m.flange_height(-88.0), m.FLANGE_Z_MIN_M)
        self.assertLess(m.flange_height(-80.0), m.FLANGE_Z_MIN_M)

    def test_pose(self):
        m = load("j1c_pose")
        self.assertTrue(m.check_pose(HOME_Q)[0])
        q = list(HOME_Q); q[0] = -math.pi
        self.assertTrue(m.check_pose(q)[0], "J0 = −180° 與 180° 等價")
        q = list(HOME_Q); q[0] = math.radians(170.0)
        self.assertFalse(m.check_pose(q)[0], "J0 不在 180° 應拒絕")
        q = list(HOME_Q); q[1] += math.radians(0.4)
        self.assertTrue(m.check_pose(q)[0])
        q = list(HOME_Q); q[1] += math.radians(0.6)
        self.assertFalse(m.check_pose(q)[0], "J1 不在 −90° 應拒絕")
        q = list(HOME_Q); q[2] = float("nan")
        self.assertFalse(m.check_pose(q)[0])
        q = list(HOME_Q); q[3] = None
        self.assertFalse(m.check_pose(q)[0])
        self.assertFalse(m.check_pose(EXT45_Q)[0])
        self.assertFalse(m.check_pose(None)[0])
        self.assertFalse(m.check_pose(HOME_Q[:5])[0])

    def test_urscript(self):
        m = load("j1c_script")
        for d, op in ((+1, "<"), (-1, ">")):
            _, _, segs = plan(m, "FULL", d)
            s = m.build_urscript(segs, d, 1, HOME_Q, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                 m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=m.END_J1_DEG)
            self.assertTrue(s.startswith("def j1_const_velocity_program():"))
            self.assertTrue(urscript_indent_ok(s))
            self.assertTrue(all(ord(ch) < 128 for ch in s), "URScript 應只含 ASCII")
            self.assertEqual(s.count("while (n < "), len(segs), "每段都要有次數上限")
            self.assertEqual(s.count(f"q_now[1] {op} "), len(segs))
            self.assertEqual(s.count("stopj("), len(segs))
            tg = movej_targets(s)
            self.assertEqual(len(tg), len(segs) + 1, "最後要多一個回原點的 movej")
            for v in tg:
                for j in (0, 2, 3, 4, 5):
                    self.assertAlmostEqual(v[j], HOME_Q[j], places=6, msg="movej 只改 J1")
                self.assertLessEqual(math.degrees(v[1]), -90.0 + 1e-4, "movej 目標不可低於診斷原點")
            for v, seg in zip(tg, segs):
                self.assertAlmostEqual(math.degrees(v[1]), seg["start_deg"], places=4)
            self.assertAlmostEqual(math.degrees(tg[-1][1]), -90.0, places=4)
            self.assertGreater(s.rfind("movej("), s.rfind("stopj("), "回原點應在最後一段 stopj 之後")
            for line in s.splitlines():
                if "speedj(" in line:
                    qd = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    self.assertEqual([j for j in range(6) if qd[j] != 0.0], [1], "speedj 只動 J1")
                    self.assertEqual(math.copysign(1.0, qd[1]), float(d))

    def test_payload(self):
        m = load("j1c_payload")
        for bad in (None, -1.0, float("nan"), "abc"):
            self.assertFalse(m.check_payload_setting(bad)[0], repr(bad))
        self.assertTrue(m.check_payload_setting(1.3)[0])
        self.assertTrue(m.check_payload_entry(1.3, " 1.30 ")[0])
        self.assertTrue(m.check_payload_entry(1.3, "1.304")[0])
        for bad in ("1.31", "2.30", "abc", ""):
            self.assertFalse(m.check_payload_entry(1.3, bad)[0], repr(bad))
        self.assertFalse(m.check_payload_entry(None, "1.3")[0])
        self.assertEqual(m.payload_tag(1.3), "p130")
        self.assertEqual(m.payload_tag(2.3), "p230")


class TestMain(Base):
    def _run(self, m, stage="QUICK", directions=(+1,), wait=None):
        """送出的腳本照收；wait 不為 None 時取代等待函式。回傳 (送出的腳本, {檔名: 內容})。"""
        sent, files = [], {}
        m.send_urscript = lambda ip, txt, *a, **k: sent.append(txt)
        if wait is not None:
            m.wait_for_motion_complete = wait
        buf = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            for d in directions:
                m.TEST_STAGE, m.DIRECTION = stage, d
                with contextlib.redirect_stdout(buf):
                    m.main()
            for root, _, names in os.walk(tmp):
                for n in names:
                    with open(os.path.join(root, n), encoding="utf-8") as f:
                        files[n] = f.read()
        self.last_out = buf.getvalue()
        return sent, files

    def test_window_unset_aborts(self):
        m = load("j1c_m1")
        m.J1_SAFE_MIN_DEG = m.J1_SAFE_MAX_DEG = None
        self.assertEqual(self._run(m)[0], [], "視窗未設定時不該送出任何軌跡")
        self.assertIn("安全檢查未通過", self.last_out)

    def test_payload_unset_aborts(self):
        m = load("j1c_m2")
        m.PAYLOAD_KG = None
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("安全檢查未通過", self.last_out)

    def test_confirm_cancel(self):
        m = load("j1c_m3")
        builtins.input = answers(confirm="j1")
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("已取消", self.last_out)

    def test_payload_mismatch_aborts(self):
        m = load("j1c_m4")
        builtins.input = answers(payload="2.30")
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("負載確認未通過", self.last_out)

    def test_rtde_failure_aborts(self):
        m = load("j1c_m5")

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("RTDE 記錄未能啟動", self.last_out)

    def test_wrong_pose_aborts(self):
        q_j1 = list(HOME_Q); q_j1[1] = math.radians(-100.0)
        for k, q in enumerate(([math.radians(170.0)] + HOME_Q[1:], q_j1, EXT45_Q)):
            m = load(f"j1c_m6{k}")
            FakeRTDE.q = list(q)
            self.assertEqual(self._run(m)[0], [], f"案例 {k}")
            self.assertIn("姿態檢查未通過", self.last_out)

    def test_far_start_aborts(self):
        m = load("j1c_m7")
        m.MAX_INITIAL_MOVE_DEG = 20.0                       # QUICK 正轉第一次 movej 為 22.4°
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("起始位置檢查未通過", self.last_out)

    def test_speed_fraction_aborts(self):
        for k, sf in enumerate((0.52, None, float("nan"))):
            m = load(f"j1c_m8{k}")
            FakeRTDE.sf = sf
            self.assertEqual(self._run(m)[0], [], f"案例 {k}")
            self.assertIn("速度滑桿檢查未通過", self.last_out)

    def test_quick_end_to_end(self):
        m = load("j1c_m9")
        sent, files = self._run(m, directions=(+1, -1),
                                wait=lambda logger, planned_s, **k: ("done", planned_s, planned_s))
        self.assertEqual(len(sent), 2)
        for txt in sent:
            self.assertTrue(txt.startswith("def j1_const_velocity_program():"))
            self.assertEqual(len(movej_targets(txt)), 2, "QUICK：1 段 + 回原點")
        self.assertIn("q_now[1] < ", sent[0])
        self.assertIn("q_now[1] > ", sent[1])
        for tag in ("pos", "neg"):
            data = [n for n in files if n.startswith(f"j1constvel_data_{tag}_p130_quick_") and n.endswith(".csv")]
            self.assertEqual(len(data), 1, tag)
            rows = list(csv.DictReader(io.StringIO(files[data[0]])))
            self.assertGreater(len(rows), 0)
            self.assertAlmostEqual(float(rows[0]["target_moment_1"]), TM_HOME[1])
            self.assertTrue(all(float(r["target_speed_fraction"]) == 1.0 for r in rows))
            info = next(v for n, v in files.items() if n.startswith(f"j1constvel_info_{tag}_p130_quick_"))
            self.assertIn("執行狀態: done", info)
            self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
            self.assertIn("J1：目前   -90.00°", info)
        self.assertIn("下一步", self.last_out)
        self.assertNotIn("這次的資料不可用", self.last_out)


def holding_torque(q_deg):
    """合成用的重力保持力矩 [N·m]（官方連桿質量 + 1.30 kg 工具的形狀：−90° 時 −23.3，−119.5° 附近為零）。"""
    return -47.38 * math.cos(math.radians(q_deg + 29.49))


class TestAnalysis(unittest.TestCase):
    def test_recovers_B_Tc_gravity_and_K(self):
        m = load("j1c_ana")
        B0, TC0, K0 = 30.0, 9.0, 11.0                       # 真值（N·m 單位）與真正的換算比例 [N·m/A]
        scale = m.KT_OUT / K0                               # 電流 × KT_OUT 相對於真值的比例
        rng = np.random.default_rng(0)
        res = {}
        for d in (+1, -1):
            ta, tb, segs = plan(m, "FULL", d)
            rows, t = [], 0.0
            for s in segs:
                v, q_deg, end = s["speed"], s["start_deg"], s["stop_cmd_deg"]
                while (q_deg - end) * d < 0:
                    q_rad = math.radians(q_deg)
                    phi = (q_rad * m.GEAR_RATIO) % (2 * math.pi)
                    ripple = 1.2 * math.cos(5 * phi + 0.7) + 0.8 * math.sin(12 * phi) + 0.4 * math.cos(2 * phi)
                    h = holding_torque(q_deg)
                    tau_true = B0 * d * v + TC0 * d + h + ripple + rng.normal(0, 0.4)
                    rows.append({"timestamp": t, "actual_q_1": q_rad, "actual_qd_1": d * v + rng.normal(0, 0.002),
                                 "actual_current_1": tau_true / K0, "target_qd_1": d * v, "target_moment_1": h})
                    q_deg += d * math.degrees(v * m.DT)
                    t += m.DT
                t += 5.0
            with tempfile.TemporaryDirectory() as tmp:
                p = os.path.join(tmp, "syn.csv")
                with open(p, "w", newline="", encoding="utf-8") as f:
                    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                    w.writeheader(); w.writerows(rows)
                with contextlib.redirect_stdout(io.StringIO()):
                    res[d] = m.analyze_angle_window(p, 1, m.KT_OUT, d, m.SPEED_LEVELS, ta, tb)
            self.assertEqual(len(res[d]["segments"]), 12)
            for sg in res[d]["segments"]:
                self.assertAlmostEqual(sg["revs"], 6.0, delta=0.05)
            self.assertAlmostEqual(res[d]["B_up"], res[d]["B_down"], delta=0.5)
        with contextlib.redirect_stdout(io.StringIO()):
            c = m.combine_directions(res[+1], res[-1])
        h_bar = float(np.mean([holding_torque(x) for x in np.linspace(ta, tb, 20001)]))
        self.assertAlmostEqual(h_bar, -8.57, delta=0.05, msg="合成的窗口平均重力保持力矩")
        self.assertAlmostEqual(c["B"], B0 * scale, delta=0.01 * B0 * scale)
        self.assertAlmostEqual(c["Tc"], TC0 * scale, delta=0.05)
        self.assertAlmostEqual(c["H"], h_bar * scale, delta=0.05)
        self.assertAlmostEqual(c["tm_win"], h_bar, delta=0.02)
        self.assertAlmostEqual(c["K_mean"], K0, delta=0.01 * K0)
        self.assertAlmostEqual(c["K_slope"], K0, delta=0.02 * K0)
        # 只看單一方向時 C 混入重力：C₊ = Tc + H̄、C₋ = Tc − H̄
        self.assertAlmostEqual(res[+1]["C"], (TC0 + h_bar) * scale, delta=0.08)
        self.assertAlmostEqual(res[-1]["C"], (TC0 - h_bar) * scale, delta=0.08)

    def test_combine_needs_regression(self):
        m = load("j1c_ana2")
        with self.assertRaises(ValueError):
            m.combine_directions({"segments": []}, {"segments": []})


if __name__ == "__main__":
    unittest.main(verbosity=2)
