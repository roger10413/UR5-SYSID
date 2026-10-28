# -*- coding: utf-8 -*-
"""
J1 靜止電流測試（ur5_j1_static_current.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式、假的 rtde_receive（手臂停在診斷原點）。

測試內容
--------
1. 語法正確。
2. 安全檢查：預設值通過；移動角度、速度、加速度超過上限或不是有限正數、靜止時間太短、姿態未知時 FAIL。
3. 姿態檢查：六軸（含 J0 = 180°）都在容差內才通過；未知姿態、讀不到、NaN 拒絕。
4. URScript：縮排規則、只含 ASCII、movej 只改 J1、正負向靠近交替、最後回到開始位置。
5. 力臂：DIAG 0.487 m、EXT45 0.787 m（與 POSES 一致），與工具軸上的位置無關。
6. 分析（合成資料）：找回中心電流、靜摩擦帶寬與靠近方向；夾砝碼前後找回設定的 K。
7. main()：負載未設定、RTDE 失敗、姿態不符、速度滑桿不是 100% 時都在送出前中止；QUICK 端對端送出一次腳本，
   產生自動命名的 CSV 與摘要檔。
"""

import os
import sys
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
import csv

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROBOT = os.path.join(ROOT, "robot") if os.path.isdir(os.path.join(ROOT, "robot")) else ROOT
J1S = os.path.join(ROBOT, "ur5_j1_static_current.py")

HOME_Q = [3.14159265, -1.57079633, 1.57079633, -1.57079633, -1.57079633, 0.0]
EXT45_Q = [HOME_Q[0]] + [math.radians(x) for x in (-45.0, 45.0, -90.0, -90.0, 0.0)]
PAYLOAD = 1.3


class FakeRTDE:
    q = list(HOME_Q)
    sf = 1.0
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getTargetMoment(self): return [0.0, -40.0, 12.0, 1.5, 0.2, 0.0]
    def getJointTemperatures(self): return [30.0] * 6
    def getTargetSpeedFraction(self): return FakeRTDE.sf
    def disconnect(self): pass


def load(name):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, J1S)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PAYLOAD_KG = PAYLOAD
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(pose="DIAG", payload=str(PAYLOAD)):
    """模擬操作者輸入：問負載質量時回答 payload，問姿態時回答 pose，其他提示回答空字串。"""
    def ans(prompt="", *a, **k):
        p = str(prompt)
        if "負載質量" in p:
            return payload
        if "輸入" in p and "繼續" in p:
            return pose
        return ""
    return ans


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


def synth_csv(path, m, q_start, n_cycle, i_center, f_half, tm, noise=0.05, seed=0):
    """
    依 plan_moves 產生 125 Hz 的合成記錄：移動時 target_qd 為 ±MOVE_VEL，靜止時為 0；
    靜止時電流 = i_center + 靠近方向 × f_half + 白雜訊。開頭有 1 s 靜止（之前沒有運動，不算）。
    """
    rng = np.random.default_rng(seed)
    dt = 1.0 / m.SAMPLE_HZ
    q1 = math.degrees(q_start[1])
    rows, t, q_now, last_dir = [], 0.0, q1, 0

    def add(n, tqd, cur, qdeg):
        nonlocal t
        for _ in range(n):
            r = {"timestamp": round(t, 6)}
            for j in range(6):
                r[f"actual_q_{j}"] = q_start[j] if j != 1 else math.radians(qdeg)
                r[f"actual_qd_{j}"] = tqd if j == 1 else 0.0
                r[f"actual_current_{j}"] = cur + rng.normal(0, noise) if j == 1 else 0.0
                r[f"target_qd_{j}"] = tqd if j == 1 else 0.0
                r[f"joint_temp_{j}"] = 30.0
                r[f"target_moment_{j}"] = tm if j == 1 else 0.0
            rows.append(r)
            t += dt

    add(int(1.0 / dt), 0.0, i_center, q1)
    for deg, is_hold, approach in m.plan_moves(q1, m.STEP_DEG, n_cycle):
        sgn = 1 if deg > q_now else -1
        n_mv = int(m.movej_time(math.radians(abs(deg - q_now)), m.MOVE_VEL, m.MOVE_ACCEL) / dt)
        add(n_mv, sgn * m.MOVE_VEL, i_center + 2 * sgn * f_half, (q_now + deg) / 2)
        q_now, last_dir = deg, sgn
        dur = m.HOLD_S if is_hold else m.PAUSE_TIME
        add(int(dur / dt), 0.0, i_center + last_dir * f_half, deg)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


class TestPieces(Base):
    def test_compile(self):
        py_compile.compile(J1S, doraise=True)

    def test_safety(self):
        m = load("s_safe")
        ok = lambda **kw: m.safety_check(**{**dict(step_deg=m.STEP_DEG, v=m.MOVE_VEL, a=m.MOVE_ACCEL,
                                                   hold_s=m.HOLD_S, skip_s=m.SKIP_S, pose="DIAG"), **kw})[0]
        self.assertTrue(ok())
        self.assertTrue(ok(pose="EXT45"))
        for kw in (dict(step_deg=0.0), dict(step_deg=1.5), dict(step_deg=None), dict(step_deg=float("nan")),
                   dict(v=0.2), dict(v=-0.05), dict(a=1.0), dict(hold_s=1.5), dict(hold_s=None),
                   dict(pose="diag"), dict(pose="FOO")):
            self.assertFalse(ok(**kw), repr(kw))

    def test_pose(self):
        m = load("s_pose")
        self.assertTrue(m.check_pose(HOME_Q, "DIAG")[0])
        self.assertTrue(m.check_pose(EXT45_Q, "EXT45")[0])
        self.assertFalse(m.check_pose(HOME_Q, "EXT45")[0])
        q = list(HOME_Q); q[0] = math.radians(170.0)
        self.assertFalse(m.check_pose(q, "DIAG")[0], "J0 不在 180° 應拒絕")
        q = list(HOME_Q); q[1] += math.radians(0.4)
        self.assertTrue(m.check_pose(q, "DIAG")[0])
        q = list(HOME_Q); q[2] = float("nan")
        self.assertFalse(m.check_pose(q, "DIAG")[0])
        self.assertFalse(m.check_pose(None, "DIAG")[0])
        self.assertFalse(m.check_pose(HOME_Q, "FOO")[0])

    def test_urscript(self):
        m = load("s_script")
        q1 = math.degrees(HOME_Q[1])
        moves = m.plan_moves(q1, m.STEP_DEG, 3)
        txt = m.build_urscript(moves, HOME_Q, 1, m.MOVE_VEL, m.MOVE_ACCEL, m.HOLD_S, m.PAUSE_TIME)
        self.assertTrue(urscript_indent_ok(txt))
        self.assertTrue(all(ord(ch) < 128 for ch in txt))
        tg = movej_targets(txt)
        self.assertEqual(len(tg), 12)
        for v in tg:
            for j in (0, 2, 3, 4, 5):
                self.assertAlmostEqual(v[j], HOME_Q[j], places=6)
            self.assertLessEqual(abs(math.degrees(v[1]) - q1), m.STEP_DEG + 1e-4)    # URScript 寫到小數 6 位
        self.assertAlmostEqual(tg[-1][1], HOME_Q[1], places=6, msg="最後回到開始位置")
        holds = [(d, a) for d, h, a in moves if h]
        self.assertEqual([a for _, a in holds], [+1, -1] * 3, "正負向靠近交替")
        self.assertTrue(all(abs(d - q1) < 1e-9 for d, _ in holds))
        self.assertEqual(txt.count(f"sleep({m.HOLD_S})"), 6)
        self.assertEqual(m.stage_cycles("QUICK"), 1)
        self.assertEqual(m.stage_cycles("FULL"), m.N_CYCLE)
        with self.assertRaises(ValueError):
            m.stage_cycles("MID")

    def test_lever(self):
        m = load("s_lever")
        for q, pose in ((HOME_Q, "DIAG"), (EXT45_Q, "EXT45")):
            lv = m.gravity_lever(q, 1)
            self.assertAlmostEqual(abs(lv), m.POSES[pose]["d_j1_m"], delta=0.001, msg=pose)
            self.assertAlmostEqual(m.gravity_lever(q, 1, 0.25), lv, places=6, msg="工具垂直時與砝碼位置無關")

    def test_analyze_and_known_mass(self):
        m = load("s_ana")
        with tempfile.TemporaryDirectory() as tmp:
            p0, p1 = os.path.join(tmp, "a.csv"), os.path.join(tmp, "b.csv")
            K_true, mass = 12.0, 1.0
            lever = m.gravity_lever(EXT45_Q, 1)
            d_tau = mass * m.G * lever
            synth_csv(p0, m, EXT45_Q, 8, i_center=-3.0, f_half=0.6, tm=-36.0, seed=1)
            synth_csv(p1, m, EXT45_Q, 8, i_center=-3.0 + d_tau / K_true, f_half=0.6, tm=-36.0 + d_tau, seed=2)
            with contextlib.redirect_stdout(io.StringIO()):
                r0 = m.analyze_static_current(p0)
                r1 = m.analyze_static_current(p1)
                c = m.compare_known_mass(r1, r0, mass)
        self.assertEqual((r0["n_pos"], r0["n_neg"]), (8, 8))
        self.assertAlmostEqual(r0["i_mid"], -3.0, delta=0.005)
        self.assertAlmostEqual(r0["i_half"], 0.6, delta=0.005)
        self.assertAlmostEqual(c["K"], K_true, delta=0.15)
        self.assertAlmostEqual(c["tm_ratio"], 1.0, places=6)
        self.assertTrue(c["same_sign"])


class TestMain(Base):
    def _run(self, m):
        sent, files = [], {}
        m.send_urscript = lambda ip, txt, *a, **k: sent.append(txt)
        buf = io.StringIO()
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            with contextlib.redirect_stdout(buf):
                m.main()
            for root, _, names in os.walk(tmp):
                for n in names:
                    with open(os.path.join(root, n), encoding="utf-8") as f:
                        files[n] = f.read()
        self.last_out = buf.getvalue()
        return sent, files

    def test_payload_unset_aborts(self):
        m = load("s_m1")
        m.PAYLOAD_KG = None
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("安全檢查未通過", self.last_out)

    def test_payload_mismatch_aborts(self):
        m = load("s_m2")
        builtins.input = answers(payload="2.30")
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("負載確認未通過", self.last_out)

    def test_pose_name_cancel(self):
        m = load("s_m3")
        builtins.input = answers(pose="no")
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("已取消", self.last_out)

    def test_rtde_failure_aborts(self):
        m = load("s_m4")

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("RTDE 記錄未能啟動", self.last_out)

    def test_wrong_pose_aborts(self):
        for k, q in enumerate((EXT45_Q, [math.radians(170.0)] + HOME_Q[1:])):
            m = load(f"s_m5{k}")
            FakeRTDE.q = list(q)
            self.assertEqual(self._run(m)[0], [], f"案例 {k}")
            self.assertIn("姿態檢查未通過", self.last_out)

    def test_speed_fraction_aborts(self):
        for k, sf in enumerate((0.52, None, float("nan"))):
            m = load(f"s_m6{k}")
            FakeRTDE.sf = sf
            self.assertEqual(self._run(m)[0], [], f"案例 {k}")
            self.assertIn("速度滑桿檢查未通過", self.last_out)

    def test_quick_end_to_end(self):
        m = load("s_m7")
        m.wait_for_motion_complete = lambda logger, planned_s, **k: ("done", planned_s, planned_s)
        sent, files = self._run(m)
        self.assertEqual(len(sent), 1)
        self.assertTrue(sent[0].startswith("def j1_static_current_program():"))
        self.assertEqual(len(movej_targets(sent[0])), 4)
        self.assertTrue(any(f.startswith("j1static_data_diag_p130_quick_") and f.endswith(".csv") for f in files))
        info = next(v for n, v in files.items() if n.startswith("j1static_info_diag_p130_quick_"))
        self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
        self.assertIn("執行狀態: done", info)


if __name__ == "__main__":
    unittest.main(verbosity=2)
