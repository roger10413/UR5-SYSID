# -*- coding: utf-8 -*-
"""
自動換姿態（robot/ur5_goto_pose.py）的離線測試。
執行：python -m unittest discover -s tests -v

不連手臂：rtde_receive 換成假物件，socket 送出的 URScript 只收集不送。
1. 姿態和各鑑別程式 POSES 的 q_deg 完全一致（UP_J3 = J3 程式的 UP）。
2. 差值換到 ±180°：讀到 −180° 不會整圈轉；未知姿態、NaN 丟錯。
3. 安全檢查：DIAG↔UP↔UP_J3 都通過且路徑最低點 > 0.40 m；轉太多、路徑太低拒絕。
4. URScript 結構：def 頂格、其餘縮排、end 頂格，movej 目標正確。
5. main()：已在目標不送程式；輸入 y 才送；取消不送；到位回 0；未知姿態回 1。
"""

import os
import sys
import math
import types
import builtins
import unittest
import importlib.util
import contextlib
import io

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROBOT = os.path.join(ROOT, "robot")


def load(name):
    spec = importlib.util.spec_from_file_location(name, os.path.join(ROBOT, name + ".py"))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


gp = load("ur5_goto_pose")


class FakeRTDE:
    q_deg = list(gp.POSES["DIAG"])
    move_to = None          # 送出程式後下一次讀取就到這個角度
    def __init__(self, ip): pass
    def getActualQ(self):
        return [math.radians(x) for x in FakeRTDE.q_deg]
    def getActualQd(self): return [0.0] * 6
    def disconnect(self): pass


class PoseConsistency(unittest.TestCase):
    def test_same_as_ident_programs(self):
        expect = {
            "ur5_j1_velocity_modulation_ident": {"UP": "UP"},
            "ur5_j1_const_velocity_ident": {"UP": "UP", "DIAG": "DIAG"},
            "ur5_j2_const_velocity_ident": {"UP": "UP", "DIAG": "DIAG"},
            "ur5_j3_const_velocity_ident": {"UP": "UP_J3", "DIAG": "DIAG"},
            "ur5_j3_velocity_modulation_ident": {"UP": "UP_J3", "DIAG": "DIAG"},
            "ur5_j4_const_velocity_ident": {"UP": "UP", "DIAG": "DIAG"},
            "ur5_j4_velocity_modulation_ident": {"UP": "UP", "DIAG": "DIAG"},
            "ur5_j5_const_velocity_ident": {"UP": "UP", "DIAG": "DIAG"},
            "ur5_j5_velocity_modulation_ident": {"UP": "UP", "DIAG": "DIAG"},
        }
        for modname, m in expect.items():
            mod = load(modname)
            for their, ours in m.items():
                self.assertEqual(tuple(mod.POSES[their]["q_deg"]), tuple(gp.POSES[ours]), f"{modname} {their}")


class Planning(unittest.TestCase):
    def test_wrap_no_full_turn(self):
        cur = [-180.0, -90.0, 90.0, -90.0, -90.0, 0.0]
        goal, delta = gp.plan_move(cur, "UP")
        self.assertAlmostEqual(goal[0], -180.0)
        self.assertAlmostEqual(delta[0], 0.0)
        self.assertAlmostEqual(delta[2], -90.0)

    def test_bad_inputs(self):
        with self.assertRaises(ValueError):
            gp.plan_move(list(gp.POSES["DIAG"]), "NOPE")
        with self.assertRaises(ValueError):
            gp.plan_move([float("nan")] + [0.0] * 5, "UP")
        with self.assertRaises(ValueError):
            gp.plan_move([0.0] * 5, "UP")

    def test_all_pairs_safe(self):
        for a in gp.POSES:
            for b in gp.POSES:
                cur = list(gp.POSES[a])
                goal, delta = gp.plan_move(cur, b)
                ok, lines = gp.safety_check(cur, goal, delta)
                self.assertTrue(ok, f"{a}->{b}: {lines}")
                self.assertGreater(gp.path_min_z(cur, goal), 0.40, f"{a}->{b}")

    def test_refuse_large_move(self):
        cur = [180.0, -90.0, -30.0, -90.0, -90.0, 0.0]
        goal, delta = gp.plan_move(cur, "DIAG")      # J2 要轉 120°
        ok, lines = gp.safety_check(cur, goal, delta)
        self.assertFalse(ok)
        self.assertTrue(any("超過上限" in s for s in lines))

    def test_refuse_low_path(self):
        cur = [180.0, -10.0, 90.0, -90.0, -90.0, 0.0]   # 大臂幾乎水平往前伸、前臂朝下
        goal, delta = gp.plan_move(cur, "DIAG")
        ok, lines = gp.safety_check(cur, goal, delta, max_move_deg=180.0)
        self.assertFalse(ok)
        self.assertTrue(any("低於" in s for s in lines))

    def test_move_time(self):
        self.assertAlmostEqual(gp.estimate_move_time(math.pi / 2), math.pi / 2 / 0.3 + 1.0)
        self.assertAlmostEqual(gp.estimate_move_time(0.27), 2 * math.sqrt(0.9))

    def test_urscript(self):
        s = gp.build_urscript(list(gp.POSES["UP"]))
        lines = s.rstrip("\n").split("\n")
        self.assertTrue(lines[0].startswith("def "))
        self.assertEqual(lines[-1], "end")
        for ln in lines[1:-1]:
            self.assertTrue(ln.startswith(" "))
        self.assertIn("movej([3.141593, -1.570796, 0.000000, -1.570796, -1.570796, 0.000000]", s)
        s.encode("ascii")


class Main(unittest.TestCase):
    def setUp(self):
        self.sent = []
        self._rtde, self._send, self._input = gp.rtde_receive, gp.send_urscript, builtins.input
        self._time = gp.time
        self.clock = [0.0]
        # 換掉 gp 自己的 time 參照（不要改到全域 time 模組，否則會影響其他測試）
        gp.time = types.SimpleNamespace(sleep=lambda s: None, time=lambda: self.clock[0])
        gp.rtde_receive = type("M", (), {"RTDEReceiveInterface": FakeRTDE})

        def fake_send(ip, text, port=30002):
            self.sent.append(text)
            if FakeRTDE.move_to is not None and "movej" in text:
                FakeRTDE.q_deg = list(FakeRTDE.move_to)
        gp.send_urscript = fake_send

    def tearDown(self):
        gp.rtde_receive, gp.send_urscript, builtins.input = self._rtde, self._send, self._input
        gp.time = self._time
        FakeRTDE.q_deg, FakeRTDE.move_to = list(gp.POSES["DIAG"]), None

    def run_main(self, argv, answer="y"):
        builtins.input = lambda prompt="": answer
        out = io.StringIO()
        with contextlib.redirect_stdout(out):
            rc = gp.main(argv)
        return rc, out.getvalue()

    def test_already_there(self):
        FakeRTDE.q_deg = [x + 0.3 for x in gp.POSES["UP"]]
        rc, out = self.run_main(["up"])
        self.assertEqual(rc, 0)
        self.assertEqual(self.sent, [])

    def test_move_diag_to_up(self):
        FakeRTDE.move_to = gp.POSES["UP"]
        rc, out = self.run_main(["UP"])
        self.assertEqual(rc, 0, out)
        self.assertEqual(len(self.sent), 1)
        self.assertIn("[完成]", out)

    def test_cancel(self):
        rc, out = self.run_main(["UP_J3"], answer="n")
        self.assertEqual(rc, 1)
        self.assertEqual(self.sent, [])

    def test_unknown(self):
        rc, out = self.run_main(["HOME"])
        self.assertEqual(rc, 1)
        self.assertEqual(self.sent, [])

    def test_timeout_sends_stop(self):
        def fake_time():
            self.clock[0] += 5.0
            return self.clock[0]
        gp.time.time = fake_time
        rc, out = self.run_main(["UP"])          # move_to=None：手臂不動
        self.assertEqual(rc, 1)
        self.assertIn("逾時", out)
        self.assertTrue(any("stopj(2.0)" in s for s in self.sent))


if __name__ == "__main__":
    unittest.main()
