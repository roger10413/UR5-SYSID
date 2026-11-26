# -*- coding: utf-8 -*-
"""
取放情境鑑別（robot/ur5_pick_place_ident.py）的離線測試。
執行：python -m unittest tests.test_pick_place_offline -v

1. 幾何：下方姿態 IK（法蘭水平位置、姿態不變，高度降 APPROACH_DZ）、預設兩站同高、工具朝下一致。
2. 規劃：動作順序、movej 定速段速度、定速窗口落在最高速度的定速段內、整數馬達圈；FULL 時間與速度上限。
3. URScript：只含 ASCII、movej/movel 數量。
4. 名目剛體模型 = 模擬用的 Newton-Euler 扣掉轉子。
5. 模擬取放（剛體 + 已知摩擦）：搬運段定速分析算回 B、Tc；全週期回歸算回 B、Tc、轉子。
"""
import math
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "robot"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "six_axis"))

import ur5_pick_place_ident as pp
import sim_six_axis as sim
import sim_pick_place as spp

OPEN = ([-999.0] * 6, [999.0] * 6)


class Geometry(unittest.TestCase):
    def test_lowered_pose(self):
        for name in ("A", "B"):
            up = pp.STATIONS[name]
            dn = pp.lowered_pose(up, pp.APPROACH_DZ)
            p0, p1 = pp.flange_pos(up), pp.flange_pos(dn)
            self.assertAlmostEqual(p0[0], p1[0], places=6)
            self.assertAlmostEqual(p0[1], p1[1], places=6)
            self.assertAlmostEqual(p0[2] - p1[2], pp.APPROACH_DZ, places=6)
            R0 = pp.six.fk_frames(np.radians(up))[6][:3, :3]
            R1 = pp.six.fk_frames(np.radians(dn))[6][:3, :3]
            self.assertLess(np.abs(R0 - R1).max(), 1e-9)

    def test_default_stations_same_height_tool_down(self):
        A, B = pp.STATIONS["A"], pp.STATIONS["B"]
        self.assertAlmostEqual(pp.flange_pos(A)[2], pp.flange_pos(B)[2], places=3)
        for q in (A, B):
            self.assertAlmostEqual(q[1] + q[2] + q[3], -90.0, places=6)
            self.assertEqual(q[4], -90.0)
            R = pp.six.fk_frames(np.radians(q))[6][:3, :3]
            self.assertAlmostEqual(R[2, 2], -1.0, places=6)          # 法蘭 z 軸朝下


class Plan(unittest.TestCase):
    def test_sequence(self):
        moves, poses = pp.plan_pick_place([0.05, 0.1])
        kinds = [m["kind"] for m in moves]
        self.assertEqual(kinds, (["movel", "dwell", "movel", "movej"] * 2) * 2)
        mj = [m for m in moves if m["kind"] == "movej"]
        self.assertEqual([m["d"] for m in mj], [1, -1, 1, -1])
        for a, b in zip(moves[:-1], moves[1:]):                       # 動作首尾相接
            self.assertTrue(np.allclose(a["target_deg"], b["start_deg"]))
        self.assertTrue(np.allclose(moves[-1]["target_deg"], poses["A_up"]))
        m = mj[0]
        lead = max(abs(x) for x in m["qd_plateau"])
        self.assertAlmostEqual(lead, 0.05)
        self.assertAlmostEqual(m["qd_plateau"][0], 0.05)              # J0 是帶頭軸

    def test_windows_inside_plateau(self):
        W = pp.transfer_windows()
        self.assertIsNone(W[3]); self.assertIsNone(W[4])               # J3、J4 量不到
        rev = 360.0 / pp.GEAR_RATIO
        A, B = pp.STATIONS["A"], pp.STATIONS["B"]
        lead = math.radians(max(abs(b - a) for a, b in zip(A, B)))
        vmax, acc = max(pp.TRANSFER_SPEEDS), pp.TRANSFER_ACCEL
        f = (vmax ** 2 / (2 * acc) + vmax * pp.SETTLE_TIME) / lead
        for j in (0, 1, 2, 5):
            lo, hi = W[j]
            n = (hi - lo) / rev
            self.assertAlmostEqual(n, round(n), places=6)
            self.assertGreaterEqual(round(n), 1)
            lam = sorted(((x - A[j]) / (B[j] - A[j]) for x in (lo, hi)))
            self.assertGreaterEqual(lam[0], f - 1e-9); self.assertLessEqual(lam[1], 1 - f + 1e-9)

    def test_full_safe_and_time(self):
        moves, poses = pp.plan_pick_place(pp.stage_speeds("FULL"))
        ok, lines, _ = pp.safety_check(moves, poses, *OPEN)
        self.assertTrue(ok, "\n".join(lines))
        self.assertLess(pp.estimate_total_time(moves), 10 * 60)

    def test_unset_window_refuses(self):
        moves, poses = pp.plan_pick_place([0.05])
        ok, lines, _ = pp.safety_check(moves, poses)
        self.assertFalse(ok)
        self.assertTrue(any("視窗未設定" in l for l in lines))

    def test_too_fast_refuses(self):
        moves, poses = pp.plan_pick_place([0.40])
        ok, _, _ = pp.safety_check(moves, poses, *OPEN)
        self.assertFalse(ok)


class Script(unittest.TestCase):
    def test_ascii_and_counts(self):
        moves, poses = pp.plan_pick_place([0.05, 0.12])
        txt = pp.build_urscript(moves, poses)
        self.assertTrue(txt.isascii())
        self.assertEqual(txt.count("movej("), 1 + 4)
        self.assertEqual(txt.count("movel("), 8)
        self.assertTrue(txt.startswith("def pick_place_ident_program():"))
        self.assertTrue(txt.rstrip().endswith("end"))


class Model(unittest.TestCase):
    def test_rigid_body_matches_sim(self):
        rng = np.random.default_rng(1)
        for _ in range(5):
            q, qd, qdd = rng.uniform(-2, 2, 6), rng.uniform(-0.3, 0.3, 6), rng.uniform(-1, 1, 6)
            a = pp.rigid_body_torque(q, qd, qdd, sim.TOOL_KG)
            b = sim.rnea(q, qd, qdd) - sim.ROTOR * qdd
            self.assertLess(np.abs(a - b).max(), 1e-9)


class SimPickPlace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        _, _, cls.path = spp.run([0.08, 0.20, 0.30], os.path.join(cls.tmp.name, "sim.csv"))
        cls.W = pp.transfer_windows()
        cls.D = pp.load_csv(cls.path)

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_transfer_recovers_friction(self):
        r = pp.analyze_transfer(self.D, self.W)
        for j in (0, 1, 2, 5):
            self.assertEqual((r[j]["n_pos"], r[j]["n_neg"]), (3, 3), f"J{j}")
            self.assertAlmostEqual(r[j]["B"] / sim.B_TRUE[j], 1.0, delta=0.01, msg=f"J{j} B")
            self.assertAlmostEqual(r[j]["Tc"] / sim.TC_TRUE[j], 1.0, delta=0.01, msg=f"J{j} Tc")
        for j in (3, 4):
            self.assertIn("skip", r[j])

    def test_full_cycle_recovers(self):
        r = pp.analyze_full_cycle(self.D, sim.TOOL_KG, step=3)
        for j in (0, 1, 2, 3, 5):
            self.assertAlmostEqual(r[j]["B"] / sim.B_TRUE[j], 1.0, delta=0.02, msg=f"J{j} B")
            self.assertAlmostEqual(r[j]["Tc"] / sim.TC_TRUE[j], 1.0, delta=0.02, msg=f"J{j} Tc")
            self.assertAlmostEqual(r[j]["J_r"] / sim.ROTOR[j], 1.0, delta=0.1, msg=f"J{j} J_r")
        self.assertIn("skip", r[4])

    def test_summary(self):
        lines = pp.summary_lines(pp.analyze_transfer(self.D, self.W), {}, self.W)
        self.assertTrue(any(l.startswith("J3") and "沒有定速窗口" in l for l in lines))


if __name__ == "__main__":
    unittest.main()
