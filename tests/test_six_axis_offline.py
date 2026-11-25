# -*- coding: utf-8 -*-
"""
六軸同時鑑別（robot/ur5_six_axis_ident.py）的離線測試。
執行：python -m unittest discover -s tests -v

1. 條件：同一條件六軸頻率互不相同、整數週期；改成重複頻率會 FAIL。
2. 定速規劃：六軸同時進出窗口、路徑是同一條直線、去回交替。
3. 速度調變規劃：調變段長度 = 巡航 × T_MOD；起點錯開。
4. 安全檢查：視窗未設定 FAIL；寬視窗 OK；法蘭高度門檻提高就 FAIL。
5. URScript：def/end 頂格、六軸都有角度保護、只含 ASCII。
6. 名目重力：DIAG 時 J1 約 −23.3 N·m、J0 為 0。
7. 模擬（剛體模型 + 已知摩擦）：MID 階段算回 B、Tc 誤差 < 1%，J 誤差 < 10%。
"""

import os
import sys
import math
import unittest
import tempfile
import importlib.util

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "robot"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "six_axis"))

import numpy as np
import ur5_six_axis_ident as six
import sim_six_axis as sim

Q0 = list(six.POSES["SIX"]["q_deg"])
WIDE = dict(safe_min=(-1e9,) * 6, safe_max=(1e9,) * 6)


class Conditions(unittest.TestCase):
    def test_default_ok(self):
        ok, lines = six.check_conditions()
        self.assertTrue(ok, lines)
        for ci in range(len(six.TIERS)):
            fs = [six.joint_conditions(j)[ci][0] for j in range(6)]
            self.assertEqual(len(set(fs)), 6)
            for f in fs:
                self.assertAlmostEqual(f * six.T_MOD, round(f * six.T_MOD))

    def test_duplicate_freq_fails(self):
        ok, lines = six.check_conditions(offsets=(0, 0, 0.25, -0.25, 0.75, 0.5))
        self.assertFalse(ok)
        self.assertTrue(any("重複" in s for s in lines))

    def test_low_speed_fails(self):
        ok, _ = six.check_conditions(v_cruise=(0.05,) * 6)
        self.assertFalse(ok)


class ConstvelPlan(unittest.TestCase):
    def setUp(self):
        self.segs = six.plan_constvel([0.02, 0.12], Q0)

    def test_alternate_and_sync(self):
        self.assertEqual([s["d"] for s in self.segs], [1, -1, 1, -1])
        for s in self.segs:
            pre = [abs(a - b) for a, b in zip(s["start_deg"], s["win_start_deg"])]
            self.assertTrue(all(abs(p - pre[0]) < 1e-9 for p in pre))          # 六軸同時進窗口
            for j in range(6):
                self.assertEqual(s["dj"][j], s["d"] * six.PATH_SIGN[j])

    def test_same_line(self):
        for s in self.segs:
            for key in ("start_deg", "win_start_deg", "win_end_deg", "stop_deg"):
                dev = [(x - c) * sg for x, c, sg in zip(s[key], Q0, six.PATH_SIGN)]
                self.assertTrue(all(abs(d - dev[0]) < 1e-9 for d in dev), key)   # 都在 q = 中心 + λ·S 上

    def test_window_width(self):
        w = six.N_MOTOR_REVS_WINDOW * 360 / six.GEAR_RATIO
        for lo, hi in six.window_bounds(Q0):
            self.assertAlmostEqual(hi - lo, w)


class VelmodPlan(unittest.TestCase):
    def test_lengths_offsets(self):
        segs = six.plan_velmod([0, 1], 3, Q0)
        self.assertEqual(len(segs), 2 * 3 * 2)
        for s in segs:
            for j in range(6):
                L = abs(s["mod_end_deg"][j] - s["mod_start_deg"][j])
                self.assertAlmostEqual(L, math.degrees(six.V_CRUISE[j] * six.T_MOD))
        offs = sorted({round(s["offset_deg"], 6) for s in segs})
        self.assertEqual(len(offs), 3)
        self.assertAlmostEqual(sum(offs), 0.0)


class Safety(unittest.TestCase):
    def test_windows_required(self):
        segs = six.plan_constvel([0.05], Q0)
        ok, lines, _ = six.safety_check(segs, "CONSTVEL", Q0, safe_min=(None,) * 6, safe_max=(None,) * 6)
        self.assertFalse(ok)
        self.assertTrue(any("視窗未設定" in s for s in lines))

    def test_full_ok_with_configured_windows(self):
        for mode in ("CONSTVEL", "VELMOD"):
            c = six.stage_config(mode, "FULL")
            segs = six.plan_constvel(c, Q0) if mode == "CONSTVEL" else six.plan_velmod(c[0], c[1], Q0)
            ok, lines, _ = six.safety_check(segs, mode, Q0, n_check=11,
                                            safe_min=six.SAFE_MIN_DEG, safe_max=six.SAFE_MAX_DEG)
            self.assertTrue(ok, (mode, lines))

    def test_full_ok_with_wide_windows(self):
        for mode in ("CONSTVEL", "VELMOD"):
            c = six.stage_config(mode, "FULL")
            segs = six.plan_constvel(c, Q0) if mode == "CONSTVEL" else six.plan_velmod(c[0], c[1], Q0)
            ok, lines, (lo, hi) = six.safety_check(segs, mode, Q0, n_check=11, **WIDE)
            self.assertTrue(ok, (mode, lines))

    def test_window_too_small(self):
        segs = six.plan_constvel([0.05], Q0)
        ok, _, _ = six.safety_check(segs, "CONSTVEL", Q0, safe_min=[q - 5 for q in Q0], safe_max=[q + 5 for q in Q0])
        self.assertFalse(ok)

    def test_flange_limit(self):
        segs = six.plan_constvel([0.05], Q0)
        ok, lines, _ = six.safety_check(segs, "CONSTVEL", Q0, flange_z_min=0.9, n_check=11, **WIDE)
        self.assertFalse(ok)

    def test_start_position(self):
        segs = six.plan_constvel([0.05], Q0)
        self.assertTrue(six.check_start_position(Q0, segs[0]["start_deg"])[0])
        self.assertFalse(six.check_start_position([q + 40 for q in Q0], segs[0]["start_deg"])[0])


class Script(unittest.TestCase):
    def check(self, txt):
        lines = txt.strip().split("\n")
        self.assertTrue(lines[0].startswith("def "))
        self.assertEqual(lines[-1], "end")
        for ln in lines[1:-1]:
            self.assertTrue(ln.startswith("  "))
        txt.encode("ascii")

    def test_constvel(self):
        txt = six.build_urscript(six.plan_constvel([0.05], Q0), "CONSTVEL", end_deg=Q0)
        self.check(txt)
        w = [ln for ln in txt.split("\n") if ln.strip().startswith("while")]
        self.assertEqual(len(w), 2)
        for ln in w:
            self.assertEqual(sum(f"q_now[{j}]" in ln for j in range(6)), 6)

    def test_velmod(self):
        txt = six.build_urscript(six.plan_velmod([0], 1, Q0), "VELMOD", end_deg=Q0)
        self.check(txt)
        qd = [ln for ln in txt.split("\n") if "qd_t = [" in ln]
        self.assertEqual(len(qd), 2)
        self.assertEqual(qd[0].count("sin("), 6)


class Gravity(unittest.TestCase):
    def test_diag(self):
        g = six.gravity_torque(np.radians([180, -90, 90, -90, -90, 0]), 1.30)
        self.assertAlmostEqual(g[0], 0.0, places=9)
        self.assertAlmostEqual(g[1], -23.3, delta=0.1)

    def test_matches_rnea(self):
        q = np.radians([30, -60, 40, -120, 80, 20]); z = np.zeros(6)
        np.testing.assert_allclose(six.gravity_torque(q, 1.30), sim.rnea(q, z, z), atol=1e-9)


class Simulation(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()

    @classmethod
    def tearDownClass(cls):
        cls.tmp.cleanup()

    def test_constvel_mid(self):
        _, res, _ = sim.run("CONSTVEL", "MID", self.tmp.name, noise_a=0.0)
        for j in range(6):
            self.assertAlmostEqual(res[j]["B"] / sim.B_TRUE[j], 1.0, delta=0.01, msg=f"J{j} B")
            self.assertAlmostEqual(res[j]["Tc"] / sim.TC_TRUE[j], 1.0, delta=0.01, msg=f"J{j} Tc")

    def test_velmod_mid(self):
        _, res, Jt = sim.run("VELMOD", "MID", self.tmp.name, noise_a=0.0)
        ci_mid, _ = six.stage_config("VELMOD", "MID")
        for j in range(6):
            self.assertEqual(sorted(res[j]["per"]), ci_mid, f"J{j} 每個條件都要有結果")
            for ci, pp in res[j]["per"].items():
                self.assertEqual(set(pp), {"pos", "neg"})
                Jm = np.mean([pp[t]["J"] for t in pp])
                self.assertAlmostEqual(Jm / Jt[j], 1.0, delta=0.10, msg=f"J{j} 條件 {ci}")


if __name__ == "__main__":
    unittest.main()
