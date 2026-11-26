# -*- coding: utf-8 -*-
"""
六軸暖機（robot/ur5_six_axis_warmup.py）的離線測試。
執行：python -m unittest discover -s tests -v

1. 一批的規劃：和鑑別程式同一個 plan_constvel，安全檢查（使用者視窗）通過，約 2.3 分鐘。
2. 模擬一批（剛體模型 + 已知摩擦）：block_result 算回 B、Tc 誤差 < 1%，溫度欄平均正確。
3. 收斂判定：最少時間、批數、只看 J0～J2、門檻、缺資料。
4. 溫度上限。
"""

import os
import sys
import math
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "robot"))
sys.path.insert(0, os.path.join(ROOT, "analysis", "six_axis"))

import ur5_six_axis_ident as six
import ur5_six_axis_warmup as wu
import sim_six_axis as sim

Q0 = list(six.POSES["SIX"]["q_deg"])


def hist(*rows):
    return [dict(B=list(r)) for r in rows]


class Plan(unittest.TestCase):
    def test_block_safe_and_short(self):
        segs = six.plan_constvel(wu.BLOCK_SPEEDS, Q0)
        self.assertEqual(len(segs), 2 * len(wu.BLOCK_SPEEDS))
        ok, lines, _ = six.safety_check(segs, "CONSTVEL", Q0)
        self.assertTrue(ok, "\n".join(lines))
        t = six.estimate_total_time(segs, Q0, "CONSTVEL", Q0)
        self.assertLess(t, 3 * 60)
        txt = six.build_urscript(segs, "CONSTVEL", end_deg=Q0)
        self.assertTrue(txt.isascii())


class SimBlock(unittest.TestCase):
    def test_recovers_friction(self):
        segs = six.plan_constvel(wu.BLOCK_SPEEDS, Q0)
        T, Q, QD, QDD = sim.trajectory(segs, "CONSTVEL", Q0)
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "block_01.csv")
            sim.write_csv(p, T, Q, QD, QDD, noise_a=0.0)
            r = wu.block_result(p)
        for j in range(6):
            self.assertAlmostEqual(r["B"][j] / sim.B_TRUE[j], 1.0, delta=0.01, msg=f"J{j} B")
            self.assertAlmostEqual(r["Tc"][j] / sim.TC_TRUE[j], 1.0, delta=0.01, msg=f"J{j} Tc")
            self.assertAlmostEqual(r["T"][j], 30.0)
            self.assertEqual(r["n_seg"][j], 2 * len(wu.BLOCK_SPEEDS))


class Converged(unittest.TestCase):
    W = (25.6, 19.5, 27.3, 2.6, 2.3, 3.0)

    def test_min_time(self):
        ok, why, _ = wu.is_converged(hist(self.W, self.W, self.W), elapsed_min=3.0)
        self.assertFalse(ok); self.assertIn("最少", why)

    def test_too_few_blocks(self):
        ok, _, _ = wu.is_converged(hist(self.W, self.W), elapsed_min=99)
        self.assertFalse(ok)

    def test_stable(self):
        rows = [(30.5, 23.4, 32.3, 3.0, 2.8, 4.0), self.W, [x * 1.01 for x in self.W], [x * 0.995 for x in self.W]]
        ok, why, sp = wu.is_converged(hist(*rows), elapsed_min=99)
        self.assertTrue(ok, why)
        self.assertLess(max(sp[:3]), wu.B_STABLE_FRAC)

    def test_big_joint_drifting(self):
        rows = [(27.4, 19.1, 27.5, 2.8, 2.4, 3.3), (25.9, 19.7, 27.2, 2.7, 2.4, 2.9), (25.6, 19.5, 27.3, 2.6, 2.3, 3.0)]
        ok, why, _ = wu.is_converged(hist(*rows), elapsed_min=99)   # 11/26 熱機三輪：J0 還在降 6.7%
        self.assertFalse(ok); self.assertIn("J0", why)

    def test_wrist_noise_ignored(self):
        rows = [self.W, (25.6, 19.5, 27.3, 2.9, 2.0, 4.0), self.W]
        ok, why, sp = wu.is_converged(hist(*rows), elapsed_min=99)
        self.assertTrue(ok, why)
        self.assertGreater(sp[5], 0.2)

    def test_missing(self):
        rows = [self.W, self.W, (float("nan"),) + self.W[1:]]
        ok, why, _ = wu.is_converged(hist(*rows), elapsed_min=99)
        self.assertFalse(ok); self.assertIn("缺資料", why)


class Temperature(unittest.TestCase):
    def test_abort(self):
        r = dict(T_max=[30.0, 41.0, float("nan"), 39.9, 40.0, 20.0])
        self.assertEqual(wu.temperature_abort(r, 40.0), [1, 4])


if __name__ == "__main__":
    unittest.main()
