# -*- coding: utf-8 -*-
"""
暖機程式（ur5_j0_warmup.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式、假的 rtde_receive（手臂停在 J0=180°）、
姿態模組設為 None。腳本位置自動判斷：有 robot/ 資料夾就用 robot/，否則用專案根目錄。

測試內容
--------
1. 語法正確。
2. 安全檢查：視窗與溫度上限都設定時通過；任一未設定、或範圍設定錯誤時 FAIL。
3. URScript：縮排規則、只含 ASCII、movej 只改 J0、迴圈有固定次數。
4. 摩擦指標與暖機完成判定（純函式）。
5. main()：溫度上限未設定、RTDE 失敗、起始位置過遠時，都在送出前中止；QUICK 端對端會送出
   一批與一次結束定位，並產生自動命名的 CSV 與摘要檔。
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

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
ROBOT = os.path.join(ROOT, "robot") if os.path.isdir(os.path.join(ROOT, "robot")) else ROOT
WARM = os.path.join(ROBOT, "ur5_j0_warmup.py")

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
    spec = importlib.util.spec_from_file_location(name, WARM)
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
        FakeRTDE.q = list(HOME_Q)

    def tearDown(self):
        socket.socket = self._sock
        builtins.input = self._input
        FakeRTDE.q = list(HOME_Q)


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


class TestPieces(Base):
    def test_compile(self):
        py_compile.compile(WARM, doraise=True)

    def test_safety(self):
        m = load("w_safe")
        args = (m.WARMUP_LO_DEG, m.WARMUP_HI_DEG, m.END_DEG, m.WARMUP_VEL, m.WARMUP_ACCEL,
                m.QD_MAX, m.QDD_CHECK_LIMIT)
        self.assertTrue(m.safety_check(*args, 65, 300, 50.0, m.SAFETY_MARGIN_DEG)[0])
        self.assertFalse(m.safety_check(*args, None, None, 50.0, m.SAFETY_MARGIN_DEG)[0], "視窗未設定應 FAIL")
        self.assertFalse(m.safety_check(*args, 65, 300, None, m.SAFETY_MARGIN_DEG)[0], "溫度上限未設定應 FAIL")
        self.assertFalse(m.safety_check(*args, 160, 300, 50.0, m.SAFETY_MARGIN_DEG)[0], "超出視窗應 FAIL")
        bad = (205.0, 155.0) + args[2:]
        self.assertFalse(m.safety_check(*bad, 65, 300, 50.0, m.SAFETY_MARGIN_DEG)[0], "範圍顛倒應 FAIL")

    def test_urscript(self):
        m = load("w_script")
        s = m.build_block_urscript(9, m.WARMUP_LO_DEG, m.WARMUP_HI_DEG, HOME_Q, 0, m.WARMUP_VEL, m.WARMUP_ACCEL)
        e = m.build_end_urscript(m.END_DEG, HOME_Q, 0, m.WARMUP_VEL, m.WARMUP_ACCEL)
        for txt in (s, e):
            self.assertTrue(urscript_indent_ok(txt))
            self.assertTrue(all(ord(ch) < 128 for ch in txt))
            for line in txt.splitlines():
                if "movej(" in line:
                    vals = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                    for j in range(1, 6):
                        self.assertAlmostEqual(vals[j], HOME_Q[j], places=5)
        self.assertEqual(s.count("movej("), 3)
        self.assertIn("while k < 9:", s)

    def test_block_metrics(self):
        m = load("w_metric")
        v = m.WARMUP_VEL
        qd = np.concatenate([np.full(100, v), np.full(50, 0.05), np.full(100, -v)])
        cur = np.concatenate([np.full(100, 0.8), np.full(50, 0.3), np.full(100, -0.9)])
        temp = np.linspace(30.0, 30.5, len(qd))
        f, r, t = m.block_metrics(qd, cur, temp, v, m.CRUISE_TOL, m.KT_OUT)
        self.assertAlmostEqual(f, 0.8 * m.KT_OUT, places=6)
        self.assertAlmostEqual(r, 0.9 * m.KT_OUT, places=6)
        self.assertAlmostEqual(t, 30.5, places=6)

    def test_is_warm(self):
        m = load("w_iswarm")
        kw = dict(min_minutes=5, temp_window_min=5, temp_stable_c=0.1, n_stable_blocks=3,
                  friction_stable_frac=0.01)
        rising = [dict(t_min=2 * k, temp=28.6 + 0.4 * k, fwd=12.0 - 0.2 * k, rev=13.0 - 0.2 * k) for k in range(1, 6)]
        self.assertFalse(m.is_warm(rising, 10, **kw)[0], "溫度仍在上升不該判定完成")
        stable = [dict(t_min=2 * k, temp=31.77, fwd=12.00 + 0.01 * (k % 2), rev=13.00, ) for k in range(1, 8)]
        self.assertTrue(m.is_warm(stable, 14, **kw)[0])
        self.assertFalse(m.is_warm(stable[:2], 4, **kw)[0], "未滿最少時間不該判定完成")
        drifting = [dict(t_min=2 * k, temp=31.77, fwd=12.0 + 0.1 * k, rev=13.0) for k in range(1, 8)]
        self.assertFalse(m.is_warm(drifting, 14, **kw)[0], "摩擦指標仍在變不該判定完成")


class TestMain(Base):
    def _run(self, m):
        sent = []
        m.send_urscript = lambda ip, txt, *a, **k: sent.append(txt)
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            with contextlib.redirect_stdout(io.StringIO()):
                m.main()
            files = [f for _, _, fs in os.walk(tmp) for f in fs]
        return sent, files

    def test_temp_limit_unset_aborts(self):
        m = load("w_m1")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG, m.TEMP_ABORT_C = 65, 300, None
        self.assertEqual(self._run(m)[0], [])

    def test_rtde_failure_aborts(self):
        m = load("w_m2")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG, m.TEMP_ABORT_C = 65, 300, 50.0

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m)[0], [])

    def test_far_start_aborts(self):
        m = load("w_m3")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG, m.TEMP_ABORT_C = 65, 300, 50.0
        FakeRTDE.q = [math.radians(100.0)] + HOME_Q[1:]
        self.assertEqual(self._run(m)[0], [])

    def test_quick_end_to_end(self):
        m = load("w_m4")
        m.J0_SAFE_MIN_DEG, m.J0_SAFE_MAX_DEG, m.TEMP_ABORT_C = 65, 300, 50.0
        m.TEST_STAGE = "QUICK"
        m.wait_for_motion_complete = lambda logger, planned_s, **k: ("done", planned_s, planned_s)
        sent, files = self._run(m)
        self.assertEqual(len(sent), 2, "QUICK 應送出 1 批 + 1 次結束定位")
        self.assertIn("while k < 2:", sent[0])
        self.assertTrue(sent[1].startswith("def j0_warmup_end():"))
        self.assertTrue(any(f.startswith("warmup_data_quick_") and f.endswith(".csv") for f in files))
        self.assertTrue(any(f.startswith("warmup_info_quick_") and f.endswith(".txt") for f in files))


if __name__ == "__main__":
    unittest.main(verbosity=2)
