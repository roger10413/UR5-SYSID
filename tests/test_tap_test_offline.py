# -*- coding: utf-8 -*-
"""
敲擊測試（robot/ur5_tap_test_record.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線）、假的 rtde_receive（手臂停在姿態不動）。

測試內容
--------
1. 語法正確；程式碼裡沒有任何送指令的東西（不 import socket，沒有 speedj、movej、stopj、servoj、30002、RTDEControl）；
   假的 RTDE 物件只被呼叫讀取方法（get 開頭）與 disconnect。
2. 設定檢查：敲的方向未知、記錄時間超出範圍、提示時間太密或太晚 → FAIL。
3. 姿態（UP、DIAG、EXT45 各自 ±0.5°）、靜止（任一軸速度超過上限或讀不到 → FAIL）、
   機器人狀態（robot_mode 須為 7、safety_mode 須為 1，讀不到 → FAIL）。
4. main()：負載未設定、敲的方向未知、姿態不符、手臂在動、機器人模式不對、安全模式不對、未輸入姿態名稱、
   負載不符、RTDE 連線失敗時都在記錄前中止，不產生資料檔。
5. main() 端對端（記錄 1 s）：資料夾 taptest_up_p130_fb_…，CSV 有六軸角度、電流、工具加速度計、模式欄位，
   摘要寫出姿態與負載；記錄中安全模式變成保護性停止時警告「這次資料不可用」；
   加速度計讀不到時該欄留空並警告。
6. 分析（合成資料：4 下敲擊，餘振 8 Hz（加速度計 x、J1 角度、J1 電流）與 13 Hz（加速度計 y、J0 角度），加雜訊）：
   找到 4 下、時間正確，峰值頻率誤差 0.3 Hz 內；沒有敲擊的資料當作背景處理、不當機；
   加速度計欄位是空的時改用電流找敲擊。
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
TT = os.path.join(ROBOT, "ur5_tap_test_record.py")

POSE_Q = {"UP": [math.radians(x) for x in (180.0, -90.0, 0.0, -90.0, -90.0, 0.0)],
          "DIAG": [math.radians(x) for x in (180.0, -90.0, 90.0, -90.0, -90.0, 0.0)],
          "EXT45": [math.radians(x) for x in (180.0, -45.0, 45.0, -90.0, -90.0, 0.0)]}
PAYLOAD = 1.3


class FakeRTDE:
    q = list(POSE_Q["UP"])
    qd = [0.0] * 6
    robot_mode = 7
    safety_mode = 1
    safety_after = None            # (列數, 新的 safety_mode)：記錄中途改變安全模式
    acc = [0.1, -0.2, 9.8]
    calls = []
    n_safety = 0

    def __init__(self, ip): pass

    def __getattribute__(self, name):
        if not name.startswith("_"):
            FakeRTDE.calls.append(name)
        return object.__getattribute__(self, name)

    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return list(FakeRTDE.qd)
    def getActualCurrent(self): return [0.1] * 6
    def getTargetQ(self): return list(FakeRTDE.q)
    def getJointTemperatures(self): return [30.0] * 6
    def getTargetMoment(self): return [0.0] * 6
    def getTargetCurrent(self): return [0.0] * 6
    def getJointControlOutput(self): return [0.0] * 6
    def getActualToolAccelerometer(self):
        if FakeRTDE.acc is None:
            raise RuntimeError("沒有加速度計")
        return list(FakeRTDE.acc)
    def getActualTCPForce(self): return [0.0] * 6
    def getRobotMode(self): return FakeRTDE.robot_mode

    def getSafetyMode(self):
        FakeRTDE.n_safety += 1
        if FakeRTDE.safety_after is not None and FakeRTDE.n_safety > FakeRTDE.safety_after[0]:
            return FakeRTDE.safety_after[1]
        return FakeRTDE.safety_mode

    def disconnect(self): pass


def reset_fake():
    FakeRTDE.q, FakeRTDE.qd = list(POSE_Q["UP"]), [0.0] * 6
    FakeRTDE.robot_mode, FakeRTDE.safety_mode, FakeRTDE.safety_after = 7, 1, None
    FakeRTDE.acc, FakeRTDE.calls, FakeRTDE.n_safety = [0.1, -0.2, 9.8], [], 0


def load(name):
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    spec = importlib.util.spec_from_file_location(name, TT)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    mod.PAYLOAD_KG = PAYLOAD
    mod.TAP_LABEL = "fb"
    mod.RECORD_S, mod.TAP_TIMES, mod.COUNTDOWN_S, mod.RECORD_S_MAX = 1.0, (), 0, 120.0
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(default="UP", payload=str(PAYLOAD)):
    return lambda prompt="", *a, **k: payload if "負載質量" in str(prompt) else default


class Base(unittest.TestCase):
    def setUp(self):
        self._sock, self._input = socket.socket, builtins.input
        socket.socket = no_network
        builtins.input = answers()
        reset_fake()

    def tearDown(self):
        socket.socket, builtins.input = self._sock, self._input
        reset_fake()
        sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE


class TestReadOnly(Base):
    def test_compile(self):
        py_compile.compile(TT, doraise=True)

    def test_no_command_code(self):
        with open(TT, encoding="utf-8") as f:
            src = f.read()
        for word in ("import socket", "speedj", "movej", "stopj", "servoj", "30002", "RTDEControl", "rtde_control",
                     "sendall", "send_urscript"):
            self.assertNotIn(word, src, word)

    def test_only_getters_called(self):
        m = load("tt_getters")
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            with contextlib.redirect_stdout(io.StringIO()):
                m.main()
        used = set(FakeRTDE.calls)
        self.assertTrue(used, "應該有讀取")
        for name in used:
            self.assertTrue(name.startswith("get") or name == "disconnect", name)


class TestChecks(Base):
    def test_record_settings(self):
        m = load("tt_set")
        self.assertTrue(m.check_record_settings("fb", 14.0, (2.0, 5.0, 8.0, 11.0))[0])
        self.assertTrue(m.check_record_settings("none", 14.0, ())[0])
        self.assertTrue(m.check_record_settings(m.TAP_LABEL, 14.0, (2.0, 5.0, 8.0, 11.0))[0])
        self.assertFalse(m.check_record_settings("xx", 14.0, (2.0,))[0], "未知方向")
        self.assertFalse(m.check_record_settings("fb", 0.5, ())[0], "記錄太短")
        self.assertFalse(m.check_record_settings("fb", 500.0, ())[0], "記錄太長")
        self.assertFalse(m.check_record_settings("fb", None, ())[0])
        self.assertFalse(m.check_record_settings("fb", 14.0, (0.5, 5.0))[0], "第一下太早")
        self.assertFalse(m.check_record_settings("fb", 14.0, (2.0, 3.0))[0], "間隔太密")
        self.assertFalse(m.check_record_settings("fb", 14.0, (2.0, 13.0))[0], "最後一下太晚")
        self.assertFalse(m.check_record_settings("fb", 14.0, (2.0, "x"))[0])

    def test_defaults(self):
        spec = importlib.util.spec_from_file_location("tt_defaults", TT)
        m = importlib.util.module_from_spec(spec)
        sys.modules["rtde_receive"] = types.ModuleType("rtde_receive")
        sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE
        spec.loader.exec_module(m)
        self.assertIsNone(m.PAYLOAD_KG)
        self.assertEqual(m.TAP_LABEL, "none")
        self.assertEqual(m.POSE, "UP")
        self.assertEqual(set(m.POSES), {"UP", "DIAG", "EXT45"})
        self.assertTrue(m.check_record_settings(m.TAP_LABEL, m.RECORD_S, m.TAP_TIMES)[0])
        self.assertTrue(m.check_record_settings("fb", m.RECORD_S, m.TAP_TIMES)[0])

    def test_check_pose(self):
        m = load("tt_pose")
        for pose, q in POSE_Q.items():
            self.assertTrue(m.check_pose(q, pose)[0], pose)
        self.assertFalse(m.check_pose(POSE_Q["DIAG"], "UP")[0])
        self.assertFalse(m.check_pose(POSE_Q["UP"][:2] + [math.radians(0.6)] + POSE_Q["UP"][3:], "UP")[0])
        self.assertFalse(m.check_pose([float("nan")] + POSE_Q["UP"][1:], "UP")[0])
        self.assertFalse(m.check_pose(None, "UP")[0])
        self.assertFalse(m.check_pose(POSE_Q["UP"], "XX")[0])

    def test_check_still_and_state(self):
        m = load("tt_state")
        self.assertTrue(m.check_still([0.0] * 6)[0])
        self.assertTrue(m.check_still([0.001, -0.001, 0, 0, 0, 0])[0])
        self.assertFalse(m.check_still([0.0, 0.01, 0, 0, 0, 0])[0], "在動")
        self.assertFalse(m.check_still(None)[0])
        self.assertFalse(m.check_still([0.0] * 5)[0])
        self.assertFalse(m.check_still([float("nan")] + [0.0] * 5)[0])
        self.assertTrue(m.check_robot_state(7, 1)[0])
        self.assertFalse(m.check_robot_state(5, 1)[0], "剎車未放開")
        self.assertFalse(m.check_robot_state(7, 3)[0], "保護性停止")
        self.assertFalse(m.check_robot_state(None, 1)[0])
        self.assertFalse(m.check_robot_state(7, None)[0])
        self.assertFalse(m.check_robot_state(7.5, 1)[0])

    def test_payload(self):
        m = load("tt_pay")
        self.assertFalse(m.check_payload_setting(None)[0])
        self.assertFalse(m.check_payload_setting(-1)[0])
        self.assertTrue(m.check_payload_entry(1.3, "1.30")[0])
        self.assertFalse(m.check_payload_entry(1.3, "2.30")[0])
        self.assertFalse(m.check_payload_entry(1.3, "")[0])
        self.assertEqual(m.payload_tag(1.3), "p130")


class TestMain(Base):
    def _run(self, m):
        out, files = io.StringIO(), {}
        with tempfile.TemporaryDirectory() as tmp:
            m.OUTPUT_DIR = tmp
            with contextlib.redirect_stdout(out):
                m.main()
            for root, _, names in os.walk(tmp):
                for n in names:
                    with open(os.path.join(root, n), encoding="utf-8") as f:
                        files[os.path.join(os.path.basename(root), n)] = f.read()
        self.last_out = out.getvalue()
        return files

    def test_aborts_before_recording(self):
        cases = [("payload_none", "負載設定未通過"), ("label", "記錄設定未通過"), ("pose_name", "姿態檢查未通過"),
                 ("pose", "姿態檢查未通過"), ("moving", "靜止檢查未通過"), ("robot_mode", "機器人狀態檢查未通過"),
                 ("safety_mode", "機器人狀態檢查未通過"), ("name", "已取消"), ("payload", "負載確認未通過"),
                 ("rtde", "RTDE 連線失敗")]
        for k, (case, reason) in enumerate(cases):
            reset_fake()
            builtins.input = answers()
            m = load(f"tt_abort{k}")
            if case == "payload_none":
                m.PAYLOAD_KG = None
            elif case == "label":
                m.TAP_LABEL = "sideways"
            elif case == "pose_name":
                m.POSE = "XX"
            elif case == "pose":
                FakeRTDE.q = list(POSE_Q["DIAG"])
            elif case == "moving":
                FakeRTDE.qd = [0.0, 0.05, 0, 0, 0, 0]
            elif case == "robot_mode":
                FakeRTDE.robot_mode = 5
            elif case == "safety_mode":
                FakeRTDE.safety_mode = 3
            elif case == "name":
                builtins.input = answers("yes")
            elif case == "payload":
                builtins.input = answers(payload="2.30")
            elif case == "rtde":
                class Failing:
                    def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
                sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
            files = self._run(m)
            self.assertEqual(files, {}, case)
            self.assertIn(reason, self.last_out, case)
            self.assertNotIn("開始記錄", self.last_out, case)
            sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE

    def test_end_to_end(self):
        for pose in ("UP", "EXT45"):
            reset_fake()
            FakeRTDE.q = list(POSE_Q[pose])
            builtins.input = answers(pose)
            m = load(f"tt_e2e_{pose}")
            m.POSE = pose
            files = self._run(m)
            data = [n for n in files if os.path.basename(n).startswith(f"taptest_data_{pose.lower()}_p130_fb_")]
            info = [n for n in files if os.path.basename(n).startswith(f"taptest_info_{pose.lower()}_p130_fb_")]
            self.assertEqual((len(data), len(info)), (1, 1), list(files))
            self.assertTrue(os.path.dirname(data[0]).startswith(f"taptest_{pose.lower()}_p130_fb_"))
            rows = list(csv.DictReader(io.StringIO(files[data[0]])))
            self.assertGreater(len(rows), 20)
            for key in ("timestamp", "actual_q_0", "actual_q_5", "actual_current_1", "target_q_2",
                        "joint_control_output_3", "tool_acc_x", "tool_acc_z", "tcp_force_rz", "robot_mode", "safety_mode"):
                self.assertIn(key, rows[0], key)
            self.assertEqual(rows[0]["tool_acc_z"], "9.8")
            self.assertEqual(rows[0]["robot_mode"], "7")
            self.assertIn(f"姿態: {pose}", files[info[0]])
            self.assertIn("負載: 1.300 kg", files[info[0]])
            self.assertNotIn("[警告] 記錄中", files[info[0]])

    def test_none_label_and_cues(self):
        m = load("tt_none")
        m.TAP_LABEL, m.RECORD_S, m.TAP_TIMES = "none", 3.2, (1.0,)
        files = self._run(m)
        self.assertTrue(any("_none_" in n for n in files))
        self.assertIn("這次不敲", self.last_out)
        self.assertNotIn("→ 敲", self.last_out)
        reset_fake()
        m = load("tt_cue")
        m.RECORD_S, m.TAP_TIMES = 3.2, (1.0,)
        self._run(m)
        self.assertEqual(self.last_out.count("→ 敲"), 1)

    def test_protective_stop_warns(self):
        m = load("tt_pstop")
        FakeRTDE.safety_after = (10, 3)
        files = self._run(m)
        self.assertIn("這次資料不可用", self.last_out)
        info = next(v for n, v in files.items() if "taptest_info_" in n)
        self.assertIn("安全模式不是正常", info)

    def test_missing_accelerometer(self):
        m = load("tt_noacc")
        FakeRTDE.acc = None
        files = self._run(m)
        rows = list(csv.DictReader(io.StringIO(next(v for n, v in files.items() if "taptest_data_" in n))))
        self.assertEqual(rows[0]["tool_acc_x"], "")
        self.assertIn("工具加速度計有", self.last_out)


# ------------------------------------------------------------
# 合成資料
# ------------------------------------------------------------

def synth_csv(path, taps=(2.0, 5.0, 8.0, 11.0), dur=14.0, dt=0.008, acc=True, seed=0):
    """敲擊後的餘振：8 Hz（加速度計 x、J1 角度、J1 電流）、13 Hz（加速度計 y、J0 角度），阻尼比 3%，加白雜訊。"""
    rng = np.random.default_rng(seed)
    n = int(round(dur / dt))
    t = np.arange(n) * dt

    def ring(f, amp):
        y = np.zeros(n)
        for tk in taps:
            m = t >= tk
            y[m] += amp * np.exp(-0.03 * 2 * math.pi * f * (t[m] - tk)) * np.sin(2 * math.pi * f * (t[m] - tk))
        return y
    sig = {"tool_acc_x": 0.1 + ring(8.0, 2.0), "tool_acc_y": -0.2 + ring(13.0, 1.0), "tool_acc_z": 9.8 + ring(8.0, 0.3)}
    for j in range(6):
        sig[f"actual_q_{j}"] = POSE_Q["UP"][j] + np.zeros(n)
        sig[f"actual_current_{j}"] = 0.1 + np.zeros(n)
    sig["actual_q_1"] = sig["actual_q_1"] + ring(8.0, 2e-4)
    sig["actual_q_0"] = sig["actual_q_0"] + ring(13.0, 1e-4)
    sig["actual_current_1"] = sig["actual_current_1"] + ring(8.0, 0.3)
    noise = {"tool_acc": 0.02, "actual_q": 5e-6, "actual_current": 0.005}
    for k in sig:
        sig[k] = sig[k] + rng.normal(0, next(v for p, v in noise.items() if k.startswith(p)), n)
    keys = ["timestamp"] + list(sig)
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f)
        w.writerow(keys)
        for i in range(n):
            w.writerow([round(t[i], 6)] + [("" if (not acc and k.startswith("tool_acc")) else sig[k][i]) for k in sig])


class TestAnalysis(unittest.TestCase):
    def _analyze(self, **kw):
        m = load("tt_ana")
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "syn.csv")
            synth_csv(p, **kw)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                r = m.analyze_tap_test(p)
        return r, out.getvalue()

    def test_finds_taps_and_frequencies(self):
        r, text = self._analyze()
        self.assertEqual(r["mode"], "tap")
        self.assertEqual(len(r["taps"]), 4)
        for got, want in zip(r["taps"], (2.0, 5.0, 8.0, 11.0)):
            self.assertAlmostEqual(got, want, delta=0.05)
        self.assertEqual(r["n_win"], 4)
        for name, f in (("tool_acc_x", 8.0), ("actual_q_1", 8.0), ("actual_current_1", 8.0),
                        ("tool_acc_y", 13.0), ("actual_q_0", 13.0)):
            self.assertAlmostEqual(r["peaks"][name][0][0], f, delta=0.3, msg=name)
        self.assertIn("工具加速度計找到 4 下", text)

    def test_no_taps_is_background(self):
        r, text = self._analyze(taps=())
        self.assertEqual(r["mode"], "background")
        self.assertEqual(r["taps"], [])
        self.assertGreater(r["n_win"], 0)
        self.assertIn("背景", text)

    def test_falls_back_to_current(self):
        r, text = self._analyze(acc=False)
        self.assertEqual(r["mode"], "tap")
        self.assertEqual(len(r["taps"]), 4)
        self.assertEqual(r["trigger"], "關節電流")
        self.assertAlmostEqual(r["peaks"]["actual_q_1"][0][0], 8.0, delta=0.3)
        self.assertNotIn("tool_acc_x", r["peaks"])

    def test_too_short(self):
        m = load("tt_short")
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "s.csv")
            synth_csv(p, taps=(), dur=0.2)
            with contextlib.redirect_stdout(io.StringIO()):
                r = m.analyze_tap_test(p)
        self.assertEqual(r["mode"], "none")
        self.assertEqual(m.tap_analysis_lines(r), ["資料太少，沒有分析。"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
