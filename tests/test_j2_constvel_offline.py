# -*- coding: utf-8 -*-
"""
J2 固定角度窗口定速（robot/ur5_j2_const_velocity_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -v

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在開始姿態不動）。

測試內容
--------
1. 語法正確；姿態表：UP（主要，J2 = 0°）為預設、DIAG（對照，診斷原點），兩者 J1 都在 −90°。
2. 重力（UR 官方連桿質量與重心 + 末端工具 1.30 kg，以官方 DH 計算）：UP 窗口中心 −5.7° 為 J2 重力零點，
   窗口內 −4.34～+4.36 N·m、平均約 0；DIAG 窗口內 +21.42～+23.41 N·m、平均 +22.68。
   對照：UP 的 J2 重力 = 同姿態的 J1 重力（上臂直立，重心在 J1 正上方），J2 = 0° 時 2.34；DIAG 23.32（同 J1）。
3. 規劃：窗口 6 個馬達圈；三個階段、兩個方向的起點與煞車後最遠點，UP 在 −24.41°～+13.01°，
   DIAG 在 52.29°～89.71°（不超過診斷原點）；正反轉以窗口中心對稱。
4. 安全檢查：兩個姿態三個階段都通過；J2 視窗未設定、速度超限、規劃範圍超出視窗、法蘭線速度超限、
   窗口移到前臂往下的一側（法蘭過低）、法蘭高度下限不是數值時 FAIL。
5. 法蘭：離 J2 軸 0.494 m（兩姿態相同，只轉 J2 時不變）；離底座安裝面 UP 0.96～1.01 m、診斷原點 0.432 m。
6. 姿態檢查：依姿態名稱檢查六軸；名稱未知、UP 與 DIAG 互換（J2 不符）、J0、J1 不符、NaN、讀不到時拒絕；±360° 等價。
7. URScript：縮排規則、只含 ASCII、每段一個有次數上限的迴圈、正轉用 "<"、反轉用 ">"、
   movej 只改 J2、speedj 只動 J2、最後回到開始角度（UP 0°、DIAG 90°）。
8. 負載設定與輸入確認。
9. main()：姿態名稱未知、J2 視窗未設定、負載未設定、輸入的姿態名稱不對、負載不符、RTDE 失敗、
   實際姿態不符（例如 POSE = UP 但 J2 還在 90°）、第一次 movej 過遠、速度滑桿不是 100% 時都在送出前中止；
   UP 兩個方向、DIAG 一個方向端對端各送出一次腳本，產生含姿態名稱的 CSV（含 target_moment）與摘要檔；
   UP 跑完提醒以教導器移回診斷原點。
10. 分析（合成資料，重力保持力矩用 2. 的官方模型）：
   - 摩擦與受力無關：已知 B、Tc 與換算比例，加上跟馬達角度鎖定的漣波與雜訊，兩個姿態都找回 B、Tc、
     窗口平均重力保持力矩；UP 以斜率找回換算比例、C₊ ≈ C₋、不以平均值估；DIAG 以平均值找回換算比例。
   - 摩擦隨減速機傳遞的力矩增加（0.1 × |力矩|）：UP 的 H̄ 仍正確；DIAG 的 H̄ 偏約 4.6%。
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
J2C = os.path.join(ROBOT, "ur5_j2_const_velocity_ident.py")

HOME_Q = [3.14159265, -1.57079633, 1.57079633, -1.57079633, -1.57079633, 0.0]   # 診斷原點（DIAG）
UP_Q = HOME_Q[:2] + [0.0] + HOME_Q[3:]                                          # 手臂朝上（UP）：J2 = 0°
EXT45_Q = [HOME_Q[0]] + [math.radians(x) for x in (-45.0, 45.0, -90.0, -90.0, 0.0)]
POSE_Q = {"UP": UP_Q, "DIAG": HOME_Q}
PAYLOAD = 1.3
SAFE = {"UP": (-28.0, 17.0),          # 測試用 J2 視窗（含 FULL 規劃範圍 −26.41°～+15.01°）
        "DIAG": (48.0, 93.0)}         # （含 FULL 規劃範圍 50.29°～92.00°）
TM = {"UP": [0.0, -2.3, 0.5, 1.0, 0.1, 0.0], "DIAG": [0.0, -23.3, -23.3, 1.0, 0.1, 0.0]}   # 假的 target_moment

# UR 官方 DH 頁面（UR5，CB3）的連桿質量 [kg] 與重心（各連桿 DH 座標系）[m]、末端工具（同 analysis/1029/known_mass_1029.py）
LINK_M = [3.7, 8.393, 2.33, 1.219, 1.219, 0.1879]
LINK_COM = [[0, -0.02561, 0.00193], [0.2125, 0, 0.11336], [0.15, 0, 0.0265],
            [0, -0.0018, 0.01634], [0, 0.0018, 0.01634], [0, 0, -0.001159]]
TOOL_KG, TOOL_COG = 1.30, np.array([0.001, 0.001, 0.074])
G = 9.81
K0 = 11.0                             # 合成資料真正的換算比例 [N·m/A]


class FakeRTDE:
    q = list(UP_Q)
    tm = list(TM["UP"])
    sf = 1.0
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getTargetQ(self): return list(FakeRTDE.q)
    def getTargetMoment(self): return list(FakeRTDE.tm)
    def getJointTemperatures(self): return [30.0] * 6
    def getTargetSpeedFraction(self): return FakeRTDE.sf
    def disconnect(self): pass


def load(name, pose=None):
    """載入程式；pose 不為 None 時改 POSE（否則用程式裡的預設），J2 視窗依姿態填 SAFE。"""
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, J2C)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if pose is not None:
        mod.POSE = pose
    mod.PAYLOAD_KG = PAYLOAD
    mod.J2_SAFE_MIN_DEG, mod.J2_SAFE_MAX_DEG = SAFE[mod.POSE]
    return mod


def arm_at(pose):
    """假的手臂停在 pose（六軸角度與 target_moment）。"""
    FakeRTDE.q, FakeRTDE.tm = list(POSE_Q[pose]), list(TM[pose])


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(confirm="UP", payload=str(PAYLOAD)):
    """模擬操作者輸入：問負載質量時回答 payload，要求輸入姿態名稱時回答 confirm，其他提示回答空字串。"""
    def ans(prompt="", *a, **k):
        p = str(prompt)
        if "負載質量" in p:
            return payload
        if "輸入" in p and "繼續" in p:
            return confirm
        return ""
    return ans


def done(logger, planned_s, **k):
    """取代等待函式：立即回報運動結束。"""
    return "done", planned_s, planned_s


def plan(m, stage, direction, pose="UP", center=None):
    levels, order, n_revs = m.stage_config(stage)
    c = m.POSES[pose]["center_deg"] if center is None else center
    ta, tb = m.window_bounds_deg(c, n_revs, m.GEAR_RATIO)
    segs = m.plan_segments(m.speed_sequence(levels, order), direction, ta, tb, m.SPEEDJ_ACCEL,
                           m.STOPJ_DECEL, m.SETTLE_TIME, m.OVERRUN_DEG, m.LOOP_TIME_FACTOR, m.DT)
    return ta, tb, segs


def check(m, segs, pose="UP", safe=None, qd_max=None, **kw):
    q = m.POSES[pose]["q_deg"]
    lo, hi = SAFE[pose] if safe is None else safe
    return m.safety_check(segs, m.QD_MAX if qd_max is None else qd_max, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL,
                          m.STOPJ_DECEL, m.MOVEJ_VEL, m.MOVEJ_ACCEL, lo, hi, m.SAFETY_MARGIN_DEG, q,
                          end_deg=q[2], **kw)


def gravity_joint(m, qj_deg, pose, j=2):
    """官方模型在 J(j) 上的重力力矩 τ_g [N·m]（正方向為正），只把 J(j) 換成 qj_deg；撐住重力的保持力矩 H = −τ_g。"""
    q = [math.radians(x) for x in m.POSES[pose]["q_deg"]]
    q[j] = math.radians(qj_deg)
    fr = m.fk_frames(q)
    o, z = fr[j][:3, 3], fr[j][:3, 2]                         # J(j) 轉軸
    w = np.array([0.0, 0.0, -G])
    tau = 0.0
    for i in range(j, 6):                                     # J(j) 以後的連桿
        c = fr[i + 1][:3, :3] @ np.array(LINK_COM[i]) + fr[i + 1][:3, 3]
        tau += np.cross(c - o, LINK_M[i] * w) @ z
    c = fr[6][:3, :3] @ TOOL_COG + fr[6][:3, 3]
    tau += np.cross(c - o, TOOL_KG * w) @ z
    return float(tau)


def holding_fn(m, pose):
    """官方模型的重力保持力矩 H(q2) = −τ_g；只轉 J2 時為 a·cos q2 + b·sin q2（以 0°、90° 兩點求係數）。"""
    a, b = -gravity_joint(m, 0.0, pose), -gravity_joint(m, 90.0, pose)
    return lambda q_deg: a * math.cos(math.radians(q_deg)) + b * math.sin(math.radians(q_deg))


class Base(unittest.TestCase):
    def setUp(self):
        self._sock = socket.socket
        self._input = builtins.input
        socket.socket = no_network
        builtins.input = answers()
        arm_at("UP")
        FakeRTDE.sf = 1.0

    def tearDown(self):
        socket.socket = self._sock
        builtins.input = self._input
        arm_at("UP")
        FakeRTDE.sf = 1.0


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


def movej_targets(text):
    return [[float(x) for x in l.split("[")[1].split("]")[0].split(",")] for l in text.splitlines() if "movej(" in l]


class TestPlan(Base):
    def test_compile(self):
        py_compile.compile(J2C, doraise=True)

    def test_poses(self):
        m = load("j2c_poses")
        self.assertEqual(m.JOINT_INDEX, 2)
        self.assertEqual(m.POSE, "UP", "預設姿態為 UP")
        self.assertEqual(set(m.POSES), {"UP", "DIAG"})
        self.assertEqual(tuple(m.POSES["UP"]["q_deg"]), (180.0, -90.0, 0.0, -90.0, -90.0, 0.0))
        self.assertEqual(tuple(m.POSES["DIAG"]["q_deg"]), (180.0, -90.0, 90.0, -90.0, -90.0, 0.0))
        self.assertEqual(m.POSES["UP"]["center_deg"], -5.7)
        self.assertEqual(m.POSES["DIAG"]["center_deg"], 71.0)
        for pose, q in POSE_Q.items():
            for a, b in zip(m.POSES[pose]["q_deg"], q):
                self.assertAlmostEqual(math.radians(a), b, places=6)

    def test_gravity_official_model(self):
        """窗口中心的依據：UP 窗口中心為 J2 重力零點、窗口內正負對稱；DIAG 窗口內重力大、都在同一邊。"""
        m = load("j2c_grav")
        rng = np.random.default_rng(1)
        for pose in ("UP", "DIAG"):
            h = holding_fn(m, pose)
            for q in rng.uniform(-180.0, 180.0, 5):
                self.assertAlmostEqual(h(q), -gravity_joint(m, q, pose), places=9, msg="只轉 J2 時為正弦")
        # 上臂直立（J1 = −90°）時，上臂的重心在 J1 正上方：J1 與 J2 的重力相同
        for pose in ("UP", "DIAG"):
            q2 = m.POSES[pose]["q_deg"][2]
            self.assertAlmostEqual(gravity_joint(m, q2, pose), gravity_joint(m, -90.0, pose, j=1), places=9, msg=pose)
        self.assertAlmostEqual(gravity_joint(m, 0.0, "UP"), 2.34, delta=0.01)      # 同 J1 定速 UP 的 −90°
        self.assertAlmostEqual(gravity_joint(m, 90.0, "DIAG"), 23.32, delta=0.01)  # 診斷原點（同 J1）
        self.assertLess(abs(gravity_joint(m, m.POSES["UP"]["center_deg"], "UP")), 0.02, "UP 窗口中心為重力零點")
        for pose, lo, hi, mean in (("UP", -4.34, 4.36, 0.01), ("DIAG", 21.42, 23.41, 22.68)):
            ta, tb = m.window_bounds_deg(m.POSES[pose]["center_deg"], m.N_MOTOR_REVS_WINDOW, m.GEAR_RATIO)
            g = [gravity_joint(m, x, pose) for x in np.linspace(ta, tb, 2001)]
            self.assertAlmostEqual(min(g), lo, delta=0.01, msg=pose)
            self.assertAlmostEqual(max(g), hi, delta=0.01, msg=pose)
            self.assertAlmostEqual(float(np.mean(g)), mean, delta=0.01, msg=pose)

    def test_window_and_range(self):
        m = load("j2c_plan")
        for pose, lo_pts, hi_pts in (("UP", -24.41, 13.01), ("DIAG", 52.29, 89.71)):
            c = m.POSES[pose]["center_deg"]
            ta, tb, sp = plan(m, "FULL", +1, pose)
            _, _, sn = plan(m, "FULL", -1, pose)
            self.assertAlmostEqual(tb - ta, 6 * 360 / 101, places=9)
            self.assertAlmostEqual((ta + tb) / 2, c, places=9)
            self.assertEqual(len(sp), 12)
            self.assertEqual([s["speed"] for s in sp], m.SPEED_LEVELS + m.SPEED_LEVELS[::-1])
            for a, b in zip(sp, sn):
                self.assertAlmostEqual(a["start_deg"] - c, -(b["start_deg"] - c), places=9)
                self.assertAlmostEqual(a["extreme_deg"] - c, -(b["extreme_deg"] - c), places=9)
                self.assertLess(a["start_deg"], ta)
                self.assertGreater(a["stop_cmd_deg"], tb)
            self.assertAlmostEqual(min(s["start_deg"] for s in sp), lo_pts, delta=0.01)
            self.assertAlmostEqual(max(s["start_deg"] for s in sn), hi_pts, delta=0.01)
            for stage in ("QUICK", "MID", "FULL"):
                for d in (+1, -1):
                    _, _, segs = plan(m, stage, d, pose)
                    pts = [s["start_deg"] for s in segs] + [s["extreme_deg"] for s in segs]
                    self.assertGreaterEqual(min(pts), lo_pts - 0.01, f"{pose} {stage} {d:+d}")
                    self.assertLessEqual(max(pts), hi_pts + 0.01, f"{pose} {stage} {d:+d}")
                    first = abs(segs[0]["start_deg"] - m.POSES[pose]["q_deg"][2])
                    self.assertLess(first, m.MAX_INITIAL_MOVE_DEG, f"{pose} {stage} {d:+d}：第一次 movej")

    def test_safety(self):
        m = load("j2c_safety")
        for pose in ("UP", "DIAG"):
            for stage in ("QUICK", "MID", "FULL"):
                for d in (+1, -1):
                    _, _, segs = plan(m, stage, d, pose)
                    ok, lines, _ = check(m, segs, pose)
                    self.assertTrue(ok, f"{pose} {stage} {d:+d}：" + "\n".join(lines))
        for pose, lo, hi, narrow in (("UP", -26.41, 15.01, (-24.0, 17.0)),
                                     ("DIAG", 50.29, 92.0, (55.0, 93.0))):
            _, _, sp = plan(m, "FULL", +1, pose)
            _, _, sn = plan(m, "FULL", -1, pose)
            self.assertAlmostEqual(check(m, sp, pose)[2][0], lo, delta=0.01, msg=pose)
            self.assertAlmostEqual(check(m, sn, pose)[2][1], hi, delta=0.01, msg=pose)
            self.assertFalse(check(m, sp, pose, safe=(None, None))[0], "J2 視窗未設定時應 FAIL")
            self.assertFalse(check(m, sp, pose, safe=(None, hi + 1.0))[0])
            self.assertFalse(check(m, sp, pose, qd_max=0.15)[0], "速度超過上限時應 FAIL")
            self.assertFalse(check(m, sp, pose, safe=narrow)[0], "規劃範圍超出 J2 視窗時應 FAIL")
            self.assertFalse(check(m, sp, pose, flange_z_min=None)[0], "法蘭高度下限不是數值時應 FAIL")
            self.assertFalse(check(m, sp, pose, flange_z_min=float("nan"))[0])
            # 法蘭線速度 = 法蘭離 J2 軸距離 0.494 m × 0.20 rad/s = 0.0988 m/s
            ok, lines, _ = check(m, sp, pose, flange_speed_max=0.09)
            self.assertFalse(ok, pose)
            k = next(i for i, l in enumerate(lines) if l.startswith("法蘭線速度"))
            self.assertIn("0.0988 m/s", lines[k])
            self.assertEqual(lines[k + 1].strip(), "[FAIL]")
            self.assertEqual(sum("[FAIL]" in l for l in lines), 1, "只有法蘭線速度 FAIL")
            self.assertTrue(check(m, sp, pose, flange_speed_max=0.10)[0])
        # 窗口設到前臂往下的一側：J2 視窗再寬，法蘭高度檢查也要擋下
        for pose, center, safe in (("DIAG", 100.0, (48.0, 125.0)), ("UP", 110.0, (-28.0, 135.0))):
            _, _, low = plan(m, "FULL", -1, pose, center=center)
            ok, lines, _ = check(m, low, pose, safe=safe)
            self.assertFalse(ok, pose)
            k = next(i for i, l in enumerate(lines) if l.startswith("法蘭離底座"))
            self.assertEqual(lines[k + 1].strip(), "[FAIL]", pose)
            self.assertEqual(sum("[FAIL]" in l for l in lines), 1, f"{pose}：其他項目應通過，只有法蘭高度 FAIL")

    def test_flange_geometry(self):
        m = load("j2c_fk")
        up, diag = m.POSES["UP"]["q_deg"], m.POSES["DIAG"]["q_deg"]
        self.assertAlmostEqual(m.flange_radius(up), 0.494, delta=0.001)
        self.assertAlmostEqual(m.flange_radius(diag), 0.494, delta=0.001)
        q = list(up); q[2] = 30.0
        self.assertAlmostEqual(m.flange_radius(q), m.flange_radius(up), places=9, msg="只轉 J2 時距離不變")
        self.assertAlmostEqual(m.flange_height(0.0, up), 1.001, delta=0.001)
        zs = [m.flange_height(a, up) for a in np.linspace(-26.41, 15.01, 200)]
        self.assertGreater(min(zs), 0.96)
        self.assertLess(max(zs), 1.01)
        self.assertAlmostEqual(m.flange_height(90.0, diag), 0.432, delta=0.001)
        self.assertAlmostEqual(m.flange_height(52.29, diag), 0.747, delta=0.001)
        self.assertGreater(m.flange_height(92.0, diag), m.FLANGE_Z_MIN_M)
        self.assertLess(m.flange_height(95.0, diag), m.FLANGE_Z_MIN_M)

    def test_pose(self):
        m = load("j2c_pose")
        self.assertTrue(m.check_pose(UP_Q, "UP")[0])
        self.assertTrue(m.check_pose(HOME_Q, "DIAG")[0])
        self.assertFalse(m.check_pose(HOME_Q, "UP")[0], "J2 還在 90° 時不可當成 UP")
        self.assertFalse(m.check_pose(UP_Q, "DIAG")[0])
        ok, lines = m.check_pose(UP_Q, "FOO")
        self.assertFalse(ok)
        self.assertIn("未知的姿態", lines[0])
        for pose, base in POSE_Q.items():
            q = list(base); q[0] = -math.pi
            self.assertTrue(m.check_pose(q, pose)[0], "J0 = −180° 與 180° 等價")
            q = list(base); q[0] = math.radians(170.0)
            self.assertFalse(m.check_pose(q, pose)[0], "J0 不在 180° 應拒絕")
            q = list(base); q[2] += math.radians(0.4)
            self.assertTrue(m.check_pose(q, pose)[0])
            q = list(base); q[2] += math.radians(0.6)
            self.assertFalse(m.check_pose(q, pose)[0], "J2 不在開始角度應拒絕")
            q = list(base); q[1] += math.radians(0.6)
            self.assertFalse(m.check_pose(q, pose)[0], "J1 不在 −90° 應拒絕")
            q = list(base); q[2] = float("nan")
            self.assertFalse(m.check_pose(q, pose)[0])
            q = list(base); q[3] = None
            self.assertFalse(m.check_pose(q, pose)[0])
            self.assertFalse(m.check_pose(EXT45_Q, pose)[0])
            self.assertFalse(m.check_pose(None, pose)[0])
            self.assertFalse(m.check_pose(base[:5], pose)[0])

    def test_urscript(self):
        m = load("j2c_script")
        for pose, lo, hi, end in (("UP", -24.41, 13.01, 0.0), ("DIAG", 52.29, 90.0, 90.0)):
            base = POSE_Q[pose]
            for d, op in ((+1, "<"), (-1, ">")):
                _, _, segs = plan(m, "FULL", d, pose)
                s = m.build_urscript(segs, d, 2, base, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                     m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=end)
                self.assertTrue(s.startswith("def j2_const_velocity_program():"))
                self.assertTrue(urscript_indent_ok(s))
                self.assertTrue(all(ord(ch) < 128 for ch in s), "URScript 應只含 ASCII")
                self.assertEqual(s.count("while (n < "), len(segs), "每段都要有次數上限")
                self.assertEqual(s.count(f"q_now[2] {op} "), len(segs))
                self.assertNotIn("q_now[1]", s)
                self.assertEqual(s.count("stopj("), len(segs))
                tg = movej_targets(s)
                self.assertEqual(len(tg), len(segs) + 1, "最後要多一個回原點的 movej")
                for v in tg:
                    for j in (0, 1, 3, 4, 5):
                        self.assertAlmostEqual(v[j], base[j], places=6, msg="movej 只改 J2")
                    self.assertGreaterEqual(math.degrees(v[2]), lo - 0.01, pose)
                    self.assertLessEqual(math.degrees(v[2]), hi + 0.01, pose)
                for v, seg in zip(tg, segs):
                    self.assertAlmostEqual(math.degrees(v[2]), seg["start_deg"], places=4)
                self.assertAlmostEqual(math.degrees(tg[-1][2]), end, places=4)
                self.assertGreater(s.rfind("movej("), s.rfind("stopj("), "回原點應在最後一段 stopj 之後")
                for line in s.splitlines():
                    if "speedj(" in line:
                        qd = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                        self.assertEqual([j for j in range(6) if qd[j] != 0.0], [2], "speedj 只動 J2")
                        self.assertEqual(math.copysign(1.0, qd[2]), float(d))

    def test_payload(self):
        m = load("j2c_payload")
        for bad in (None, -1.0, float("nan"), "abc"):
            self.assertFalse(m.check_payload_setting(bad)[0], repr(bad))
        self.assertTrue(m.check_payload_setting(1.3)[0])
        self.assertTrue(m.check_payload_entry(1.3, " 1.30 ")[0])
        self.assertTrue(m.check_payload_entry(1.3, "1.304")[0])
        for bad in ("1.31", "2.30", "abc", ""):
            self.assertFalse(m.check_payload_entry(1.3, bad)[0], repr(bad))
        self.assertFalse(m.check_payload_entry(None, "1.3")[0])
        self.assertEqual(m.payload_tag(1.3), "p130")


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

    def test_unknown_pose_aborts(self):
        m = load("j2c_m0")
        m.POSE = "FOO"
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("未知的姿態", self.last_out)

    def test_window_unset_aborts(self):
        m = load("j2c_m1")
        m.J2_SAFE_MIN_DEG = m.J2_SAFE_MAX_DEG = None
        self.assertEqual(self._run(m)[0], [], "視窗未設定時不該送出任何軌跡")
        self.assertIn("安全檢查未通過", self.last_out)

    def test_payload_unset_aborts(self):
        m = load("j2c_m2")
        m.PAYLOAD_KG = None
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("安全檢查未通過", self.last_out)

    def test_confirm_cancel(self):
        for k, typed in enumerate(("DIAG", "up", "J2", "")):
            m = load(f"j2c_m3{k}")
            builtins.input = answers(confirm=typed)
            self.assertEqual(self._run(m)[0], [], repr(typed))
            self.assertIn("已取消", self.last_out)

    def test_payload_mismatch_aborts(self):
        m = load("j2c_m4")
        builtins.input = answers(payload="2.30")
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("負載確認未通過", self.last_out)

    def test_rtde_failure_aborts(self):
        m = load("j2c_m5")

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")
        sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("RTDE 記錄未能啟動", self.last_out)

    def test_wrong_pose_aborts(self):
        q_j2 = list(UP_Q); q_j2[2] = math.radians(10.0)
        q_j1 = list(UP_Q); q_j1[1] = math.radians(-100.0)
        for k, q in enumerate((HOME_Q, [math.radians(170.0)] + UP_Q[1:], q_j1, q_j2, EXT45_Q)):
            m = load(f"j2c_m6{k}")
            FakeRTDE.q = list(q)
            self.assertEqual(self._run(m)[0], [], f"案例 {k}")
            self.assertIn("姿態檢查未通過", self.last_out)
        m = load("j2c_m6d", "DIAG")
        builtins.input = answers(confirm="DIAG")
        self.assertEqual(self._run(m)[0], [], "POSE = DIAG 但手臂在 UP")
        self.assertIn("姿態檢查未通過", self.last_out)

    def test_far_start_aborts(self):
        m = load("j2c_m7")
        m.MAX_INITIAL_MOVE_DEG = 3.0                        # UP QUICK 正轉第一次 movej 為 9.06°
        self.assertEqual(self._run(m)[0], [])
        self.assertIn("起始位置檢查未通過", self.last_out)

    def test_speed_fraction_aborts(self):
        for k, sf in enumerate((0.52, None, float("nan"))):
            m = load(f"j2c_m8{k}")
            FakeRTDE.sf = sf
            self.assertEqual(self._run(m)[0], [], f"案例 {k}")
            self.assertIn("速度滑桿檢查未通過", self.last_out)

    def test_quick_end_to_end_up(self):
        m = load("j2c_m9")
        sent, files = self._run(m, directions=(+1, -1), wait=done)
        self.assertEqual(len(sent), 2)
        for txt in sent:
            self.assertTrue(txt.startswith("def j2_const_velocity_program():"))
            tg = movej_targets(txt)
            self.assertEqual(len(tg), 2, "QUICK：1 段 + 回原點")
            for v in tg:
                self.assertAlmostEqual(math.degrees(v[1]), -90.0, places=4, msg="J1 維持 −90°")
            self.assertAlmostEqual(math.degrees(tg[-1][2]), 0.0, places=4)
        self.assertIn("q_now[2] < ", sent[0])
        self.assertIn("q_now[2] > ", sent[1])
        self.assertAlmostEqual(math.degrees(movej_targets(sent[0])[0][2]), -9.06, delta=0.01)
        self.assertAlmostEqual(math.degrees(movej_targets(sent[1])[0][2]), -2.34, delta=0.01)
        for tag in ("pos", "neg"):
            data = [n for n in files if n.startswith(f"j2constvel_data_up_{tag}_p130_quick_") and n.endswith(".csv")]
            self.assertEqual(len(data), 1, tag)
            rows = list(csv.DictReader(io.StringIO(files[data[0]])))
            self.assertGreater(len(rows), 0)
            self.assertAlmostEqual(float(rows[0]["target_moment_2"]), TM["UP"][2])
            self.assertTrue(all(float(r["target_speed_fraction"]) == 1.0 for r in rows))
            info = next(v for n, v in files.items() if n.startswith(f"j2constvel_info_up_{tag}_p130_quick_"))
            self.assertIn("姿態: UP", info)
            self.assertIn("執行狀態: done", info)
            self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
            self.assertIn("J1：目前   -90.00°", info)
            self.assertIn("J2：目前     0.00°", info)
        self.assertIn("下一步", self.last_out)
        self.assertIn("[提醒]", self.last_out)
        self.assertNotIn("這次的資料不可用", self.last_out)

    def test_quick_end_to_end_diag(self):
        m = load("j2c_m10", "DIAG")
        arm_at("DIAG")
        builtins.input = answers(confirm="DIAG")
        sent, files = self._run(m, wait=done)
        self.assertEqual(len(sent), 1)
        tg = movej_targets(sent[0])
        self.assertAlmostEqual(math.degrees(tg[0][2]), 67.64, delta=0.01)
        self.assertAlmostEqual(math.degrees(tg[-1][2]), 90.0, places=4)
        for v in tg:
            self.assertAlmostEqual(math.degrees(v[1]), -90.0, places=4)
        data = [n for n in files if n.startswith("j2constvel_data_diag_pos_p130_quick_") and n.endswith(".csv")]
        self.assertEqual(len(data), 1)
        rows = list(csv.DictReader(io.StringIO(files[data[0]])))
        self.assertAlmostEqual(float(rows[0]["target_moment_2"]), TM["DIAG"][2])
        info = next(v for n, v in files.items() if n.startswith("j2constvel_info_diag_pos_p130_quick_"))
        self.assertIn("姿態: DIAG", info)
        self.assertIn("執行狀態: done", info)
        self.assertNotIn("[提醒]", self.last_out)


def synth(m, pose, stage, d, tau_true, noise=0.0, ripple=0.0, seed=0):
    """
    合成一個方向的記錄並分析：J2 依規劃以定速走完每一段，tau_true(d, v, h) 為真正的力矩 [N·m]
    （h：官方模型的重力保持力矩），另加跟馬達角度鎖定的漣波（ripple 倍）與雜訊（標準差 noise）；
    電流 = 力矩 ÷ K0，target_moment = h。回傳 (THETA_A, THETA_B, analyze_angle_window 的結果)。
    """
    rng = np.random.default_rng(seed)
    hold = holding_fn(m, pose)
    ta, tb, segs = plan(m, stage, d, pose)
    rows, t = [], 0.0
    for s in segs:
        v, q_deg, end = s["speed"], s["start_deg"], s["stop_cmd_deg"]
        while (q_deg - end) * d < 0:
            q_rad = math.radians(q_deg)
            phi = (q_rad * m.GEAR_RATIO) % (2 * math.pi)
            rip = ripple * (1.2 * math.cos(5 * phi + 0.7) + 0.8 * math.sin(12 * phi) + 0.4 * math.cos(2 * phi))
            h = hold(q_deg)
            tau = tau_true(d, v, h) + rip + (rng.normal(0, noise) if noise else 0.0)
            rows.append({"timestamp": t, "actual_q_2": q_rad,
                         "actual_qd_2": d * v + (rng.normal(0, 0.002) if noise else 0.0),
                         "actual_current_2": tau / K0, "target_qd_2": d * v, "target_moment_2": h})
            q_deg += d * math.degrees(v * m.DT)
            t += m.DT
        t += 5.0
    with tempfile.TemporaryDirectory() as tmp:
        p = os.path.join(tmp, "syn.csv")
        with open(p, "w", newline="", encoding="utf-8") as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)
        with contextlib.redirect_stdout(io.StringIO()):
            res = m.analyze_angle_window(p, 2, m.KT_OUT, d, m.SPEED_LEVELS, ta, tb)
    return ta, tb, res


def window_mean(f, ta, tb):
    return float(np.mean([f(x) for x in np.linspace(ta, tb, 20001)]))


class TestAnalysis(unittest.TestCase):
    B0, TC0 = 25.0, 8.0                                     # 真值 [N·m·s/rad]、[N·m]

    def test_linear_friction_both_poses(self):
        m = load("j2c_ana")
        B0, TC0 = self.B0, self.TC0
        scale = m.KT_OUT / K0                               # 電流 × KT_OUT 相對於真值的比例
        for pose, h_expect in (("UP", -0.01), ("DIAG", -22.68)):
            res = {}
            for d in (+1, -1):
                ta, tb, res[d] = synth(m, pose, "FULL", d, lambda d_, v, h: B0 * d_ * v + TC0 * d_ + h,
                                       noise=0.4, ripple=1.0, seed=int(d > 0))
                self.assertEqual(len(res[d]["segments"]), 12)
                for sg in res[d]["segments"]:
                    self.assertAlmostEqual(sg["revs"], 6.0, delta=0.05)
                self.assertAlmostEqual(res[d]["B_up"], res[d]["B_down"], delta=0.5)
            with contextlib.redirect_stdout(io.StringIO()):
                c = m.combine_directions(res[+1], res[-1])
            h_bar = window_mean(holding_fn(m, pose), ta, tb)
            self.assertAlmostEqual(h_bar, h_expect, delta=0.01, msg=f"{pose} 窗口平均重力保持力矩")
            self.assertAlmostEqual(c["B"], B0 * scale, delta=0.01 * B0 * scale, msg=pose)
            self.assertAlmostEqual(c["Tc"], TC0 * scale, delta=0.05, msg=pose)
            self.assertAlmostEqual(c["H"], h_bar * scale, delta=0.05, msg=pose)
            self.assertAlmostEqual(c["tm_win"], h_bar, delta=0.02, msg=pose)
            # 單一方向時 C 混入重力：C₊ = Tc + H̄、C₋ = Tc − H̄
            self.assertAlmostEqual(res[+1]["C"], (TC0 + h_bar) * scale, delta=0.08, msg=pose)
            self.assertAlmostEqual(res[-1]["C"], (TC0 - h_bar) * scale, delta=0.08, msg=pose)
            if pose == "DIAG":
                self.assertAlmostEqual(c["K_mean"], K0, delta=0.01 * K0)
            else:
                self.assertAlmostEqual(c["K_slope"], K0, delta=0.03 * K0, msg="UP 以重力斜率估換算比例")
                self.assertTrue(math.isnan(c["K_mean"]), "UP 窗口平均重力約為 0，不以平均值估換算比例")
                self.assertAlmostEqual(res[+1]["C"], res[-1]["C"], delta=0.1, msg="UP：正、反轉的 C 相同")

    def test_load_dependent_friction(self):
        """
        摩擦隨減速機傳遞的力矩增加：τ = h + d·(F0 + β·|τ|)，F0 = B0·v + TC0
        → τ = (d·F0 + h)/(1 − d·β·sign(d·F0 + h))。UP 窗口內重力正負對稱，正反轉的受力相同，H̄ 仍正確；
        DIAG 往上抬與往下放的受力不同，H̄ 偏 3% 以上（約 4.6%）。
        """
        m = load("j2c_ana3")
        B0, TC0, beta = self.B0, self.TC0, 0.1
        scale = m.KT_OUT / K0

        def tau_true(d, v, h):
            x = d * (B0 * v + TC0) + h
            return x / (1.0 - d * beta * math.copysign(1.0, x))

        for pose in ("UP", "DIAG"):
            res = {}
            for d in (+1, -1):
                ta, tb, res[d] = synth(m, pose, "MID", d, tau_true)
            with contextlib.redirect_stdout(io.StringIO()):
                c = m.combine_directions(res[+1], res[-1])
            h_true = window_mean(holding_fn(m, pose), ta, tb) * scale
            if pose == "UP":
                self.assertLess(abs(c["H"] - h_true), 0.02, "UP：摩擦隨受力改變時 H̄ 仍正確")
                self.assertAlmostEqual(res[+1]["C"], res[-1]["C"], delta=0.04)
            else:
                self.assertGreater(abs(c["H"] / h_true - 1.0), 0.03, "DIAG：H̄ 偏 3% 以上（約 4.6%）")

    def test_combine_needs_regression(self):
        m = load("j2c_ana2")
        with self.assertRaises(ValueError):
            m.combine_directions({"segments": []}, {"segments": []})


if __name__ == "__main__":
    unittest.main(verbosity=2)
