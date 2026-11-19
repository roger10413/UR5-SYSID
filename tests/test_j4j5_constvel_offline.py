# -*- coding: utf-8 -*-
"""
J4、J5 固定角度窗口定速（robot/ur5_j4_const_velocity_ident.py、ur5_j5_const_velocity_ident.py）的離線測試。

執行：  python -m unittest discover -s tests -p "test_*.py"

安全設計同其他離線測試：socket 換成一定丟例外的函式（不可能連線、不可能送出 URScript）、
假的 rtde_receive（手臂停在開始姿態不動）。兩支程式由 J3 定速改寫，測試內容相同，以 AXES 設定區分。

測試內容（J4、J5 各一份）
--------
1. 語法正確；姿態表：J4 預設 UP（J2 = 0°、J3 = −90°，同 J1 速度調變的 UP）、J5 預設 DIAG（診斷原點）；窗口中心 J4 −90°、J5 0°。
2. 重力（UR 官方連桿質量與重心 + 末端工具 1.30 kg，以官方 DH 計算）：J4 UP 轉軸垂直、重力恆為 0；J4 DIAG 像單擺，
   重力零點 −89.7°、最大約 2.34 N·m，窗口（中心 −90°）內 −0.42～+0.45 N·m、平均約 0；J5 兩個姿態都 ≤ 0.02 N·m。
3. 規劃：窗口 6 個馬達圈；三個階段、兩個方向的起點與煞車後最遠點，J4 在 −108.71°～−71.29°，J5 在 ±18.71°；
   正反轉以窗口中心對稱；第一次 movej 不超過上限。
4. 安全檢查：兩個姿態三個階段都通過；視窗未設定、速度超限、規劃範圍超出視窗、法蘭高度下限不是數值或太高、
   J4 法蘭線速度超限時 FAIL。
5. 法蘭：J4 離轉軸 0.082 m、J5 為 0（只轉該軸時不變）；UP 時法蘭高度不變（1.001 m），J4 DIAG 工具朝下時最低（0.432 m）。
6. 姿態檢查：依姿態名稱檢查六軸；名稱未知、兩個姿態互換、J3 的 UP（J3 = −130°）、任一軸超差、NaN、讀不到時拒絕；±360° 等價。
7. URScript：縮排規則、只含 ASCII、每段一個有次數上限的迴圈、正轉用 "<"、反轉用 ">"、movej 與 speedj 只動該軸、最後回到開始角度。
8. 負載設定與輸入確認。
9. main()：各種錯誤都在送出前中止；兩個姿態、兩個方向端對端各送出一次腳本，產生含姿態名稱的 CSV 與摘要檔；UP 跑完提醒移回診斷原點。
10. 分析（合成資料，重力保持力矩用 2. 的官方模型）：兩個姿態都找回 B、Tc、窗口平均重力保持力矩；
    J4 DIAG 以重力斜率估換算比例（參考）；J4 UP、J5 的重力太小，K_mean、K_slope 都不算（NaN）。
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

UP_DEG = (180.0, -90.0, 0.0, -90.0, -90.0, 0.0)       # 手臂朝上、工具朝前水平（同 J1 速度調變的 UP）
DIAG_DEG = (180.0, -90.0, 90.0, -90.0, -90.0, 0.0)    # 診斷原點
J3UP_DEG = (180.0, -90.0, 0.0, -130.0, -90.0, 0.0)    # J3 程式的 UP（J3 = −130°），不是 J4、J5 的姿態
EXT45_DEG = (180.0, -45.0, 45.0, -90.0, -90.0, 0.0)
RAD = lambda q: [math.radians(x) for x in q]
POSE_Q = {"UP": RAD(UP_DEG), "DIAG": RAD(DIAG_DEG)}
PAYLOAD = 1.3

AXES = {
    4: dict(path=os.path.join(ROBOT, "ur5_j4_const_velocity_ident.py"), default="UP", center=-90.0,
            safe=(-112.0, -68.0), pts=(-108.71, -71.29), margin=(-110.71, -69.29), radius=0.0823,
            quick_first={+1: -93.36, -1: -86.64}),
    5: dict(path=os.path.join(ROBOT, "ur5_j5_const_velocity_ident.py"), default="DIAG", center=0.0,
            safe=(-22.0, 22.0), pts=(-18.71, 18.71), margin=(-20.71, 20.71), radius=0.0,
            quick_first={+1: -3.36, -1: 3.36}),
}
TM = [0.0, -2.3, 0.5, -0.07, 0.013, 0.001]             # 假的 target_moment

# UR 官方 DH 頁面（UR5，CB3）的連桿質量 [kg] 與重心（各連桿 DH 座標系）[m]、末端工具（同 test_j3_constvel_offline.py）
LINK_M = [3.7, 8.393, 2.33, 1.219, 1.219, 0.1879]
LINK_COM = [[0, -0.02561, 0.00193], [0.2125, 0, 0.11336], [0.15, 0, 0.0265],
            [0, -0.0018, 0.01634], [0, 0.0018, 0.01634], [0, 0, -0.001159]]
TOOL_KG, TOOL_COG = 1.30, np.array([0.001, 0.001, 0.074])
G = 9.81
K0 = 9.0                              # 合成資料真正的換算比例 [N·m/A]


class FakeRTDE:
    q = list(POSE_Q["UP"])
    sf = 1.0
    def __init__(self, ip): pass
    def getTimestamp(self): return time.time()
    def getActualQ(self): return list(FakeRTDE.q)
    def getActualQd(self): return [0.0] * 6
    def getActualCurrent(self): return [0.0] * 6
    def getTargetQd(self): return [0.0] * 6
    def getTargetQ(self): return list(FakeRTDE.q)
    def getTargetMoment(self): return list(TM)
    def getJointTemperatures(self): return [30.0] * 6
    def getTargetSpeedFraction(self): return FakeRTDE.sf
    def disconnect(self): pass


def load(n, name, pose=None):
    """載入 J(n) 程式；pose 不為 None 時改 POSE，視窗填 AXES[n]["safe"]、負載 1.30。"""
    fake = types.ModuleType("rtde_receive")
    fake.RTDEReceiveInterface = FakeRTDE
    sys.modules["rtde_receive"] = fake
    if ROBOT not in sys.path:
        sys.path.insert(0, ROBOT)
    spec = importlib.util.spec_from_file_location(name, AXES[n]["path"])
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if pose is not None:
        mod.POSE = pose
    mod.PAYLOAD_KG = PAYLOAD
    setattr(mod, f"J{n}_SAFE_MIN_DEG", AXES[n]["safe"][0])
    setattr(mod, f"J{n}_SAFE_MAX_DEG", AXES[n]["safe"][1])
    return mod


def no_network(*a, **k):
    raise RuntimeError("測試中禁止任何網路連線")


def answers(confirm, payload=str(PAYLOAD)):
    def ans(prompt="", *a, **k):
        p = str(prompt)
        if "負載質量" in p:
            return payload
        if "輸入" in p and "繼續" in p:
            return confirm
        return ""
    return ans


def done(logger, planned_s, **k):
    return "done", planned_s, planned_s


def plan(m, stage, direction, pose, center=None):
    levels, order, n_revs = m.stage_config(stage)
    c = m.POSES[pose]["center_deg"] if center is None else center
    ta, tb = m.window_bounds_deg(c, n_revs, m.GEAR_RATIO)
    segs = m.plan_segments(m.speed_sequence(levels, order), direction, ta, tb, m.SPEEDJ_ACCEL,
                           m.STOPJ_DECEL, m.SETTLE_TIME, m.OVERRUN_DEG, m.LOOP_TIME_FACTOR, m.DT)
    return ta, tb, segs


def gravity_joint(m, qj_deg, pose, j):
    """官方模型在 J(j) 上的重力力矩 τ_g [N·m]，只把 J(j) 換成 qj_deg；撐住重力的保持力矩 H = −τ_g。"""
    q = [math.radians(x) for x in m.POSES[pose]["q_deg"]]
    q[j] = math.radians(qj_deg)
    fr = m.fk_frames(q)
    o, z = fr[j][:3, 3], fr[j][:3, 2]
    w = np.array([0.0, 0.0, -G])
    tau = 0.0
    for i in range(j, 6):
        c = fr[i + 1][:3, :3] @ np.array(LINK_COM[i]) + fr[i + 1][:3, 3]
        tau += np.cross(c - o, LINK_M[i] * w) @ z
    c = fr[6][:3, :3] @ TOOL_COG + fr[6][:3, 3]
    tau += np.cross(c - o, TOOL_KG * w) @ z
    return float(tau)


def holding_fn(m, pose, j):
    """只轉 J(j) 時 H(q) = −τ_g = a·cos q + b·sin q（以 0°、90° 兩點求係數）。"""
    a, b = -gravity_joint(m, 0.0, pose, j), -gravity_joint(m, 90.0, pose, j)
    return lambda q_deg: a * math.cos(math.radians(q_deg)) + b * math.sin(math.radians(q_deg))


def urscript_indent_ok(text):
    lines = text.splitlines()
    return lines[0].startswith("def ") and all(l.startswith(" ") for l in lines[1:-1]) and lines[-1] == "end"


def movej_targets(text):
    return [[float(x) for x in l.split("[")[1].split("]")[0].split(",")] for l in text.splitlines() if "movej(" in l]


class Base(unittest.TestCase):
    N = None

    @property
    def cfg(self):
        return AXES[self.N]

    def setUp(self):
        self._sock, self._input = socket.socket, builtins.input
        socket.socket = no_network
        builtins.input = answers(self.cfg["default"])
        FakeRTDE.q, FakeRTDE.sf = list(POSE_Q[self.cfg["default"]]), 1.0

    def tearDown(self):
        socket.socket, builtins.input = self._sock, self._input
        FakeRTDE.q, FakeRTDE.sf = list(POSE_Q["UP"]), 1.0
        sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE

    def m(self, tag, pose=None):
        return load(self.N, f"j{self.N}c_{tag}", pose)

    def check(self, m, segs, pose, safe=None, qd_max=None, **kw):
        q = m.POSES[pose]["q_deg"]
        lo, hi = self.cfg["safe"] if safe is None else safe
        return m.safety_check(segs, m.QD_MAX if qd_max is None else qd_max, m.QDD_CHECK_LIMIT, m.SPEEDJ_ACCEL,
                              m.STOPJ_DECEL, m.MOVEJ_VEL, m.MOVEJ_ACCEL, lo, hi, m.SAFETY_MARGIN_DEG, q,
                              end_deg=q[self.N], **kw)


class PlanCases:
    def test_compile(self):
        py_compile.compile(self.cfg["path"], doraise=True)

    def test_poses(self):
        m = self.m("poses")
        n = self.N
        self.assertEqual(m.JOINT_INDEX, n)
        self.assertEqual(m.POSE, self.cfg["default"])
        self.assertEqual(set(m.POSES), {"UP", "DIAG"})
        self.assertEqual(tuple(m.POSES["UP"]["q_deg"]), UP_DEG)
        self.assertEqual(tuple(m.POSES["DIAG"]["q_deg"]), DIAG_DEG)
        for pose in ("UP", "DIAG"):
            self.assertEqual(m.POSES[pose]["center_deg"], self.cfg["center"])
        self.assertEqual(m.SPEED_LEVELS, [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16, 0.20])
        self.assertTrue(hasattr(m, f"J{n}_SAFE_MIN_DEG") and not hasattr(m, "J3_SAFE_MIN_DEG"))

    def test_gravity_official_model(self):
        m = self.m("grav")
        n = self.N
        rng = np.random.default_rng(1)
        for pose in ("UP", "DIAG"):
            h = holding_fn(m, pose, n)
            for q in rng.uniform(-180.0, 180.0, 5):
                self.assertAlmostEqual(h(q), -gravity_joint(m, q, pose, n), places=9, msg="只轉一軸時為正弦")
        g = {p: [gravity_joint(m, x, p, n) for x in np.linspace(-180.0, 180.0, 721)] for p in ("UP", "DIAG")}
        ta, tb = m.window_bounds_deg(self.cfg["center"], m.N_MOTOR_REVS_WINDOW, m.GEAR_RATIO)
        if n == 4:
            self.assertLess(max(map(abs, g["UP"])), 1e-9, "J4 UP：轉軸垂直，重力恆為 0")
            self.assertAlmostEqual(max(map(abs, g["DIAG"])), 2.34, delta=0.01, msg="J4 DIAG 重力最大約 2.34 N·m")
            self.assertLess(abs(gravity_joint(m, -89.7, "DIAG", 4)), 0.005, "J4 DIAG 重力零點 −89.7°")
            self.assertLess(abs(gravity_joint(m, -90.0, "DIAG", 4)), 0.015, "窗口中心 −90° 的重力可忽略")
            gw = [gravity_joint(m, x, "DIAG", 4) for x in np.linspace(ta, tb, 2001)]
            self.assertAlmostEqual(min(gw), -0.42, delta=0.01)       # 中心 −90° 比重力零點偏 0.3°
            self.assertAlmostEqual(max(gw), 0.45, delta=0.01)
            self.assertLess(abs(float(np.mean(gw))), 0.02)
        else:
            for pose in ("UP", "DIAG"):
                self.assertLess(max(map(abs, g[pose])), 0.02, f"J5 {pose}：重力可忽略")
            self.assertLess(max(map(abs, g["DIAG"])), 1e-9, "J5 DIAG：轉軸垂直，重力恆為 0")

    def test_window_and_range(self):
        m = self.m("plan")
        n, (lo_pts, hi_pts) = self.N, self.cfg["pts"]
        for pose in ("UP", "DIAG"):
            c = m.POSES[pose]["center_deg"]
            ta, tb, sp = plan(m, "FULL", +1, pose)
            _, _, sn = plan(m, "FULL", -1, pose)
            self.assertAlmostEqual(tb - ta, 6 * 360 / 101, places=9)
            self.assertAlmostEqual((ta + tb) / 2, c, places=9)
            self.assertEqual(len(sp), 2 * len(m.SPEED_LEVELS))
            for a, b in zip(sp, sn):
                self.assertAlmostEqual(a["start_deg"] - c, -(b["start_deg"] - c), places=9)
                self.assertAlmostEqual(a["extreme_deg"] - c, -(b["extreme_deg"] - c), places=9)
            self.assertAlmostEqual(min(s["start_deg"] for s in sp), lo_pts, delta=0.01)
            self.assertAlmostEqual(max(s["start_deg"] for s in sn), hi_pts, delta=0.01)
            for stage in ("QUICK", "MID", "FULL"):
                for d in (+1, -1):
                    _, _, segs = plan(m, stage, d, pose)
                    pts = [s["start_deg"] for s in segs] + [s["extreme_deg"] for s in segs]
                    self.assertGreaterEqual(min(pts), lo_pts - 0.01, f"{pose} {stage} {d:+d}")
                    self.assertLessEqual(max(pts), hi_pts + 0.01, f"{pose} {stage} {d:+d}")
                    first = abs(segs[0]["start_deg"] - m.POSES[pose]["q_deg"][n])
                    self.assertLess(first, 13.0, f"{pose} {stage} {d:+d}：第一次 movej")

    def test_safety(self):
        m = self.m("safety")
        lo, hi = self.cfg["margin"]
        for pose in ("UP", "DIAG"):
            for stage in ("QUICK", "MID", "FULL"):
                for d in (+1, -1):
                    _, _, segs = plan(m, stage, d, pose)
                    ok, lines, _ = self.check(m, segs, pose)
                    self.assertTrue(ok, f"{pose} {stage} {d:+d}：" + "\n".join(lines))
            _, _, sp = plan(m, "FULL", +1, pose)
            _, _, sn = plan(m, "FULL", -1, pose)
            self.assertAlmostEqual(self.check(m, sp, pose)[2][0], lo, delta=0.01)
            self.assertAlmostEqual(self.check(m, sn, pose)[2][1], hi, delta=0.01)
            self.assertFalse(self.check(m, sp, pose, safe=(None, None))[0], "視窗未設定時應 FAIL")
            self.assertFalse(self.check(m, sp, pose, safe=(None, hi + 1.0))[0])
            self.assertFalse(self.check(m, sp, pose, qd_max=0.15)[0], "速度超過上限時應 FAIL")
            self.assertFalse(self.check(m, sp, pose, safe=(lo + 3.0, hi))[0], "規劃範圍超出視窗時應 FAIL")
            self.assertFalse(self.check(m, sp, pose, flange_z_min=None)[0])
            self.assertFalse(self.check(m, sp, pose, flange_z_min=float("nan"))[0])
            ok, lines, _ = self.check(m, sp, pose, flange_z_min=1.2)          # 法蘭最高 1.001 m
            self.assertFalse(ok)
            k = next(i for i, l in enumerate(lines) if l.startswith("法蘭離底座"))
            self.assertEqual(lines[k + 1].strip(), "[FAIL]")
            self.assertEqual(sum("[FAIL]" in l for l in lines), 1, "只有法蘭高度 FAIL")
            ok, lines, _ = self.check(m, sp, pose)
            k = next(i for i, l in enumerate(lines) if l.startswith("法蘭線速度"))
            v_fl = self.cfg["radius"] * 0.20
            self.assertIn(f"{v_fl:.4f} m/s", lines[k])
            if self.N == 4:                                                  # 0.0823 m × 0.20 rad/s = 0.0165 m/s
                ok, lines, _ = self.check(m, sp, pose, flange_speed_max=0.01)
                self.assertFalse(ok)
                self.assertEqual(sum("[FAIL]" in l for l in lines), 1, "只有法蘭線速度 FAIL")

    def test_flange_geometry(self):
        m = self.m("fk")
        n = self.N
        up, diag = m.POSES["UP"]["q_deg"], m.POSES["DIAG"]["q_deg"]
        for pose_q in (up, diag):
            self.assertAlmostEqual(m.flange_radius(pose_q), self.cfg["radius"], delta=0.0005)
            q = list(pose_q); q[n] += 15.0
            self.assertAlmostEqual(m.flange_radius(q), m.flange_radius(pose_q), places=9, msg="只轉一軸時距離不變")
        lo, hi = self.cfg["margin"]
        zs = [m.flange_height(a, up) for a in np.linspace(lo, hi, 100)]
        self.assertAlmostEqual(min(zs), 1.001, delta=0.001)
        self.assertAlmostEqual(max(zs), 1.001, delta=0.001)
        zd = [m.flange_height(a, diag) for a in np.linspace(lo, hi, 100)]
        self.assertAlmostEqual(m.flange_height(diag[n], diag), 0.432, delta=0.001)
        self.assertGreaterEqual(min(zd), 0.432 - 0.001)
        if n == 4:
            self.assertGreater(m.flange_height(-110.0, diag), m.flange_height(-90.0, diag), "往側邊擺法蘭變高")
            self.assertLess(max(zd), 0.44)
        else:
            self.assertAlmostEqual(max(zd), 0.432, delta=0.001)

    def test_pose(self):
        m = self.m("pose")
        n = self.N
        self.assertTrue(m.check_pose(POSE_Q["UP"], "UP")[0])
        self.assertTrue(m.check_pose(POSE_Q["DIAG"], "DIAG")[0])
        self.assertFalse(m.check_pose(POSE_Q["DIAG"], "UP")[0])
        self.assertFalse(m.check_pose(POSE_Q["UP"], "DIAG")[0])
        self.assertFalse(m.check_pose(RAD(J3UP_DEG), "UP")[0], "J3 程式的 UP（J3 = −130°）不是這裡的 UP")
        ok, lines = m.check_pose(POSE_Q["UP"], "FOO")
        self.assertFalse(ok)
        self.assertIn("未知的姿態", lines[0])
        for pose, base in POSE_Q.items():
            q = list(base); q[0] = -math.pi
            self.assertTrue(m.check_pose(q, pose)[0], "J0 = −180° 與 180° 等價")
            for j in range(6):
                q = list(base); q[j] += math.radians(0.4)
                self.assertTrue(m.check_pose(q, pose)[0])
                q = list(base); q[j] += math.radians(0.6)
                self.assertFalse(m.check_pose(q, pose)[0], f"J{j} 超差 0.6° 應拒絕")
            q = list(base); q[n] = float("nan")
            self.assertFalse(m.check_pose(q, pose)[0])
            q = list(base); q[2] = None
            self.assertFalse(m.check_pose(q, pose)[0])
            self.assertFalse(m.check_pose(RAD(EXT45_DEG), pose)[0])
            self.assertFalse(m.check_pose(None, pose)[0])
            self.assertFalse(m.check_pose(base[:5], pose)[0])

    def test_urscript(self):
        m = self.m("script")
        n = self.N
        lo, hi = self.cfg["pts"]
        for pose in ("UP", "DIAG"):
            base, end = POSE_Q[pose], m.POSES[pose]["q_deg"][n]
            for d, op in ((+1, "<"), (-1, ">")):
                _, _, segs = plan(m, "FULL", d, pose)
                s = m.build_urscript(segs, d, n, base, m.SPEEDJ_ACCEL, m.STOPJ_DECEL,
                                     m.MOVEJ_ACCEL, m.MOVEJ_VEL, m.PAUSE_TIME, m.DT, end_deg=end)
                self.assertTrue(s.startswith(f"def j{n}_const_velocity_program():"))
                self.assertTrue(urscript_indent_ok(s))
                self.assertTrue(all(ord(ch) < 128 for ch in s), "URScript 應只含 ASCII")
                self.assertEqual(s.count("while (n < "), len(segs), "每段都要有次數上限")
                self.assertEqual(s.count(f"q_now[{n}] {op} "), len(segs))
                for j in range(6):
                    if j != n:
                        self.assertNotIn(f"q_now[{j}]", s)
                self.assertEqual(s.count("stopj("), len(segs))
                tg = movej_targets(s)
                self.assertEqual(len(tg), len(segs) + 1, "最後要多一個回原點的 movej")
                for v in tg:
                    for j in range(6):
                        if j != n:
                            self.assertAlmostEqual(v[j], base[j], places=6, msg=f"movej 只改 J{n}")
                    self.assertGreaterEqual(math.degrees(v[n]), lo - 0.01)
                    self.assertLessEqual(math.degrees(v[n]), hi + 0.01)
                for v, seg in zip(tg, segs):
                    self.assertAlmostEqual(math.degrees(v[n]), seg["start_deg"], places=4)
                self.assertAlmostEqual(math.degrees(tg[-1][n]), end, places=4)
                for line in s.splitlines():
                    if "speedj(" in line:
                        qd = [float(x) for x in line.split("[")[1].split("]")[0].split(",")]
                        self.assertEqual([j for j in range(6) if qd[j] != 0.0], [n], f"speedj 只動 J{n}")
                        self.assertEqual(math.copysign(1.0, qd[n]), float(d))

    def test_payload(self):
        m = self.m("payload")
        for bad in (None, -1.0, float("nan"), "abc"):
            self.assertFalse(m.check_payload_setting(bad)[0], repr(bad))
        self.assertTrue(m.check_payload_entry(1.3, " 1.30 ")[0])
        for bad in ("1.31", "2.30", "abc", ""):
            self.assertFalse(m.check_payload_entry(1.3, bad)[0], repr(bad))
        self.assertEqual(m.payload_tag(1.3), "p130")


class MainCases:
    def _run(self, m, stage="QUICK", directions=(+1,), wait=None):
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
                for nm in names:
                    with open(os.path.join(root, nm), encoding="utf-8") as f:
                        files[nm] = f.read()
        self.last_out = buf.getvalue()
        return sent, files

    def test_aborts_before_sending(self):
        n, default = self.N, self.cfg["default"]
        other = "DIAG" if default == "UP" else "UP"

        class Failing:
            def __init__(self, ip): raise ConnectionError("模擬 RTDE 連線失敗")

        cases = [("pose_name", "未知的姿態"), ("window", "安全檢查未通過"), ("payload_unset", "安全檢查未通過"),
                 ("confirm", "已取消"), ("payload", "負載確認未通過"), ("rtde", "RTDE 記錄未能啟動"),
                 ("wrong_pose", "姿態檢查未通過"), ("j3up", "姿態檢查未通過"), ("far", "起始位置檢查未通過"),
                 ("slider", "速度滑桿檢查未通過"), ("slider_none", "速度滑桿檢查未通過")]
        for k, (case, reason) in enumerate(cases):
            m = self.m(f"abort{k}")
            if case == "pose_name":
                m.POSE = "FOO"
            elif case == "window":
                setattr(m, f"J{n}_SAFE_MIN_DEG", None)
            elif case == "payload_unset":
                m.PAYLOAD_KG = None
            elif case == "confirm":
                builtins.input = answers(other)
            elif case == "payload":
                builtins.input = answers(default, payload="2.30")
            elif case == "rtde":
                sys.modules["rtde_receive"].RTDEReceiveInterface = Failing
            elif case == "wrong_pose":
                FakeRTDE.q = list(POSE_Q[other])
            elif case == "j3up":
                FakeRTDE.q = RAD(J3UP_DEG)
            elif case == "far":
                m.MAX_INITIAL_MOVE_DEG = 3.0                # QUICK 第一次 movej 為 3.36°
            elif case == "slider":
                FakeRTDE.sf = 0.52
            elif case == "slider_none":
                FakeRTDE.sf = None
            sent, _ = self._run(m)
            self.assertEqual(sent, [], case)
            self.assertIn(reason, self.last_out, case)
            FakeRTDE.q, FakeRTDE.sf = list(POSE_Q[default]), 1.0
            builtins.input = answers(default)
            sys.modules["rtde_receive"].RTDEReceiveInterface = FakeRTDE

    def test_quick_end_to_end(self):
        n = self.N
        for pose in ("UP", "DIAG"):
            m = self.m(f"e2e_{pose}", pose)
            FakeRTDE.q = list(POSE_Q[pose])
            builtins.input = answers(pose)
            sent, files = self._run(m, directions=(+1, -1), wait=done)
            self.assertEqual(len(sent), 2, pose)
            for d, txt in zip((+1, -1), sent):
                self.assertTrue(txt.startswith(f"def j{n}_const_velocity_program():"))
                tg = movej_targets(txt)
                self.assertEqual(len(tg), 2, "QUICK：1 段 + 回原點")
                self.assertAlmostEqual(math.degrees(tg[0][n]), self.cfg["quick_first"][d], delta=0.01)
                self.assertAlmostEqual(math.degrees(tg[-1][n]), m.POSES[pose]["q_deg"][n], places=4)
                self.assertIn(f"q_now[{n}] {'<' if d > 0 else '>'} ", txt)
            for tag in ("pos", "neg"):
                pre = f"j{n}constvel_data_{pose.lower()}_{tag}_p130_quick_"
                data = [nm for nm in files if nm.startswith(pre) and nm.endswith(".csv")]
                self.assertEqual(len(data), 1, (pose, tag, list(files)))
                rows = list(csv.DictReader(io.StringIO(files[data[0]])))
                self.assertGreater(len(rows), 0)
                self.assertAlmostEqual(float(rows[0][f"target_moment_{n}"]), TM[n])
                self.assertIn(f"joint_temp_{n}", rows[0])
                info = next(v for nm, v in files.items() if nm.startswith(f"j{n}constvel_info_{pose.lower()}_{tag}_"))
                self.assertIn(f"姿態: {pose}", info)
                self.assertIn("執行狀態: done", info)
                self.assertIn("教導器速度滑桿（記錄中最低）: 100.0%", info)
                self.assertIn(f"J{n}：目前 {m.POSES[pose]['q_deg'][n]:8.2f}°", info)
            self.assertIn("下一步", self.last_out)
            if pose == "UP":
                self.assertIn("J2 轉回 90°", self.last_out)
            else:
                self.assertNotIn("[提醒]", self.last_out)
            self.assertNotIn("這次的資料不可用", self.last_out)


class AnalysisCases:
    B0, TC0 = 6.0, 2.5                                      # 真值 [N·m·s/rad]、[N·m]（手腕較小，推測）

    def synth(self, m, pose, stage, d, noise=0.0, ripple=0.0, seed=0):
        n = self.N
        rng = np.random.default_rng(seed)
        hold = holding_fn(m, pose, n)
        ta, tb, segs = plan(m, stage, d, pose)
        rows, t = [], 0.0
        for s in segs:
            v, q_deg, end = s["speed"], s["start_deg"], s["stop_cmd_deg"]
            while (q_deg - end) * d < 0:
                q_rad = math.radians(q_deg)
                phi = (q_rad * m.GEAR_RATIO) % (2 * math.pi)
                rip = ripple * (0.6 * math.cos(5 * phi + 0.7) + 0.4 * math.sin(12 * phi) + 0.2 * math.cos(2 * phi))
                h = hold(q_deg)
                tau = self.B0 * d * v + self.TC0 * d + h + rip + (rng.normal(0, noise) if noise else 0.0)
                rows.append({"timestamp": t, f"actual_q_{n}": q_rad,
                             f"actual_qd_{n}": d * v + (rng.normal(0, 0.002) if noise else 0.0),
                             f"actual_current_{n}": tau / K0, f"target_qd_{n}": d * v, f"target_moment_{n}": h})
                q_deg += d * math.degrees(v * m.DT)
                t += m.DT
            t += 5.0
        with tempfile.TemporaryDirectory() as tmp:
            p = os.path.join(tmp, "syn.csv")
            with open(p, "w", newline="", encoding="utf-8") as f:
                w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
                w.writeheader(); w.writerows(rows)
            with contextlib.redirect_stdout(io.StringIO()):
                res = m.analyze_angle_window(p, n, m.KT_OUT, d, m.SPEED_LEVELS, ta, tb)
        return ta, tb, res

    def test_linear_friction_both_poses(self):
        m = self.m("ana")
        n = self.N
        scale = m.KT_OUT / K0
        for pose in ("UP", "DIAG"):
            res = {}
            for d in (+1, -1):
                ta, tb, res[d] = self.synth(m, pose, "FULL", d, noise=0.2, ripple=1.0, seed=int(d > 0))
                self.assertEqual(len(res[d]["segments"]), 2 * len(m.SPEED_LEVELS))
                for sg in res[d]["segments"]:
                    self.assertAlmostEqual(sg["revs"], 6.0, delta=0.05)
            with contextlib.redirect_stdout(io.StringIO()):
                c = m.combine_directions(res[+1], res[-1])
            hold = holding_fn(m, pose, n)
            h_bar = float(np.mean([hold(x) for x in np.linspace(ta, tb, 20001)]))
            self.assertLess(abs(h_bar), 0.02, f"{pose}：窗口平均重力約 0")
            msg = f"J{n} {pose}"
            self.assertAlmostEqual(c["B"], self.B0 * scale, delta=0.02 * self.B0 * scale, msg=msg)
            self.assertAlmostEqual(c["Tc"], self.TC0 * scale, delta=0.05, msg=msg)
            self.assertAlmostEqual(c["H"], h_bar * scale, delta=0.05, msg=msg)
            self.assertTrue(math.isnan(c["K_mean"]), f"{msg}：窗口平均重力太小，不以平均值估換算比例")
            if (n, pose) == (4, "DIAG"):
                self.assertGreater(abs(c["slope_tm"]), 0.02)
                self.assertAlmostEqual(c["K_slope"], K0, delta=0.2 * K0, msg="J4 DIAG 以重力斜率估換算比例（參考）")
            else:
                self.assertLess(abs(c["slope_tm"]), 0.02)
                self.assertTrue(math.isnan(c["K_slope"]), f"{msg}：重力太小，不以斜率估換算比例")

    def test_combine_needs_regression(self):
        m = self.m("ana2")
        with self.assertRaises(ValueError):
            m.combine_directions({"segments": []}, {"segments": []})


class TestJ4Plan(PlanCases, Base):
    N = 4


class TestJ5Plan(PlanCases, Base):
    N = 5


class TestJ4Main(MainCases, Base):
    N = 4


class TestJ5Main(MainCases, Base):
    N = 5


class TestJ4Analysis(AnalysisCases, Base):
    N = 4


class TestJ5Analysis(AnalysisCases, Base):
    N = 5


if __name__ == "__main__":
    unittest.main(verbosity=2)
