"""
UR5 J1 靜止電流測試 —— 以已知質量直接校正 J1 的力矩換電流比例（Kt × 減速比）
================================================================

★★★ 為什麼要這個實驗 ★★★

J0、B、Tc 的力矩都是「電流 × Kt × 101」換算而來，Kt = 0.135 借自文獻，未在這台手臂上驗證：
  - 控制器內部對 J0 用 12.44 N·m/A、J1 用 12.05 N·m/A（10/22 資料，target_current 對 target_moment 回歸），
    比本研究的 13.64 低約 10%。
  - 10/22 速度調變的兩姿態差 ΔJ 為控制器模型的 1.25～1.40 倍。
J0 是垂直軸，重力不作用在 J0 上；J1（肩關節）是水平軸，夾爪夾砝碼時 J1 的重力力矩增加
  Δτ = m·g·d（d = 砝碼重心離 J1 軸的水平距離，由 UR5 官方 DH 與實際關節角計算），
與連桿參數、控制器模型都無關。同姿態「夾砝碼 − 不夾」的電流差 Δi 即得 K = Δτ ÷ Δi。
J1 與 J0 同尺寸（UR5 前三軸額定 150 N·m，推測為同型馬達與減速機），K 可作為 J0 換算比例的參考。

做法：
  A. J0~J5 不動、停在姿態 POSE（DIAG 或 EXT45，與速度調變相同），J1 每次循環：
       movej 到 q1 − STEP_DEG → movej 回 q1（最後以正向靠近）→ 靜止 HOLD_S 秒；
       movej 到 q1 + STEP_DEG → movej 回 q1（最後以負向靠近）→ 靜止 HOLD_S 秒。
     其他五軸維持執行開始時 RTDE 量到的角度；結束時停在 q1（與開始時相同）。
  B. 靜止時減速機的靜摩擦可吸收任意大小（在 ±Tc 內）的力矩，只量一次靜止電流會混入未知的摩擦；
     從兩個方向靠近，摩擦方向相反，平均即消去（假設兩方向摩擦大小相同）。兩者差的一半為靜摩擦帶寬。
  C. 分析只用每次靜止開始 SKIP_S 秒之後到 HOLD_S 為止的資料；記錄 target_moment（控制器的重力模型，
     含教導器負載設定），可檢查教導器負載是否已更新：夾砝碼前後 target_moment 的差應為 m·g·d。
  D. 時間：QUICK 1 次循環約 10 s，FULL N_CYCLE 次約 1.5 分鐘。
  E. J1 轉 STEP_DEG 時，末端上下移動約 d × STEP_DEG（DIAG 約 ±4 mm，EXT45 約 ±7 mm）。

離線驗證（合成資料：中心電流、兩方向靜摩擦 ±0.6 A、白雜訊 0.05 A）
-------------------------------------------------------------------------------------------
  兩方向判定、中心電流與靜摩擦帶寬都正確找回；夾砝碼前後的 K 找回設定值（見 tests/test_j1_static_offline.py）。

★★★ 使用前必讀 ★★★
1. 【審查狀態】新程式，尚未審查；依 TEST_STAGE = QUICK → FULL 階梯執行。
2. 【只動 J1】每次 movej 只改 J1（±STEP_DEG，上限 STEP_MAX_DEG），速度 MOVE_VEL、加速度 MOVE_ACCEL，
   都有上限檢查；J0、J2~J5 的目標角度等於開始時量到的角度。
3. 【姿態】程式不會擺姿態：請先以教導器擺好（J0 = 180°），連上 RTDE 後以實際角度檢查六軸（±POSE_TOL_DEG），
   不符或姿態名稱未知即中止。須輸入姿態名稱，確認工具（含砝碼）下方淨空。
4. 【負載 PAYLOAD_KG】同速度調變：填教導器上設定的負載質量（只有末端工具時 1.30 kg；夾砝碼時為 1.30 + 砝碼
   質量），未設定一律中止；執行前須輸入教導器上的負載質量，不符即中止。
5. 【教導器速度滑桿】讀 target_speed_fraction，不是 100% 或讀不到一律中止。
6. 【IP 需要你自己再次確認】
"""

import csv
import math
import os
import socket
import statistics
import threading
import time
from datetime import datetime

import numpy as np

try:
    import rtde_receive
except ImportError:
    rtde_receive = None


# ============================================================
# 使用者設定
# ============================================================

ROBOT_IP = "192.168.50.114"   # 執行前請再次確認
JOINT_INDEX = 1               # J1（肩關節）

# QUICK: 1 次循環（正向靠近、負向靠近各 1 次），確認動作與記錄
# FULL : N_CYCLE 次循環
TEST_STAGE = "QUICK"
N_CYCLE = 8

STEP_DEG = 0.5                # 每次離開 q1 的角度 [deg]；馬達側 50°，足以讓摩擦完全反向
MOVE_VEL = 0.05               # movej 速度 [rad/s]
MOVE_ACCEL = 0.3              # movej 加速度 [rad/s^2]
HOLD_S = 4.0                  # 回到 q1 後的靜止時間 [s]
SKIP_S = 1.0                  # 分析時略過每次靜止開始的時間 [s]
PAUSE_TIME = 0.3              # 到達 q1 ± STEP_DEG 後的停頓 [s]

STEP_MAX_DEG = 1.0            # STEP_DEG 上限
MOVE_VEL_MAX = 0.10           # MOVE_VEL 上限 [rad/s]
MOVE_ACCEL_MAX = 0.5          # MOVE_ACCEL 上限 [rad/s^2]

SAMPLE_HZ = 125.0
GEAR_RATIO = 101
KT_OUT = GEAR_RATIO * 0.1350  # 本研究目前的換算比例 [N·m/A]（只用於列印對照）
G = 9.81

# 姿態：六軸期望角度 [deg]，J0 = 180°（速度調變結束時的位置）。
# d_j1_m：法蘭中心離 J1 軸的水平距離 [m]（UR5 官方 DH 計算；工具垂直朝下，砝碼在工具軸上任何位置都相同）。
POSE = "DIAG"
POSES = {
    "DIAG": dict(q_deg=(180.0, -90.0, 90.0, -90.0, -90.0, 0.0), d_j1_m=0.487),    # 診斷原點
    "EXT45": dict(q_deg=(180.0, -45.0, 45.0, -90.0, -90.0, 0.0), d_j1_m=0.787),   # 手臂外伸，前臂水平
}
POSE_TOL_DEG = 0.5

# 法蘭負載：填教導器上設定的負載質量 [kg]（只有末端工具 FT300、夾爪、3D 相機時 1.30；夾砝碼時為 1.30 + 砝碼質量）。
# 只用於確認與記錄（寫進輸出資料夾名稱），不影響動作；未設定則程式中止 —— 不提供預設值。
PAYLOAD_KG = None
PAYLOAD_TOL_KG = 0.005

# 教導器速度滑桿（RTDE target_speed_fraction，0～1）的下限。
SPEED_FRACTION_MIN = 0.999

# 額外記錄的 RTDE 輸出（只記錄，不影響動作）：(欄名, ur_rtde 方法)。控制器或 ur_rtde 不支援時該欄留空。
EXTRA_SIGNALS = (("target_q", "getTargetQ"), ("target_moment", "getTargetMoment"),
                 ("target_current", "getTargetCurrent"), ("joint_control_output", "getJointControlOutput"))

# UR5（CB3）官方 DH 參數
DH_D = (0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823)
DH_A = (0.0, -0.425, -0.39225, 0.0, 0.0, 0.0)
DH_ALPHA = (math.pi / 2, 0.0, 0.0, math.pi / 2, -math.pi / 2, 0.0)

OUTPUT_DIR = "."


# ============================================================
# 第一部分：規劃與安全檢查
# ============================================================

def stage_cycles(stage):
    """回傳循環次數。"""
    if stage == "QUICK":
        return 1
    if stage == "FULL":
        return N_CYCLE
    raise ValueError(f"未知的 TEST_STAGE: {stage}")


def plan_moves(q1_deg, step_deg, n_cycle):
    """
    J1 的 movej 目標序列 [(角度 deg, 是否靜止點, 靠近方向)]：每次循環
    q1 − step → q1（正向靠近，靜止）→ q1 + step → q1（負向靠近，靜止）。
    """
    moves = []
    for _ in range(n_cycle):
        moves += [(q1_deg - step_deg, False, 0), (q1_deg, True, +1),
                  (q1_deg + step_deg, False, 0), (q1_deg, True, -1)]
    return moves


def movej_time(dist_rad, v, a):
    """單軸 movej 的時間（梯形或三角形速度）。"""
    d = abs(dist_rad)
    if d < 1e-12:
        return 0.0
    if d >= v * v / a:
        return d / v + v / a
    return 2.0 * math.sqrt(d / a)


def estimate_total_time(n_cycle, step_deg, v, a, hold_s, pause_s):
    t_move = movej_time(math.radians(step_deg), v, a)
    return n_cycle * (4 * t_move + 2 * hold_s + 2 * pause_s)


def safety_check(step_deg, v, a, hold_s, skip_s, pose, poses=POSES):
    """
    預設拒絕：任一參數不是有限正數、超過上限、姿態未知，一律 FAIL。回傳 (ok, lines)。
    """
    lines, ok = [], True

    def chk(name, val, lim, unit):
        nonlocal ok
        x = _finite_float(val)
        passed = x is not None and 0.0 < x <= lim
        ok = ok and passed
        lines.append(f"{name}: {val} {unit}  (限 0 < x ≤ {lim})  {'[OK]' if passed else '[FAIL]'}")

    chk("J1 單次移動角度 STEP_DEG", step_deg, STEP_MAX_DEG, "deg")
    chk("movej 速度 MOVE_VEL", v, MOVE_VEL_MAX, "rad/s")
    chk("movej 加速度 MOVE_ACCEL", a, MOVE_ACCEL_MAX, "rad/s^2")
    h, s = _finite_float(hold_s), _finite_float(skip_s)
    if h is None or s is None or s < 0 or h < s + 1.0:
        lines.append(f"靜止時間 HOLD_S = {hold_s} s 須比 SKIP_S = {skip_s} s 長至少 1 s  [FAIL]")
        ok = False
    else:
        lines.append(f"靜止時間 {h} s，分析用 {h - s:.1f} s  [OK]")
    if pose not in poses:
        lines.append(f"未知的姿態 POSE = {pose!r}（可用：{', '.join(poses)}）  [FAIL]")
        ok = False
    elif ok:
        d = poses[pose]["d_j1_m"]
        lines.append(f"姿態 {pose}：末端離 J1 軸水平約 {d:.3f} m，J1 轉 ±{step_deg}° 時末端上下移動約 "
                     f"±{1000 * d * math.radians(step_deg):.1f} mm，速度約 {d * v:.3f} m/s")
    return ok, lines


def check_pose(q_full_rad, pose, poses=POSES, tol_deg=POSE_TOL_DEG):
    """
    以 RTDE 量到的六軸角度 [rad] 檢查 J0~J5 是否在 poses[pose] ± tol_deg 內（差值以 ±360° 等價處理）。
    預設拒絕：pose 不在 poses、讀不到六軸角度、角度不是有限值、任一軸超差，一律回傳 False。回傳 (ok, lines)。
    """
    if pose not in poses:
        return False, [f"未知的姿態 POSE = {pose!r}（可用：{', '.join(poses)}），拒絕執行"]
    if q_full_rad is None or len(q_full_rad) != 6:
        return False, ["讀不到六軸角度，拒絕執行"]
    lines, ok = [f"姿態 {pose}：J0~J5 期望 {list(poses[pose]['q_deg'])}°，容差 ±{tol_deg}°"], True
    for j, tgt in enumerate(poses[pose]["q_deg"]):
        cur = math.degrees(q_full_rad[j])
        d = (cur - tgt + 180.0) % 360.0 - 180.0
        passed = abs(d) <= tol_deg                                 # NaN 時為 False
        ok = ok and passed
        lines.append(f"  J{j}：目前 {cur:8.2f}°，期望 {tgt:8.2f}°，差 {d:+7.2f}°  {'[OK]' if passed else '[FAIL]'}")
    return ok, lines


def _finite_float(x):
    """轉成有限浮點數；None、無法轉換、NaN、無窮大時回傳 None。"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def check_payload_setting(payload_kg):
    """PAYLOAD_KG 必須是 ≥ 0 的有限數值；未設定（None）或其他值一律拒絕。回傳 (ok, 說明)。"""
    p = _finite_float(payload_kg)
    if p is None or p < 0:
        return False, f"PAYLOAD_KG = {payload_kg!r}：尚未設定或不是 ≥ 0 的數值，拒絕執行"
    return True, f"負載設定 PAYLOAD_KG = {p:.3f} kg"


def check_payload_entry(payload_kg, entered, tol_kg=PAYLOAD_TOL_KG):
    """操作者輸入的教導器負載質量須與 PAYLOAD_KG 相差 ≤ tol_kg；預設拒絕。回傳 (ok, 說明)。"""
    ok, msg = check_payload_setting(payload_kg)
    if not ok:
        return False, msg
    e = _finite_float(entered.strip() if isinstance(entered, str) else entered)
    if e is None:
        return False, f"輸入的負載質量 {entered!r} 讀不懂，拒絕執行"
    p = float(payload_kg)
    if abs(e - p) > tol_kg:
        return False, f"教導器負載 {e:.3f} kg 與 PAYLOAD_KG = {p:.3f} kg 不符，拒絕執行（請先確認哪一個才對）"
    return True, f"教導器負載 {e:.3f} kg = PAYLOAD_KG -> OK"


def payload_tag(payload_kg):
    """輸出資料夾名稱中的負載標記：1.30 kg → "p130"（單位 0.01 kg）。"""
    return f"p{int(round(float(payload_kg) * 100)):03d}"


def check_speed_fraction(fraction, min_fraction=SPEED_FRACTION_MIN):
    """教導器速度滑桿：讀不到、不是有限數值、低於 min_fraction，一律拒絕。回傳 (ok, 說明)。"""
    x = _finite_float(fraction)
    if x is None:
        return False, f"讀不到教導器速度滑桿（target_speed_fraction = {fraction!r}），拒絕執行"
    msg = f"教導器速度滑桿 {100 * x:.1f}%（須 ≥ {100 * min_fraction:.1f}%）"
    if x < min_fraction:
        return False, msg + " -> 請把滑桿調回 100% 後重新執行"
    return True, msg + " -> OK"


# ============================================================
# 第二部分：URScript 產生與送出
# ============================================================

def build_urscript(moves, q_start, joint_index, v, a, hold_s, pause_s):
    """
    q_start：六軸絕對角度 [rad]（執行開始時 RTDE 量到的值）。每個 movej 只替換 joint_index 的角度，
    其他五軸維持原值。靜止點之後 sleep(hold_s)，其他點之後 sleep(pause_s)。最後一個目標為 q_start。
    """
    lines = ["def j1_static_current_program():"]
    for k, (deg, is_hold, approach) in enumerate(moves):
        q_t = list(q_start)
        q_t[joint_index] = math.radians(deg)
        if is_hold:
            lines.append(f"  # hold {k // 2 + 1}: approach {'+' if approach > 0 else '-'}")   # URScript 內只放 ASCII
        lines.append(f"  movej([{', '.join(f'{x:.6f}' for x in q_t)}], a={a}, v={v})")
        lines.append(f"  sleep({hold_s if is_hold else pause_s})")
    lines.append("end")
    return "\n".join(lines) + "\n"


def send_urscript(ip, script_text, port=30002):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(5.0)
    s.connect((ip, port))
    s.sendall(script_text.encode("utf-8"))
    s.close()


def send_abort(ip):
    """送一段只含 stopj 的程式到 30002，取代控制器上執行中的程式。"""
    try:
        send_urscript(ip, "def abort_prog():\n  stopj(2.0)\nend\n")
        print("[中止] 已送出 stopj。")
    except Exception as e:
        print(f"[中止] stopj 送出失敗：{e} —— 請立即按下緊急停止。")


def wait_for_motion_complete(logger, planned_s, quiet_s, qd_eps=2e-3, no_motion_s=10.0, timeout_factor=3.0):
    """
    J1 實際速度連續 quiet_s 秒低於門檻才視為結束。程式中每次靜止 HOLD_S 秒，quiet_s 須大於 HOLD_S
    （main 用 HOLD_S + 2 s），才不會在中途的靜止誤判結束。
    """
    t_start = time.time()
    hard_timeout = planned_s * timeout_factor + quiet_s + 5.0
    moved, quiet_since, t_first, t_last = False, None, None, None
    while True:
        now = time.time()
        elapsed = now - t_start
        qd = logger.last_qd0
        if qd is not None and abs(qd) > qd_eps:
            if not moved:
                t_first = now
            moved, t_last, quiet_since = True, now, None
        elif moved:
            if quiet_since is None:
                quiet_since = now
            elif now - quiet_since >= quiet_s:
                return "done", elapsed, t_last - t_first
        if not moved and elapsed > no_motion_s:
            return "no_motion", elapsed, 0.0
        if elapsed > hard_timeout:
            return "timeout", elapsed, (t_last - t_first) if moved else 0.0
        time.sleep(0.02)


# ============================================================
# 第三部分：RTDE 記錄（同 ur5_velocity_modulation_ident.py；last_q0 / last_qd0 為 joint_index 的值）
# ============================================================

class RtdeLogger:
    def __init__(self, robot_ip, sample_hz, csv_path, joint_index):
        self.robot_ip = robot_ip
        self.period = 1.0 / sample_hz
        self.csv_path = csv_path
        self.joint_index = joint_index
        self._stop_flag = threading.Event()
        self._thread = None
        self.rows = self.loops = self.gaps = 0
        self.ctrl_times = []
        self.last_q0 = self.last_qd0 = self.last_q_full = None
        self.q0_min = self.q0_max = None
        self.extra_reads = {name: 0 for name, _ in EXTRA_SIGNALS}
        self.extra_missing = []
        self.last_speed_fraction = self.min_speed_fraction = None
        self.speed_fraction_missing = 0
        self.started_ok = threading.Event()
        self.error = None

    def start(self):
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        self._stop_flag.set()
        if self._thread is not None:
            self._thread.join(timeout=10.0)

    def _run(self):
        try:
            rtde_r = rtde_receive.RTDEReceiveInterface(self.robot_ip)
            csv_file = open(self.csv_path, "w", newline="", encoding="utf-8")
        except Exception as e:
            self.error = e
            return
        get_ts = getattr(rtde_r, "getTimestamp", None)
        extras = [(name, getattr(rtde_r, meth, None)) for name, meth in EXTRA_SIGNALS]
        self.extra_missing = [name for name, getter in extras if getter is None]
        get_sf = getattr(rtde_r, "getTargetSpeedFraction", None)
        writer, last_ctrl_ts, t0_ctrl = None, None, None
        poll_interval = self.period / 5.0
        try:
            while not self._stop_flag.is_set():
                self.loops += 1
                ctrl_ts = get_ts() if get_ts is not None else None
                if ctrl_ts is not None:
                    if ctrl_ts == last_ctrl_ts:
                        time.sleep(poll_interval)
                        continue
                    if last_ctrl_ts is not None and (ctrl_ts - last_ctrl_ts) > self.period * 1.5:
                        self.gaps += 1
                    last_ctrl_ts = ctrl_ts
                    if t0_ctrl is None:
                        t0_ctrl = ctrl_ts
                    t_main = ctrl_ts - t0_ctrl
                else:
                    t_main = time.perf_counter()
                q = rtde_r.getActualQ(); qd = rtde_r.getActualQd()
                i = rtde_r.getActualCurrent(); tqd = rtde_r.getTargetQd()
                try:
                    temps = rtde_r.getJointTemperatures()
                except Exception:
                    temps = [None] * 6
                row = {"timestamp": round(t_main, 6)}
                for j in range(6):
                    row[f"actual_q_{j}"] = q[j]
                    row[f"actual_qd_{j}"] = qd[j]
                    row[f"actual_current_{j}"] = i[j]
                    row[f"target_qd_{j}"] = tqd[j]
                    row[f"joint_temp_{j}"] = temps[j]
                for name, getter in extras:
                    vals = None
                    if getter is not None:
                        try:
                            vals = list(getter())
                        except Exception:
                            vals = None
                    if vals is None or len(vals) != 6:
                        vals = [None] * 6
                    else:
                        self.extra_reads[name] += 1
                    for j in range(6):
                        row[f"{name}_{j}"] = vals[j]
                sf = None
                if get_sf is not None:
                    try:
                        sf = _finite_float(get_sf())
                    except Exception:
                        sf = None
                row["target_speed_fraction"] = sf
                if writer is None:
                    writer = csv.DictWriter(csv_file, fieldnames=list(row.keys()))
                    writer.writeheader()
                writer.writerow(row)
                self.rows += 1
                self.ctrl_times.append(t_main)
                j = self.joint_index
                self.last_q0, self.last_qd0, self.last_q_full = q[j], qd[j], list(q)
                self.q0_min = q[j] if self.q0_min is None else min(self.q0_min, q[j])
                self.q0_max = q[j] if self.q0_max is None else max(self.q0_max, q[j])
                self.last_speed_fraction = sf
                if sf is None:
                    self.speed_fraction_missing += 1
                else:
                    self.min_speed_fraction = sf if self.min_speed_fraction is None else min(self.min_speed_fraction, sf)
                self.started_ok.set()
                if self.rows % 250 == 0:
                    csv_file.flush()
                time.sleep(poll_interval)
        finally:
            csv_file.flush()
            csv_file.close()
            rtde_r.disconnect()

    def summary(self):
        lines = [f"迴圈次數 : {self.loops}, 寫入列數 : {self.rows}", f"疑似漏拍 : {self.gaps} 次"]
        if len(self.ctrl_times) >= 3:
            d = [b - a for a, b in zip(self.ctrl_times[:-1], self.ctrl_times[1:])]
            lines.append(f"平均間隔 : {statistics.fmean(d)*1000:.4f} ms, std {statistics.pstdev(d)*1000:.4f} ms")
        lines.append("額外訊號（成功讀取列數）: " + "，".join(
            f"{name} {n}" + ("（ur_rtde 無此方法）" if name in self.extra_missing else "")
            for name, n in self.extra_reads.items()))
        if self.min_speed_fraction is None:
            lines.append("教導器速度滑桿：整段讀不到")
        else:
            lines.append(f"教導器速度滑桿（記錄中最低）: {100 * self.min_speed_fraction:.1f}%"
                         + (f"，讀不到 {self.speed_fraction_missing} 列" if self.speed_fraction_missing else ""))
        return "\n".join(lines)

    def speed_fraction_held(self, min_fraction=SPEED_FRACTION_MIN):
        """整段記錄中滑桿都讀得到，且不低於 min_fraction。"""
        return (self.rows > 0 and self.speed_fraction_missing == 0
                and self.min_speed_fraction is not None and self.min_speed_fraction >= min_fraction)


# ============================================================
# 第四部分：離線分析
# ============================================================

def fk_frames(q):
    """UR5 官方 DH 的各座標系（4×4），frames[k] 為第 k 軸之後的座標系；frames[j] 的 z 軸為 J(j) 的轉軸。"""
    T = np.eye(4)
    frames = [T.copy()]
    for th, d, a, al in zip(q, DH_D, DH_A, DH_ALPHA):
        ct, st, ca, sa = math.cos(th), math.sin(th), math.cos(al), math.sin(al)
        T = T @ np.array([[ct, -st * ca, st * sa, a * ct],
                          [st, ct * ca, -ct * sa, a * st],
                          [0.0, sa, ca, d],
                          [0.0, 0.0, 0.0, 1.0]])
        frames.append(T.copy())
    return frames


def gravity_lever(q, joint_index=1, tool_offset_m=0.0):
    """
    沿工具 z 軸距法蘭 tool_offset_m 的點質量，每 1 N 重力對 J(joint_index) 轉軸的力矩 [m]
    （正負號依 DH 的轉軸方向）。工具垂直朝下時與 tool_offset_m 無關。
    """
    fr = fk_frames(q)
    p = fr[6][:3, 3] + fr[6][:3, 2] * tool_offset_m
    o, z = fr[joint_index][:3, 3], fr[joint_index][:3, 2]
    return float(z @ np.cross(p - o, np.array([0.0, 0.0, -1.0])))


def find_holds(t, tqd, min_s, qd_eps=1e-6):
    """
    target_qd ≡ 0 且持續 ≥ min_s 秒的區段，以及靠近方向（區段之前最後一個非零 target_qd 的正負號）。
    記錄開始時的靜止（之前沒有運動）不算。回傳 [(i0, i1, 方向)]，i1 不含。
    """
    still = np.abs(tqd) <= qd_eps
    out = []
    idx = np.where(still)[0]
    if not len(idx):
        return out
    for run in np.split(idx, np.where(np.diff(idx) != 1)[0] + 1):
        i0, i1 = run[0], run[-1] + 1
        if t[run[-1]] - t[i0] < min_s or i0 == 0:
            continue
        out.append((i0, i1, int(np.sign(tqd[i0 - 1]))))
    return out


def analyze_static_current(csv_path, joint_index=JOINT_INDEX, hold_s=HOLD_S, skip_s=SKIP_S, kt_out=KT_OUT):
    """
    每次靜止取 [開始 + skip_s, 開始 + hold_s) 的平均：電流、target_moment、六軸角度、溫度。
      i_mid  = (正向靠近平均 + 負向靠近平均) / 2   —— 兩方向靜摩擦相反，平均後消去
      i_half = (正向靠近平均 − 負向靠近平均) / 2   —— 靜摩擦帶寬（電流）
    sd_mid 為各次循環（相鄰的一正一負）中心電流的標準差。
    """
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    def col(key):
        if key not in rows[0]:
            return np.full(len(rows), np.nan)
        return np.array([float(r[key]) if r.get(key) not in (None, "", "None") else np.nan for r in rows])

    t = col("timestamp")
    tqd = col(f"target_qd_{joint_index}")
    cur = col(f"actual_current_{joint_index}")
    tm = col(f"target_moment_{joint_index}")
    temp = col(f"joint_temp_{joint_index}")
    q_all = np.column_stack([col(f"actual_q_{j}") for j in range(6)])

    holds = []
    for i0, i1, sgn in find_holds(t, tqd, 0.8 * hold_s):
        m = (t[i0:i1] >= t[i0] + skip_s) & (t[i0:i1] < t[i0] + hold_s)
        sl = np.arange(i0, i1)[m]
        if len(sl) < 10 or sgn == 0:
            continue
        holds.append(dict(dir=sgn, n=len(sl), i=float(np.mean(cur[sl])), i_sd=float(np.std(cur[sl])),
                          tm=float(np.nanmean(tm[sl])) if np.any(np.isfinite(tm[sl])) else float("nan"),
                          temp=float(np.nanmean(temp[sl])) if np.any(np.isfinite(temp[sl])) else float("nan"),
                          q=np.mean(q_all[sl], axis=0)))
    pos = [h for h in holds if h["dir"] > 0]
    neg = [h for h in holds if h["dir"] < 0]
    res = dict(holds=holds, n_pos=len(pos), n_neg=len(neg))
    if not pos or not neg:
        print(f"[警告] 正向靠近 {len(pos)} 次、負向靠近 {len(neg)} 次，無法消去靜摩擦。")
        return res
    i_pos, i_neg = np.mean([h["i"] for h in pos]), np.mean([h["i"] for h in neg])
    mids = [(p["i"] + n["i"]) / 2 for p, n in zip(pos, neg)]
    res.update(i_pos=float(i_pos), i_neg=float(i_neg), i_mid=float((i_pos + i_neg) / 2),
               i_half=float((i_pos - i_neg) / 2), sd_mid=float(np.std(mids, ddof=1)) if len(mids) > 1 else float("nan"),
               n_pair=len(mids), tm_mid=float(np.mean([h["tm"] for h in holds])),
               q_mean=np.mean([h["q"] for h in holds], axis=0), temp=float(np.nanmean([h["temp"] for h in holds])))

    print(f"[J{joint_index} 靜止電流] 正向靠近 {len(pos)} 次、負向靠近 {len(neg)} 次")
    print(f"  正向靠近平均 {i_pos:+.4f} A，負向靠近平均 {i_neg:+.4f} A")
    print(f"  中心電流 i_mid = {res['i_mid']:+.4f} A（各次循環的標準差 {res['sd_mid']:.4f} A，{len(mids)} 次）")
    print(f"  靜摩擦帶寬 ±{abs(res['i_half']):.4f} A（× {kt_out:.2f} = ±{abs(res['i_half']) * kt_out:.2f} N·m）")
    if math.isfinite(res["tm_mid"]) and abs(res["i_mid"]) > 1e-6:
        print(f"  控制器重力模型 target_moment = {res['tm_mid']:+.3f} N·m；÷ 中心電流 = "
              f"{res['tm_mid'] / res['i_mid']:.3f} N·m/A（只在控制器的連桿與負載模型正確時等於換算比例）")
    return res


def compare_known_mass(res_with, res_without, mass_kg, joint_index=JOINT_INDEX, tool_offset_m=0.0, g=G):
    """
    同姿態「夾砝碼 − 不夾」：Δτ = m·g·d（d 由兩次的平均關節角以 UR5 DH 計算），K = |Δτ| ÷ |Δi|。
    target_moment 的差 ÷ Δτ 應接近 1（教導器負載已更新、d 正確）；Δi 與 Δtarget_moment 應同號。
    """
    q = (np.asarray(res_with["q_mean"]) + np.asarray(res_without["q_mean"])) / 2
    lever = gravity_lever(q, joint_index, tool_offset_m)
    d_tau = mass_kg * g * lever
    d_i = res_with["i_mid"] - res_without["i_mid"]
    sd = math.sqrt(res_with["sd_mid"] ** 2 / res_with["n_pair"] + res_without["sd_mid"] ** 2 / res_without["n_pair"])
    K = abs(d_tau) / abs(d_i)
    d_tm = res_with["tm_mid"] - res_without["tm_mid"]
    out = dict(lever=lever, d_tau=d_tau, d_i=d_i, sd_d_i=sd, K=K, sd_K=K * sd / abs(d_i), Kt=K / GEAR_RATIO,
               d_tm=d_tm, tm_ratio=abs(d_tm) / abs(d_tau) if abs(d_tau) > 0 else float("nan"),
               same_sign=bool(np.sign(d_i) == np.sign(d_tm)) if math.isfinite(d_tm) else None)
    print(f"[已知質量] m = {mass_kg:.4f} kg，離 J{joint_index} 軸水平 {abs(lever):.4f} m，Δτ = {abs(d_tau):.3f} N·m")
    print(f"  Δi = {d_i:+.4f} A（標準差 {sd:.4f} A）→ K = {K:.3f} ± {out['sd_K']:.3f} N·m/A"
          f"（Kt = {out['Kt']:.4f} N·m/A，以減速比 {GEAR_RATIO} 換算）")
    if math.isfinite(d_tm):
        print(f"  控制器 target_moment 的差 {d_tm:+.3f} N·m，為 m·g·d 的 {out['tm_ratio']:.3f} 倍"
              + ("" if out["same_sign"] else "  [警告] 與 Δi 反號，請檢查資料"))
    return out


# ============================================================
# 主程式
# ============================================================

def main():
    print("=" * 70)
    print(" J1 靜止電流測試（兩方向靠近）：安全檢查")
    print("=" * 70)
    print(f"姿態: {POSE}（J0~J5 期望 {list(POSES[POSE]['q_deg']) if POSE in POSES else '未知'}°）")
    print(f"負載: PAYLOAD_KG = {PAYLOAD_KG}")

    n_cycle = stage_cycles(TEST_STAGE)
    print(f"模式: {TEST_STAGE}，{n_cycle} 次循環（正向靠近、負向靠近各 {n_cycle} 次），J1 ±{STEP_DEG}°")

    ok, report = safety_check(STEP_DEG, MOVE_VEL, MOVE_ACCEL, HOLD_S, SKIP_S, POSE)
    print("\n".join(report))
    payload_ok, payload_msg = check_payload_setting(PAYLOAD_KG)
    print(payload_msg + ("  [OK]" if payload_ok else "  [FAIL]"))
    if not (ok and payload_ok):
        print("\n[中止] 安全檢查未通過。")
        return

    print("\n[通過安全檢查]")
    if ROBOT_IP is None:
        print("\n[中止] ROBOT_IP 未設定。")
        return
    if rtde_receive is None:
        print("\n[中止] 找不到 rtde_receive，僅完成離線設計檢查。")
        return

    d = POSES[POSE]["d_j1_m"]
    print(f"\n[注意] 姿態 {POSE}：六軸須先以教導器擺到 {list(POSES[POSE]['q_deg'])}°（程式不會擺姿態）。")
    print(f"       J1 會來回轉 ±{STEP_DEG}°，末端（含砝碼）上下移動約 ±{1000 * d * math.radians(STEP_DEG):.0f} mm。")
    print("       請確認工具與砝碼下方淨空、J1 附近線材不會被拉扯；連上 RTDE 後會再以實際角度檢查六軸。")
    if input(f"已確認，輸入 {POSE} 繼續：").strip() != POSE:
        print("已取消。")
        return

    print("\n[負載確認] 請看教導器上設定的負載質量（夾上或放下砝碼後必須先在教導器更新質量與重心）。")
    entry_ok, entry_msg = check_payload_entry(PAYLOAD_KG, input("輸入教導器上的負載質量 [kg]："))
    print(entry_msg)
    if not entry_ok:
        print("\n[中止] 負載確認未通過，未送出任何軌跡。")
        return

    est_total = estimate_total_time(n_cycle, STEP_DEG, MOVE_VEL, MOVE_ACCEL, HOLD_S, PAUSE_TIME)
    print(f"\n即將對 IP={ROBOT_IP} 送出軌跡，只動 J{JOINT_INDEX}，姿態 {POSE}，負載 {PAYLOAD_KG} kg，"
          f"預估 {est_total:.0f} 秒。")
    print("*** 送出後 Ctrl+C 會送出 stopj，但最可靠的仍是緊急停止按鈕 ***")
    input("按 Enter 繼續，或 Ctrl+C 取消...")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"{POSE.lower()}_{payload_tag(PAYLOAD_KG)}_{TEST_STAGE.lower()}_{ts}"
    session_dir = os.path.join(OUTPUT_DIR, f"j1static_{tag}")
    os.makedirs(session_dir, exist_ok=True)
    csv_path = os.path.join(session_dir, f"j1static_data_{tag}.csv")
    info_path = os.path.join(session_dir, f"j1static_info_{tag}.txt")

    logger = RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, JOINT_INDEX)
    logger.start()
    if not logger.started_ok.wait(5.0):
        logger.stop()
        print(f"\n[中止] RTDE 記錄未能啟動：{logger.error}。未送出任何軌跡，手臂未動作。")
        return

    pose_ok, pose_lines = check_pose(logger.last_q_full, POSE)
    print("\n[姿態檢查（RTDE 實際角度）]\n" + "\n".join(pose_lines))
    if not pose_ok:
        logger.stop()
        print("\n[中止] 姿態檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    sf_ok, sf_msg = check_speed_fraction(logger.last_speed_fraction)
    print(f"\n[速度滑桿檢查] {sf_msg}")
    if not sf_ok:
        logger.stop()
        print("\n[中止] 速度滑桿檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    q_start = list(logger.last_q_full)
    q1_deg = math.degrees(q_start[JOINT_INDEX])
    moves = plan_moves(q1_deg, STEP_DEG, n_cycle)
    script = build_urscript(moves, q_start, JOINT_INDEX, MOVE_VEL, MOVE_ACCEL, HOLD_S, PAUSE_TIME)
    status, elapsed, span = "unknown", 0.0, 0.0
    try:
        print(f"[資訊] 送出 URScript（{len(script.splitlines())} 行，預估 {est_total:.1f} 秒）...")
        send_urscript(ROBOT_IP, script)
        status, elapsed, span = wait_for_motion_complete(logger, est_total, quiet_s=HOLD_S + 2.0)
        if status == "no_motion":
            print(f"\n[警告] 送出後 {elapsed:.1f}s 內未偵測到 J{JOINT_INDEX} 運動。")
            print("       控制器可能未接受腳本（不在 remote control 模式，或腳本解析失敗）。")
        elif status == "timeout":
            print(f"\n[警告] 超過硬性逾時（實際 {elapsed:.1f}s / 預估 {est_total:.1f}s），主動中止。")
            send_abort(ROBOT_IP)
        else:
            print(f"\n[資訊] 運動結束，實際 {span:.1f}s / 預估 {est_total:.1f}s")
    except KeyboardInterrupt:
        print("\n[中止] 收到 Ctrl+C。")
        send_abort(ROBOT_IP)
        raise
    except BaseException:
        send_abort(ROBOT_IP)
        raise
    finally:
        logger.stop()

    measured = []
    if logger.q0_min is not None and status == "done":
        q_lo, q_hi = math.degrees(logger.q0_min), math.degrees(logger.q0_max)
        q_end = math.degrees(logger.last_q0)
        measured = [f"實測 J{JOINT_INDEX} 範圍 : {q_lo:.3f}° 至 {q_hi:.3f}°（規劃 {q1_deg - STEP_DEG:.3f}° 至 "
                    f"{q1_deg + STEP_DEG:.3f}°）",
                    f"結束位置 : {q_end:.3f}°（開始 {q1_deg:.3f}°）"]
        if q_lo < q1_deg - STEP_DEG - 0.2 or q_hi > q1_deg + STEP_DEG + 0.2:
            measured.append("[警告] 實測範圍超出規劃範圍 0.2° 以上，請先查明原因再跑下一階段。")
        if abs(q_end - q1_deg) > 0.05:
            measured.append("[警告] 結束位置與開始位置相差超過 0.05°。")
        print("\n" + "\n".join(measured))
    if not logger.speed_fraction_held():
        print("\n[警告] 執行中教導器速度滑桿曾低於 100% 或讀不到（見摘要），這次的資料不可用。")

    with open(info_path, "w", encoding="utf-8") as f:
        f.write(f"模式: {TEST_STAGE}\n循環次數: {n_cycle}\n負載 PAYLOAD_KG: {PAYLOAD_KG} kg\n")
        f.write("\n".join(pose_lines) + "\n")
        f.write(f"J1 移動: ±{STEP_DEG}°，速度 {MOVE_VEL} rad/s，加速度 {MOVE_ACCEL} rad/s^2\n")
        f.write(f"靜止時間: {HOLD_S}s（分析略過前 {SKIP_S}s）\n預估總時長: {est_total:.1f}s\n執行狀態: {status}\n")
        f.write(logger.summary() + "\n")
        if measured:
            f.write("\n" + "\n".join(measured) + "\n")
    print(f"\n[完成] 資料存至 {csv_path}")
    print(f"[完成] 摘要存至 {info_path}")
    if status == "done":
        print()
        analyze_static_current(csv_path)
    print("\n下一步：同姿態夾砝碼與不夾各跑一次後，"
          "compare_known_mass(analyze_static_current(有砝碼的 csv), analyze_static_current(沒有的 csv), 砝碼質量)")


if __name__ == "__main__":
    main()
