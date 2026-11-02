"""
UR5 J1 固定角度窗口定速實驗 —— 鑑別 J1 的 B、Tc（正反轉經過同一段角度，消去重力）
================================================================

★★★ 跟 J0 定速（ur5_const_velocity_angle_ident.py）的差異 ★★★

J0 是垂直軸，重力不作用在 J0 上；J1（肩關節）是水平軸，定速時馬達除了摩擦，還要撐住手臂的重力：
  tau = B·qd + Tc·d + H(q1)      d：方向（+1、−1）；H(q1)：撐住重力的力矩，隨 J1 角度改變
（tau 一律為 電流 × KT_OUT。）沿用 J0 的固定角度窗口：所有速度、正反轉都以定速掃過**同一段 J1 角度窗口**
（整數個馬達圈），各段取窗口內的平均。H 平均後為窗口平均 H̄，每個速度、兩個方向都相同，
因此每個方向分開回歸 tau = B·qd + C·d：
  正轉 C₊ = Tc + H̄，反轉 C₋ = Tc − H̄  →  Tc = (C₊ + C₋)/2，H̄ = (C₊ − C₋)/2
不需要連桿質量模型（combine_directions）。

另記錄控制器的 target_moment（重力模型：UR 官方連桿質量 + 教導器負載；只有 J1 定速時即為 H）。
量到的 H̄ 與 target_moment 窗口平均之比、窗口內兩者對角度的斜率之比，可作為 J1 換算比例的參考
（推測用：前提是控制器的重力模型正確）。

窗口與動作（J0 = 180°，J2~J5 為診斷原點；J1 角度減少時手臂往回抬，增加時手臂下降）：
  - 窗口 6 個馬達圈（21.39°），中心 −109°：[−119.69°, −98.31°]。
  - 每一段同 J0：movej 到起點 → speedj 加速並穩定 → 定速掃過窗口，控制器每週期讀 J1 角度，超過終點才 stopj。
  - 所有段的起點、煞車後最遠點都在 −127.71°～−90.29°（−127.71° 為 FULL 正轉 0.20 rad/s 的起點）：
    手臂只在診斷原點往回抬的一側活動，不會低於診斷原點。全部段跑完後 movej 回到 J1 = −90°（診斷原點）。
  - 法蘭離底座安裝面的高度：−90° 時 0.43 m，−127.7° 時 0.66 m；手肘最多移到底座後方約 0.25 m（高度約 0.43 m）。
  - 重力力矩（官方連桿質量 + 末端工具 1.30 kg）在 −119.5° 附近為零，窗口內 −0.2～+17.1 N·m（平均 +8.6）。
    正轉（手臂下降）時馬達力矩在窗口內由正變負（減速機受力面換邊）；反轉不換邊。
  - FULL 每個方向 12 段，約 2.4 分鐘。

★★★ 使用前必讀 ★★★
1. 【新程式，尚未經學長審查，也尚未上機】依 TEST_STAGE = QUICK → MID → FULL 階梯執行。
2. 【兩個方向分開執行】DIRECTION = +1 跑完、人工確認後再改 −1。
3. 【J1 視窗未設定一律中止】J1_SAFE_MIN_DEG / J1_SAFE_MAX_DEG 由操作者依現場淨空填寫；所有段的起點、
   煞車後最遠點、結束點各外推 SAFETY_MARGIN_DEG 後都必須在視窗內（FULL 為 −129.71°～−88.00°）。
4. 【法蘭高度】以 UR5 官方 DH 計算規劃範圍內（含餘裕）法蘭離底座安裝面的最低高度，低於 FLANGE_Z_MIN_M 即中止
   （J1 角度增加時手臂下降；防止窗口設到下降側）。
5. 【姿態】程式不會擺姿態：須先以教導器擺到診斷原點（J0 = 180°）；連上 RTDE 後以實際角度檢查六軸（±0.5°），
   不符即中止。執行前須輸入 J1，確認手臂往回抬時上方與底座後方淨空、線材不會被拉扯。
6. 【負載 PAYLOAD_KG】同速度調變：填教導器上設定的負載質量（只有末端工具時 1.30 kg），未設定一律中止；
   執行前須輸入教導器上的負載質量，不符即中止。
7. 【第一次 movej 的距離有上限】J1 與第一段起點相差超過 MAX_INITIAL_MOVE_DEG 就中止（各階段最大為 MID 正轉的 31.3°）。
8. 【教導器速度滑桿】讀 target_speed_fraction，不是 100% 或讀不到一律中止；執行中的最低值寫進摘要。
9. 【IP 需要你自己再次確認】
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

DIRECTION = +1        # +1 = J1 角度增加（手臂下降），-1 = J1 角度減少（手臂往上抬）。兩個方向分開執行，不要自動連續做
# QUICK: 單一檔 0.05 rad/s、窗口 1 個馬達圈（3.56°），驗證方向、角度停止與 movej 定位
# MID  : 三檔 0.05 / 0.12 / 0.20 rad/s、完整窗口、只升速一趟，驗證多段與較長腳本
# FULL : 六檔、完整窗口、先升後降（12 段）
TEST_STAGE = "QUICK"

SPEED_LEVELS = [0.02, 0.05, 0.08, 0.12, 0.16, 0.20]   # rad/s（同 J0）
SPEED_ORDER = "UP_DOWN"                              # FULL 用；"UP" 或 "UP_DOWN"

GEAR_RATIO = 101
N_MOTOR_REVS_WINDOW = 6         # 窗口長度：6 個馬達圈 = 6 × 360/101 = 21.39°
WINDOW_CENTER_DEG = -109.0      # 窗口中心（絕對角度）；FULL 反轉 0.20 rad/s 的起點為 −90.29°，不低於診斷原點
END_J1_DEG = -90.0              # 全部段跑完後 movej 回到的 J1 角度（= 診斷原點）
SETTLE_TIME = 0.5               # 進入窗口前，在目標速度至少穩定的時間 [s]
OVERRUN_DEG = 0.5               # 越過窗口終點後再走多少才 stopj [deg]

SPEEDJ_ACCEL = 0.5              # speedj 加速度 [rad/s^2]
STOPJ_DECEL = 1.0               # stopj 減速度 [rad/s^2]
MOVEJ_ACCEL = 0.5               # 段間 movej 加速度 [rad/s^2]
MOVEJ_VEL = 0.15                # 段間 movej 速度 [rad/s]
PAUSE_TIME = 0.5                # 每段 stopj 後、movej 後的靜止時間 [s]（須小於結束判定的 1 s 靜止窗）

QD_MAX = 0.30                   # 速度上限（同 J0：每圈 12 次的漣波約 0.32 rad/s 以上在 125 Hz 下混疊）
QDD_CHECK_LIMIT = 3.0           # 安全檢查用加速度上限
LOOP_TIME_FACTOR = 1.5          # URScript 迴圈次數上限 = 預估時間 × 此係數 / dt
SAFETY_MARGIN_DEG = 2.0         # 絕對角度檢查時，起點與最遠點各再外推的餘裕 [deg]
MAX_INITIAL_MOVE_DEG = 40.0     # 目前 J1 到第一段起點的最大允許距離 [deg]
FLANGE_Z_MIN_M = 0.40           # 規劃範圍（含餘裕）內法蘭離底座安裝面的最低高度 [m]（診斷原點為 0.43 m）
SPEED_FRACTION_MIN = 0.999      # 教導器速度滑桿（RTDE target_speed_fraction，0～1）的下限

SAMPLE_HZ = 125.0
DT = 1.0 / SAMPLE_HZ
KT_OUT = GEAR_RATIO * 0.1350    # 本研究目前的換算比例 [N·m/A]（同 J0）

# J1 允許的絕對角度活動視窗。必須由操作者依現場線材與淨空狀況填寫，
# 未填寫則程式中止 —— 不提供預設值。
J1_SAFE_MIN_DEG = None
J1_SAFE_MAX_DEG = None

# 開始姿態：診斷原點，六軸 [deg]（J0 = 180°）。程式不會擺姿態，連上 RTDE 後以實際角度檢查。
START_Q_DEG = (180.0, -90.0, 90.0, -90.0, -90.0, 0.0)
POSE_TOL_DEG = 0.5

# 法蘭負載：填教導器上設定的負載質量 [kg]（只有末端工具 FT300、夾爪、3D 相機時 1.30）。
# 只用於確認與記錄（寫進輸出資料夾名稱），不影響動作；未設定則程式中止 —— 不提供預設值。
PAYLOAD_KG = None
PAYLOAD_TOL_KG = 0.005

# 額外記錄的 RTDE 輸出（只記錄，不影響動作）：(欄名, ur_rtde 方法)。控制器或 ur_rtde 不支援時該欄留空。
EXTRA_SIGNALS = (("target_q", "getTargetQ"), ("target_moment", "getTargetMoment"),
                 ("target_current", "getTargetCurrent"), ("joint_control_output", "getJointControlOutput"))

# UR5（CB3）官方 DH 參數（只用於法蘭高度檢查）
DH_D = (0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823)
DH_A = (0.0, -0.425, -0.39225, 0.0, 0.0, 0.0)
DH_ALPHA = (math.pi / 2, 0.0, 0.0, math.pi / 2, -math.pi / 2, 0.0)

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "outputs")   # 存到專案根目錄的 outputs/，不和程式混在一起


# ============================================================
# 第一部分：角度窗口、每段規劃與安全檢查
# ============================================================

def window_bounds_deg(center_deg, n_revs, gear_ratio):
    """回傳窗口 (THETA_A, THETA_B)，THETA_A < THETA_B，長度 = n_revs 個馬達圈。"""
    width = n_revs * 360.0 / gear_ratio
    return center_deg - width / 2.0, center_deg + width / 2.0


def speed_sequence(levels, order):
    if order == "UP":
        return list(levels)
    if order == "UP_DOWN":
        return list(levels) + list(reversed(levels))
    raise ValueError(f"未知的 SPEED_ORDER: {order}")


def plan_segments(speeds, direction, theta_a_deg, theta_b_deg, accel, decel,
                  settle_time, overrun_deg, loop_time_factor, dt):
    """
    每一段的絕對角度規劃（單位：度；同 J0 定速）。direction=+1 從 THETA_A 掃到 THETA_B，-1 反之。
      start_deg    : movej 目標（窗口起點往回退「加速距離 + 穩定距離」）
      win_start_deg/win_end_deg : 窗口（分析只用這段）
      stop_cmd_deg : 越過此角度就 stopj
      extreme_deg  : 煞車後理論最遠點
      n_max        : URScript 迴圈次數上限
    """
    win_start = theta_a_deg if direction > 0 else theta_b_deg
    win_end = theta_b_deg if direction > 0 else theta_a_deg
    segs = []
    for k, v in enumerate(speeds):
        d_ramp = v * v / (2.0 * accel)                   # rad
        d_settle = v * settle_time                       # rad
        pre_deg = math.degrees(d_ramp + d_settle)
        start = win_start - direction * pre_deg
        stop_cmd = win_end + direction * overrun_deg
        extreme = stop_cmd + direction * math.degrees(v * v / (2.0 * decel))
        travel_rad = math.radians(abs(stop_cmd - start))
        t_expect = travel_rad / v + v / accel            # 加速段比定速多花 v/(2a)，取 v/a 保守
        n_max = int(math.ceil(loop_time_factor * t_expect / dt)) + 10
        segs.append(dict(index=k, speed=v, start_deg=start, win_start_deg=win_start,
                         win_end_deg=win_end, stop_cmd_deg=stop_cmd, extreme_deg=extreme,
                         t_expect=t_expect, n_max=n_max))
    return segs


def estimate_movej_time(dist_rad, v_peak, accel):
    """movej 時間估計（僅用於預估總時長，不是安全門檻）。"""
    if dist_rad <= 0:
        return 0.0
    t_up = v_peak / accel
    d_up = 0.5 * accel * t_up ** 2
    if dist_rad <= 2 * d_up:
        return 2.0 * math.sqrt(dist_rad / accel)
    return 2.0 * t_up + (dist_rad - 2 * d_up) / v_peak


def estimate_total_time(segs, q0_deg, decel, end_deg=None):
    total, pos = 0.0, q0_deg
    for s in segs:
        total += estimate_movej_time(math.radians(abs(s["start_deg"] - pos)), MOVEJ_VEL, MOVEJ_ACCEL)
        total += PAUSE_TIME + s["t_expect"] + s["speed"] / decel + PAUSE_TIME
        pos = s["extreme_deg"]
    if end_deg is not None:                                   # 結束後回原點
        total += estimate_movej_time(math.radians(abs(end_deg - pos)), MOVEJ_VEL, MOVEJ_ACCEL) + PAUSE_TIME
    return total


def fk_frames(q):
    """UR5 官方 DH 的各座標系（4×4），frames[k] 為第 k 軸之後的座標系；frames[6] 為法蘭。"""
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


def flange_height(q1_deg, q_deg=START_Q_DEG, joint_index=JOINT_INDEX):
    """法蘭中心離底座安裝面的高度 [m]：六軸角度 q_deg [deg]，其中 joint_index 換成 q1_deg。"""
    q = [math.radians(x) for x in q_deg]
    q[joint_index] = math.radians(q1_deg)
    return float(fk_frames(q)[6][2, 3])


def safety_check(segs, qd_max, qdd_max, accel, decel, movej_vel, movej_accel,
                 safe_min_deg, safe_max_deg, margin_deg, end_deg=None,
                 q_pose_deg=START_Q_DEG, flange_z_min=FLANGE_Z_MIN_M):
    """
    回傳 (ok, lines, (planned_min_deg, planned_max_deg))。
    1. 速度、加速度上限（同 J0）。
    2. 絕對角度：所有段的起點、煞車後最遠點與結束點，各外推 margin_deg，都必須在 J1 視窗內。
       J1 視窗未設定（None）一律 FAIL（預設拒絕）。
    3. 法蘭高度：規劃範圍（含餘裕）內每 0.25°，以 UR5 官方 DH（其他五軸 = q_pose_deg）計算法蘭離底座安裝面的高度，
       最低值須 ≥ flange_z_min；flange_z_min 不是有限數值時 FAIL。
    """
    lines, ok = [], True

    def chk(name, val, lim, unit):
        nonlocal ok
        lines.append(f"{name} : {val:.4f} {unit}  (限 {lim})")
        if val > lim:
            lines.append("  [FAIL]"); ok = False
        else:
            lines.append("  [OK]")

    chk("最高速度", max(s["speed"] for s in segs), qd_max, "rad/s")
    chk("speedj 加速度", accel, qdd_max, "rad/s^2")
    chk("stopj 減速度", decel, qdd_max, "rad/s^2")
    chk("movej 速度", movej_vel, qd_max, "rad/s")
    chk("movej 加速度", movej_accel, qdd_max, "rad/s^2")

    pts = [s["start_deg"] for s in segs] + [s["extreme_deg"] for s in segs]
    if end_deg is not None:
        pts.append(end_deg)
    lo_plan, hi_plan = min(pts) - margin_deg, max(pts) + margin_deg
    lines.append(f"規劃 J1 絕對角度範圍（含餘裕 {margin_deg:.1f}°） : {lo_plan:.2f}° 至 {hi_plan:.2f}°")
    if safe_min_deg is None or safe_max_deg is None:
        lines.append("  [FAIL] J1_SAFE_MIN_DEG / J1_SAFE_MAX_DEG 尚未設定，拒絕執行")
        ok = False
    else:
        lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
        if lo <= lo_plan and hi_plan <= hi:
            lines.append(f"  [OK] 在 J1 視窗 [{lo:.1f}°, {hi:.1f}°] 內")
        else:
            lines.append(f"  [FAIL] 超出 J1 視窗 [{lo:.1f}°, {hi:.1f}°]")
            ok = False

    grid = list(np.arange(lo_plan, hi_plan, 0.25)) + [hi_plan]
    z_min = min(flange_height(a, q_pose_deg) for a in grid)
    z_lim = _finite_float(flange_z_min)
    lines.append(f"法蘭離底座安裝面最低 : {z_min:.3f} m（規劃範圍含餘裕；其他五軸為開始姿態）  (限 ≥ {flange_z_min} m)")
    if z_lim is None or z_min < z_lim:
        lines.append("  [FAIL]"); ok = False
    else:
        lines.append("  [OK]")
    return ok, lines, (lo_plan, hi_plan)


def check_start_position(q1_deg, first_start_deg, max_move_deg, safe_min_deg, safe_max_deg):
    """手臂目前位置：J1 必須在 J1 視窗內，且離第一段起點不超過 max_move_deg。"""
    if safe_min_deg is None or safe_max_deg is None:
        return False, "J1_SAFE_MIN_DEG / J1_SAFE_MAX_DEG 尚未設定，拒絕執行"
    lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
    move = abs(first_start_deg - q1_deg)
    msg = (f"J1 目前 {q1_deg:.2f}°，第一段起點 {first_start_deg:.2f}°（需移動 {move:.2f}°，"
           f"上限 {max_move_deg:.1f}°），J1 視窗 [{lo:.1f}°, {hi:.1f}°]")
    if not (lo <= q1_deg <= hi):
        return False, msg + " -> 目前角度在視窗外"
    if move > max_move_deg:
        return False, msg + " -> 第一次 movej 距離過大"
    return True, msg + " -> OK"


def check_pose(q_full_rad, target_deg=START_Q_DEG, tol_deg=POSE_TOL_DEG):
    """
    以 RTDE 量到的六軸角度 [rad] 檢查 J0~J5 是否在 target_deg ± tol_deg 內（差值以 ±360° 等價處理）。
    預設拒絕：讀不到六軸角度、角度不是有限值、任一軸超差，一律回傳 False。回傳 (ok, lines)。
    """
    if q_full_rad is None or len(q_full_rad) != 6:
        return False, ["讀不到六軸角度，拒絕執行"]
    lines, ok = [f"開始姿態（診斷原點）：J0~J5 期望 {list(target_deg)}°，容差 ±{tol_deg}°"], True
    for j, tgt in enumerate(target_deg):
        x = _finite_float(q_full_rad[j])
        cur = math.degrees(x) if x is not None else float("nan")
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

def build_urscript(segs, direction, joint_index, q_others, accel, decel,
                   movej_accel, movej_vel, pause_time, dt, end_deg=None):
    """
    同 J0 定速。q_others：六軸絕對角度 [rad]（取自執行開始時 RTDE 量到的值），只替換 joint_index 當作各段
    movej 目標，其他五軸維持原值不動。
    每段：movej → 停 → 迴圈（每週期 speedj 一次，讀 J1 實際角度，越過 stop 角度或達次數上限即離開）
          → stopj → 停。
    end_deg 不為 None 時，全部段結束後再 movej 回 J1 = end_deg（其他五軸不變）。
    """
    lines = ["def j1_const_velocity_program():"]
    cmp_op = "<" if direction > 0 else ">"
    for s in segs:
        q_t = list(q_others)
        q_t[joint_index] = math.radians(s["start_deg"])
        q_t_str = "[" + ", ".join(f"{v:.6f}" for v in q_t) + "]"
        qd = ["0.0"] * 6
        qd[joint_index] = f"{direction * s['speed']:.6f}"
        qd_str = "[" + ", ".join(qd) + "]"
        stop_rad = math.radians(s["stop_cmd_deg"])
        lines.append(f"  # segment {s['index'] + 1}: {s['speed']:.2f} rad/s")   # URScript 內只放 ASCII
        lines.append(f"  movej({q_t_str}, a={movej_accel}, v={movej_vel})")
        lines.append(f"  sleep({pause_time})")
        lines.append("  q_now = get_actual_joint_positions()")
        lines.append("  n = 0")
        lines.append(f"  while (n < {s['n_max']}) and (q_now[{joint_index}] {cmp_op} {stop_rad:.6f}):")
        lines.append(f"    speedj({qd_str}, a={accel}, t={dt})")
        lines.append("    q_now = get_actual_joint_positions()")
        lines.append("    n = n + 1")
        lines.append("  end")
        lines.append(f"  stopj({decel})")
        lines.append(f"  sleep({pause_time})")
    if end_deg is not None:
        q_e = list(q_others)
        q_e[joint_index] = math.radians(end_deg)
        lines.append("  # return to home J1")
        lines.append(f"  movej([{', '.join(f'{v:.6f}' for v in q_e)}], a={movej_accel}, v={movej_vel})")
        lines.append(f"  sleep({pause_time})")
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


def wait_for_motion_complete(logger, planned_s, qd_eps=2e-3,
                             quiet_s=1.0, no_motion_s=10.0, timeout_factor=3.0):
    """同 J0 定速：J1 實際速度連續 quiet_s 秒低於門檻才視為結束（每段的靜止 PAUSE_TIME 須小於 quiet_s）。"""
    t_start = time.time()
    hard_timeout = planned_s * timeout_factor + 5.0
    moved = False
    quiet_since = None
    t_first_move = None
    t_last_move = None
    while True:
        now = time.time()
        elapsed = now - t_start
        qd0 = logger.last_qd0
        if qd0 is not None and abs(qd0) > qd_eps:
            if not moved:
                t_first_move = now
            moved = True
            t_last_move = now
            quiet_since = None
        elif moved:
            if quiet_since is None:
                quiet_since = now
            elif now - quiet_since >= quiet_s:
                return "done", elapsed, t_last_move - t_first_move
        if not moved and elapsed > no_motion_s:
            return "no_motion", elapsed, 0.0
        if elapsed > hard_timeout:
            span = (t_last_move - t_first_move) if moved else 0.0
            return "timeout", elapsed, span
        time.sleep(0.02)


# ============================================================
# 第三部分：RTDE 記錄（同 ur5_j1_static_current.py；last_q0 / last_qd0 為 joint_index 的值）
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

def analyze_angle_window(csv_path, joint_index, kt_out, direction, speed_levels,
                         theta_a_deg, theta_b_deg, gear_ratio=GEAR_RATIO):
    """
    同 J0 定速：每一段取 target_qd 等於該檔速度、且實際角度落在窗口 [THETA_A, THETA_B] 內的樣本
    （窗口為整數個馬達圈，漣波在平均中抵銷）。以各段平均（實測 qd、tau）回歸
      tau = B·qd + C·direction        C = Tc + direction·H̄（H̄：窗口內撐住重力的平均力矩）
    ± 為以段為點的單輪回歸標準差。各段另記 target_moment 的平均（控制器的重力模型；未記錄時為 NaN），
    以及 tau、target_moment 對角度的斜率 [N·m/°]（窗口內重力的變化）。
    回傳 dict：segments、B、C、sd_B、sd_C、升降速趟分開的 B、C，與各段斜率的平均 slope_tau、slope_tm。
    """
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    def col(key):
        if not rows or key not in rows[0]:
            return np.full(len(rows), np.nan)
        return np.array([float(r[key]) if r.get(key) not in (None, "", "None") else np.nan for r in rows])

    q = np.degrees(col(f"actual_q_{joint_index}"))
    qd = col(f"actual_qd_{joint_index}")
    tqd = col(f"target_qd_{joint_index}")
    tau = col(f"actual_current_{joint_index}") * kt_out
    tm = col(f"target_moment_{joint_index}")
    t = col("timestamp")
    d = float(direction)
    in_win = (q >= theta_a_deg) & (q <= theta_b_deg)
    rev_deg = 360.0 / gear_ratio

    segments = []
    for v in speed_levels:
        on = np.abs(tqd - d * v) < 1e-4
        idx = np.where(on & in_win)[0]
        if len(idx) == 0:
            continue
        for r in np.split(idx, np.where(np.diff(idx) != 1)[0] + 1):
            if len(r) < 10:
                continue
            has_tm = bool(np.all(np.isfinite(tm[r])))
            segments.append(dict(speed=v, idx=r, t0=float(t[r[0]]), n=len(r),
                                 revs=float((q[r].max() - q[r].min()) / rev_deg),
                                 qd_mean=float(np.mean(qd[r])), tau_mean=float(np.mean(tau[r])),
                                 tm_mean=float(np.mean(tm[r])) if has_tm else float("nan"),
                                 slope_tau=float(np.polyfit(q[r], tau[r], 1)[0]),
                                 slope_tm=float(np.polyfit(q[r], tm[r], 1)[0]) if has_tm else float("nan")))
    segments.sort(key=lambda s: s["t0"])

    print(f"{'段':>3}{'速度':>7}{'樣本':>6}{'覆蓋馬達圈':>10}{'平均qd':>10}{'平均tau':>10}{'平均tm':>10}")
    for k, s in enumerate(segments):
        print(f"{k+1:3d}{s['speed']:7.2f}{s['n']:6d}{s['revs']:10.2f}{s['qd_mean']:10.4f}"
              f"{s['tau_mean']:10.3f}{s['tm_mean']:10.3f}")
    out = dict(segments=segments, q_deg=q, tau=tau, tm=tm, t=t)
    if len(segments) < 2:
        print("[提示] 有效段數少於 2（QUICK 只有 1 段），不做回歸。")
        return out

    x = np.array([s["qd_mean"] for s in segments])
    y = np.array([s["tau_mean"] for s in segments])
    X = np.column_stack([x, np.full(len(x), d)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    dof = max(len(y) - 2, 1)
    sd = np.sqrt(np.diag(res @ res / dof * np.linalg.inv(X.T @ X)))
    out.update(B=float(coef[0]), C=float(coef[1]), sd_B=float(sd[0]), sd_C=float(sd[1]))
    print(f"B = {coef[0]:.3f} ± {sd[0]:.3f}，C = {coef[1]:.3f} ± {sd[1]:.3f}"
          f"（{len(y)} 段，± 為單輪回歸的標準差；C = Tc + 方向 × H̄）")
    # 升速趟、降速趟分開回歸（速度先升後降時）：兩者差異 = 時間或順序效應
    sp = [seg["speed"] for seg in segments]
    k = int(np.argmax(sp))
    if k + 1 < len(sp) and sp[k + 1] == sp[k]:
        k += 1                                        # 最高速有兩段時，第二段起算降速趟
    if k >= 2 and len(sp) - k >= 2:
        for name, sl in (("up", slice(0, k)), ("down", slice(k, None))):
            c, *_ = np.linalg.lstsq(X[sl], y[sl], rcond=None)
            out[f"B_{name}"], out[f"C_{name}"] = float(c[0]), float(c[1])
        print(f"升速趟 B = {out['B_up']:.3f}、C = {out['C_up']:.3f}；"
              f"降速趟 B = {out['B_down']:.3f}、C = {out['C_down']:.3f}")
    # tau、target_moment 對角度的斜率（各段直線擬合後平均）：窗口內重力的變化
    out["slope_tau"] = float(np.mean([s["slope_tau"] for s in segments]))
    out["slope_tm"] = float(np.mean([s["slope_tm"] for s in segments]))
    print(f"對角度斜率（各段平均）：tau {out['slope_tau']:+.4f}、target_moment {out['slope_tm']:+.4f} N·m/°")
    return out


def combine_directions(res_pos, res_neg, kt_out=KT_OUT):
    """
    正轉（+1）、反轉（−1）的 analyze_angle_window 結果合併（兩者須為同一個窗口）：
      Tc = (C₊ + C₋)/2，H̄ = (C₊ − C₋)/2，B = (B₊ + B₋)/2；± 由兩次單輪回歸的標準差合成。
    H̄：量到的窗口平均重力保持力矩（電流 × kt_out）。tm_win：控制器 target_moment 的窗口平均（兩方向所有段）。
      K_mean  = kt_out × tm_win ÷ H̄
      K_slope = kt_out × （target_moment 對角度的斜率）÷（tau 對角度的斜率）—— 只用窗口內重力的變化，與摩擦無關
    兩者都假設控制器的重力模型正確，只作 J1 換算比例的參考（推測），不是校正。
    """
    for r, name in ((res_pos, "正轉"), (res_neg, "反轉")):
        if "C" not in r:
            raise ValueError(f"{name}的結果沒有回歸係數（有效段數少於 2）")
    cp, cn = res_pos["C"], res_neg["C"]
    sd_c = 0.5 * math.hypot(res_pos["sd_C"], res_neg["sd_C"])
    out = dict(B=0.5 * (res_pos["B"] + res_neg["B"]), sd_B=0.5 * math.hypot(res_pos["sd_B"], res_neg["sd_B"]),
               B_pos=res_pos["B"], B_neg=res_neg["B"],
               Tc=0.5 * (cp + cn), H=0.5 * (cp - cn), sd_Tc=sd_c, sd_H=sd_c)
    segs = res_pos["segments"] + res_neg["segments"]
    tm_win = float(np.mean([s["tm_mean"] for s in segs]))
    s_tau = float(np.mean([s["slope_tau"] for s in segs]))
    s_tm = float(np.mean([s["slope_tm"] for s in segs]))
    out.update(tm_win=tm_win, slope_tau=s_tau, slope_tm=s_tm,
               K_mean=kt_out * tm_win / out["H"] if math.isfinite(tm_win) and abs(out["H"]) > 1e-9 else float("nan"),
               K_slope=kt_out * s_tm / s_tau if math.isfinite(s_tm) and abs(s_tau) > 1e-9 else float("nan"))
    print(f"[兩方向合併] B = {out['B']:.3f} ± {out['sd_B']:.3f}（正轉 {out['B_pos']:.3f}、反轉 {out['B_neg']:.3f}）")
    print(f"  Tc = {out['Tc']:.3f} ± {sd_c:.3f}，窗口平均重力保持力矩 H̄ = {out['H']:+.3f} ± {sd_c:.3f} N·m"
          f"（電流 × {kt_out:.2f}；± 為單輪回歸的標準差合成）")
    if math.isfinite(tm_win):
        print(f"  控制器 target_moment 窗口平均 {tm_win:+.3f} N·m → K_mean = {out['K_mean']:.2f} N·m/A；"
              f"對角度斜率 tau {s_tau:+.4f}、target_moment {s_tm:+.4f} N·m/° → K_slope = {out['K_slope']:.2f} N·m/A")
        print("  （推測用：假設控制器的重力模型正確）")
    return out


# ============================================================
# 主流程
# ============================================================

def stage_config(stage):
    if stage == "QUICK":
        return [0.05], "UP", 1
    if stage == "MID":
        return [0.05, 0.12, 0.20], "UP", N_MOTOR_REVS_WINDOW
    if stage == "FULL":
        return list(SPEED_LEVELS), SPEED_ORDER, N_MOTOR_REVS_WINDOW
    raise ValueError(f"未知的 TEST_STAGE: {stage}")


def main():
    print("=" * 70)
    print(" J1 固定角度窗口定速實驗：安全檢查")
    print("=" * 70)
    print(f"方向: {'+1（J1 角度增加，手臂下降）' if DIRECTION > 0 else '-1（J1 角度減少，手臂往上抬）'}")
    print(f"負載: PAYLOAD_KG = {PAYLOAD_KG}")

    levels, order, n_revs = stage_config(TEST_STAGE)
    speeds = speed_sequence(levels, order)
    theta_a, theta_b = window_bounds_deg(WINDOW_CENTER_DEG, n_revs, GEAR_RATIO)
    segs = plan_segments(speeds, DIRECTION, theta_a, theta_b, SPEEDJ_ACCEL, STOPJ_DECEL,
                         SETTLE_TIME, OVERRUN_DEG, LOOP_TIME_FACTOR, DT)
    print(f"模式: {TEST_STAGE}，速度順序: {[f'{v:.2f}' for v in speeds]}")
    print(f"角度窗口: {theta_a:.3f}° 至 {theta_b:.3f}°（{n_revs} 個馬達圈 = {theta_b - theta_a:.3f}°）")

    ok, report, (lo_plan, hi_plan) = safety_check(
        segs, QD_MAX, QDD_CHECK_LIMIT, SPEEDJ_ACCEL, STOPJ_DECEL, MOVEJ_VEL, MOVEJ_ACCEL,
        J1_SAFE_MIN_DEG, J1_SAFE_MAX_DEG, SAFETY_MARGIN_DEG, end_deg=END_J1_DEG)
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

    pts = [s["start_deg"] for s in segs] + [s["extreme_deg"] for s in segs] + [END_J1_DEG]
    print(f"\n[注意] 六軸須先以教導器擺到診斷原點 {list(START_Q_DEG)}°（程式不會擺姿態）。")
    print(f"       J1 會在 {min(pts):.1f}° 至 {max(pts):.1f}° 之間來回（角度減少時手臂往回抬），"
          f"法蘭離底座安裝面最高約 {flange_height(min(pts)):.2f} m，結束時回到 {END_J1_DEG:.0f}°。")
    print("       請確認手臂上方與底座後方淨空、J1 附近線材不會被拉扯；連上 RTDE 後會再以實際角度檢查六軸。")
    if input("已確認，輸入 J1 繼續：").strip() != "J1":
        print("已取消。")
        return

    print("\n[負載確認] 請看教導器上設定的負載質量。")
    entry_ok, entry_msg = check_payload_entry(PAYLOAD_KG, input("輸入教導器上的負載質量 [kg]："))
    print(entry_msg)
    if not entry_ok:
        print("\n[中止] 負載確認未通過，未送出任何軌跡。")
        return

    print(f"\n即將對 IP={ROBOT_IP} 送出軌跡，關節 J{JOINT_INDEX}，方向 {DIRECTION:+d}，共 {len(segs)} 段，"
          f"負載 {PAYLOAD_KG} kg。")
    print("*** 送出後 Ctrl+C 會送出 stopj，但最可靠的仍是緊急停止按鈕 ***")
    input("按 Enter 繼續，或 Ctrl+C 取消...")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"{'pos' if DIRECTION > 0 else 'neg'}_{payload_tag(PAYLOAD_KG)}_{TEST_STAGE.lower()}_{ts}"
    session_dir = os.path.join(OUTPUT_DIR, f"j1constvel_{tag}")
    os.makedirs(session_dir, exist_ok=True)
    csv_path = os.path.join(session_dir, f"j1constvel_data_{tag}.csv")
    info_path = os.path.join(session_dir, f"j1constvel_info_{tag}.txt")

    logger = RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, JOINT_INDEX)
    logger.start()
    if not logger.started_ok.wait(5.0):
        logger.stop()
        print(f"\n[中止] RTDE 記錄未能啟動：{logger.error}")
        print("       未送出任何軌跡，手臂未動作。")
        return

    pose_ok, pose_lines = check_pose(logger.last_q_full)
    print("\n[姿態檢查（RTDE 實際角度）]\n" + "\n".join(pose_lines))
    if not pose_ok:
        logger.stop()
        print("\n[中止] 姿態檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    q1_deg = math.degrees(logger.last_q0)
    pos_ok, pos_msg = check_start_position(q1_deg, segs[0]["start_deg"], MAX_INITIAL_MOVE_DEG,
                                           J1_SAFE_MIN_DEG, J1_SAFE_MAX_DEG)
    print(f"\n[起始位置檢查] {pos_msg}")
    if not pos_ok:
        logger.stop()
        print("\n[中止] 起始位置檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    sf_ok, sf_msg = check_speed_fraction(logger.last_speed_fraction)
    print(f"\n[速度滑桿檢查] {sf_msg}")
    if not sf_ok:
        logger.stop()
        print("\n[中止] 速度滑桿檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    script_text = build_urscript(segs, DIRECTION, JOINT_INDEX, logger.last_q_full,
                                 SPEEDJ_ACCEL, STOPJ_DECEL, MOVEJ_ACCEL, MOVEJ_VEL, PAUSE_TIME, DT,
                                 end_deg=END_J1_DEG)
    est_total = estimate_total_time(segs, q1_deg, STOPJ_DECEL, end_deg=END_J1_DEG)

    status, elapsed, motion_span = "unknown", 0.0, 0.0
    try:
        print(f"[資訊] 送出 URScript（預估總時長 {est_total:.1f} 秒）...")
        send_urscript(ROBOT_IP, script_text)
        status, elapsed, motion_span = wait_for_motion_complete(logger, est_total)
        if status == "no_motion":
            print(f"\n[警告] 送出後 {elapsed:.1f}s 內未偵測到 J{JOINT_INDEX} 運動。")
            print("       控制器可能未接受腳本（不在 remote control 模式，或腳本解析失敗）。")
        elif status == "timeout":
            print(f"\n[警告] 超過硬性逾時（實際 {elapsed:.1f}s / 預估 {est_total:.1f}s），主動中止。")
            send_abort(ROBOT_IP)
        else:
            print(f"\n[資訊] 運動結束，實際運動 {motion_span:.1f}s / 預估 {est_total:.1f}s")
    except KeyboardInterrupt:
        print("\n[中止] 收到 Ctrl+C。")
        send_abort(ROBOT_IP)
        raise
    except BaseException:
        send_abort(ROBOT_IP)
        raise
    finally:
        logger.stop()

    measured_lines = []
    if logger.q0_min is not None and status == "done":
        q_lo, q_hi = math.degrees(logger.q0_min), math.degrees(logger.q0_max)
        measured_lines = [
            f"實測 J1 範圍 : {q_lo:.3f}° 至 {q_hi:.3f}°（規劃含餘裕 {lo_plan:.3f}° 至 {hi_plan:.3f}°）",
            f"結束位置     : {math.degrees(logger.last_q0):.3f}°",
            f"實際運動時長 / 預估 : {motion_span:.2f}s / {est_total:.2f}s",
        ]
        if q_lo < lo_plan or q_hi > hi_plan:
            measured_lines.append("[警告] 實測範圍超出規劃範圍，請先查明原因再跑下一階段。")
        print("\n" + "\n".join(measured_lines))
    if not logger.speed_fraction_held():
        print("\n[警告] 執行中教導器速度滑桿曾低於 100% 或讀不到（見摘要），這次的資料不可用。")

    with open(info_path, "w", encoding="utf-8") as f:
        f.write(f"方向: {DIRECTION:+d}\n模式: {TEST_STAGE}\n負載 PAYLOAD_KG: {PAYLOAD_KG} kg\n")
        f.write(f"速度順序: {speeds}\n")
        f.write(f"角度窗口: {theta_a:.4f} 至 {theta_b:.4f} deg（{n_revs} 個馬達圈）\n")
        f.write(f"規劃範圍（含餘裕）: {lo_plan:.3f} 至 {hi_plan:.3f} deg\n")
        f.write("\n".join(pose_lines) + "\n")
        f.write(f"預估總時長: {est_total:.1f}s\n執行狀態: {status}\n")
        f.write(logger.summary() + "\n")
        if measured_lines:
            f.write("\n" + "\n".join(measured_lines) + "\n")

    print(f"\n[完成] 資料存至 {csv_path}")
    print(f"[完成] 摘要存至 {info_path}")
    if status == "done":
        print()
        try:
            analyze_angle_window(csv_path, JOINT_INDEX, KT_OUT, DIRECTION, SPEED_LEVELS, theta_a, theta_b)
        except Exception as e:
            print(f"[警告] 離線分析失敗（資料已存檔，不影響實驗）：{e}")
    if len(segs) >= 2:
        print(f"\n下一步：正、反轉都跑完後 combine_directions("
              f"analyze_angle_window(正轉 csv, {JOINT_INDEX}, KT_OUT, +1, SPEED_LEVELS, {theta_a:.4f}, {theta_b:.4f}), "
              f"analyze_angle_window(反轉 csv, {JOINT_INDEX}, KT_OUT, -1, SPEED_LEVELS, {theta_a:.4f}, {theta_b:.4f}))")
    else:
        print("\n下一步：QUICK 只確認動作（1 段，不回歸）；正、反轉都確認後再跑 MID。")


if __name__ == "__main__":
    main()
