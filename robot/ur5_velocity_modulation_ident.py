"""
UR5 J0 速度調變實驗 —— 鑑別 J（等效慣量），並檢查 J 是否隨頻率、振幅改變
================================================================

★★★ 為什麼要這個實驗 ★★★

定加速度 v3（10/13、10/15，三角形與梯形）的 J 會隨模型寫法改變（1.3～1.9）：
  (1) 每次循環都從靜止起步，加速度只出現在起步後的前幾度行程與減速段；起步後的摩擦變化
      （起步後行程效應）、減速段多出的摩擦，都和 J·a 綁在一起。
  (2) 10/15 殘差分解中，只改起步後行程項的節點，J 就從 1.47 變為 1.61。
這是資料本身分不開，需要換激勵。

這支的做法：
  A. 每一段先以 RAMP_ACCEL 加速到巡航速度 V_CRUISE，定速 SETTLE_TIME 後，在巡航速度上疊加
     正弦速度  v(t) = V_CRUISE + A·sin(2π f t)，持續 T_MOD 秒（整數個週期），之後 stopj。
     速度永遠為正（不歸零、不反向），庫倫摩擦方向固定。
  B. 只分析調變段（起步後行程約 6.3° 之後）：加速度是正弦，與起步、行程無關。
     力矩中與加速度同相的成分 → J；與速度同相的成分 → 局部摩擦斜率 B_local。
  C. 六個條件 (f, A) = (1 Hz, 0.06)、(2 Hz, 0.06)、(3 Hz, 0.06)、(2 Hz, 0.03)、(1 Hz, 0.03)、(3 Hz, 0.04) rad/s，
     加速度振幅 2πf·A = 0.38、0.75、1.13、0.38、0.19、0.75 rad/s²（前四個為 10/22 上機的條件）：
       - 各條件的 J 一致，且不隨模型寫法改變 → J 的方法可以定案。
       - 加速度振幅相同、頻率不同：(1 Hz, 0.06) 與 (2 Hz, 0.03)，(2 Hz, 0.06) 與 (3 Hz, 0.04)；
         (2 Hz, 0.06) 與 (2 Hz, 0.03) 頻率相同、振幅不同。若 J 隨振幅變化（J ≈ J∞ + k/A），表示有與加減速相關的
         摩擦差（v3 減速段多出的摩擦即屬此類）；若只隨頻率變化，則屬動態效應。
       - 10/22 實測 J 隨加速度振幅增加（0.38 → 1.13 rad/s² 時增加 13～16%），以 J∞ + k/A 與直線 a + b·A
         外推的結果不同；(1 Hz, 0.03) 的 0.19 rad/s² 可分辨兩種寫法。(3 Hz, 0.04) 檢查 3 Hz 的 J 偏高
         來自頻率還是振幅（10/22 的 3 Hz 同時改變了兩者）。
       - 電流相對編碼器若有延遲 Δ，B_local 會增加 J·Δ·ω²、J 會減少 B·Δ；
         由 B_local 隨頻率的變化可檢查此延遲。
  D. 調變段中心固定在 WINDOW_CENTER_DEG（180°，同定速實驗窗口中心），長 V_CRUISE·T_MOD = 34.4°，
     正轉與反轉掃過同一段角度。FULL 的第 k 次重複，起點錯開 (k − 3)/7 個馬達圈：
     馬達角漣波（10/15 實測以每圈 2、4、5、10、12 次為主）與調變的相位關係在 7 次重複中平均分布，
     7 的倍數次諧波以外都會平均掉（10/15 實測 7、14、21 次振幅皆 ≤ 0.03 N·m）。
     條件依重複輪流執行（k = 0 的各條件 → k = 1 的各條件 → …），溫度與時間漂移平均分配到各條件。
  E. 控制器迴圈：每個控制週期（8 ms）以 speedj 送一次速度，正弦由 URScript 內的計數器計算；
     同時讀 J0 實際角度，越過「調變終點 + OVERRUN_DEG」即離開迴圈（角度保護），另有次數上限。
     9/11 掃頻測試（J0 正弦運動，0.2～2 Hz）中，實際速度跟隨目標速度的增益約 1.00、延遲 2～3.5 ms。
  F. 時間：QUICK 1 段約 13 s，MID 6 段約 75 s，FULL 42 段約 8.6 分鐘（每個方向；10/22 的四個條件 FULL 為 5.7 分鐘）。
  G. 姿態 POSE：預設 DIAG（診斷原點，J1~J5 = −90°, 90°, −90°, −90°, 0°，PHM 基準）。
     EXT45 為手臂外伸（J1、J2 = −45°、45°，其餘同診斷原點，前臂仍水平），J0 的連桿慣量由約 0.9
     增為約 2.8 kg·m²（以 UR5 官方 DH、質量、質心與文獻連桿慣量張量計算）。兩個姿態的 J 差只來自連桿：
     J 對連桿慣量的斜率可檢查 Kt 的比例，截距為馬達轉子與減速機部分。程式不移動 J1~J5。
  H. 另記錄控制器的 target_qdd、target_moment、target_current、joint_control_output（只記錄，不影響動作）。
     target_moment 對 target_qdd 回歸得到控制器模型在此姿態的 J0 慣量 J_model，可與電流算出的 J 比較
     （Clochiatti et al. 2024 將 UR 的 target_moment 視為逆動力學前饋）。
  I. 已知質量測試：法蘭加上秤過的質量（裝在法蘭軸線上），同姿態加質量前後 J 的差應為 m·r²
     （r = 法蘭離 J0 軸距離，DIAG 0.50 m、EXT45 0.79 m），可直接檢查實測 J 的比例；控制器的 J_model
     會隨教導器的負載設定改變，可檢查設定是否正確。負載質量 PAYLOAD_KG 寫進輸出資料夾名稱。

離線驗證（合成資料：J = 1.5、10/15 實測的漣波與起步後行程曲線、白雜訊 0.4 N·m 加慢變雜訊 0.2 N·m）
-------------------------------------------------------------------------------------------
  各條件 J 的偏差 ≤ 0.02（1 Hz 約 −1%，見 analyze_velocity_modulation 的說明），標準差 0.01～0.04；
  減速時摩擦多 0.3 N·m：各條件 J 為 0.96～1.34（隨振幅變化），J 對 1/A 回歸找回 1.51 ± 0.02、減速摩擦 0.29；
  電流延遲 2 ms：各條件 J 為 1.41～1.47（≈ J − B_local·Δ = 1.43），由 B_local 隨頻率的變化估得延遲 2.0 ± 0.8 ms；
  模型不放漣波項時，起點錯開的 J 仍無偏差；不錯開則 2 Hz 的 J 偏低 0.09～0.17。

★★★ 使用前必讀 ★★★
1. 【審查狀態】commit 72591a7 的版本已通過學長審查，10/22 上機（前四個條件）。之後新增條件 5、6、負載確認、
   教導器速度滑桿檢查，尚未審查；依 TEST_STAGE = QUICK → MID → FULL 階梯執行。
2. 【兩個方向分開執行】DIRECTION=+1 跑完、人工確認後再改 −1。
3. 【新的 URScript 結構】調變段迴圈每週期計算 sin 後 speedj（加速度上限 MOD_ACCEL_LIMIT），
   這是新結構，需學長確認。迴圈仍遵守「第一行 def 頂格、其餘行至少縮排一格」。
4. 【J0 視窗預設拒絕】J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 未設定一律中止；所有段的起點、
   角度保護觸發後的最遠點（以最高速度煞車）、回原點角度 ± SAFETY_MARGIN_DEG 都須在視窗內。
   規劃範圍在固定角度窗口定速實驗 FULL 的規劃範圍（152.2°～207.8°，含餘裕）之內。
5. 【第一次 movej 的距離有上限】目前 J0 與第一段起點相差超過 MAX_INITIAL_MOVE_DEG 就中止。
6. 【結束時 movej 回 J0 = END_J0_DEG（180°）】不需再手動回原點。
7. 【IP 需要你自己再次確認】
8. 【姿態 POSE】程式不會移動 J1~J5，請先以教導器擺好姿態；連上 RTDE 後以實際角度檢查 J1~J5
   （±POSE_TOL_DEG），不符即中止，未知的姿態名稱也一律中止。EXT45 時法蘭離 J0 軸約 0.80 m，須輸入
   姿態名稱，確認 J0 掃過的扇形範圍（含工具）淨空；結束後 J0 回 180°，J1~J5 仍在外伸姿態，請以教導器移回。
9. 【負載 PAYLOAD_KG】填教導器上設定的負載質量（只有夾爪時 1.30 kg；加已知質量時為夾爪 + 質量，教導器的
   負載質量與重心也要一併改）。未設定一律中止；執行前須輸入教導器上的負載質量，與 PAYLOAD_KG 不符即中止。
10.【教導器速度滑桿】連上 RTDE 後讀 target_speed_fraction，不是 100% 或讀不到一律中止（10/22 冷機定速因
   滑桿停在 52% 而作廢）；執行中的最低值寫進摘要，低於 100% 時會警告這次資料不可用。
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

try:
    import ur5_home_pose as home_pose
except ImportError:
    home_pose = None


# ============================================================
# 使用者設定
# ============================================================

ROBOT_IP = "192.168.50.114"   # 執行前請再次確認
JOINT_INDEX = 0

DIRECTION = +1        # +1 = 正轉，-1 = 反轉。兩個方向分開執行，不要自動連續做
# QUICK: 第一個條件（1 Hz, ±0.06）× 1 次，確認正弦迴圈、角度保護與 movej 歸位
# MID  : 六個條件 × 各 1 次（含最大加速度 1.13 rad/s² 的 3 Hz），確認多段與較長腳本
# FULL : 六個條件 × 各 N_REPEAT 次，起點錯開、條件輪流
TEST_STAGE = "QUICK"

V_CRUISE = 0.12                  # 巡航速度 [rad/s]
# (頻率 Hz, 速度振幅 rad/s)。前四個為 10/22 上機的條件；後兩個：加速度振幅 0.19（低振幅）、
# 0.75（3 Hz，與 2 Hz × 0.06 同振幅）。分析以頻率與振幅辨認條件，10/22 的資料也可用這份清單分析。
CONDITIONS = [(1.0, 0.06), (2.0, 0.06), (3.0, 0.06), (2.0, 0.03), (1.0, 0.03), (3.0, 0.04)]
T_MOD = 5.0                      # 調變時間 [s]：六個條件都是整數個週期（5、10、15、10、5、15）
N_REPEAT = 7                     # FULL 每個條件的次數；第 k 次起點錯開 (k − 3)/7 個馬達圈

GEAR_RATIO = 101
WINDOW_CENTER_DEG = 180.0        # 調變段中心（絕對角度，同定速實驗窗口中心）
END_J0_DEG = WINDOW_CENTER_DEG   # 全部段跑完後 movej 回到的 J0 角度（= 診斷原點的 180°）
SETTLE_TIME = 0.8                # 到達巡航速度後、開始調變前的定速時間 [s]；調變從起步後約 6.3° 開始
RAMP_ACCEL = 0.5                 # 加速到巡航速度的 speedj 加速度 [rad/s^2]（同定速實驗）
MOD_ACCEL_LIMIT = 1.5            # 調變段 speedj 的加速度上限 [rad/s^2]，須大於最大調變加速度 2πf·A（1.13）
STOPJ_DECEL = 1.0                # stopj 減速度 [rad/s^2]（同定速實驗）
MOVEJ_ACCEL = 0.5                # 段間 movej 加速度 [rad/s^2]
MOVEJ_VEL = 0.15                 # 段間 movej 速度 [rad/s]
PAUSE_TIME = 0.5                 # 每段 stopj 後、movej 後的靜止時間 [s]（須小於結束判定的 1 s 靜止窗）
OVERRUN_DEG = 1.0                # 角度保護：越過調變終點此角度仍未結束即離開迴圈 [deg]
V_MIN_ALLOWED = 0.03             # 調變中最低速度的下限 [rad/s]：速度不歸零、不反向

QD_MAX = 0.30                    # 速度上限；J0 需低於約 0.32 rad/s，否則每圈 12 次的漣波在 125 Hz 下混疊
QDD_CHECK_LIMIT = 3.0            # 安全檢查用加速度上限
SAFETY_MARGIN_DEG = 2.0          # 絕對角度檢查時，起點與最遠點各再外推的餘裕 [deg]
MAX_INITIAL_MOVE_DEG = 40.0      # 目前 J0 到第一段起點的最大允許距離 [deg]

SAMPLE_HZ = 125.0
DT = 1.0 / SAMPLE_HZ
KT_OUT = 101 * 0.1350

# J0 允許的絕對角度活動視窗。必須由操作者依現場線纜與淨空狀況填寫，
# 未填寫則程式中止 —— 不提供預設值。
J0_SAFE_MIN_DEG = None
J0_SAFE_MAX_DEG = None

# 姿態：J1~J5 期望角度 [deg]，以及法蘭中心離 J0 軸的水平距離 [m]（UR5 官方 DH 計算，不含工具）。
# 程式不會移動 J1~J5，只檢查；DIAG 以外的姿態都須在現場確認淨空。
POSE = "DIAG"
POSES = {
    "DIAG": dict(q_deg=(-90.0, 90.0, -90.0, -90.0, 0.0), r_flange_m=0.50),    # 診斷原點（PHM 基準）
    "EXT45": dict(q_deg=(-45.0, 45.0, -90.0, -90.0, 0.0), r_flange_m=0.80),   # 手臂外伸，前臂水平
}
POSE_TOL_DEG = 0.5               # J1~J5 與期望角度的最大允許差 [deg]（同 ur5_home_pose）
FLANGE_SPEED_MAX = 0.25          # 法蘭水平線速度上限 [m/s]：(巡航 + 最大振幅) × 法蘭離 J0 軸距離

# 法蘭負載：填教導器上設定的負載質量 [kg]（只有夾爪時 1.30；加已知質量時為夾爪 + 質量）。
# 只用於確認與記錄（寫進輸出資料夾名稱），不影響動作；未設定則程式中止 —— 不提供預設值。
PAYLOAD_KG = None
PAYLOAD_TOL_KG = 0.005           # 執行前輸入的教導器負載質量與 PAYLOAD_KG 的最大允許差 [kg]

# 教導器速度滑桿（RTDE target_speed_fraction，0～1）的下限：滑桿會等比例縮小速度與加速度。
SPEED_FRACTION_MIN = 0.999

# 額外記錄的 RTDE 輸出（只記錄，不影響動作）：(欄名, ur_rtde 方法)。控制器或 ur_rtde 不支援時該欄留空。
EXTRA_SIGNALS = (("target_qdd", "getTargetQdd"), ("target_moment", "getTargetMoment"),
                 ("target_current", "getTargetCurrent"), ("joint_control_output", "getJointControlOutput"))

OUTPUT_DIR = "."


# ============================================================
# 第一部分：每段規劃與安全檢查
# ============================================================

def stage_config(stage):
    """回傳 (條件清單, 每個條件的次數)。"""
    if stage == "QUICK":
        return list(CONDITIONS[:1]), 1
    if stage == "MID":
        return list(CONDITIONS), 1
    if stage == "FULL":
        return list(CONDITIONS), N_REPEAT
    raise ValueError(f"未知的 TEST_STAGE: {stage}")


def start_offsets_deg(n_repeat, gear_ratio=GEAR_RATIO):
    """第 k 次重複的起點錯開量（k = 0..n_repeat−1）：(k − (n−1)/2)/n 個馬達圈，以 0 為中心。"""
    rev = 360.0 / gear_ratio
    return [(k - (n_repeat - 1) / 2.0) / n_repeat * rev for k in range(n_repeat)]


def plan_segments(conditions, n_repeat, direction, v_cruise, t_mod, settle_time, ramp_accel, decel,
                  overrun_deg, center_deg, gear_ratio, dt):
    """
    每一段的絕對角度規劃（單位：度）。條件依重複輪流：k = 0 的各條件 → k = 1 的各條件 → …
      start_deg     : movej 目標（調變起點往回退「加速距離 + 定速距離」）
      mod_start_deg / mod_end_deg : 調變段（分析只用這段）
      guard_deg     : 角度保護，越過即離開迴圈
      end_deg       : 正常結束（計數器到）時的停止點（以巡航速度煞車）
      extreme_deg   : 最壞情況的最遠點（角度保護觸發時以最高速度 V_CRUISE + A 煞車）
      n_pre / n_mod : 加速＋定速、調變兩個迴圈的次數（每次一個控制週期）
    """
    d = direction
    mod_len = math.degrees(v_cruise * t_mod)                       # 整數週期：正弦的位移積分為 0
    pre_deg = math.degrees(v_cruise ** 2 / (2.0 * ramp_accel) + v_cruise * settle_time)
    n_pre = int(round((v_cruise / ramp_accel + settle_time) / dt))
    n_mod = int(round(t_mod / dt))
    segs = []
    for k, off in enumerate(start_offsets_deg(n_repeat, gear_ratio)):
        for ci, (f, amp) in enumerate(conditions):
            mod_start = center_deg - d * mod_len / 2.0 + d * off
            mod_end = mod_start + d * mod_len
            guard = mod_end + d * overrun_deg
            segs.append(dict(
                index=len(segs), cond=ci, freq=f, amp=amp, repeat=k, offset_deg=off,
                start_deg=mod_start - d * pre_deg, mod_start_deg=mod_start, mod_end_deg=mod_end,
                guard_deg=guard,
                end_deg=mod_end + d * math.degrees(v_cruise ** 2 / (2.0 * decel)),
                extreme_deg=guard + d * math.degrees((v_cruise + amp) ** 2 / (2.0 * decel)),
                n_pre=n_pre, n_mod=n_mod, w_dt=2.0 * math.pi * f * dt,
                t_active=(n_pre + n_mod) * dt))
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
        total += PAUSE_TIME + s["t_active"] + V_CRUISE / decel + PAUSE_TIME
        pos = s["end_deg"]
    if end_deg is not None:                                   # 結束後回原點
        total += estimate_movej_time(math.radians(abs(end_deg - pos)), MOVEJ_VEL, MOVEJ_ACCEL) + PAUSE_TIME
    return total


def safety_check(segs, v_cruise, t_mod, qd_max, qdd_max, v_min_allowed, ramp_accel, mod_accel, decel,
                 movej_vel, movej_accel, safe_min_deg, safe_max_deg, margin_deg, end_deg=None):
    """
    回傳 (ok, lines, (planned_min_deg, planned_max_deg))。
    速度、加速度上限；調變中最低速度（不歸零、不反向）；speedj 加速度上限須大於調變所需（否則正弦被削平）；
    每個條件在 T_MOD 內為整數個週期；絕對角度：所有段的起點、最壞情況最遠點、回原點角度，
    各外推 margin_deg，都必須在 J0 視窗內。J0 視窗未設定（None）一律 FAIL（預設拒絕）。
    """
    lines, ok = [], True

    def chk(name, val, lim, unit):
        nonlocal ok
        lines.append(f"{name} : {val:.4f} {unit}  (限 {lim})")
        if val > lim:
            lines.append("  [FAIL]"); ok = False
        else:
            lines.append("  [OK]")

    amp_max = max(s["amp"] for s in segs)
    acc_mod = max(2.0 * math.pi * s["freq"] * s["amp"] for s in segs)
    chk("最高速度（巡航 + 振幅）", v_cruise + amp_max, qd_max, "rad/s")
    v_low = v_cruise - amp_max
    lines.append(f"最低速度（巡航 − 振幅） : {v_low:.4f} rad/s  (下限 {v_min_allowed})")
    if v_low < v_min_allowed:
        lines.append("  [FAIL] 速度可能過低、歸零或反向"); ok = False
    else:
        lines.append("  [OK]")
    chk("最大調變加速度 2πf·A", acc_mod, mod_accel, "rad/s^2")
    chk("調變段 speedj 加速度上限", mod_accel, qdd_max, "rad/s^2")
    chk("加速到巡航的 speedj 加速度", ramp_accel, qdd_max, "rad/s^2")
    chk("stopj 減速度", decel, qdd_max, "rad/s^2")
    chk("movej 速度", movej_vel, qd_max, "rad/s")
    chk("movej 加速度", movej_accel, qdd_max, "rad/s^2")
    for f in sorted({s["freq"] for s in segs}):
        n_cyc = f * t_mod
        lines.append(f"{f:.2f} Hz 在 {t_mod:.2f} s 內的週期數 : {n_cyc:.4f}")
        if abs(n_cyc - round(n_cyc)) > 1e-6:
            lines.append("  [FAIL] 不是整數個週期（調變結束時速度不等於巡航速度）"); ok = False
        else:
            lines.append("  [OK]")

    pts = [s["start_deg"] for s in segs] + [s["extreme_deg"] for s in segs]
    if end_deg is not None:
        pts.append(end_deg)
    lo_plan, hi_plan = min(pts) - margin_deg, max(pts) + margin_deg
    lines.append(f"規劃絕對角度範圍（含餘裕 {margin_deg:.1f}°） : {lo_plan:.2f}° 至 {hi_plan:.2f}°")
    if safe_min_deg is None or safe_max_deg is None:
        lines.append("  [FAIL] J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 尚未設定，拒絕執行")
        ok = False
    else:
        lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
        if lo <= lo_plan and hi_plan <= hi:
            lines.append(f"  [OK] 在 J0 視窗 [{lo:.1f}°, {hi:.1f}°] 內")
        else:
            lines.append(f"  [FAIL] 超出 J0 視窗 [{lo:.1f}°, {hi:.1f}°]")
            ok = False
    return ok, lines, (lo_plan, hi_plan)


def check_start_position(q0_deg, first_start_deg, max_move_deg, safe_min_deg, safe_max_deg):
    """手臂目前位置：必須在 J0 視窗內，且離第一段起點不超過 max_move_deg。"""
    if safe_min_deg is None or safe_max_deg is None:
        return False, "J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 尚未設定，拒絕執行"
    lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
    move = abs(first_start_deg - q0_deg)
    msg = (f"J0 目前 {q0_deg:.2f}°，第一段起點 {first_start_deg:.2f}°（需移動 {move:.2f}°，"
           f"上限 {max_move_deg:.1f}°），J0 視窗 [{lo:.1f}°, {hi:.1f}°]")
    if not (lo <= q0_deg <= hi):
        return False, msg + " -> 目前角度在視窗外"
    if move > max_move_deg:
        return False, msg + " -> 第一次 movej 距離過大"
    return True, msg + " -> OK"


def check_pose(q_full_rad, pose, poses=POSES, tol_deg=POSE_TOL_DEG):
    """
    以 RTDE 量到的六軸角度 [rad] 檢查 J1~J5 是否在 poses[pose] 的期望角度 ± tol_deg 內
    （差值以 ±360° 等價處理，同 ur5_home_pose；J0 由 J0 視窗與起始位置檢查負責）。
    預設拒絕：pose 不在 poses、讀不到六軸角度、角度不是有限值、任一軸超差，一律回傳 False。
    回傳 (ok, lines)。
    """
    if pose not in poses:
        return False, [f"未知的姿態 POSE = {pose!r}（可用：{', '.join(poses)}），拒絕執行"]
    if q_full_rad is None or len(q_full_rad) != 6:
        return False, ["讀不到六軸角度，拒絕執行"]
    lines, ok = [f"姿態 {pose}：J1~J5 期望 {list(poses[pose]['q_deg'])}°，容差 ±{tol_deg}°"], True
    for j, tgt in zip(range(1, 6), poses[pose]["q_deg"]):
        cur = math.degrees(q_full_rad[j])
        d = (cur - tgt + 180.0) % 360.0 - 180.0
        passed = abs(d) <= tol_deg                                 # NaN 時為 False
        ok = ok and passed
        lines.append(f"  J{j}：目前 {cur:8.2f}°，期望 {tgt:8.2f}°，差 {d:+7.2f}°  {'[OK]' if passed else '[FAIL]'}")
    return ok, lines


def check_flange_speed(pose, conditions, v_cruise, poses=POSES, v_max=FLANGE_SPEED_MAX):
    """法蘭最高水平線速度 = (巡航 + 最大振幅) × 法蘭離 J0 軸距離，須 ≤ v_max。pose 未知時拒絕。回傳 (ok, 說明)。"""
    if pose not in poses:
        return False, f"未知的姿態 POSE = {pose!r}（可用：{', '.join(poses)}），拒絕執行"
    r = poses[pose]["r_flange_m"]
    v = (v_cruise + max(amp for _, amp in conditions)) * r
    return v <= v_max, (f"法蘭最高水平線速度（姿態 {pose}，法蘭離 J0 軸 {r:.2f} m）: {v:.4f} m/s  (限 {v_max})"
                        + ("  [OK]" if v <= v_max else "  [FAIL]"))


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
    """
    操作者輸入的教導器負載質量 entered（字串）須與 PAYLOAD_KG 相差 ≤ tol_kg。
    預設拒絕：PAYLOAD_KG 不合格、輸入讀不懂、相差超過 tol_kg，一律回傳 False。回傳 (ok, 說明)。
    """
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
    """
    教導器速度滑桿（RTDE target_speed_fraction，0～1）：滑桿會等比例縮小 URScript 的速度與加速度
    （10/22 冷機定速因滑桿停在 52% 而作廢）。預設拒絕：讀不到（None、ur_rtde 沒有此方法）、
    不是有限數值、低於 min_fraction，一律回傳 False。回傳 (ok, 說明)。
    """
    x = _finite_float(fraction)
    if x is None:
        return False, f"讀不到教導器速度滑桿（target_speed_fraction = {fraction!r}），拒絕執行"
    msg = f"教導器速度滑桿 {100 * x:.1f}%（須 ≥ {100 * min_fraction:.1f}%）"
    if x < min_fraction:
        return False, msg + " -> 請把滑桿調回 100% 後重新執行"
    return True, msg + " -> OK"


# ============================================================
# 第二部分：URScript 產生
# ============================================================

def build_urscript(segs, direction, joint_index, q_others, v_cruise, ramp_accel, mod_accel, decel,
                   movej_accel, movej_vel, pause_time, dt, end_deg=None):
    """
    q_others：六軸絕對角度 [rad]（取自執行開始時 RTDE 量到的值），只替換 J0 當作各段 movej 目標，
    其他五軸維持原值不動。
    每段：movej → 停 →
          迴圈 1（n_pre 次）：speedj(巡航速度, RAMP_ACCEL) —— 加速到巡航後定速
          迴圈 2（n_mod 次）：qd_t = 巡航速度 + A·sin(ω·dt·n)，speedj(qd_t, MOD_ACCEL_LIMIT)
          兩個迴圈每週期都讀 J0 實際角度，越過 guard 角度即離開 → stopj → 停。
    end_deg 不為 None 時，全部段結束後再 movej 回 J0 = end_deg（其他五軸不變）。
    """
    lines = ["def velocity_modulation_program():"]
    cmp_op = "<" if direction > 0 else ">"
    pm = "+" if direction > 0 else "-"
    qd_c = ["0.0"] * 6
    qd_c[joint_index] = f"{direction * v_cruise:.6f}"
    qd_c_str = "[" + ", ".join(qd_c) + "]"
    qd_m = ["0.0"] * 6
    qd_m[joint_index] = "qd_t"
    qd_m_str = "[" + ", ".join(qd_m) + "]"
    for s in segs:
        q_t = list(q_others)
        q_t[joint_index] = math.radians(s["start_deg"])
        q_t_str = "[" + ", ".join(f"{v:.6f}" for v in q_t) + "]"
        guard = math.radians(s["guard_deg"])
        cond = f"(q_now[{joint_index}] {cmp_op} {guard:.6f})"
        lines.append(f"  # segment {s['index'] + 1}: f = {s['freq']:.2f} Hz, A = {s['amp']:.3f} rad/s, "
                     f"offset {s['offset_deg']:+.4f} deg")      # URScript 內只放 ASCII
        lines.append(f"  movej({q_t_str}, a={movej_accel}, v={movej_vel})")
        lines.append(f"  sleep({pause_time})")
        lines.append("  q_now = get_actual_joint_positions()")
        lines.append("  n = 0")
        lines.append(f"  while (n < {s['n_pre']}) and {cond}:")
        lines.append(f"    speedj({qd_c_str}, a={ramp_accel}, t={dt})")
        lines.append("    q_now = get_actual_joint_positions()")
        lines.append("    n = n + 1")
        lines.append("  end")
        lines.append("  n = 0")
        lines.append(f"  qd_t = {direction * v_cruise:.6f}")
        lines.append(f"  while (n < {s['n_mod']}) and {cond}:")
        lines.append(f"    qd_t = {direction * v_cruise:.6f} {pm} {s['amp']:.6f} * sin({s['w_dt']:.9f} * n)")
        lines.append(f"    speedj({qd_m_str}, a={mod_accel}, t={dt})")
        lines.append("    q_now = get_actual_joint_positions()")
        lines.append("    n = n + 1")
        lines.append("  end")
        lines.append(f"  stopj({decel})")
        lines.append(f"  sleep({pause_time})")
    if end_deg is not None:
        q_e = list(q_others)
        q_e[joint_index] = math.radians(end_deg)
        lines.append("  # return to home J0")
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
    """同 ur5_const_velocity_angle_ident.py：J0 實際速度連續 quiet_s 秒低於門檻才視為結束。"""
    t_start = time.time()
    hard_timeout = planned_s * timeout_factor + 5.0
    moved, quiet_since, t_first, t_last = False, None, None, None
    while True:
        now = time.time()
        elapsed = now - t_start
        qd0 = logger.last_qd0
        if qd0 is not None and abs(qd0) > qd_eps:
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
# 第三部分：RTDE 記錄（同 ur5_const_velocity_angle_ident.py，另記錄 EXTRA_SIGNALS）
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
        self.extra_reads = {name: 0 for name, _ in EXTRA_SIGNALS}    # 各額外訊號成功讀取的列數
        self.extra_missing = []                                      # ur_rtde 沒有對應方法的訊號
        self.last_speed_fraction = self.min_speed_fraction = None    # 教導器速度滑桿（target_speed_fraction）
        self.speed_fraction_missing = 0                              # 讀不到滑桿的列數
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
                for name, getter in extras:                    # 讀不到時該列留空，不影響主要欄位
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

def _runs(mask):
    idx = np.where(mask)[0]
    return np.split(idx, np.where(np.diff(idx) != 1)[0] + 1) if len(idx) else []


def _fit_sine(t, y, f, tc):
    """固定頻率 f：y ≈ c0 + c1·cos(ω(t − tc)) + c2·sin(ω(t − tc))。回傳 (係數, 殘差 RMS)。"""
    x = t - tc
    w = 2.0 * math.pi * f
    A = np.column_stack([np.ones_like(x), np.cos(w * x), np.sin(w * x)])
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    return c, float(np.sqrt(np.mean((y - A @ c) ** 2)))


def _refine_freq(t, y, f0, tc, rel=0.1, n_grid=201):
    """在 f0·(1 ± rel) 內找殘差最小的頻率（格點 + 對殘差平方做拋物線內插）。"""
    fs = np.linspace(f0 * (1 - rel), f0 * (1 + rel), n_grid)
    e2 = np.array([_fit_sine(t, y, f, tc)[1] ** 2 for f in fs])
    i = int(np.argmin(e2))
    if 0 < i < n_grid - 1:
        den = e2[i - 1] - 2 * e2[i] + e2[i + 1]
        if den > 0:
            return float(fs[i] + 0.5 * (e2[i - 1] - e2[i + 1]) / den * (fs[1] - fs[0]))
    return float(fs[i])


def find_modulation_segments(t, tq, v_cruise, conditions, dt, skip_s=0.06, trim_s=0.06, tol=3e-3):
    """
    tq：運動方向的目標速度。每次「目標速度 > 巡航速度的 1/10」的連續區段中（調變時最低速度
    V_CRUISE − A 不會低於此門檻，一段不會被切開；movej 回程方向相反，不會被選入）：
      調變起點 = 到達巡航速度後，目標速度第一次偏離巡航速度的前一個樣本；
      以起點後 1 s 擬合各條件的頻率與振幅，取吻合者（頻率差 < 10%、振幅差 < 25%）；
      調變終點 = 目標速度第一次偏離擬合正弦超過 tol 的樣本（stopj 開始，或角度保護觸發）。
    分析區段去掉開頭 skip_s、結尾 trim_s 秒。回傳 (段清單, 無法配對的區段數)。
    """
    segs, skipped = [], 0
    n_skip, n_trim, n_1s = int(round(skip_s / dt)), int(round(trim_s / dt)), int(round(1.0 / dt))
    for r in _runs(tq > 0.1 * v_cruise):
        x, tr = tq[r], t[r]
        cr = np.flatnonzero(np.abs(x - v_cruise) < 1e-6)
        if len(cr) == 0:
            continue
        dv = np.flatnonzero(np.abs(x[cr[0]:] - v_cruise) > 1e-4)
        if len(dv) == 0:
            continue                                            # 只有定速、沒有調變
        i_s = cr[0] + dv[0] - 1
        j0 = i_s + 1                                            # 擬合從調變後第一個樣本起，不受起點轉折影響
        n1 = min(len(r) - j0, n_1s)
        best = None
        for ci, (f, amp) in enumerate(conditions):
            fr = _refine_freq(tr[j0:j0 + n1], x[j0:j0 + n1], f, tr[i_s])
            c, rms = _fit_sine(tr[j0:j0 + n1], x[j0:j0 + n1], fr, tr[i_s])
            a_fit = math.hypot(c[1], c[2])
            if abs(fr / f - 1) < 0.1 and abs(a_fit / amp - 1) < 0.25 and (best is None or rms < best[2]):
                best = (ci, fr, rms, c)
        if best is None:
            skipped += 1
            continue
        ci, fr, _, c = best
        w = 2.0 * math.pi * fr
        xx = tr[j0:] - tr[i_s]
        model = c[0] + c[1] * np.cos(w * xx) + c[2] * np.sin(w * xx)
        bad = np.flatnonzero(np.abs(x[j0:] - model) > tol)
        i_e = j0 + (int(bad[0]) if len(bad) else len(r) - j0)
        fr = _refine_freq(tr[j0:i_e], x[j0:i_e], fr, tr[i_s], rel=0.01)
        c, rms = _fit_sine(tr[j0:i_e], x[j0:i_e], fr, tr[i_s])
        k0, k1 = i_s + n_skip, i_e - n_trim
        if k1 - k0 < n_1s:
            skipped += 1
            continue
        segs.append(dict(cond=ci, idx=r[k0:k1], run_start=int(r[0]), freq=fr, amp=math.hypot(c[1], c[2]),
                         n_cycles=(tr[i_e - 1] - tr[i_s]) * fr, target_rms=rms))
    return segs, skipped


def _fit_kinematics(t, y, f, n_harm=3, poly=3):
    """
    以「三次多項式 + 調變頻率的 1～n_harm 次諧波」擬合角度 y(t)，解析微分得到速度、加速度。
    角度量測雜訊約 2e-5 rad，直接數值二次微分的雜訊會使回歸的 J 偏低（迴歸變數含誤差）；
    擬合後的加速度幾乎沒有雜訊。回傳 (v, a, 調變頻率的速度相量, 加速度振幅, 擬合殘差 RMS)。
    """
    x = t - 0.5 * (t[0] + t[-1])
    w = 2.0 * math.pi * f
    z = np.zeros_like(x)
    cols = [x ** p for p in range(poly + 1)]
    d1 = [p * x ** (p - 1) if p >= 1 else z for p in range(poly + 1)]
    d2 = [p * (p - 1) * x ** (p - 2) if p >= 2 else z for p in range(poly + 1)]
    for k in range(1, n_harm + 1):
        wk = k * w
        cs, sn = np.cos(wk * x), np.sin(wk * x)
        cols += [cs, sn]
        d1 += [-wk * sn, wk * cs]
        d2 += [-wk ** 2 * cs, -wk ** 2 * sn]
    A = np.column_stack(cols)
    c, *_ = np.linalg.lstsq(A, y, rcond=None)
    cc, cs1 = c[poly + 1], c[poly + 2]                           # 角度的 cos、sin 係數（基頻）
    v_ph = w * cs1 + 1j * w * cc                                  # 速度相量（y = Re(Y·e^{jωx})）
    return (np.column_stack(d1) @ c, np.column_stack(d2) @ c, v_ph, w ** 2 * math.hypot(cc, cs1),
            float(np.sqrt(np.mean((y - A @ c) ** 2))))


def _hat_basis(x, knots):
    """分段線性（帽形）基底，去掉第一個節點：各段常數已涵蓋常數項。"""
    x = np.clip(x, knots[0], knots[-1])
    Bm = np.zeros((len(x), len(knots)))
    for j in range(len(knots)):
        if j > 0:
            m = (x >= knots[j - 1]) & (x <= knots[j])
            Bm[m, j] = (x[m] - knots[j - 1]) / (knots[j] - knots[j - 1])
        if j < len(knots) - 1:
            m = (x >= knots[j]) & (x < knots[j + 1])
            Bm[m, j] = (knots[j + 1] - x[m]) / (knots[j + 1] - knots[j])
    return Bm[:, 1:]


def _ripple_names(k_ripple, only=None):
    ks = range(1, k_ripple + 1) if only is None else only
    return [n for k in ks for n in (f"c{k}", f"s{k}")]


def _solve_normal(XtX, Xty, sel):
    """XtX、Xty：全部欄位的 XᵀX、Xᵀy；sel：使用的欄位索引。"""
    sel = np.asarray(sel)
    return np.linalg.solve(XtX[np.ix_(sel, sel)], Xty[sel])


def _hac_cov(Xs, es, XtX, lag):
    """係數的共變異數，考慮段內殘差自相關（Bartlett 權重、最大延遲 lag 個樣本，段與段之間視為獨立）：
    Var = (XᵀX)⁻¹ S (XᵀX)⁻¹，S = Σ_段 Σ_|l|≤lag (1 − |l|/(lag+1)) Σ_t z_t z_{t−l}ᵀ，z_t = x_t·e_t。"""
    S = np.zeros_like(XtX)
    for X, e in zip(Xs, es):
        Z = X * e[:, None]
        S += Z.T @ Z
        for l in range(1, min(lag, len(e) - 1) + 1):
            G = Z[l:].T @ Z[:-l]
            S += (1.0 - l / (lag + 1.0)) * (G + G.T)
    Ai = np.linalg.inv(XtX)
    return Ai @ S @ Ai


def _wls(X, y, se):
    W = 1.0 / np.asarray(se) ** 2
    cov = np.linalg.inv((X.T * W) @ X)
    return cov @ ((X.T * W) @ y), np.sqrt(np.diag(cov))


def _fit_target_moment(cs, tm, tqdd, tq):
    """
    控制器模型（運動方向）：target_moment = J_model·target_qdd + B_model·target_qd + 各段常數。
    target_moment 由控制器以目標軌跡算出、沒有量測雜訊，不另估標準差；改列殘差 RMS，
    殘差大表示 target_moment 還含有此式以外的成分。有效樣本不足時回傳 J_model = None。
    """
    Xs, ys = [], []
    for g, sg in enumerate(cs):
        ii = sg["idx"][np.isfinite(tm[sg["idx"]]) & np.isfinite(tqdd[sg["idx"]])]
        if len(ii) < 10:
            continue
        D = np.zeros((len(ii), len(cs)))
        D[:, g] = 1.0
        Xs.append(np.column_stack([tqdd[ii], tq[ii], D]))
        ys.append(tm[ii])
    if not Xs:
        return dict(J_model=None, B_model=None, model_resid_rms=None, model_n=0)
    X, y = np.vstack(Xs), np.concatenate(ys)
    X = X[:, np.r_[True, True, np.any(X[:, 2:] != 0, axis=0)]]   # 去掉沒有資料的段
    p, *_ = np.linalg.lstsq(X, y, rcond=None)
    return dict(J_model=float(p[0]), B_model=float(p[1]),
                model_resid_rms=float(np.sqrt(np.mean((y - X @ p) ** 2))), model_n=int(len(y)))


def analyze_velocity_modulation(csv_path, joint_index, kt_out, direction, conditions, v_cruise,
                                gear_ratio=GEAR_RATIO, dt=DT, k_ripple=24, n_knots=9, skip_s=0.06, trim_s=0.06,
                                hac_s=0.25):
    """
    每段調變區段（見 find_modulation_segments）：角度以多項式 + 諧波擬合得到速度 v、加速度 a。
    每個條件分開回歸（運動方向的力矩）：
        tau = J·a + B·v + C·v² + 漣波(馬達角 1～k_ripple 次) + 各段常數 + 各段行程斜率
    各段常數與斜率吸收溫度、起步後行程等慢變化；調變段內的振盪成分決定 J 與 B_local = B + 2C·v̄。
    ± 為考慮段內殘差自相關（最大延遲 hac_s 秒）的回歸係數標準差。每個條件只有 7 段，以段重抽樣的
    標準差不穩定（合成資料中低估 20～25%）；此估計在合成資料中與重複產生資料的實際散布相符（1 Hz 偏大）。
    並以不同模型寫法重算 J（無漣波項、只取主要諧波、無 v²、行程項改為各段共用曲線、無行程項），
    檢查 J 是否隨寫法改變。
    注意：調變時行程也跟著振盪（1 Hz、±0.06 rad/s 時 ±0.55°，與加速度反相），起步後行程效應的彎曲部分
    會被當成慣性：合成資料（10/15 實測的行程曲線）中 1 Hz 的 J 偏低約 1%，2、3 Hz 可忽略。
    改用共用曲線偏差較小，但 1 Hz 的 J 標準差約加倍（曲線越細越不穩定），因此只列為對照。
    跨條件：
        J 對 1/A 回歸（J = J∞ + k/A）：k ≈ 0 表示 J 不隨振幅變化；減速時摩擦多一個常數 c 時 k = −2c/π。
        B_local 對 ω² 回歸：電流相對編碼器延遲 Δ 時，斜率 = J·Δ（也可能來自摩擦的動態特性）。
    控制器模型（CSV 有 target_moment、target_qdd 欄位時；每個條件只要 1 段即可，QUICK、MID 也算）：
        target_moment 對 target_qdd、target_qd 回歸（見 _fit_target_moment），J_model 為控制器模型在此姿態的
        J0 慣量。J ÷ J_model 不是 1 時，可能是 Kt 不準，也可能是控制器模型與實際不同，須比較多個姿態才能分開。
    """
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    def col(key):
        if key not in rows[0]:
            return np.full(len(rows), np.nan)
        return np.array([float(r[key]) if r.get(key) not in (None, "", "None") else np.nan for r in rows])

    s = float(direction)
    t = np.array([float(r["timestamp"]) for r in rows])
    q = np.array([float(r[f"actual_q_{joint_index}"]) for r in rows])
    tq = s * np.array([float(r[f"target_qd_{joint_index}"]) for r in rows])
    tau = s * np.array([float(r[f"actual_current_{joint_index}"]) for r in rows]) * kt_out
    temp = col(f"joint_temp_{joint_index}")
    tm, tqdd = s * col(f"target_moment_{joint_index}"), s * col(f"target_qdd_{joint_index}")

    segs, skipped = find_modulation_segments(t, tq, v_cruise, conditions, dt, skip_s, trim_s)
    if not segs:
        print("[警告] 沒有找到調變區段。")
        return dict(segments=[], skipped=skipped, conditions={}, model={})

    for sg in segs:
        ii = sg["idx"]
        v, a, v_ph, a1, res = _fit_kinematics(t[ii], s * q[ii], sg["freq"])
        c_t, _ = _fit_sine(t[ii], tq[ii], sg["freq"], 0.5 * (t[ii[0]] + t[ii[-1]]))
        t_ph = c_t[1] - 1j * c_t[2]                               # 目標速度相量
        h = v_ph / t_ph
        sg.update(v=v, a=a, a_amp=a1, q_fit_rms=res, track_gain=abs(h),
                  track_lag_ms=-np.angle(h) / (2 * math.pi * sg["freq"]) * 1000.0,
                  temp=float(np.nanmean(temp[ii])) if np.any(np.isfinite(temp[ii])) else float("nan"))

    rip_all = _ripple_names(k_ripple)
    rip_main = _ripple_names(k_ripple, only=[k for k in (2, 4, 5, 10, 12) if k <= k_ripple])
    res_c = {}
    for ci, (f, amp) in enumerate(conditions):
        cs = [sg for sg in segs if sg["cond"] == ci]
        if len(cs) < 2:
            continue
        seg_n = [f"seg{g}" for g in range(len(cs))]
        trv_n = [f"tr{g}" for g in range(len(cs))]
        trav = [s * np.degrees(q[sg["idx"]] - q[sg["run_start"]]) for sg in cs]      # 起步後行程 [deg]
        knots = np.linspace(min(x.min() for x in trav), max(x.max() for x in trav), n_knots)
        g_n = [f"g{j}" for j in range(1, n_knots)]
        cols_all = ["J", "B", "C"] + rip_all + seg_n + trv_n + g_n
        Xs, ys = [], []
        for g, sg in enumerate(cs):
            ii = sg["idx"]
            phi = gear_ratio * q[ii]
            C = {"J": sg["a"], "B": sg["v"], "C": sg["v"] ** 2}
            for k in range(1, k_ripple + 1):
                C[f"c{k}"], C[f"s{k}"] = np.cos(k * phi), np.sin(k * phi)
            one, trv, zero = np.ones(len(ii)), s * (q[ii] - q[ii[0]]), np.zeros(len(ii))
            for h in range(len(cs)):
                C[f"seg{h}"] = one if h == g else zero
                C[f"tr{h}"] = trv if h == g else zero
            Hb = _hat_basis(trav[g], knots)
            for j, n in enumerate(g_n):
                C[n] = Hb[:, j]
            Xs.append(np.column_stack([C[n] for n in cols_all]))
            ys.append(tau[ii])
        XtX = sum(X.T @ X for X in Xs)
        Xty = sum(X.T @ y for X, y in zip(Xs, ys))
        pos = {n: i for i, n in enumerate(cols_all)}
        full = ["J", "B", "C"] + rip_all + seg_n + trv_n
        variants = {
            "完整": full,
            "無漣波項": ["J", "B", "C"] + seg_n + trv_n,
            "漣波只取 2、4、5、10、12 次": ["J", "B", "C"] + rip_main + seg_n + trv_n,
            "無 v² 項": ["J", "B"] + rip_all + seg_n + trv_n,
            f"行程項改為共用曲線（{n_knots} 節點）": ["J", "B", "C"] + rip_all + seg_n + g_n,
            "無行程項（各段只用常數）": ["J", "B", "C"] + rip_all + seg_n,
        }
        lag = int(round(hac_s / dt))

        def fit(names):
            si = [pos[n] for n in names]
            pp = _solve_normal(XtX, Xty, si)
            ee = [y - X[:, si] @ pp for X, y in zip(Xs, ys)]
            return pp, ee, _hac_cov([X[:, si] for X in Xs], ee, XtX[np.ix_(si, si)], lag)

        v_bar = float(np.mean(np.concatenate([sg["v"] for sg in cs])))
        p, es, cov = fit(full)
        J, B_loc = p[0], p[1] + 2 * p[2] * v_bar
        e = np.concatenate(es)
        J_var = {}
        for k, names in variants.items():
            pv, _, cv = (p, es, cov) if k == "完整" else fit(names)
            J_var[k] = (float(pv[0]), float(np.sqrt(cv[0, 0])))
        res_c[ci] = dict(
            freq=f, amp=amp, n_seg=len(cs), J=float(J), se_J=float(np.sqrt(cov[0, 0])),
            B_local=float(B_loc), v_bar=v_bar,
            se_B_local=float(np.sqrt(cov[1, 1] + 4 * v_bar ** 2 * cov[2, 2] + 4 * v_bar * cov[1, 2])),
            a_amp=float(np.mean([sg["a_amp"] for sg in cs])),
            freq_meas=float(np.mean([sg["freq"] for sg in cs])),
            track_gain=float(np.median([sg["track_gain"] for sg in cs])),
            track_lag_ms=float(np.median([sg["track_lag_ms"] for sg in cs])),
            temp=float(np.nanmean([sg["temp"] for sg in cs])) if any(np.isfinite(sg["temp"]) for sg in cs) else float("nan"),
            resid_rms=float(np.sqrt(np.mean(e ** 2))),
            J_variants={k: v[0] for k, v in J_var.items()}, se_J_variants={k: v[1] for k, v in J_var.items()})
        res_c[ci]["J_spread"] = float(np.ptp(list(res_c[ci]["J_variants"].values())))

    model = {}
    for ci in sorted({sg["cond"] for sg in segs}):
        mr = _fit_target_moment([sg for sg in segs if sg["cond"] == ci], tm, tqdd, tq)
        if mr["J_model"] is not None:
            model[ci] = mr
    out = dict(segments=segs, skipped=skipped, conditions=res_c, model=model)
    if model:
        out["J_model_mean"] = float(np.mean([mr["J_model"] for mr in model.values()]))

    def print_model():
        if not model:
            return
        print("控制器模型：target_moment = J_model·target_qdd + B_model·target_qd + 各段常數（運動方向，目標值無量測雜訊）")
        for ci, mr in model.items():
            ratio = (f"；J ÷ J_model = {res_c[ci]['J'] / mr['J_model']:.3f}"
                     if ci in res_c and abs(mr["J_model"]) > 1e-9 else "")
            print(f"  條件 {ci + 1}：J_model = {mr['J_model']:.4f}，B_model = {mr['B_model']:.4f}，"
                  f"殘差 RMS {mr['model_resid_rms']:.4f} N·m（{mr['model_n']} 點）{ratio}")

    used = sorted(res_c)
    print(f"調變區段：{len(segs)} 段（無法配對 {skipped} 段）；每段去掉開頭 {skip_s * 1000:.0f} ms、"
          f"結尾 {trim_s * 1000:.0f} ms；角度擬合殘差中位數 {np.median([sg['q_fit_rms'] for sg in segs]):.2e} rad")
    if not used:
        print("[提示] 每個條件至少要 2 段才做回歸：QUICK、MID 只用來確認動作，J 請用 FULL 的資料。各段：")
        for g, sg in enumerate(segs):
            print(f"  第 {g + 1} 段：條件 {sg['cond'] + 1}，f = {sg['freq']:.3f} Hz，A = {sg['amp']:.4f} rad/s，"
                  f"{sg['n_cycles']:.2f} 個週期，實測加速度振幅 {sg['a_amp']:.3f} rad/s²，"
                  f"追隨增益 {sg['track_gain']:.3f}，延遲 {sg['track_lag_ms']:.2f} ms")
        print_model()
        return out
    print(f"{'條件':>4}{'f[Hz]':>7}{'A':>7}{'段數':>5}{'實測a振幅':>10}{'J':>9}{'±':>7}{'B_local':>9}{'±':>7}"
          f"{'追隨增益':>9}{'延遲ms':>8}{'溫度':>7}")
    for ci in used:
        r = res_c[ci]
        print(f"{ci + 1:4d}{r['freq_meas']:7.3f}{r['amp']:7.3f}{r['n_seg']:5d}{r['a_amp']:10.3f}{r['J']:9.4f}"
              f"{r['se_J']:7.4f}{r['B_local']:9.3f}{r['se_B_local']:7.3f}{r['track_gain']:9.3f}"
              f"{r['track_lag_ms']:8.2f}{r['temp']:7.2f}")
    print("J 隨模型寫法（各條件，J ± 標準差；最後為最大差）：")
    for ci in used:
        r = res_c[ci]
        print(f"  條件 {ci + 1}：" + "，".join(f"{k} {v:.4f}±{r['se_J_variants'][k]:.4f}"
                                            for k, v in r["J_variants"].items())
              + f"；最大差 {r['J_spread']:.4f}")
    print_model()
    if len(used) >= 2:
        Jc = np.array([res_c[c]["J"] for c in used]); sJ = np.array([res_c[c]["se_J"] for c in used])
        p_m, s_m = _wls(np.ones((len(used), 1)), Jc, sJ)
        out["J_mean"], out["se_J_mean"] = float(p_m[0]), float(s_m[0])
        chi2 = float(np.sum(((Jc - p_m[0]) / sJ) ** 2))
        out["J_chi2"], out["J_dof"] = chi2, len(used) - 1
        print(f"J（各條件加權平均）= {p_m[0]:.4f} ± {s_m[0]:.4f}；各條件與平均差的卡方 = {chi2:.1f}"
              f"（自由度 {len(used) - 1}，遠大於自由度表示各條件的 J 不一致）")
        A = np.array([res_c[c]["a_amp"] for c in used])
        if np.ptp(A) > 0.1 * np.max(A):
            p_a, s_a = _wls(np.column_stack([np.ones(len(used)), 1.0 / A]), Jc, sJ)
            out.update(J_inf=float(p_a[0]), se_J_inf=float(s_a[0]), k_invA=float(p_a[1]), se_k_invA=float(s_a[1]))
            print(f"J 對 1/A 回歸：J∞ = {p_a[0]:.4f} ± {s_a[0]:.4f}，k = {p_a[1]:+.4f} ± {s_a[1]:.4f} N·m"
                  f"（k ≈ 0：J 不隨振幅變化；若 k 顯著為負，相當於減速時摩擦多 {-math.pi * p_a[1] / 2:+.3f} N·m）")
        wv = np.array([2 * math.pi * res_c[c]["freq_meas"] for c in used])
        if np.ptp(wv) > 0.1 * np.max(wv):
            Bc = np.array([res_c[c]["B_local"] for c in used]); sB = np.array([res_c[c]["se_B_local"] for c in used])
            p_b, s_b = _wls(np.column_stack([np.ones(len(used)), wv ** 2]), Bc, sB)
            out.update(B_slope_w2=float(p_b[1]), se_B_slope_w2=float(s_b[1]),
                       delay_ms=float(p_b[1] / p_m[0] * 1000), se_delay_ms=float(s_b[1] / p_m[0] * 1000))
            print(f"B_local 對 ω² 回歸：斜率 {p_b[1]:+.5f} ± {s_b[1]:.5f}；若全部來自電流相對編碼器的延遲，"
                  f"延遲 = {out['delay_ms']:+.2f} ± {out['se_delay_ms']:.2f} ms（也可能來自摩擦的動態特性）")
    return out


# ============================================================
# 主流程
# ============================================================

def main():
    print("=" * 70)
    print(" 速度調變實驗（巡航速度上疊加正弦）：安全檢查")
    print("=" * 70)
    print(f"方向: {'正轉 (+1)' if DIRECTION > 0 else '反轉 (-1)'}")
    print(f"姿態: {POSE}（J1~J5 期望 {list(POSES[POSE]['q_deg']) if POSE in POSES else '未知'}°）")
    print(f"負載: PAYLOAD_KG = {PAYLOAD_KG}")

    conds, n_repeat = stage_config(TEST_STAGE)
    segs = plan_segments(conds, n_repeat, DIRECTION, V_CRUISE, T_MOD, SETTLE_TIME, RAMP_ACCEL, STOPJ_DECEL,
                         OVERRUN_DEG, WINDOW_CENTER_DEG, GEAR_RATIO, DT)
    print(f"模式: {TEST_STAGE}，條件 (Hz, rad/s): {conds} × 各 {n_repeat} 次，共 {len(segs)} 段")
    print(f"巡航 {V_CRUISE} rad/s，調變 {T_MOD} s；調變段 {segs[0]['mod_start_deg']:.3f}° 至 "
          f"{segs[0]['mod_end_deg']:.3f}°（第一段），起點錯開 {[round(o, 4) for o in start_offsets_deg(n_repeat)]}°")

    ok, report, (lo_plan, hi_plan) = safety_check(
        segs, V_CRUISE, T_MOD, QD_MAX, QDD_CHECK_LIMIT, V_MIN_ALLOWED, RAMP_ACCEL, MOD_ACCEL_LIMIT,
        STOPJ_DECEL, MOVEJ_VEL, MOVEJ_ACCEL, J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG, SAFETY_MARGIN_DEG,
        end_deg=END_J0_DEG)
    print("\n".join(report))
    speed_ok, speed_msg = check_flange_speed(POSE, conds, V_CRUISE)
    print(speed_msg)
    payload_ok, payload_msg = check_payload_setting(PAYLOAD_KG)
    print(payload_msg + ("  [OK]" if payload_ok else "  [FAIL]"))
    if not (ok and speed_ok and payload_ok):
        print("\n[中止] 安全檢查未通過。")
        return

    print("\n[通過安全檢查]")
    if ROBOT_IP is None:
        print("\n[中止] ROBOT_IP 未設定。")
        return
    if rtde_receive is None:
        print("\n[中止] 找不到 rtde_receive，僅完成離線設計檢查。")
        return

    if POSE == "DIAG":
        print("\n執行前姿態確認（J1~J5 診斷原點）")
        if home_pose is None:
            print("[警告] 找不到 ur5_home_pose 模組，請手動確認 J1~J5 為 -90°, 90°, -90°, -90°, 0°")
            if input("已確認，輸入 yes 繼續：").strip().lower() != "yes":
                print("已取消。")
                return
        elif not home_pose.ensure_home_pose(ROBOT_IP, mode="CHECK"):
            print("\n[中止] 姿態不符合診斷原點。")
            return
    else:
        r = POSES[POSE]["r_flange_m"]
        print(f"\n[注意] 姿態 {POSE}：J1~J5 須先以教導器擺到 {list(POSES[POSE]['q_deg'])}°（程式不會移動這五軸）。")
        print(f"       法蘭離 J0 軸約 {r:.2f} m（不含工具），J0 會掃過 {lo_plan:.1f}° 至 {hi_plan:.1f}°（含餘裕）。")
        print("       請確認此扇形範圍（加上工具與線材）沒有人員與物品；連上 RTDE 後會再以實際角度檢查 J1~J5。")
        if input(f"已確認，輸入 {POSE} 繼續：").strip() != POSE:
            print("已取消。")
            return

    print("\n[負載確認] 請看教導器上設定的負載質量（加減質量後必須先在教導器更新質量與重心）。")
    entry_ok, entry_msg = check_payload_entry(PAYLOAD_KG, input("輸入教導器上的負載質量 [kg]："))
    print(entry_msg)
    if not entry_ok:
        print("\n[中止] 負載確認未通過，未送出任何軌跡。")
        return

    print(f"\n即將對 IP={ROBOT_IP} 送出軌跡，關節 J{JOINT_INDEX}，方向 {DIRECTION:+d}，姿態 {POSE}，"
          f"負載 {PAYLOAD_KG} kg，共 {len(segs)} 段。")
    print("*** 送出後 Ctrl+C 會送出 stopj，但最可靠的仍是緊急停止按鈕 ***")
    input("按 Enter 繼續，或 Ctrl+C 取消...")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_label = "pos" if DIRECTION > 0 else "neg"
    tag = f"{POSE.lower()}_{payload_tag(PAYLOAD_KG)}_{dir_label}_{TEST_STAGE.lower()}_{ts}"
    session_dir = os.path.join(OUTPUT_DIR, f"velmod_{tag}")
    os.makedirs(session_dir, exist_ok=True)
    csv_path = os.path.join(session_dir, f"velmod_data_{tag}.csv")
    info_path = os.path.join(session_dir, f"velmod_info_{tag}.txt")

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

    q0_deg = math.degrees(logger.last_q0)
    pos_ok, pos_msg = check_start_position(q0_deg, segs[0]["start_deg"], MAX_INITIAL_MOVE_DEG,
                                           J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG)
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

    script = build_urscript(segs, DIRECTION, JOINT_INDEX, logger.last_q_full, V_CRUISE, RAMP_ACCEL,
                            MOD_ACCEL_LIMIT, STOPJ_DECEL, MOVEJ_ACCEL, MOVEJ_VEL, PAUSE_TIME, DT,
                            end_deg=END_J0_DEG)
    est_total = estimate_total_time(segs, q0_deg, STOPJ_DECEL, end_deg=END_J0_DEG)
    status, elapsed, span = "unknown", 0.0, 0.0
    try:
        print(f"[資訊] 送出 URScript（{len(script.splitlines())} 行，預估 {est_total:.1f} 秒）...")
        send_urscript(ROBOT_IP, script)
        status, elapsed, span = wait_for_motion_complete(logger, est_total)
        if status == "no_motion":
            print(f"\n[警告] 送出後 {elapsed:.1f}s 內未偵測到 J0 運動。")
            print("       控制器可能未接受腳本（不在 remote control 模式，或腳本解析失敗）。")
        elif status == "timeout":
            print(f"\n[警告] 超過硬性逾時（實際 {elapsed:.1f}s / 預估 {est_total:.1f}s），主動中止。")
            send_abort(ROBOT_IP)
        else:
            print(f"\n[資訊] 運動結束，實際運動 {span:.1f}s / 預估 {est_total:.1f}s")
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
        measured = [
            f"實測 J0 範圍 : {q_lo:.3f}° 至 {q_hi:.3f}°（規劃含餘裕 {lo_plan:.3f}° 至 {hi_plan:.3f}°）",
            f"結束位置     : {math.degrees(logger.last_q0):.3f}°",
            f"實際運動時長 / 預估 : {span:.2f}s / {est_total:.2f}s",
        ]
        if q_lo < lo_plan or q_hi > hi_plan:
            measured.append("[警告] 實測範圍超出規劃範圍，請先查明原因再跑下一階段。")
        print("\n" + "\n".join(measured))
    if not logger.speed_fraction_held():
        print("\n[警告] 執行中教導器速度滑桿曾低於 100% 或讀不到（見摘要），這次的資料不可用。")

    with open(info_path, "w", encoding="utf-8") as f:
        f.write(f"方向: {DIRECTION:+d}\n模式: {TEST_STAGE}\n負載 PAYLOAD_KG: {PAYLOAD_KG} kg\n")
        f.write("\n".join(pose_lines) + "\n")
        f.write(f"條件 (Hz, rad/s): {conds}\n每個條件次數: {n_repeat}\n")
        f.write(f"巡航速度: {V_CRUISE}\n調變時間: {T_MOD}s\n定速時間: {SETTLE_TIME}s\n")
        f.write(f"起點錯開 [deg]: {[round(o, 6) for o in start_offsets_deg(n_repeat)]}\n")
        f.write("段順序（條件編號）: " + ", ".join(str(sg["cond"] + 1) for sg in segs) + "\n")
        f.write(f"預估總時長: {est_total:.1f}s\n執行狀態: {status}\n")
        f.write(logger.summary() + "\n")
        if measured:
            f.write("\n" + "\n".join(measured) + "\n")
    print(f"\n[完成] 資料存至 {csv_path}")
    print(f"[完成] 摘要存至 {info_path}")
    if POSE != "DIAG":
        print(f"\n[提醒] J0 已回 {END_J0_DEG:.0f}°，J1~J5 仍在 {POSE} 姿態；請以教導器移回診斷原點。")
    print("\n下一步：analyze_velocity_modulation(csv_path, JOINT_INDEX, KT_OUT, DIRECTION, CONDITIONS, V_CRUISE)")


if __name__ == "__main__":
    main()
