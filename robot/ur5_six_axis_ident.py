"""
UR5 六軸同時鑑別（雛形）—— 一次動作量六軸的 B、Tc（定速）與 J（速度調變）
================================================================

★★★ 想法 ★★★

單軸程式一次只動一軸，六軸 × 定速、速度調變 × 正反轉要上機十幾次。這支讓六軸**沿同一條關節空間直線同步移動**：
  q(λ) = 中心 + λ · S · 半行程（S：各軸的方向，見 PATH_SIGN）
1. 定速（MODE = "CONSTVEL"）：六軸同一個速度 v，同時掃過各自的角度窗口（6 個馬達圈），去程、回程走同一條路。
   重力 g(q) 只跟位置有關；科氏力、離心力和速度平方成正比 —— 兩者正反轉都相同。
   摩擦和運動方向同號，所以每一軸照單軸定速的方法分正反轉回歸、合併（combine_directions），
   重力、科氏力在 Tc = (C₊ + C₋)/2 的相消中一起消掉。去程、回程在同一次執行裡交替，不必分兩次。
2. 速度調變（MODE = "VELMOD"）：六軸在各自的巡航速度上疊加正弦，**每軸用不同頻率**（FREQ_OFFSET，0.25 Hz 間隔），
   調變 T_MOD = 4 s 是每個頻率的整數週期，不同軸的加速度在調變段內互相正交，
   各軸照單軸速度調變的回歸（tau = J·a + B·v + …）只會抓到自己頻率的成分，得到對角項 J（含轉子）。
   其他軸的耦合（M_ij·a_j）落在別的頻率，正交後只剩修邊造成的少量洩漏（tests 以剛體模型模擬檢查）。

分析沿用單軸程式的函式（ur5_j3_const_velocity_ident.analyze_angle_window、combine_directions，
ur5_j3_velocity_modulation_ident.analyze_j3_velocity_modulation；兩者都以 joint_index 參數化），每軸各跑一次。

★★★ 姿態 SIX ★★★

(180, −90, 70, −70, −90, 0)：從診斷原點把 J2 轉到 70°（前臂往上 20°）、J3 轉到 −70°（工具仍垂直朝下）。
方向 S = (+1, +1, −1, −1, +1, +1)：J1、J2 反向，前臂幾乎平移、法蘭高度變化小。
以 UR5 官方 DH 算，整條路徑（含餘裕）法蘭離底座安裝面 ≥ 0.58 m（診斷原點 0.43 m）。
可以用 ur5_goto_pose.py SIX 自動轉到這個姿態。

★★★ 使用前必讀 ★★★
1. 【新程式，尚未上機】依 TEST_STAGE = QUICK → MID → FULL 階梯執行，MODE 兩種各自從 QUICK 開始。
2. 【六軸視窗未設定一律中止】SAFE_MIN_DEG / SAFE_MAX_DEG（六軸各一個）由操作者依現場線材與淨空填寫；
   所有段的起點、煞車後最遠點、結束點各外推 SAFETY_MARGIN_DEG 後都必須在視窗內。
3. 【法蘭高度、工具尖端高度與線速度】沿整條規劃路徑以官方 DH 計算：法蘭高度 ≥ FLANGE_Z_MIN_M，
   工具尖端（法蘭 z 軸往外 TOOL_LENGTH_M）高度 ≥ TOOL_TIP_Z_MIN_M，法蘭線速度 ≤ FLANGE_SPEED_MAX。
4. 【姿態】程式不會擺姿態：連上 RTDE 後以實際角度檢查六軸（±0.5°）。
5. 【負載 PAYLOAD_KG】填教導器上設定的負載質量，未設定一律中止；執行前須輸入，不符即中止。
6. 【第一次 movej 的距離有上限】任何一軸與第一段起點相差超過 MAX_INITIAL_MOVE_DEG 就中止。
7. 【教導器速度滑桿】不是 100% 或讀不到一律中止。
8. 【IP 需要你自己再次確認】（第二台 UR5 只改 ROBOT_IP；裸機 PAYLOAD_KG = 0）
"""

import math
import os
import sys
import time
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ur5_j3_const_velocity_ident as cv      # 定速分析（以 joint_index 參數化）
import ur5_j3_velocity_modulation_ident as vm  # 速度調變分析、RTDE 記錄、URScript 傳送、檢查函式

rtde_receive = vm.rtde_receive


# ============================================================
# 使用者設定
# ============================================================

ROBOT_IP = "192.168.50.114"   # 執行前請再次確認
MODE = "CONSTVEL"             # "CONSTVEL"（B、Tc）或 "VELMOD"（J）
# QUICK: CONSTVEL 單一速度 0.05 去回各一段；VELMOD 第一個條件去回各一段 —— 驗證六軸同步、角度停止與 movej 定位
#        （VELMOD QUICK 每條件只有 1 段，不算 J）
# MID  : CONSTVEL 三檔速度；VELMOD 條件 2、4（2 Hz 大／小振幅）各兩次
# FULL : CONSTVEL 九檔先升後降；VELMOD 六個條件 × N_REPEAT 次（起點錯開）
TEST_STAGE = "QUICK"

POSE = "SIX"
POSES = {"SIX": dict(q_deg=(180.0, -90.0, 70.0, -70.0, -90.0, 0.0))}
PATH_SIGN = (+1, +1, -1, -1, +1, +1)   # 「去程」時各軸角度變化方向；回程全部相反
POSE_TOL_DEG = 0.5

GEAR_RATIO = 101
SAMPLE_HZ = 125.0
DT = 1.0 / SAMPLE_HZ
KT_OUT = GEAR_RATIO * 0.1350           # 文獻值 13.64（同單軸程式）；分析結果再乘 K_AXIS / KT_OUT
K_AXIS = (10.8, 12.03, 12.10, 9.31, 9.31, 9.31)   # 各軸換算比例：J0 11/4 砝碼、J1 J2 11/12 砝碼、J3～J5 暫定（101 × 0.0922）

# ---------- 定速 ----------
SPEED_LEVELS = [0.005, 0.01, 0.02, 0.03, 0.05, 0.08, 0.12, 0.16, 0.20]   # rad/s（同單軸；六軸同時動時 0.20 的法蘭線速度約 0.15 m/s）
N_MOTOR_REVS_WINDOW = 6        # 各軸窗口 6 個馬達圈 = 21.39°，中心 = 姿態角度
SETTLE_TIME = 0.5              # 進入窗口前穩定時間 [s]
OVERRUN_DEG = 0.5              # 越過窗口終點後再走多少才 stopj [deg]
SPEEDJ_ACCEL = 0.5             # 定速 speedj 加速度 [rad/s^2]

# ---------- 速度調變 ----------
T_MOD = 4.0                    # 調變時間 [s]：所有頻率都是 0.25 Hz 的倍數 → 整數週期，各軸加速度互相正交
TIERS = [(1.0, "L"), (2.0, "L"), (3.0, "L"), (2.0, "S"), (1.0, "S"), (3.0, "S")]   # 六個條件：基準頻率 × 振幅檔（同單軸順序）
FREQ_OFFSET = (-0.50, 0.00, +0.25, -0.25, +0.75, +0.50)   # 各軸頻率 = 基準 + 偏移；J1 維持 1、2、3 Hz 方便和單軸比
V_CRUISE = (0.10, 0.10, 0.10, 0.15, 0.15, 0.15)           # 巡航速度 [rad/s]；手腕提高，最低速度離開低速區（11/19 手腕 J 為負）
AMP = {"L": (0.04, 0.04, 0.04, 0.06, 0.06, 0.06),          # 速度振幅 [rad/s]
       "S": (0.02, 0.02, 0.02, 0.03, 0.03, 0.03)}
N_REPEAT = 3                   # FULL 每個條件的次數；第 k 次起點錯開 (k − 1)/3 個馬達圈
MOD_SETTLE_TIME = 0.8          # 到達巡航後、開始調變前的定速時間 [s]
RAMP_ACCEL = 0.5               # 加速到巡航的 speedj 加速度（最快的一軸）[rad/s^2]
MOD_ACCEL_LIMIT = 2.0          # 調變段 speedj 加速度上限 [rad/s^2]
V_MIN_ALLOWED = 0.03           # 調變中最低速度的下限 [rad/s]

# ---------- 共用 ----------
STOPJ_DECEL = 1.0              # stopj 減速度（最快的一軸）[rad/s^2]
MOVEJ_ACCEL = 0.5
MOVEJ_VEL = 0.15
PAUSE_TIME = 0.5
QD_MAX = 0.30                  # 單軸速度上限（漣波在 125 Hz 下混疊）
LOOP_TIME_FACTOR = 1.5
SAFETY_MARGIN_DEG = 2.0
MAX_INITIAL_MOVE_DEG = 30.0
FLANGE_Z_MIN_M = 0.40
TOOL_LENGTH_M = 0.20           # 末端工具長度（法蘭 z 軸方向），只用於工具尖端高度檢查
TOOL_TIP_Z_MIN_M = 0.30
FLANGE_SPEED_MAX = 0.25        # 法蘭線速度上限 [m/s]
SPEED_FRACTION_MIN = 0.999

# 六軸允許的絕對角度視窗（未填一律中止）。11/25 使用者採用：FULL 規劃範圍（含 2° 餘裕）再往外約 5°
SAFE_MIN_DEG = (155.0, -115.0, 45.0, -105.0, -125.0, -35.0)
SAFE_MAX_DEG = (205.0, -65.0, 95.0, -35.0, -55.0, 35.0)

PAYLOAD_KG = 1.30             # 末端工具（FT300、夾爪、3D 相機）；教導器負載須相同，執行時還要再輸入一次
PAYLOAD_TOL_KG = 0.005

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


# ============================================================
# 第一部分：條件、規劃
# ============================================================

def joint_conditions(j, tiers=TIERS, offsets=FREQ_OFFSET, amp=AMP):
    """第 j 軸在各條件的 (頻率, 振幅)。"""
    return [(f + offsets[j], amp[a][j]) for f, a in tiers]


def check_conditions(tiers=TIERS, offsets=FREQ_OFFSET, amp=AMP, v_cruise=V_CRUISE, t_mod=T_MOD,
                     v_min=V_MIN_ALLOWED, mod_accel=MOD_ACCEL_LIMIT, qd_max=QD_MAX):
    """回傳 (ok, lines)：同一條件六軸頻率互不相同、整數週期、最低速度、加速度、速度上限。"""
    lines, ok = [], True
    for ci, (f0, a) in enumerate(tiers):
        fs = [f0 + offsets[j] for j in range(6)]
        if len({round(f, 6) for f in fs}) < 6 or min(fs) <= 0:
            ok = False; lines.append(f"[FAIL] 條件 {ci + 1}：六軸頻率有重複或 ≤ 0：{fs}")
        for j in range(6):
            f, A = fs[j], amp[a][j]
            cyc = f * t_mod
            if abs(cyc - round(cyc)) > 1e-6:
                ok = False; lines.append(f"[FAIL] J{j} {f} Hz：T_MOD 內不是整數週期（{cyc:.3f}）")
            if v_cruise[j] - A < v_min:
                ok = False; lines.append(f"[FAIL] J{j}：最低速度 {v_cruise[j] - A:.3f} < {v_min}")
            if 2 * math.pi * f * A > 0.9 * mod_accel:
                ok = False; lines.append(f"[FAIL] J{j} {f} Hz：調變加速度 {2*math.pi*f*A:.2f} > 0.9 × {mod_accel}")
            if v_cruise[j] + A > qd_max:
                ok = False; lines.append(f"[FAIL] J{j}：最高速度 {v_cruise[j] + A:.3f} > {qd_max}")
    if ok:
        lines.append(f"[OK] 六個條件、六軸頻率互不相同、整數週期；最低速度 ≥ {v_min}、調變加速度 ≤ 0.9 × {mod_accel}")
    return ok, lines


def stage_config(mode, stage):
    """CONSTVEL → 速度清單；VELMOD → (條件索引清單, 重複次數)。"""
    if mode == "CONSTVEL":
        if stage == "QUICK":
            return [0.05]
        if stage == "MID":
            return [0.02, 0.05, 0.12]
        if stage == "FULL":
            return list(SPEED_LEVELS) + list(reversed(SPEED_LEVELS))
    elif mode == "VELMOD":
        if stage == "QUICK":
            return [0], 1
        if stage == "MID":
            return [1, 3], 2               # 分析每個條件每個方向至少要 2 段才算 J
        if stage == "FULL":
            return list(range(len(TIERS))), N_REPEAT
    raise ValueError(f"未知的 MODE / TEST_STAGE: {mode} / {stage}")


def window_bounds(center_deg, n_revs=N_MOTOR_REVS_WINDOW, gear_ratio=GEAR_RATIO):
    w = n_revs * 360.0 / gear_ratio
    return [(c - w / 2.0, c + w / 2.0) for c in center_deg]


def plan_constvel(speeds, center_deg, sign=PATH_SIGN, accel=SPEEDJ_ACCEL, decel=STOPJ_DECEL,
                  settle_time=SETTLE_TIME, overrun_deg=OVERRUN_DEG, loop_time_factor=LOOP_TIME_FACTOR, dt=DT):
    """
    每個速度依序做去程（d = +1）、回程（d = −1）。各軸方向 dj = d·S_j，窗口、加速距離、停止點都以度為單位（六軸各一）。
    六軸速度相同、加速度相同 → 同時進出窗口。
    """
    win = window_bounds(center_deg)
    segs = []
    for v in speeds:
        pre = math.degrees(v * v / (2.0 * accel) + v * settle_time)
        brake = math.degrees(v * v / (2.0 * decel))
        w = win[0][1] - win[0][0]
        t_expect = math.radians(w + pre + overrun_deg) / v + v / accel
        for d in (+1, -1):
            dj = [d * s for s in sign]
            ws = [lo if x > 0 else hi for x, (lo, hi) in zip(dj, win)]
            we = [hi if x > 0 else lo for x, (lo, hi) in zip(dj, win)]
            start = [a - x * pre for a, x in zip(ws, dj)]
            stop = [b + x * overrun_deg for b, x in zip(we, dj)]
            segs.append(dict(index=len(segs), speed=v, d=d, dj=dj, start_deg=start, win_start_deg=ws, win_end_deg=we,
                             stop_deg=stop, extreme_deg=[s + x * brake for s, x in zip(stop, dj)],
                             qd=[x * v for x in dj], t_expect=t_expect,
                             n_max=int(math.ceil(loop_time_factor * t_expect / dt)) + 10))
    return segs


def start_offsets_deg(n_repeat, gear_ratio=GEAR_RATIO):
    rev = 360.0 / gear_ratio
    return [(k - (n_repeat - 1) / 2.0) / n_repeat * rev for k in range(n_repeat)]


def plan_velmod(cond_idx, n_repeat, center_deg, sign=PATH_SIGN, v_cruise=V_CRUISE, t_mod=T_MOD,
                settle_time=MOD_SETTLE_TIME, ramp_accel=RAMP_ACCEL, decel=STOPJ_DECEL, overrun_deg=1.0, dt=DT):
    """
    每個條件做去程、回程各一段。speedj／stopj 的加速度以最快的一軸為準、其他軸等比例（六軸同時到達巡航、同時停）：
      加速時間 t_r = max(vc)/ramp_accel，第 j 軸加速＋定速距離 = vc_j·(t_r/2 + settle)；
      最壞煞車時間 t_s = max(vc + A)/decel，第 j 軸煞車距離 ≤ (vc_j + A_j)·t_s/2。
    """
    conds = [joint_conditions(j) for j in range(6)]
    t_r = max(v_cruise) / ramp_accel
    n_pre = int(round((t_r + settle_time) / dt))
    n_mod = int(round(t_mod / dt))
    segs = []
    for k, off in enumerate(start_offsets_deg(n_repeat)):
        for ci in cond_idx:
            fa = [conds[j][ci] for j in range(6)]
            t_s = max(v_cruise[j] + fa[j][1] for j in range(6)) / decel
            for d in (+1, -1):
                dj = [d * s for s in sign]
                L = [math.degrees(v_cruise[j] * t_mod) for j in range(6)]
                ms = [center_deg[j] - dj[j] * L[j] / 2.0 + dj[j] * off for j in range(6)]
                me = [ms[j] + dj[j] * L[j] for j in range(6)]
                pre = [math.degrees(v_cruise[j] * (t_r / 2.0 + settle_time)) for j in range(6)]
                guard = [me[j] + dj[j] * overrun_deg for j in range(6)]
                segs.append(dict(index=len(segs), cond=ci, repeat=k, offset_deg=off, d=d, dj=dj, fa=fa,
                                 start_deg=[ms[j] - dj[j] * pre[j] for j in range(6)], mod_start_deg=ms, mod_end_deg=me,
                                 guard_deg=guard,
                                 extreme_deg=[guard[j] + dj[j] * math.degrees((v_cruise[j] + fa[j][1]) * t_s / 2.0)
                                              for j in range(6)],
                                 n_pre=n_pre, n_mod=n_mod, t_active=(n_pre + n_mod) * dt))
    return segs


def estimate_movej_time(dist_rad, v_peak=MOVEJ_VEL, accel=MOVEJ_ACCEL):
    return vm.estimate_movej_time(dist_rad, v_peak, accel)


def estimate_total_time(segs, q0_deg, mode, end_deg=None):
    total, pos = 0.0, list(q0_deg)
    for s in segs:
        total += estimate_movej_time(math.radians(max(abs(a - b) for a, b in zip(s["start_deg"], pos))))
        if mode == "CONSTVEL":
            total += PAUSE_TIME + s["t_expect"] + s["speed"] / STOPJ_DECEL + PAUSE_TIME
        else:
            total += PAUSE_TIME + s["t_active"] + 0.3 + PAUSE_TIME
        pos = s["extreme_deg"]
    if end_deg is not None:
        total += estimate_movej_time(math.radians(max(abs(a - b) for a, b in zip(end_deg, pos)))) + PAUSE_TIME
    return total


# ============================================================
# 第二部分：安全檢查
# ============================================================

def fk_frames(q):
    return cv.fk_frames(q)


def flange_and_tip_z(q_deg, tool_len=TOOL_LENGTH_M):
    F = fk_frames([math.radians(x) for x in q_deg])[6]
    return float(F[2, 3]), float((F[:3, :3] @ np.array([0.0, 0.0, tool_len]) + F[:3, 3])[2])


def flange_speed(q_deg, qd):
    """法蘭線速度 [m/s]（數值微分的位置 Jacobian）。"""
    q = np.radians(q_deg)
    p0 = fk_frames(q)[6][:3, 3]
    v = np.zeros(3)
    for j in range(6):
        dq = np.zeros(6); dq[j] = 1e-6
        v += (fk_frames(q + dq)[6][:3, 3] - p0) / 1e-6 * qd[j]
    return float(np.linalg.norm(v))


def planned_range(segs, end_deg, margin_deg):
    pts = [s[k] for s in segs for k in ("start_deg", "extreme_deg")] + ([list(end_deg)] if end_deg else [])
    lo = [min(p[j] for p in pts) - margin_deg for j in range(6)]
    hi = [max(p[j] for p in pts) + margin_deg for j in range(6)]
    return lo, hi


def safety_check(segs, mode, q_pose_deg, safe_min=SAFE_MIN_DEG, safe_max=SAFE_MAX_DEG, margin_deg=SAFETY_MARGIN_DEG,
                 flange_z_min=FLANGE_Z_MIN_M, tip_z_min=TOOL_TIP_Z_MIN_M, flange_speed_max=FLANGE_SPEED_MAX,
                 qd_max=QD_MAX, n_check=41):
    """
    回傳 (ok, lines, (lo, hi))。
    1. 各軸速度上限；調變條件（VELMOD）。
    2. 六軸絕對角度：所有段起點、煞車後最遠點、結束點各外推 margin_deg，須在 SAFE_MIN/MAX 內；未設定一律 FAIL。
    3. 每一段從起點到最遠點的直線上（含 margin，沿路徑方向外推）檢查法蘭高度、工具尖端高度，
       以該段最高速度算法蘭線速度。
    """
    lines, ok = [], True
    if mode == "VELMOD":
        c_ok, c_lines = check_conditions()
        ok &= c_ok; lines += c_lines
        vmax = [max(V_CRUISE[j] + max(a[j] for a in AMP.values()), MOVEJ_VEL) for j in range(6)]
    else:
        vmax = [max(max(s["speed"] for s in segs), MOVEJ_VEL)] * 6
    if max(vmax) > qd_max:
        ok = False; lines.append(f"[FAIL] 速度 {max(vmax):.3f} > {qd_max}")
    lo, hi = planned_range(segs, q_pose_deg, margin_deg)
    for j in range(6):
        if safe_min[j] is None or safe_max[j] is None:
            ok = False; lines.append(f"[FAIL] J{j} 視窗未設定（SAFE_MIN_DEG / SAFE_MAX_DEG）")
        elif lo[j] < safe_min[j] or hi[j] > safe_max[j]:
            ok = False; lines.append(f"[FAIL] J{j} 規劃範圍 {lo[j]:.2f}°～{hi[j]:.2f}° 超出視窗 {safe_min[j]}°～{safe_max[j]}°")
        else:
            lines.append(f"[OK] J{j} 規劃範圍 {lo[j]:.2f}°～{hi[j]:.2f}°（視窗 {safe_min[j]}°～{safe_max[j]}°）")
    zmin, tmin, vfl = 9.0, 9.0, 0.0
    for s in segs:
        a, b = np.array(s["start_deg"]), np.array(s["extreme_deg"])
        u = np.sign(b - a)
        a, b = a - u * margin_deg, b + u * margin_deg
        qd = np.array(s["qd"]) if mode == "CONSTVEL" else np.array(
            [s["dj"][j] * (V_CRUISE[j] + s["fa"][j][1]) for j in range(6)])
        for lam in np.linspace(0.0, 1.0, n_check):
            qq = a + lam * (b - a)
            z, tp = flange_and_tip_z(qq)
            zmin, tmin = min(zmin, z), min(tmin, tp)
            vfl = max(vfl, flange_speed(qq, qd))
    for name, val, lim, cmp in (("法蘭高度", zmin, flange_z_min, ">="), ("工具尖端高度", tmin, tip_z_min, ">="),
                                ("法蘭線速度", vfl, flange_speed_max, "<=")):
        good = (val >= lim) if cmp == ">=" else (val <= lim)
        if not (np.isfinite(lim) and good):
            ok = False
        lines.append(f"[{'OK' if good else 'FAIL'}] {name} {val:.3f}（須 {cmp} {lim}）")
    return ok, lines, (lo, hi)


def check_start_position(q_now_deg, first_start_deg, max_move_deg=MAX_INITIAL_MOVE_DEG):
    d = max(abs(a - b) for a, b in zip(q_now_deg, first_start_deg))
    ok = d <= max_move_deg
    return ok, f"[{'OK' if ok else 'FAIL'}] 第一次 movej 最大單軸距離 {d:.2f}°（上限 {max_move_deg}°）"


# ============================================================
# 第三部分：URScript
# ============================================================

def _vec(v):
    return "[" + ", ".join(f"{x:.6f}" for x in v) + "]"


def _guard(limits_rad, dj):
    """任何一軸越過自己的極限就離開迴圈（六軸各自比較，全部未越過才繼續）。"""
    return " and ".join(f"(q_now[{j}] {'<' if dj[j] > 0 else '>'} {limits_rad[j]:.6f})" for j in range(6))


def build_urscript(segs, mode, end_deg=None, dt=DT):
    lines = ["def six_axis_ident_program():"]
    for s in segs:
        start = [math.radians(x) for x in s["start_deg"]]
        lines.append(f"  # segment {s['index'] + 1}: d = {s['d']:+d}" +
                     (f", v = {s['speed']:.3f} rad/s" if mode == "CONSTVEL" else f", condition {s['cond'] + 1}, repeat {s['repeat'] + 1}"))
        lines.append(f"  movej({_vec(start)}, a={MOVEJ_ACCEL}, v={MOVEJ_VEL})")
        lines.append(f"  sleep({PAUSE_TIME})")
        lines.append("  q_now = get_actual_joint_positions()")
        lines.append("  n = 0")
        if mode == "CONSTVEL":
            cond = _guard([math.radians(x) for x in s["stop_deg"]], s["dj"])
            lines.append(f"  while (n < {s['n_max']}) and {cond}:")
            lines.append(f"    speedj({_vec(s['qd'])}, a={SPEEDJ_ACCEL}, t={dt})")
            lines.append("    q_now = get_actual_joint_positions()")
            lines.append("    n = n + 1")
            lines.append("  end")
        else:
            cond = _guard([math.radians(x) for x in s["guard_deg"]], s["dj"])
            qc = [s["dj"][j] * V_CRUISE[j] for j in range(6)]
            lines.append(f"  while (n < {s['n_pre']}) and {cond}:")
            lines.append(f"    speedj({_vec(qc)}, a={RAMP_ACCEL}, t={dt})")
            lines.append("    q_now = get_actual_joint_positions()")
            lines.append("    n = n + 1")
            lines.append("  end")
            lines.append("  n = 0")
            terms = [f"{s['dj'][j] * V_CRUISE[j]:.6f} + {s['dj'][j] * s['fa'][j][1]:.6f} * sin({2 * math.pi * s['fa'][j][0] * dt:.9f} * n)"
                     for j in range(6)]
            lines.append(f"  while (n < {s['n_mod']}) and {cond}:")
            lines.append("    qd_t = [" + ", ".join(terms) + "]")
            lines.append(f"    speedj(qd_t, a={MOD_ACCEL_LIMIT}, t={dt})")
            lines.append("    q_now = get_actual_joint_positions()")
            lines.append("    n = n + 1")
            lines.append("  end")
        lines.append(f"  stopj({STOPJ_DECEL})")
        lines.append(f"  sleep({PAUSE_TIME})")
    if end_deg is not None:
        lines.append("  # return to the start pose")
        lines.append(f"  movej({_vec([math.radians(x) for x in end_deg])}, a={MOVEJ_ACCEL}, v={MOVEJ_VEL})")
        lines.append(f"  sleep({PAUSE_TIME})")
    lines.append("end")
    return "\n".join(lines) + "\n"


# ============================================================
# 第四部分：離線分析（每軸沿用單軸函式）
# ============================================================

def _quiet(fn, *a, **k):
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def analyze_constvel(csv_path, center_deg=POSES["SIX"]["q_deg"], speeds=SPEED_LEVELS, k_axis=K_AXIS, joints=range(6)):
    """每軸：窗口內、target_qd = ±v 的段落分正反轉回歸，再合併。回傳 {j: dict(B, Tc, H, …)}（已乘 K_AXIS / KT_OUT）。"""
    win = window_bounds(center_deg)
    out = {}
    for j in joints:
        lo, hi = win[j]
        rp = _quiet(cv.analyze_angle_window, csv_path, j, KT_OUT, +1, speeds, lo, hi)
        rn = _quiet(cv.analyze_angle_window, csv_path, j, KT_OUT, -1, speeds, lo, hi)
        if "B" not in rp or "B" not in rn:
            out[j] = dict(n_pos=len(rp["segments"]), n_neg=len(rn["segments"]))
            continue
        c = _quiet(cv.combine_directions, rp, rn)
        s = k_axis[j] / KT_OUT
        out[j] = dict(B=c["B"] * s, sd_B=c["sd_B"] * s, Tc=c["Tc"] * s, sd_Tc=c["sd_Tc"] * s, H=c["H"] * s,
                      B_pos=c["B_pos"] * s, B_neg=c["B_neg"] * s, n_pos=len(rp["segments"]), n_neg=len(rn["segments"]),
                      K_used=k_axis[j])
    return out


LINK_M = (3.7, 8.393, 2.33, 1.219, 1.219, 0.1879)        # UR 官方 DH 頁面（UR5，CB3）連桿質量 [kg]
LINK_COM = ((0, -0.02561, 0.00193), (0.2125, 0, 0.11336), (0.15, 0, 0.0265),
            (0, -0.0018, 0.01634), (0, 0.0018, 0.01634), (0, 0, -0.001159))   # 重心（各連桿 DH 座標系）[m]
TOOL_COG = (0.001, 0.001, 0.074)                           # 末端工具重心（法蘭座標，教導器設定）[m]


def gravity_torque(q, payload_kg, g=9.81):
    """撐住重力的六軸力矩 [N·m]（名目模型：官方連桿點質量 + 負載點質量；正負號同 target_moment）。"""
    fr = fk_frames(q)
    pts = [(LINK_M[i], fr[i + 1][:3, :3] @ np.array(LINK_COM[i]) + fr[i + 1][:3, 3], i) for i in range(6)]
    pts.append((payload_kg, fr[6][:3, :3] @ np.array(TOOL_COG) + fr[6][:3, 3], 5))
    up = np.array([0.0, 0.0, g])
    return np.array([sum(float(np.cross(c - fr[j][:3, 3], m * up) @ fr[j][:3, 2]) for m, c, i in pts if i >= j)
                     for j in range(6)])


def gravity_compensated_csv(csv_path, payload_kg, k_axis=K_AXIS):
    """
    把名目模型的重力（沿實際六軸角度）從電流扣掉，另存 *_gcomp.csv（電流欄 = 電流 − g_j(q)/K_j）。
    六軸同時動時，第 j 軸的重力隨所有軸改變；單軸回歸的重力項（sin、cos q_j）只描述自己的角度，
    調變造成的角度振盪會把 ∂g/∂q 的誤差帶進 J（偏差約 −∂g/∂q/ω²，低頻最大）。扣掉名目重力後只剩模型誤差。
    """
    import csv as _csv
    out_path = os.path.splitext(csv_path)[0] + "_gcomp.csv"
    with open(csv_path, encoding="utf-8") as fi, open(out_path, "w", newline="", encoding="utf-8") as fo:
        rd = _csv.DictReader(fi)
        wr = _csv.DictWriter(fo, fieldnames=rd.fieldnames)
        wr.writeheader()
        for r in rd:
            q = [float(r[f"actual_q_{j}"]) for j in range(6)]
            gq = gravity_torque(q, payload_kg)
            for j in range(6):
                r[f"actual_current_{j}"] = float(r[f"actual_current_{j}"]) - gq[j] / k_axis[j]
            wr.writerow(r)
    return out_path


def analyze_velmod(csv_path, k_axis=K_AXIS, joints=range(6), variant="完整", payload_kg=None):
    """
    每軸：自己的條件（頻率 = 基準 + 偏移）、自己的巡航速度，正反轉各跑一次單軸速度調變回歸。
    payload_kg 不為 None 時先扣掉名目重力（gravity_compensated_csv）。
    """
    if payload_kg is not None:
        csv_path = gravity_compensated_csv(csv_path, payload_kg, k_axis)
    out = {}
    for j in joints:
        conds = joint_conditions(j)
        s = k_axis[j] / KT_OUT
        per = {}
        for d, tag in ((+1, "pos"), (-1, "neg")):
            r = _quiet(vm.analyze_j3_velocity_modulation, csv_path, j, KT_OUT, d, conds, V_CRUISE[j])
            for ci, c in r["conditions"].items():
                J = c["J_variants"].get(variant, c["J"])
                per.setdefault(ci, {})[tag] = dict(J=J * s, se=c["se_J"] * s, B=c["B_local"] * s,
                                                   n_seg=c["n_seg"], gain=c["track_gain"],
                                                   Jm=r["model"].get(ci, {}).get("J_model"))
        out[j] = dict(conds=conds, per=per, K_used=k_axis[j])
    return out


def summary_lines(mode, res):
    lines = []
    if mode == "CONSTVEL":
        lines.append(f"{'軸':>3}{'B':>10}{'±':>7}{'Tc':>9}{'±':>7}{'H̄':>9}{'段數(正/反)':>13}  K")
        for j, r in res.items():
            if "B" not in r:
                lines.append(f"J{j}  段數不足（正 {r['n_pos']}、反 {r['n_neg']}）")
                continue
            lines.append(f"J{j} {r['B']:9.3f}{r['sd_B']:7.3f}{r['Tc']:9.3f}{r['sd_Tc']:7.3f}{r['H']:+9.3f}"
                         f"{r['n_pos']:>7}/{r['n_neg']:<5}  {r['K_used']}")
    else:
        for j, r in res.items():
            lines.append(f"J{j}（K {r['K_used']}）：")
            if not r["per"]:
                lines.append("   沒有結果（每個條件每個方向至少要 2 段；QUICK 只驗證動作，不算 J）")
            for ci in sorted(r["per"]):
                f, A = r["conds"][ci]
                p, n = r["per"][ci].get("pos", {}), r["per"][ci].get("neg", {})
                jm = p.get("Jm") or n.get("Jm")
                lines.append(f"   {f:.2f} Hz A={A:.3f}: J+ {p.get('J', float('nan')):.4f}±{p.get('se', float('nan')):.4f}"
                             f"  J- {n.get('J', float('nan')):.4f}±{n.get('se', float('nan')):.4f}"
                             f"  J_model {jm if jm is None else round(jm, 4)}")
    return lines


# ============================================================
# 主程式
# ============================================================

def main():
    print("=" * 70)
    print(f" UR5 六軸同時鑑別（雛形）：MODE = {MODE}，TEST_STAGE = {TEST_STAGE}")
    print("=" * 70)
    if POSE not in POSES:
        print(f"[中止] 未知的姿態 {POSE!r}"); return
    pose_q = list(POSES[POSE]["q_deg"])
    if MODE == "CONSTVEL":
        segs = plan_constvel(stage_config(MODE, TEST_STAGE), pose_q)
    elif MODE == "VELMOD":
        ci, nr = stage_config(MODE, TEST_STAGE)
        segs = plan_velmod(ci, nr, pose_q)
        for j in range(6):
            print(f"  J{j}: 巡航 {V_CRUISE[j]} rad/s，條件 {joint_conditions(j)}")
    else:
        print(f"[中止] 未知的 MODE {MODE!r}"); return
    print(f"共 {len(segs)} 段，預估 {estimate_total_time(segs, pose_q, MODE, pose_q) / 60:.1f} 分鐘")

    ok, report, (lo, hi) = safety_check(segs, MODE, pose_q, safe_min=SAFE_MIN_DEG, safe_max=SAFE_MAX_DEG)
    print("\n".join(report))
    p_ok, p_msg = vm.check_payload_setting(PAYLOAD_KG)
    print(p_msg + ("  [OK]" if p_ok else "  [FAIL]"))
    if not (ok and p_ok):
        print("\n[中止] 安全檢查未通過。"); return
    if ROBOT_IP is None or rtde_receive is None:
        print("\n[中止] ROBOT_IP 未設定或找不到 rtde_receive，僅完成離線設計檢查。"); return

    print(f"\n[注意] 六軸須先擺到 {pose_q}°（可用 ur5_goto_pose.py SIX）；六軸會同時在下列範圍內移動：")
    for j in range(6):
        print(f"   J{j}: {lo[j]:.1f}° ～ {hi[j]:.1f}°")
    print("       請確認手臂周圍淨空、手腕線材與相機線不會被拉扯；連上 RTDE 後會再以實際角度檢查六軸。")
    if input(f"已確認，輸入 {POSE} 繼續：").strip() != POSE:
        print("已取消。"); return
    e_ok, e_msg = vm.check_payload_entry(PAYLOAD_KG, input("輸入教導器上的負載質量 [kg]："))
    print(e_msg)
    if not e_ok:
        print("[中止]"); return

    tag = f"six_{MODE.lower()}_{POSE.lower()}_{vm.payload_tag(PAYLOAD_KG)}_{TEST_STAGE.lower()}_{datetime.now():%Y%m%d_%H%M%S}"
    out_dir = os.path.join(OUTPUT_DIR, tag)
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, tag + ".csv")
    logger = vm.RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, 0)
    logger.start()
    if not logger.started_ok.wait(timeout=10.0):
        print(f"[中止] RTDE 沒有資料：{logger.error}"); logger.stop(); return
    try:
        p_ok2, p_lines = vm.check_pose(logger.last_q_full, POSE, POSES, POSE_TOL_DEG)
        print("\n".join(p_lines))
        q_now = [math.degrees(x) for x in logger.last_q_full]
        s_ok, s_msg = check_start_position(q_now, segs[0]["start_deg"])
        print(s_msg)
        f_ok, f_msg = vm.check_speed_fraction(logger.last_speed_fraction)
        print(f_msg)
        if not (p_ok2 and s_ok and f_ok):
            print("[中止] 開始前檢查未通過。"); return
        script = build_urscript(segs, MODE, end_deg=pose_q)
        with open(os.path.join(out_dir, "program.script"), "w", encoding="utf-8") as fh:
            fh.write(script)
        planned = estimate_total_time(segs, q_now, MODE, pose_q)
        vm.send_urscript(ROBOT_IP, script)
        print(f"已送出程式，預估 {planned / 60:.1f} 分鐘；Ctrl+C 中止（會送 stopj）。")
        status, elapsed, moving = vm.wait_for_motion_complete(logger, planned)
        print(f"結束狀態：{status}，{elapsed:.1f} s（動作 {moving:.1f} s）")
        if status != "done":
            vm.send_abort(ROBOT_IP)
    except KeyboardInterrupt:
        vm.send_abort(ROBOT_IP)
        print("[中止] 使用者中斷。")
    finally:
        time.sleep(0.5)
        logger.stop()
        print(logger.summary())
    if not logger.speed_fraction_held():
        print("[警告] 執行中教導器速度滑桿低於 100% 或讀不到，速度、加速度可能不是程式設定值。")

    print("\n===== 離線分析")
    res = analyze_constvel(csv_path) if MODE == "CONSTVEL" else analyze_velmod(csv_path, payload_kg=float(PAYLOAD_KG))
    lines = summary_lines(MODE, res)
    print("\n".join(lines))
    with open(os.path.join(out_dir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n存檔：{out_dir}")


if __name__ == "__main__":
    main()
