"""
UR5 取放情境鑑別（雛形）—— 照工業取放的動作走，從動作中的定速段算各軸 B、Tc
================================================================

★★★ 想法（學長提議：最貼近工業使用情境）★★★

手臂反覆做「取 → 搬 → 放 → 搬回」：
  A 上方 →（movel 下降）A 下方 → 停（夾取）→（movel 上升）A 上方 →（movej 搬運）B 上方
  →（movel 下降）B 下方 → 停（放置）→（movel 上升）B 上方 →（movej 搬回）A 上方
每一個循環用不同的搬運速度（帶頭軸的 movej 速度 v），速度先升後降。

1. 搬運段（movej）是**關節空間的直線**：各軸同時起動、同時到達，速度是梯形（加速、定速、減速），
   定速段第 j 軸的速度 = v × |Δq_j| / |Δq_帶頭軸|。去程 A→B、回程 B→A 走同一條路、方向相反。
   → 和六軸定速一樣：每軸取同一個角度窗口（整數個馬達圈，漣波抵銷）內的定速樣本，
     重力只跟位置有關、科氏力和速度平方成正比，正反方向都相同；摩擦和方向同號 →
     正反轉分開回歸 tau = B·qd + C·方向，再合併（Tc = (C₊ + C₋)/2）。
   一個循環兩段（去、回），每段每軸一個點；FULL 有 8 檔速度 × 先升後降 → 每軸每方向 16 個點。
2. 取放段（movel 下降、上升）走工具空間直線，各軸速度不是定值 → 不進定速分析，
   只進「全週期回歸」（參考）：扣掉名目剛體模型後 tau = J_r·qdd + B·qd + Tc·sign(qd) + 重力修正。

★★★ 哪些軸看得到 ★★★
預設的 A、B（見 STATIONS，兩站同高、水平距離不同）：J0 轉 40°（帶頭軸）、J1 −19°、J2 +17°、J5 轉 30°（放置時轉向）。
工具全程垂直朝下 → J4 不動；兩站同高時 J3 只差 1.6°（J3 = −90° − J1 − J2，J1、J2 幾乎互相抵銷）。
→ 搬運段量得到 J0、J1、J2、J5；**J3、J4 量不到**（J3 只在 movel 取放段動，進全週期回歸）。
各軸的定速速度 = 帶頭軸速度 × |Δq_j| / 40°，所以 J1、J2 只到帶頭軸速度的約 0.45 倍。
定速窗口的寬度由最高速度時的加減速距離決定（transfer_windows），窗口不到 1 個馬達圈的軸不做定速分析。

★★★ 使用前必讀 ★★★
1. 【新程式，尚未上機】依 TEST_STAGE = QUICK → MID → FULL 執行。
2. 【取放位置】STATIONS 的 A、B 上方姿態是**預設值，不是現場教導的位置**：請在教導器上確認或重新教導，
   兩者都要工具朝下（J4 相同、J1+J2+J3 相同），下方姿態由程式沿垂直線往下 APPROACH_DZ 算出（IK）。
3. 【六軸視窗未設定一律中止】SAFE_MIN_DEG / SAFE_MAX_DEG 由操作者依現場填寫；
   所有動作（含下方姿態與 movel 路徑）外推 SAFETY_MARGIN_DEG 後都必須在視窗內。離線執行會印出規劃範圍。
4. 【法蘭高度、工具尖端高度與線速度】沿每一段路徑以官方 DH 檢查（同六軸程式的門檻）。
5. 【沒有真的夾東西】取、放只停 DWELL_TIME 秒，負載全程不變（PAYLOAD_KG）；
   真的夾放時去回兩段的負載不同，正反轉相消就不成立（要改成來回都夾著同一個物件）。
6. 【開始位置】第一個動作是慢速 movej 到 A 上方；任何一軸距離超過 MAX_INITIAL_MOVE_DEG 就中止。
7. 【教導器速度滑桿】不是 100% 或讀不到一律中止。負載質量執行時要再輸入一次。
8. 【IP 需要你自己再次確認】（第二台 UR5 只改 ROBOT_IP；裸機 PAYLOAD_KG = 0）
"""

import math
import os
import sys
import time
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ur5_six_axis_ident as six             # 安全檢查門檻、FK、重力模型、K_AXIS
vm = six.vm                                   # RTDE 記錄、URScript 傳送、檢查函式


# ============================================================
# 使用者設定
# ============================================================

ROBOT_IP = "192.168.50.114"   # 執行前請再次確認
# QUICK: 一個速度（0.05）一個循環 —— 驗證取放動作、IK、movel、時間
# MID  : 三個速度（0.05、0.12、0.20）各一個循環
# FULL : 八個速度先升後降（每個速度兩個循環）
TEST_STAGE = "QUICK"

# 取放位置：上方姿態 [deg]（工具朝下：J4 = −90、J1 + J2 + J3 = −90，同 SIX 姿態）
#   預設：兩站法蘭高度都是 0.60 m，A 水平距離 0.50 m、B 0.36 m（取料、放料在同一個工作台高度）
STATIONS = {"A": (160.0, -84.63, 58.94, -64.31, -90.0, -15.0),
            "B": (200.0, -103.62, 76.29, -62.67, -90.0, 15.0)}
APPROACH_DZ = 0.08            # 下降距離 [m]（垂直）
APPROACH_V = 0.05             # movel 工具速度 [m/s]
APPROACH_A = 0.30             # movel 工具加速度 [m/s²]
DWELL_TIME = 0.5              # 取、放停留時間 [s]（代替夾爪動作）
PAUSE_TIME = 0.3              # 每個動作之間的停頓 [s]

TRANSFER_SPEEDS = [0.03, 0.05, 0.08, 0.12, 0.16, 0.20, 0.25, 0.30]   # 搬運 movej 帶頭軸速度 [rad/s]
TRANSFER_ACCEL = 0.8          # 搬運 movej 帶頭軸加速度 [rad/s²]
SETTLE_TIME = 0.25            # 加速結束後到進窗口的穩定時間 [s]
V_MIN_PLATEAU = 0.003         # 定速段最低速度 [rad/s]（比這低的軸、段不用）

START_MOVEJ_VEL = 0.15        # 開始時慢速 movej 到 A 上方
START_MOVEJ_ACCEL = 0.5

GEAR_RATIO = six.GEAR_RATIO
SAMPLE_HZ = six.SAMPLE_HZ
DT = six.DT
K_AXIS = six.K_AXIS           # 電流 → 力矩換算比例（同六軸程式）
QD_MAX = six.QD_MAX
SAFETY_MARGIN_DEG = six.SAFETY_MARGIN_DEG
MAX_INITIAL_MOVE_DEG = six.MAX_INITIAL_MOVE_DEG
FLANGE_Z_MIN_M = six.FLANGE_Z_MIN_M
TOOL_TIP_Z_MIN_M = six.TOOL_TIP_Z_MIN_M
FLANGE_SPEED_MAX = six.FLANGE_SPEED_MAX

SAFE_MIN_DEG = (None,) * 6    # 六軸允許的絕對角度視窗：未填一律中止（離線執行會印出規劃範圍供參考）
SAFE_MAX_DEG = (None,) * 6

PAYLOAD_KG = six.PAYLOAD_KG   # 1.30（FT300、夾爪、3D 相機）；教導器負載須相同
FULL_CYCLE_VMIN = 0.01        # 全週期回歸只用 |qd| 大於此值的樣本 [rad/s]

OUTPUT_DIR = six.OUTPUT_DIR


# ============================================================
# 第一部分：幾何（IK）與規劃
# ============================================================

def flange_pos(q_deg):
    return six.fk_frames([math.radians(x) for x in q_deg])[6][:3, 3].copy()


def lowered_pose(q_up_deg, dz, tol=1e-9, max_iter=50):
    """
    上方姿態沿垂直線往下 dz [m] 的姿態：J0、J4、J5 不變，J3 = (J1+J2+J3)_up − J1 − J2（法蘭姿態不變），
    解 J1、J2 使法蘭水平距離不變、高度降 dz（數值 Newton）。回傳 [deg]。
    """
    q = np.radians(np.array(q_up_deg, float))
    s0 = q[1] + q[2] + q[3]
    p0 = six.fk_frames(q)[6][:3, 3]
    tgt = np.array([math.hypot(p0[0], p0[1]), p0[2] - dz])

    def f(x):
        qq = q.copy(); qq[1], qq[2] = x; qq[3] = s0 - x[0] - x[1]
        p = six.fk_frames(qq)[6][:3, 3]
        return np.array([math.hypot(p[0], p[1]), p[2]]) - tgt, qq

    x = q[1:3].copy()
    for _ in range(max_iter):
        r, qq = f(x)
        if np.max(np.abs(r)) < tol:
            return [math.degrees(v) for v in qq]
        Jm = np.column_stack([(f(x + h)[0] - r) / 1e-7 for h in np.eye(2) * 1e-7])
        x = x - np.linalg.solve(Jm, r)
    raise RuntimeError("lowered_pose：IK 未收斂（下降距離太大或接近奇異姿態）")


def stage_speeds(stage):
    if stage == "QUICK":
        return [0.05]
    if stage == "MID":
        return [0.05, 0.12, 0.20]
    if stage == "FULL":
        return list(TRANSFER_SPEEDS) + list(reversed(TRANSFER_SPEEDS))
    raise ValueError(f"未知的 TEST_STAGE: {stage}")


def movej_profile(dist_rad, v, a):
    """梯形（或三角形）速度：回傳 (總時間, 加速時間, 實際最高速度)。"""
    if dist_rad <= 0:
        return 0.0, 0.0, 0.0
    if v * v / a >= dist_rad:                 # 達不到 v：三角形
        vp = math.sqrt(dist_rad * a)
        return 2 * vp / a, vp / a, vp
    return dist_rad / v + v / a, v / a, v


def movel_time(dist_m, v=APPROACH_V, a=APPROACH_A):
    return movej_profile(dist_m, v, a)[0]


def plan_pick_place(speeds, stations=STATIONS, dz=APPROACH_DZ, accel=TRANSFER_ACCEL):
    """
    動作清單（依執行順序）。每個循環從 A 上方開始：A 取（下、停、上）→ 搬 A→B → B 放（下、停、上）→ 搬 B→A。
    movej 段記錄帶頭軸速度 v、方向 d（A→B 為 +1）與每軸定速段速度 qd_plateau。
    """
    A, B = list(stations["A"]), list(stations["B"])
    Ad, Bd = lowered_pose(A, dz), lowered_pose(B, dz)
    moves = []

    def add(kind, start, target, **kw):
        moves.append(dict(index=len(moves), kind=kind, start_deg=list(start), target_deg=list(target), **kw))

    def station(name, up, down):
        add("movel", up, down, label=f"{name} down", v=APPROACH_V, a=APPROACH_A, t=movel_time(dz))
        add("dwell", down, down, label=f"{name} {'pick' if name == 'A' else 'place'}", t=DWELL_TIME)
        add("movel", down, up, label=f"{name} up", v=APPROACH_V, a=APPROACH_A, t=movel_time(dz))

    for k, v in enumerate(speeds):
        for d, (s_name, s_up, s_dn, e_up) in ((+1, ("A", A, Ad, B)), (-1, ("B", B, Bd, A))):
            station(s_name, s_up, s_dn)
            dq = [math.radians(e - s) for s, e in zip(s_up, e_up)]
            lead = max(abs(x) for x in dq)
            t_tot, t_acc, vp = movej_profile(lead, v, accel)
            add("movej", s_up, e_up, label=f"cycle {k + 1} {'A->B' if d > 0 else 'B->A'}", cycle=k, d=d, v=v, a=accel,
                v_peak=vp, t=t_tot, t_acc=t_acc, qd_plateau=[vp * x / lead for x in dq])
    return moves, dict(A_up=A, A_down=Ad, B_up=B, B_down=Bd)


def transfer_windows(stations=STATIONS, speeds=TRANSFER_SPEEDS, accel=TRANSFER_ACCEL, settle=SETTLE_TIME,
                     gear_ratio=GEAR_RATIO, v_min=V_MIN_PLATEAU):
    """
    每軸的定速窗口 (lo, hi) [deg]，或 None（量不到）。
    movej 各軸同步：路徑參數 λ ∈ [0, 1]，最高速度 v_max 時加速＋穩定佔 f = (v²/2a + v·settle)/|Δq_帶頭|，
    所有速度的定速段都涵蓋 λ ∈ [f, 1 − f]。窗口 = 中心 ± 整數個馬達圈的一半，須落在這段內。
    最低速度時該軸定速段速度 < v_min 也照樣給窗口（分析時該段會被略過）。
    """
    A, B = stations["A"], stations["B"]
    dq = [b - a for a, b in zip(A, B)]
    lead = math.radians(max(abs(x) for x in dq))
    vmax = max(speeds)
    _, t_acc, vp = movej_profile(lead, vmax, accel)
    f = (vp * vp / (2 * accel) + vp * settle) / lead
    rev = 360.0 / gear_ratio
    out = []
    for j in range(6):
        usable = abs(dq[j]) * (1 - 2 * f)
        n = int(math.floor(usable / rev + 1e-9))
        if n < 1 or vmax * abs(dq[j]) / math.degrees(lead) < v_min:
            out.append(None); continue
        c = 0.5 * (A[j] + B[j])
        out.append((c - n * rev / 2.0, c + n * rev / 2.0))
    return out


def estimate_total_time(moves, q0_deg=None):
    total = 0.0
    if q0_deg is not None:
        total += vm.estimate_movej_time(math.radians(max(abs(a - b) for a, b in zip(q0_deg, moves[0]["start_deg"]))),
                                        START_MOVEJ_VEL, START_MOVEJ_ACCEL) + PAUSE_TIME
    for m in moves:
        total += m["t"] + (PAUSE_TIME if m["kind"] != "dwell" else 0.0)
    return total


# ============================================================
# 第二部分：安全檢查
# ============================================================

def path_samples(m, n=21):
    """沿一個動作的路徑取樣 [deg]：movej 關節空間直線；movel 工具垂直線（以 IK 沿線求解）。"""
    a, b = np.array(m["start_deg"]), np.array(m["target_deg"])
    if m["kind"] != "movel":
        return [a + lam * (b - a) for lam in np.linspace(0.0, 1.0, n)]
    up, down = (a, b) if flange_pos(a)[2] > flange_pos(b)[2] else (b, a)
    dz = flange_pos(up)[2] - flange_pos(down)[2]
    pts = [np.array(lowered_pose(up, lam * dz)) if lam > 0 else up for lam in np.linspace(0.0, 1.0, n)]
    return pts if (a == up).all() else pts[::-1]


def max_joint_speed_movel(m, v=APPROACH_V, n=21):
    """movel 以工具速度 v 走時，各軸最大速度 [rad/s]（沿線數值微分）。"""
    pts = np.radians(np.array(path_samples(m, n)))
    dz = abs(flange_pos(m["start_deg"])[2] - flange_pos(m["target_deg"])[2])
    dqds = np.abs(np.diff(pts, axis=0)) / (dz / (n - 1))
    return dqds.max(axis=0) * v


def safety_check(moves, poses, safe_min=SAFE_MIN_DEG, safe_max=SAFE_MAX_DEG, margin_deg=SAFETY_MARGIN_DEG,
                 flange_z_min=FLANGE_Z_MIN_M, tip_z_min=TOOL_TIP_Z_MIN_M, flange_speed_max=FLANGE_SPEED_MAX,
                 qd_max=QD_MAX, start_vel=START_MOVEJ_VEL):
    """
    回傳 (ok, lines, (lo, hi))。
    1. 各軸速度：movej 定速段、movel（工具速度換算）、開始的慢速 movej，都 ≤ qd_max。
    2. 六軸絕對角度：所有路徑取樣點外推 margin_deg，須在 SAFE_MIN/MAX 內；未設定一律 FAIL。
    3. 每段路徑上的法蘭高度、工具尖端高度；movej 以定速段速度、movel 以工具速度算法蘭線速度。
    4. 取放兩站的上方姿態工具朝下一致（J4 相同、J1+J2+J3 相同），下方姿態 IK 有解。
    """
    lines, ok = [], True
    A, B = poses["A_up"], poses["B_up"]
    if abs(A[4] - B[4]) > 1e-6 or abs((A[1] + A[2] + A[3]) - (B[1] + B[2] + B[3])) > 1e-6:
        ok = False; lines.append("[FAIL] A、B 上方姿態的工具方向不同（J4 或 J1+J2+J3 不同）")
    vj = max([max(abs(x) for x in m["qd_plateau"]) for m in moves if m["kind"] == "movej"] + [start_vel])
    vl = max([float(max_joint_speed_movel(m).max()) for m in moves if m["kind"] == "movel"] + [0.0])
    good = max(vj, vl) <= qd_max
    ok &= good
    lines.append(f"[{'OK' if good else 'FAIL'}] 最高關節速度：movej {vj:.3f}、movel {vl:.3f} rad/s（上限 {qd_max}）")
    pts, zmin, tmin, vfl = [], 9.0, 9.0, 0.0
    for m in moves:
        if m["kind"] == "dwell":
            continue
        ps = path_samples(m)
        pts += ps
        qd = (np.array(m["qd_plateau"]) if m["kind"] == "movej" else None)
        for qq in ps:
            z, tp = six.flange_and_tip_z(qq)
            zmin, tmin = min(zmin, z), min(tmin, tp)
            if qd is not None:
                vfl = max(vfl, six.flange_speed(qq, qd))
    vfl = max(vfl, APPROACH_V)
    lo = [min(p[j] for p in pts) - margin_deg for j in range(6)]
    hi = [max(p[j] for p in pts) + margin_deg for j in range(6)]
    for j in range(6):
        if safe_min[j] is None or safe_max[j] is None:
            ok = False; lines.append(f"[FAIL] J{j} 視窗未設定（規劃範圍含餘裕 {lo[j]:.2f}°～{hi[j]:.2f}°）")
        elif lo[j] < safe_min[j] or hi[j] > safe_max[j]:
            ok = False; lines.append(f"[FAIL] J{j} 規劃範圍 {lo[j]:.2f}°～{hi[j]:.2f}° 超出視窗 {safe_min[j]}°～{safe_max[j]}°")
        else:
            lines.append(f"[OK] J{j} 規劃範圍 {lo[j]:.2f}°～{hi[j]:.2f}°（視窗 {safe_min[j]}°～{safe_max[j]}°）")
    for name, val, lim, cmp in (("法蘭高度", zmin, flange_z_min, ">="), ("工具尖端高度", tmin, tip_z_min, ">="),
                                ("法蘭線速度", vfl, flange_speed_max, "<=")):
        good = (val >= lim) if cmp == ">=" else (val <= lim)
        ok &= bool(good)
        lines.append(f"[{'OK' if good else 'FAIL'}] {name} {val:.3f}（須 {cmp} {lim}）")
    return ok, lines, (lo, hi)


# ============================================================
# 第三部分：URScript
# ============================================================

def build_urscript(moves, poses, start_vel=START_MOVEJ_VEL, start_accel=START_MOVEJ_ACCEL, pause=PAUSE_TIME):
    """只含 ASCII。movel 以關節角度當目標（控制器用 FK 換成工具位姿，沿工具空間直線走）。無 blend（r = 0）。"""
    v6 = lambda q: six._vec([math.radians(x) for x in q])
    lines = ["def pick_place_ident_program():",
             "  # go to A (above)",
             f"  movej({v6(poses['A_up'])}, a={start_accel}, v={start_vel})",
             f"  sleep({pause})"]
    for m in moves:
        lines.append(f"  # move {m['index'] + 1}: {m['label']}")
        if m["kind"] == "movel":
            lines.append(f"  movel({v6(m['target_deg'])}, a={m['a']}, v={m['v']})")
            lines.append(f"  sleep({pause})")
        elif m["kind"] == "dwell":
            lines.append(f"  sleep({m['t']})")
        else:
            lines.append(f"  movej({v6(m['target_deg'])}, a={m['a']}, v={m['v']:.4f})")
            lines.append(f"  sleep({pause})")
    lines.append("end")
    return "\n".join(lines) + "\n"


# ============================================================
# 第四部分：離線分析
# ============================================================

def load_csv(csv_path):
    import pandas as pd
    return pd.read_csv(csv_path)


def plateau_segments(D, j, window, v_min=V_MIN_PLATEAU, tol=1e-4, min_n=10):
    """
    第 j 軸在窗口 (lo, hi) [deg] 內的定速段：target_qd 連續保持同一個值（|差| < tol）且 |target_qd| ≥ v_min，
    而且這一段的實際角度要完整涵蓋窗口（窗口是整數個馬達圈，漣波才會在平均中抵銷）。
    回傳 [dict(v=目標速度, d=方向, qd_mean, tau_mean, n, t0)]，tau = 電流 × K_AXIS[j]。
    """
    lo, hi = window
    tqd = D[f"target_qd_{j}"].to_numpy()
    q = np.degrees(D[f"actual_q_{j}"].to_numpy())
    qd = D[f"actual_qd_{j}"].to_numpy()
    tau = D[f"actual_current_{j}"].to_numpy() * K_AXIS[j]
    t = D["timestamp"].to_numpy()
    flat = np.r_[False, np.abs(np.diff(tqd)) < tol] & (np.abs(tqd) >= v_min)
    segs = []
    idx = np.flatnonzero(flat)
    if len(idx) == 0:
        return segs
    for r in np.split(idx, np.flatnonzero(np.diff(idx) != 1) + 1):
        if len(r) < min_n:
            continue
        qq = q[r]
        if qq.min() > lo or qq.max() < hi:            # 沒有完整涵蓋窗口
            continue
        w = r[(qq >= lo) & (qq <= hi)]
        if len(w) < min_n:
            continue
        segs.append(dict(v=float(np.median(tqd[r])), d=int(np.sign(np.median(tqd[r]))), qd_mean=float(qd[w].mean()),
                         tau_mean=float(tau[w].mean()), n=len(w), t0=float(t[w[0]])))
    return segs


def _fit_direction(segs, d):
    x = np.array([s["qd_mean"] for s in segs]); y = np.array([s["tau_mean"] for s in segs])
    X = np.column_stack([x, np.full(len(x), float(d))])
    c, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ c
    dof = max(len(y) - 2, 1)
    sd = np.sqrt(np.diag(res @ res / dof * np.linalg.inv(X.T @ X)))
    return float(c[0]), float(c[1]), float(sd[0]), float(sd[1]), float(np.sqrt(res @ res / dof))


def analyze_transfer(csv_path_or_df, windows, joints=range(6), min_segs=3):
    """
    搬運段定速分析：每軸正反方向各自回歸 tau = B·qd + C·方向，再合併
      B = (B₊ + B₋)/2，Tc = (C₊ + C₋)/2，H̄ = (C₊ − C₋)/2（窗口平均重力保持力矩）。
    回傳 {j: dict(B, sd_B, Tc, sd_Tc, H, B_pos, B_neg, n_pos, n_neg, v_range, rms)} 或 {j: dict(skip=原因)}。
    """
    D = load_csv(csv_path_or_df) if isinstance(csv_path_or_df, str) else csv_path_or_df
    out = {}
    for j in joints:
        if windows[j] is None:
            out[j] = dict(skip="沒有定速窗口（取放中幾乎不動，或行程不到 1 個馬達圈）"); continue
        segs = plateau_segments(D, j, windows[j])
        sp, sn = [s for s in segs if s["d"] > 0], [s for s in segs if s["d"] < 0]
        if len(sp) < min_segs or len(sn) < min_segs:
            out[j] = dict(skip=f"定速段不足（正 {len(sp)}、反 {len(sn)}，各需 ≥ {min_segs}）", n_pos=len(sp), n_neg=len(sn))
            continue
        bp, cp, sbp, scp, rp = _fit_direction(sp, +1)
        bn, cn, sbn, scn, rn = _fit_direction(sn, -1)
        vs = [abs(s["qd_mean"]) for s in segs]
        out[j] = dict(B=0.5 * (bp + bn), sd_B=0.5 * math.hypot(sbp, sbn), Tc=0.5 * (cp + cn), sd_Tc=0.5 * math.hypot(scp, scn),
                      H=0.5 * (cp - cn), B_pos=bp, B_neg=bn, n_pos=len(sp), n_neg=len(sn), v_range=(min(vs), max(vs)),
                      rms=0.5 * (rp + rn), segments=segs)
    return out


LINK_M, LINK_COM, TOOL_COG = six.LINK_M, six.LINK_COM, six.TOOL_COG


def rigid_body_torque(q, qd, qdd, payload_kg, gravity=True, g=9.81):
    """名目剛體逆動力學（Newton–Euler；官方連桿質量、重心當點質量 + 負載點質量；不含轉子、摩擦）[N·m]。"""
    fr = six.fk_frames(q)
    z = [f[:3, 2] for f in fr]; p = [f[:3, 3] for f in fr]
    w, dw = np.zeros(3), np.zeros(3)
    a = np.array([0.0, 0.0, g if gravity else 0.0])
    Fl, Ml = [], []
    for i in range(6):
        w_new = w + z[i] * qd[i]
        dw = dw + z[i] * qdd[i] + np.cross(w, z[i] * qd[i])
        w = w_new
        R = fr[i + 1][:3, :3]
        bodies = [(LINK_M[i], np.array(LINK_COM[i]))] + ([(payload_kg, np.array(TOOL_COG))] if i == 5 else [])
        F, M = np.zeros(3), np.zeros(3)
        for m, com in bodies:
            r = R @ com + p[i + 1] - p[i]
            f = m * (a + np.cross(dw, r) + np.cross(w, np.cross(w, r)))
            F += f; M += np.cross(r, f)
        Fl.append(F); Ml.append(M)
        dd = p[i + 1] - p[i]
        a = a + np.cross(dw, dd) + np.cross(w, np.cross(w, dd))
    tau, f, n = np.zeros(6), np.zeros(3), np.zeros(3)
    for i in range(5, -1, -1):
        n = Ml[i] + n + np.cross(p[i + 1] - p[i], f)
        f = Fl[i] + f
        tau[i] = n @ z[i]
    return tau


def analyze_full_cycle(csv_path_or_df, payload_kg, joints=range(6), v_min=FULL_CYCLE_VMIN, step=1):
    """
    全週期回歸（參考，含 movel 取放段與加減速）：
      y_j = 電流 × K_j − 名目剛體力矩_j(q, qd, target_qdd)
      y_j = J_r·qdd_j + B·qd_j + Tc·sign(qd_j) + c₀ + c₁·g_j(q)
    J_r：轉子＋對角慣量誤差；c₀、c₁·g_j：重力模型的偏差。只用 |qd_j| > v_min 的樣本（離開低速、靜止區）。
    假設摩擦對速度是直線 —— 手腕在低速是彎的（11/26），結果只當參考，和定速分析比。
    """
    D = load_csv(csv_path_or_df) if isinstance(csv_path_or_df, str) else csv_path_or_df
    D = D.iloc[::step].reset_index(drop=True)
    Q = D[[f"actual_q_{j}" for j in range(6)]].to_numpy()
    QD = D[[f"actual_qd_{j}" for j in range(6)]].to_numpy()
    if all(f"target_qdd_{j}" in D for j in range(6)) and D[[f"target_qdd_{j}" for j in range(6)]].notna().all().all():
        QDD = D[[f"target_qdd_{j}" for j in range(6)]].to_numpy()
    else:
        QDD = np.gradient(D[[f"target_qd_{j}" for j in range(6)]].to_numpy(), D["timestamp"].to_numpy(), axis=0)
    moving = np.any(np.abs(QD) > v_min, axis=1)
    TRB = np.zeros_like(Q); G = np.zeros_like(Q)
    for k in np.flatnonzero(moving):
        TRB[k] = rigid_body_torque(Q[k], QD[k], QDD[k], payload_kg)
        G[k] = six.gravity_torque(Q[k], payload_kg)
    out = {}
    for j in joints:
        m = np.abs(QD[:, j]) > v_min
        if m.sum() < 50:
            out[j] = dict(skip=f"動作樣本不足（{int(m.sum())}）"); continue
        y = D[f"actual_current_{j}"].to_numpy()[m] * K_AXIS[j] - TRB[m, j]
        X = np.column_stack([QDD[m, j], QD[m, j], np.sign(QD[m, j]), np.ones(m.sum()), G[m, j]])
        c, *_ = np.linalg.lstsq(X, y, rcond=None)
        res = y - X @ c
        dof = max(len(y) - X.shape[1], 1)
        sd = np.sqrt(np.diag(res @ res / dof * np.linalg.pinv(X.T @ X)))
        out[j] = dict(J_r=float(c[0]), sd_J_r=float(sd[0]), B=float(c[1]), sd_B=float(sd[1]), Tc=float(c[2]),
                      sd_Tc=float(sd[2]), c0=float(c[3]), c1=float(c[4]), n=int(m.sum()), rms=float(np.sqrt(res @ res / dof)))
    return out


def summary_lines(res_t, res_f, windows):
    lines = ["搬運段定速分析（主要結果；± 為以段為點的回歸標準差）",
             f"{'軸':>3}{'窗口 [deg]':>20}{'B':>9}{'±':>7}{'Tc':>8}{'±':>7}{'H̄':>8}{'段數(正/反)':>12}{'速度範圍 [rad/s]':>18}"]
    for j in range(6):
        r = res_t.get(j, {})
        w = windows[j]
        ws = f"{w[0]:.1f}～{w[1]:.1f}" if w else "—"
        if "B" not in r:
            lines.append(f"J{j} {ws:>19}  {r.get('skip', '未分析')}")
            continue
        lines.append(f"J{j} {ws:>19}{r['B']:9.3f}{r['sd_B']:7.3f}{r['Tc']:8.3f}{r['sd_Tc']:7.3f}{r['H']:+8.3f}"
                     f"{r['n_pos']:>7}/{r['n_neg']:<4}{r['v_range'][0]:>10.3f}～{r['v_range'][1]:.3f}")
    if res_f:
        lines += ["", "全週期回歸（參考；扣名目剛體模型，tau = J_r·qdd + B·qd + Tc·sign(qd) + 重力修正）",
                  f"{'軸':>3}{'J_r':>9}{'±':>7}{'B':>9}{'±':>7}{'Tc':>8}{'±':>7}{'樣本':>8}{'rms':>8}"]
        for j in range(6):
            r = res_f.get(j, {})
            if "B" not in r:
                lines.append(f"J{j}  {r.get('skip', '未分析')}"); continue
            lines.append(f"J{j} {r['J_r']:9.3f}{r['sd_J_r']:7.3f}{r['B']:9.3f}{r['sd_B']:7.3f}{r['Tc']:8.3f}{r['sd_Tc']:7.3f}"
                         f"{r['n']:8d}{r['rms']:8.3f}")
    return lines


# ============================================================
# 主程式
# ============================================================

def main():
    print("=" * 70)
    print(f" UR5 取放情境鑑別（雛形）：TEST_STAGE = {TEST_STAGE}")
    print("=" * 70)
    moves, poses = plan_pick_place(stage_speeds(TEST_STAGE))
    windows = transfer_windows()
    for k in ("A_up", "A_down", "B_up", "B_down"):
        p = flange_pos(poses[k])
        print(f"  {k:7s} {[round(x, 2) for x in poses[k]]}  法蘭 ({p[0]:+.3f}, {p[1]:+.3f}, {p[2]:.3f}) m")
    dq = [b - a for a, b in zip(poses["A_up"], poses["B_up"])]
    print("  搬運各軸行程 [deg]：" + "  ".join(f"J{j} {dq[j]:+.1f}" for j in range(6)))
    print("  定速窗口：" + "  ".join(f"J{j} " + (f"{w[1] - w[0]:.2f}°（{round((w[1] - w[0]) / (360 / GEAR_RATIO))} 圈）"
                                              if w else "—") for j, w in enumerate(windows)))
    n_mv = sum(m["kind"] == "movej" for m in moves)
    print(f"共 {len(moves)} 個動作（搬運 {n_mv} 段），預估 {estimate_total_time(moves) / 60:.1f} 分鐘")

    ok, report, (lo, hi) = safety_check(moves, poses)
    print("\n".join(report))
    p_ok, p_msg = vm.check_payload_setting(PAYLOAD_KG)
    print(p_msg + ("  [OK]" if p_ok else "  [FAIL]"))
    if not (ok and p_ok):
        print("\n[中止] 安全檢查未通過。"); return
    if ROBOT_IP is None or vm.rtde_receive is None:
        print("\n[中止] ROBOT_IP 未設定或找不到 rtde_receive，僅完成離線設計檢查。"); return

    print("\n[注意] 手臂會在 A、B 兩站之間反覆取放，六軸範圍：")
    for j in range(6):
        print(f"   J{j}: {lo[j]:.1f}° ～ {hi[j]:.1f}°")
    print(f"       下方姿態比上方低 {APPROACH_DZ * 100:.0f} cm：請確認兩站下方淨空（不會撞到桌面、治具），線材不會被拉扯。")
    if input("已確認，輸入 PICK 繼續：").strip() != "PICK":
        print("已取消。"); return
    e_ok, e_msg = vm.check_payload_entry(PAYLOAD_KG, input("輸入教導器上的負載質量 [kg]："))
    print(e_msg)
    if not e_ok:
        print("[中止]"); return

    tag = f"pick_place_{vm.payload_tag(PAYLOAD_KG)}_{TEST_STAGE.lower()}_{datetime.now():%Y%m%d_%H%M%S}"
    out_dir = os.path.join(OUTPUT_DIR, tag)
    os.makedirs(out_dir, exist_ok=True)
    csv_path = os.path.join(out_dir, tag + ".csv")
    logger = vm.RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, 1)    # 以 J1 判斷結束：取放、搬運時 J1 都在動
    logger.start()
    if not logger.started_ok.wait(timeout=10.0):
        print(f"[中止] RTDE 沒有資料：{logger.error}"); logger.stop(); return
    try:
        q_now = [math.degrees(x) for x in logger.last_q_full]
        s_ok, s_msg = six.check_start_position(q_now, poses["A_up"])
        print(s_msg)
        f_ok, f_msg = vm.check_speed_fraction(logger.last_speed_fraction)
        print(f_msg)
        if not (s_ok and f_ok):
            print("[中止] 開始前檢查未通過。"); return
        script = build_urscript(moves, poses)
        with open(os.path.join(out_dir, "program.script"), "w", encoding="utf-8") as fh:
            fh.write(script)
        planned = estimate_total_time(moves, q_now)
        vm.send_urscript(ROBOT_IP, script)
        print(f"已送出程式，預估 {planned / 60:.1f} 分鐘；Ctrl+C 中止（會送 stopj）。")
        status, elapsed, moving = vm.wait_for_motion_complete(logger, planned, quiet_s=max(DWELL_TIME, PAUSE_TIME) + 1.5)
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
    D = load_csv(csv_path)
    res_t = analyze_transfer(D, windows)
    res_f = analyze_full_cycle(D, float(PAYLOAD_KG), step=2) if TEST_STAGE != "QUICK" else {}
    lines = summary_lines(res_t, res_f, windows)
    print("\n".join(lines))
    with open(os.path.join(out_dir, "summary.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\n存檔：{out_dir}")


if __name__ == "__main__":
    main()
