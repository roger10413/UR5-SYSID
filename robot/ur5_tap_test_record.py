"""
UR5 敲擊測試 —— 手臂靜止時輕敲工具，記錄餘振，找出手臂（含底座）的固有頻率
================================================================

★★★ 這支程式只記錄，不送任何指令給手臂 ★★★

只用 RTDE 讀取（rtde_receive）。程式裡沒有任何運動指令，也不開送指令用的連線，手臂不會因為這支程式而動。
記錄時手臂停在指定姿態（沒有程式在跑，伺服保持位置），操作者用手掌或橡膠槌輕敲工具。

目的
----
速度調變量到的 J 隨頻率上升（J0 各姿態與 J1 都有，11/10 分析）。控制器回報的角度是減速機輸出側，
所以該關節自己的減速機彈性不會改變「力矩 ÷ 輸出側加速度」；剩下的可能之一是輸出側編碼器之後還有彈性
（連桿、外側關節、底座或桌子），手臂末端晃得比編碼器看到的多。要解釋 3 Hz 多 12～19%，
固有頻率要在 7～9 Hz 左右。敲擊後的餘振可以直接讀出固有頻率：有沒有 7～9 Hz 的峰。

做法
----
1. 以教導器擺好姿態（POSE），確認沒有程式在跑、手臂靜止；教導器負載照實際設定。
2. 執行程式：檢查姿態（±0.5°）、手臂靜止、機器人在一般模式（剎車放開、伺服保持）、安全狀態正常。
3. 倒數後開始記錄 RECORD_S 秒。畫面在 TAP_TIMES 的時間點提示「敲」，每次提示敲一下，敲完手離開：
     TAP_LABEL = "fb"：前後敲（J1 轉動時工具移動的方向）
                 "lr"：左右敲（J0 轉動時工具移動的方向）
                 "ud"：上下敲
                 "none"：不敲，量背景（和有敲的比較）
   敲的位置：手腕 3 的金屬外殼或工具法蘭側面。不要敲 FT300、夾爪手指、相機。
   力道：手臂看得出輕微晃動即可。太大力會觸發保護性停止（教導器跳出視窗），那次資料不可用，按掉後輕一點重來。
4. 記錄：六軸角度、速度、電流、目標角度、控制器的目標力矩、目標電流、輸出電流、關節溫度，
   工具端加速度計（actual_tool_accelerometer，工具座標 x、y、z）、TCP 力估計值，以及機器人模式、安全模式。
5. 結束後自動分析（analyze_tap_test）：以工具加速度計找出每一下敲擊，取敲擊後的餘振做頻譜，
   各下平均後列出各訊號 3～45 Hz 的峰值頻率。

輸出：outputs/taptest_{姿態}_p{負載}_{敲的方向}_{日期}_{時間}/ 內的 CSV 與摘要。

★★★ 使用前必讀 ★★★
1. 【新程式，尚未上機】第一次先用 TAP_LABEL = "none" 跑一次（不敲），確認記錄正常、加速度計讀得到。
2. 【沒有程式在跑】執行時教導器上不能有程式在執行；手臂任一軸速度超過 STILL_QD_MAX 即中止。
3. 【姿態】程式不會擺姿態：先以教導器擺好 POSE，連上 RTDE 後以實際角度檢查六軸（±0.5°），不符即中止；
   執行前須輸入姿態名稱。
4. 【負載 PAYLOAD_KG】同其他程式：填教導器上設定的負載質量（只有末端工具時 1.30 kg），未設定一律中止；
   執行前須輸入教導器上的負載質量，不符即中止。負載會改變固有頻率，所以寫進資料夾名稱。
5. 【機器人模式】須為一般模式（robot_mode = 7）、安全模式正常（safety_mode = 1），讀不到一律中止。
   剎車鎖住時量到的不是伺服保持的狀態。
6. 【人站的位置】敲的時候站在手臂側面，手敲完就離開工具；另一人在緊急停止按鈕旁。
7. 【IP 需要你自己再次確認】
"""

import csv
import math
import os
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

# 敲的方向（只用於記錄與資料夾名稱，程式不會動手臂）
TAP_LABEL = "none"
TAP_LABELS = {
    "fb": "前後敲（J1 轉動時工具移動的方向）",
    "lr": "左右敲（J0 轉動時工具移動的方向）",
    "ud": "上下敲",
    "none": "不敲，量背景",
}

RECORD_S = 14.0                    # 記錄時間 [s]
TAP_TIMES = (2.0, 5.0, 8.0, 11.0)  # 畫面提示「敲」的時間點 [s]：每下間隔 3 s，留給餘振
COUNTDOWN_S = 3                    # 開始記錄前的倒數 [s]
RECORD_S_MAX = 120.0

SAMPLE_HZ = 125.0
DT = 1.0 / SAMPLE_HZ

# 姿態：六軸角度 [deg]（J0 = 180°）。程式不動手臂，只檢查實際角度是否在這個姿態（±POSE_TOL_DEG）。
POSE = "UP"
POSES = {
    "UP": dict(q_deg=(180.0, -90.0, 0.0, -90.0, -90.0, 0.0)),      # 手臂朝上
    "DIAG": dict(q_deg=(180.0, -90.0, 90.0, -90.0, -90.0, 0.0)),   # 診斷原點
    "EXT45": dict(q_deg=(180.0, -45.0, 45.0, -90.0, -90.0, 0.0)),  # 手臂外伸
}
POSE_TOL_DEG = 0.5

STILL_QD_MAX = 2e-3                # 開始前各軸速度的上限 [rad/s]：超過表示手臂在動（有程式在跑）
MOVE_WARN_DEG = 0.5                # 記錄中任一軸角度變化超過此值就警告（敲擊只會造成很小的晃動）
ROBOT_MODE_RUNNING = 7             # RTDE robot_mode：7 = 一般模式（剎車放開、伺服保持）
SAFETY_MODE_NORMAL = 1             # RTDE safety_mode：1 = 正常

# 法蘭負載：填教導器上設定的負載質量 [kg]（只有末端工具 FT300、夾爪、3D 相機時 1.30）。
# 只用於確認與記錄（寫進輸出資料夾名稱）；未設定則程式中止 —— 不提供預設值。
PAYLOAD_KG = None
PAYLOAD_TOL_KG = 0.005             # 執行前輸入的教導器負載質量與 PAYLOAD_KG 的最大允許差 [kg]

# 記錄的 RTDE 輸出：(欄名, ur_rtde 方法)。控制器或 ur_rtde 不支援時該欄留空。
SIGNALS_JOINT = (("actual_q", "getActualQ"), ("actual_qd", "getActualQd"), ("actual_current", "getActualCurrent"),
                 ("target_q", "getTargetQ"), ("joint_temp", "getJointTemperatures"),
                 ("target_moment", "getTargetMoment"), ("target_current", "getTargetCurrent"),
                 ("joint_control_output", "getJointControlOutput"))
SIGNALS_VECTOR = (("tool_acc", "getActualToolAccelerometer", ("x", "y", "z")),
                  ("tcp_force", "getActualTCPForce", ("x", "y", "z", "rx", "ry", "rz")))
SIGNALS_SCALAR = (("robot_mode", "getRobotMode"), ("safety_mode", "getSafetyMode"))

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")   # 存到這支程式所在資料夾裡的 outputs/


# ============================================================
# 第一部分：執行前檢查
# ============================================================

def _finite_float(x):
    """轉成有限浮點數；None、無法轉換、NaN、無窮大時回傳 None。"""
    try:
        v = float(x)
    except (TypeError, ValueError):
        return None
    return v if math.isfinite(v) else None


def check_pose(q_full_rad, pose, poses=POSES, tol_deg=POSE_TOL_DEG):
    """
    以 RTDE 量到的六軸角度 [rad] 檢查 J0~J5 是否在 poses[pose] ± tol_deg 內（差值以 ±360° 等價處理）。
    預設拒絕：pose 不在 poses、讀不到六軸角度、角度不是有限值、任一軸超差，一律回傳 False。回傳 (ok, lines)。
    """
    if pose not in poses:
        return False, [f"未知的姿態 POSE = {pose!r}（可用：{', '.join(poses)}），拒絕執行"]
    if q_full_rad is None or len(q_full_rad) != 6:
        return False, ["讀不到六軸角度，拒絕執行"]
    target_deg = poses[pose]["q_deg"]
    lines, ok = [f"姿態 {pose}：J0~J5 期望 {list(target_deg)}°，容差 ±{tol_deg}°"], True
    for j, tgt in enumerate(target_deg):
        x = _finite_float(q_full_rad[j])
        cur = math.degrees(x) if x is not None else float("nan")
        d = (cur - tgt + 180.0) % 360.0 - 180.0
        passed = abs(d) <= tol_deg                                 # NaN 時為 False
        ok = ok and passed
        lines.append(f"  J{j}：目前 {cur:8.2f}°，期望 {tgt:8.2f}°，差 {d:+7.2f}°  {'[OK]' if passed else '[FAIL]'}")
    return ok, lines


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


def check_record_settings(label, record_s, tap_times, labels=TAP_LABELS, record_s_max=RECORD_S_MAX):
    """
    敲的方向須在 labels 內；記錄時間須在 1～record_s_max 秒；提示時間須遞增、間隔 ≥ 2 s，
    第一下 ≥ 1 s、最後一下後至少留 2 s。回傳 (ok, lines)。
    """
    lines, ok = [], True
    if label not in labels:
        return False, [f"未知的 TAP_LABEL = {label!r}（可用：{', '.join(labels)}），拒絕執行"]
    lines.append(f"敲的方向 {label}：{labels[label]}")
    r = _finite_float(record_s)
    if r is None or not (1.0 <= r <= record_s_max):
        return False, lines + [f"RECORD_S = {record_s!r} 須在 1～{record_s_max:.0f} 秒內，拒絕執行"]
    ts = [_finite_float(x) for x in tap_times]
    if any(x is None for x in ts):
        return False, lines + ["TAP_TIMES 含讀不懂的數值，拒絕執行"]
    if ts:
        gaps_ok = all(b - a >= 2.0 for a, b in zip(ts[:-1], ts[1:]))
        if not (ts[0] >= 1.0 and gaps_ok and r - ts[-1] >= 2.0):
            ok = False
            lines.append(f"TAP_TIMES = {list(ts)}：第一下須 ≥ 1 s、間隔須 ≥ 2 s、最後一下後須留 2 s（RECORD_S = {r:.1f}）  [FAIL]")
    if ok:
        lines.append(f"記錄 {r:.1f} s，提示時間 {list(ts)} s  [OK]")
    return ok, lines


def check_robot_state(robot_mode, safety_mode, mode_ok=ROBOT_MODE_RUNNING, safety_ok=SAFETY_MODE_NORMAL):
    """
    機器人須在一般模式（robot_mode = 7：剎車放開、伺服保持）且安全模式正常（safety_mode = 1）。
    預設拒絕：讀不到、不是整數、不符，一律回傳 False。回傳 (ok, lines)。
    """
    lines, ok = [], True
    for name, val, want in (("robot_mode", robot_mode, mode_ok), ("safety_mode", safety_mode, safety_ok)):
        x = _finite_float(val)
        passed = x is not None and int(x) == want and x == int(x)
        ok = ok and passed
        lines.append(f"{name} = {val!r}（須為 {want}）  {'[OK]' if passed else '[FAIL]'}")
    if not ok:
        lines.append("手臂須在一般模式（剎車放開、沒有保護性停止）才量得到伺服保持時的餘振")
    return ok, lines


def check_still(qd_full, qd_max=STILL_QD_MAX):
    """六軸實際速度都須低於 qd_max（手臂靜止、沒有程式在跑）。讀不到一律拒絕。回傳 (ok, 說明)。"""
    if qd_full is None or len(qd_full) != 6:
        return False, "讀不到六軸速度，拒絕執行"
    v = [_finite_float(x) for x in qd_full]
    if any(x is None for x in v):
        return False, "六軸速度含讀不懂的數值，拒絕執行"
    worst = max(abs(x) for x in v)
    if worst > qd_max:
        return False, f"手臂在動（最大速度 {worst:.4f} rad/s > {qd_max} rad/s）：請先停止教導器上的程式，拒絕執行"
    return True, f"手臂靜止（最大速度 {worst:.5f} rad/s ≤ {qd_max} rad/s） -> OK"


# ============================================================
# 第二部分：RTDE 記錄（只讀取）
# ============================================================

def _call(getter):
    """呼叫 ur_rtde 的讀取方法；沒有這個方法或呼叫失敗時回傳 None。"""
    if getter is None:
        return None
    try:
        return getter()
    except Exception:
        return None


def record(rtde_r, csv_path, duration_s, tap_times=(), period=DT, out=print):
    """
    記錄 duration_s 秒（以控制器時間計）。只呼叫 rtde_r 的讀取方法。
    到達 tap_times 的時間點時以 out 提示「敲」。回傳統計 dict。
    """
    get_ts = getattr(rtde_r, "getTimestamp", None)
    joint = [(name, getattr(rtde_r, meth, None)) for name, meth in SIGNALS_JOINT]
    vector = [(name, getattr(rtde_r, meth, None), axes) for name, meth, axes in SIGNALS_VECTOR]
    scalar = [(name, getattr(rtde_r, meth, None)) for name, meth in SIGNALS_SCALAR]
    reads = {name: 0 for name, _ in joint}
    reads.update({name: 0 for name, _, _ in vector})
    reads.update({name: 0 for name, _ in scalar})
    stats = dict(rows=0, gaps=0, reads=reads, bad_robot_mode=0, bad_safety_mode=0,
                 q_min=[None] * 6, q_max=[None] * 6, duration=0.0)
    cues = list(tap_times)
    last_ts, t0, writer = None, None, None
    wall0 = time.time()
    poll = period / 5.0
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        while True:
            if time.time() - wall0 > duration_s + 5.0:             # 控制器時間不走時的保險
                break
            ts = _call(get_ts)
            if ts is not None:
                if ts == last_ts:
                    time.sleep(poll)
                    continue
                if last_ts is not None and (ts - last_ts) > period * 1.5:
                    stats["gaps"] += 1
                last_ts = ts
                if t0 is None:
                    t0 = ts
                t = ts - t0
            else:
                t = time.time() - wall0
            if t > duration_s:
                break
            row = {"timestamp": round(t, 6)}
            for name, getter in joint:
                vals = _call(getter)
                vals = list(vals) if vals is not None and len(vals) == 6 else None
                if vals is not None:
                    reads[name] += 1
                for j in range(6):
                    row[f"{name}_{j}"] = vals[j] if vals is not None else None
                if name == "actual_q" and vals is not None:
                    for j in range(6):
                        stats["q_min"][j] = vals[j] if stats["q_min"][j] is None else min(stats["q_min"][j], vals[j])
                        stats["q_max"][j] = vals[j] if stats["q_max"][j] is None else max(stats["q_max"][j], vals[j])
            for name, getter, axes in vector:
                vals = _call(getter)
                vals = list(vals) if vals is not None and len(vals) == len(axes) else None
                if vals is not None:
                    reads[name] += 1
                for k, ax in enumerate(axes):
                    row[f"{name}_{ax}"] = vals[k] if vals is not None else None
            for name, getter in scalar:
                val = _finite_float(_call(getter))
                if val is not None:
                    reads[name] += 1
                row[name] = int(val) if val is not None else None
            if row["robot_mode"] != ROBOT_MODE_RUNNING:
                stats["bad_robot_mode"] += 1
            if row["safety_mode"] != SAFETY_MODE_NORMAL:
                stats["bad_safety_mode"] += 1
            if writer is None:
                writer = csv.DictWriter(f, fieldnames=list(row.keys()))
                writer.writeheader()
            writer.writerow(row)
            stats["rows"] += 1
            stats["duration"] = t
            while cues and t >= cues[0]:
                out(f"  {cues.pop(0):5.1f} s  → 敲！")
            time.sleep(poll)
    return stats


def record_summary(stats):
    """記錄統計的文字摘要與警告。回傳 (lines, warnings)。"""
    lines = [f"寫入列數 : {stats['rows']}（{stats['duration']:.2f} s），疑似漏拍 : {stats['gaps']} 次",
             "各訊號成功讀取列數 : " + "，".join(f"{k} {v}" for k, v in stats["reads"].items())]
    warns = []
    if stats["rows"] == 0:
        warns.append("沒有記錄到任何資料")
        return lines, warns
    span = [math.degrees(b - a) if a is not None and b is not None else float("nan")
            for a, b in zip(stats["q_min"], stats["q_max"])]
    lines.append("記錄中各軸角度變化 [deg] : " + "，".join(f"J{j} {s:.4f}" for j, s in enumerate(span)))
    if any(not (s <= MOVE_WARN_DEG) for s in span):
        warns.append(f"記錄中有軸的角度變化超過 {MOVE_WARN_DEG}°（或讀不到）：手臂不是靜止的，這次資料不可用")
    if stats["bad_safety_mode"]:
        warns.append(f"記錄中有 {stats['bad_safety_mode']} 列安全模式不是正常（可能敲太大力觸發保護性停止）：這次資料不可用")
    if stats["bad_robot_mode"]:
        warns.append(f"記錄中有 {stats['bad_robot_mode']} 列機器人模式不是一般模式：這次資料不可用")
    if stats["reads"]["tool_acc"] < stats["rows"]:
        warns.append(f"工具加速度計有 {stats['rows'] - stats['reads']['tool_acc']} 列讀不到：分析改用關節電流找敲擊")
    return lines, warns


# ============================================================
# 第三部分：離線分析
# ============================================================

def _load_csv(csv_path):
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    def col(key):
        if not rows or key not in rows[0]:
            return np.full(len(rows), np.nan)
        return np.array([float(r[key]) if r.get(key) not in (None, "", "None") else np.nan for r in rows])
    return rows, col


def find_taps(t, x, dt=DT, k_sigma=8.0, floor=0.0, min_gap_s=1.0):
    """
    在訊號 x 裡找敲擊：x 減去 0.2 s 移動中位數後取絕對值，超過 max(k_sigma × 穩健標準差, floor) 的第一個樣本為敲擊，
    之後 min_gap_s 秒內不再找。回傳樣本編號 list。
    """
    n = len(x)
    if n < 20 or not np.all(np.isfinite(x)):
        return []
    h = max(2, int(round(0.1 / dt)))
    med = np.array([np.median(x[max(0, i - h):i + h + 1]) for i in range(n)])
    r = np.abs(x - med)
    sigma = 1.4826 * np.median(r)
    thr = max(k_sigma * sigma, floor)
    if not thr > 0:
        return []
    taps, i, gap = [], 0, int(round(min_gap_s / dt))
    while i < n:
        if r[i] > thr:
            taps.append(i)
            i += gap
        else:
            i += 1
    return taps


def _spectrum_peaks(freq, power, f_lo, f_hi, n_peaks=3):
    """f_lo～f_hi 內的局部極大，依大小排序取前 n_peaks 個。回傳 [(頻率, 相對最大峰的比例)]。"""
    m = (freq >= f_lo) & (freq <= f_hi)
    idx = np.flatnonzero(m)
    pk = [k for k in idx[1:-1] if power[k] > power[k - 1] and power[k] >= power[k + 1]]
    if not pk:
        return []
    pk.sort(key=lambda k: -power[k])
    top = power[pk[0]]
    return [(float(freq[k]), float(power[k] / top)) for k in pk[:n_peaks]]


def analyze_tap_test(csv_path, dt=DT, skip_s=0.04, win_s=2.5, min_win_s=0.8, f_lo=3.0, f_hi=45.0, n_fft=4096,
                     quiet=False):
    """
    1. 找敲擊：工具加速度計三軸各自減去中位數後的合成量（讀不到時改用六軸電流變化的合成量）。
    2. 每一下取敲擊後 skip_s 起、最長 win_s 秒（到下一下敲擊或資料結束前為止，不足 min_win_s 秒的不用）。
    3. 各訊號（加速度計三軸、六軸角度、六軸電流）去掉直線趨勢後做頻譜（補零到 n_fft 點），各下的功率譜平均。
    4. 列出 f_lo～f_hi 內的峰值頻率。頻率解析度約 1 ÷ 視窗長度（2.5 s 視窗約 0.4 Hz）。
    沒有敲擊（TAP_LABEL = "none"）時，把整段切成 win_s 秒的區塊做同樣的頻譜，當作背景。
    回傳 dict(taps=敲擊時間, n_win=用到的視窗數, mode="tap"|"background", peaks={訊號: [(頻率, 比例)]}, freq, power)。
    """
    rows, col = _load_csv(csv_path)
    if len(rows) < 50:
        if not quiet:
            print("[警告] 資料太少，無法分析。")
        return dict(taps=[], n_win=0, mode="none", peaks={}, freq=None, power={})
    t = col("timestamp")
    acc = np.column_stack([col(f"tool_acc_{a}") for a in "xyz"])
    cur = np.column_stack([col(f"actual_current_{j}") for j in range(6)])
    if np.all(np.isfinite(acc)):
        trig, src, floor = np.linalg.norm(acc - np.median(acc, axis=0), axis=1), "工具加速度計", 0.2
    else:
        trig, src, floor = np.linalg.norm(cur - np.median(cur, axis=0), axis=1), "關節電流", 0.02
    taps = find_taps(t, trig, dt, floor=floor)
    n, n_skip, n_win, n_min = len(t), int(round(skip_s / dt)), int(round(win_s / dt)), int(round(min_win_s / dt))
    wins = []
    for a, k in enumerate(taps):
        end = min(n, k + n_skip + n_win, taps[a + 1] - n_skip if a + 1 < len(taps) else n)
        if end - (k + n_skip) >= n_min:
            wins.append((k + n_skip, end))
    mode = "tap"
    if not wins:
        mode = "background"
        wins = [(s, s + n_win) for s in range(0, n - n_win + 1, n_win)]
    signals = {f"tool_acc_{a}": acc[:, k] for k, a in enumerate("xyz")}
    signals.update({f"actual_q_{j}": col(f"actual_q_{j}") for j in range(6)})
    signals.update({f"actual_current_{j}": cur[:, j] for j in range(6)})
    freq = np.fft.rfftfreq(n_fft, dt)
    power, peaks = {}, {}
    for name, y in signals.items():
        acc_p, used = np.zeros(len(freq)), 0
        for s, e in wins:
            seg = y[s:e]
            if len(seg) < n_min or not np.all(np.isfinite(seg)) or np.ptp(seg) < 1e-12:   # 讀不到或完全沒有變化的不算
                continue
            x = np.arange(len(seg))
            seg = seg - np.polyval(np.polyfit(x, seg, 1), x)
            acc_p += np.abs(np.fft.rfft(seg, n_fft)) ** 2
            used += 1
        if used:
            power[name] = acc_p / used
            peaks[name] = _spectrum_peaks(freq, power[name], f_lo, f_hi)
    res = dict(taps=[float(t[k]) for k in taps], n_win=len(wins), mode=mode, peaks=peaks, freq=freq, power=power,
               trigger=src)
    if not quiet:
        print_tap_analysis(res)
    return res


def tap_analysis_lines(res):
    lines = []
    if res["mode"] == "none":
        return ["資料太少，沒有分析。"]
    if res["mode"] == "tap":
        lines.append(f"敲擊：以{res['trigger']}找到 {len(res['taps'])} 下，時間 "
                     + "、".join(f"{x:.2f}" for x in res["taps"]) + f" s；用到 {res['n_win']} 個餘振視窗")
    else:
        lines.append(f"沒有找到敲擊：把整段切成 {res['n_win']} 個區塊，當作背景頻譜")
    lines.append("各訊號的峰值頻率 [Hz]（括號內為相對該訊號最大峰的功率比例）：")
    for name, pk in res["peaks"].items():
        lines.append(f"  {name:18s} " + ("，".join(f"{f:5.2f} ({r:.2f})" for f, r in pk) if pk else "（沒有峰）"))
    return lines


def print_tap_analysis(res):
    for line in tap_analysis_lines(res):
        print(line)


# ============================================================
# 第四部分：主程式
# ============================================================

def main():
    print("=" * 64)
    print("UR5 敲擊測試（只記錄，不送任何指令給手臂）")
    print("=" * 64)

    ok, msg = check_payload_setting(PAYLOAD_KG)
    print(msg)
    if not ok:
        print("[中止] 負載設定未通過。")
        return
    ok, lines = check_record_settings(TAP_LABEL, RECORD_S, TAP_TIMES)
    for line in lines:
        print(line)
    if not ok:
        print("[中止] 記錄設定未通過。")
        return
    if POSE not in POSES:
        print(f"未知的姿態 POSE = {POSE!r}（可用：{', '.join(POSES)}）")
        print("[中止] 姿態檢查未通過。")
        return
    if rtde_receive is None:
        print("[中止] 找不到 rtde_receive（ur_rtde），無法記錄。")
        return

    try:
        rtde_r = rtde_receive.RTDEReceiveInterface(ROBOT_IP)
    except Exception as e:
        print(f"[中止] RTDE 連線失敗：{e}")
        return
    try:
        ok, lines = check_pose(_call(getattr(rtde_r, "getActualQ", None)), POSE)
        for line in lines:
            print(line)
        if not ok:
            print("[中止] 姿態檢查未通過。")
            return
        ok, msg = check_still(_call(getattr(rtde_r, "getActualQd", None)))
        print(msg)
        if not ok:
            print("[中止] 靜止檢查未通過。")
            return
        ok, lines = check_robot_state(_call(getattr(rtde_r, "getRobotMode", None)),
                                      _call(getattr(rtde_r, "getSafetyMode", None)))
        for line in lines:
            print(line)
        if not ok:
            print("[中止] 機器人狀態檢查未通過。")
            return

        print()
        print(f"請確認：教導器上沒有程式在跑、手臂停在 {POSE}、另一人在緊急停止按鈕旁。")
        print(f"這次是「{TAP_LABELS[TAP_LABEL]}」。敲手腕 3 的金屬外殼或工具法蘭側面，不要敲 FT300、夾爪手指、相機。")
        if input(f"已確認，輸入 {POSE} 繼續：").strip() != POSE:
            print("[中止] 已取消。")
            return
        ok, msg = check_payload_entry(PAYLOAD_KG, input("輸入教導器上的負載質量 [kg]："))
        print(msg)
        if not ok:
            print("[中止] 負載確認未通過。")
            return
        input("按 Enter 開始倒數，或 Ctrl+C 取消：")

        tag = f"{POSE.lower()}_{payload_tag(PAYLOAD_KG)}_{TAP_LABEL}_{datetime.now().strftime('%Y%m%d_%H%M%S')}"
        session_dir = os.path.join(OUTPUT_DIR, f"taptest_{tag}")
        os.makedirs(session_dir, exist_ok=True)
        csv_path = os.path.join(session_dir, f"taptest_data_{tag}.csv")
        info_path = os.path.join(session_dir, f"taptest_info_{tag}.txt")

        for k in range(int(COUNTDOWN_S), 0, -1):
            print(f"  {k} ...")
            time.sleep(1.0)
        cues = () if TAP_LABEL == "none" else TAP_TIMES
        print(f"開始記錄 {RECORD_S:.0f} s" + ("（這次不敲）" if not cues else "：看到「敲」就敲一下，敲完手離開"))
        stats = record(rtde_r, csv_path, RECORD_S, cues)
        print("記錄結束。")
    finally:
        try:
            rtde_r.disconnect()
        except Exception:
            pass

    lines, warns = record_summary(stats)
    report = [f"姿態: {POSE}", f"負載: {float(PAYLOAD_KG):.3f} kg", f"敲的方向: {TAP_LABEL}（{TAP_LABELS[TAP_LABEL]}）",
              f"資料: {csv_path}"] + lines
    for w in warns:
        report.append(f"[警告] {w}")
    try:
        res = analyze_tap_test(csv_path, quiet=True)
        report += tap_analysis_lines(res)
    except Exception as e:                                         # 分析失敗不影響已存的資料
        report.append(f"[警告] 分析失敗：{e}")
    print()
    for line in report:
        print(line)
    with open(info_path, "w", encoding="utf-8") as f:
        f.write("\n".join(report) + "\n")
    print(f"\n摘要已存：{info_path}")


if __name__ == "__main__":
    main()
