"""
UR5 J0 定加速度實驗 v3 —— 鑑別 J，並由同一份資料交叉驗證 B、Tc
================================================================

★★★ v3 跟 v2（ur5_const_accel_ident_v2_movej.py）的差異 ★★★

v2 的問題（2026/09/30 真機資料）：
  (1) 每次爬升都從同一個絕對角度出發、爬升段不到 1 個馬達圈，跟馬達角度鎖定的電流漣波
      在 40 次重複中每次都一樣，平均不掉（40 次電流波形相關係數 0.95～0.97）。
  (2) 只有加速段，摩擦必須借用另一次實驗的 B、Tc 扣除；日間 Tc 差 0.3～0.4 N·m，
      除以小加速度後 J 的誤差可達 1 kg·m² 以上。
  (3) a = 2.0 時爬升只有約 5 個樣本，實際加速度只達命令的 91～95%（伺服跟不上）。
  (4) 爬升峰值 0.15 rad/s，速度範圍太窄，三參數同時回歸時 B 不可靠。

v3 的做法：
  A. 每次循環走「梯形速度」：以加速度 a 爬升到 V_PEAK → 定速 HOLD（整數個馬達圈）→
     以**同樣大小的 a** 線性減速到 0（speedj 至零速，不再用 stopj 的黑盒煞車曲線）。
     同一方向、同一速度 v 下：
         加速段  tau_acc(v) =  J·a + f(v)
         減速段  tau_dec(v) = −J·a + f(v)      f(v) = B·v + Tc（摩擦）
     → J    = [tau_acc(v) − tau_dec(v)] / (a_acc + a_dec)   摩擦完全相消，不需要借用 B、Tc
     → f(v) = [tau_acc(v) + tau_dec(v)] / 2                 慣性項相消，得到摩擦曲線 → B、Tc
     定速段的力矩 tau_hold 另外與 f(V_PEAK) 比對，檢查一致性。
  B. 第 k 次循環的起始角度往運動方向錯開 k/N_REPEAT 個馬達圈（N_REPEAT 次剛好涵蓋一整圈），
     漣波在重複之間平均掉；回程與段間移動仍用 movej 絕對定位（v2 已驗證不累積）。
  C. 加速度檔改為 0.3、0.6、0.9、1.2 rad/s²（移除伺服跟不上的 2.0），V_PEAK 提高到
     0.25 rad/s（仍低於 0.30 上限，也低於每圈 12 次漣波在 125 Hz 下會混疊的約 0.32 rad/s），
     分析用速度範圍 0.05～0.2375 rad/s。
  D. 暖機：請先執行 ur5_j0_warmup.py。

★★★ 使用前必讀 ★★★
1. 【新腳本，尚未經學長審查，也尚未上機】依 TEST_STAGE = QUICK → MID → FULL 階梯執行。
2. 【方向分開兩次執行】DIRECTION=+1 先做，人工確認後才改 −1。
3. 【URScript 為逐次展開】每次循環的 movej 目標角度不同（起始角錯開），程式把每次循環都
   明確寫出來，不在 URScript 內做清單運算；FULL 約 80 次循環、約 500 行。
4. 【J0 視窗預設拒絕】J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 未設定一律中止。
5. 【IP 需要你自己再次確認】
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

DIRECTION = +1        # +1 = 正轉，-1 = 反轉。兩個方向分開執行
# QUICK: a = 1.2 × 2 次，確認梯形速度、起始角錯開、movej 歸位
# MID  : a = 0.3、0.6 × 各 4 次，確認較長腳本與行程
# FULL : 四檔 × 各 N_REPEAT 次
TEST_STAGE = "QUICK"

ACCEL_LEVELS = [0.3, 0.6, 0.9, 1.2]   # rad/s^2
V_PEAK = 0.25                          # rad/s
N_REPEAT = 20                          # 每檔次數；起始角錯開 1/N_REPEAT 個馬達圈
MID_N_REPEAT = 4
QUICK_N_REPEAT = 2
GEAR_RATIO = 101
HOLD_MOTOR_REVS = 4                    # 定速段長度：4 個馬達圈（0.25 rad/s 時約 1.0 s）
PAUSE_TIME = 0.3                       # 每次循環前後靜止時間 [s]（須小於結束判定的 1 s 靜止窗）

QD_MAX = 0.30
QDD_CHECK_LIMIT = 3.0
MAX_EXCURSION_DEG = 90.0               # 單次循環離起點最遠距離的警戒線
MOVEJ_RETURN_ACCEL = 1.0               # 同 v2
MOVEJ_RETURN_VEL = 0.15                # 同 v2
EXCURSION_MARGIN_FACTOR = 1.5          # 控制器軌跡產生器的不確定係數；QUICK 實測後回填

SAMPLE_HZ = 125.0
DT = 1.0 / SAMPLE_HZ
KT_OUT = 101 * 0.1350

J0_SAFE_MIN_DEG = None
J0_SAFE_MAX_DEG = None

OUTPUT_DIR = "."


# ============================================================
# 第一部分：規劃與安全檢查
# ============================================================

def motor_rev_rad(gear_ratio=GEAR_RATIO):
    """馬達一圈對應的關節角 [rad]。"""
    return 2.0 * math.pi / gear_ratio


def hold_time(v_peak, hold_revs, gear_ratio=GEAR_RATIO):
    """定速段時間：剛好 hold_revs 個馬達圈。"""
    return hold_revs * motor_rev_rad(gear_ratio) / v_peak


def cycle_excursion_rad(a, v_peak, t_hold):
    """單次梯形（加速 + 定速 + 同 a 減速）的理論行程 [rad]。"""
    return v_peak ** 2 / a + v_peak * t_hold


def start_offsets_rad(n_repeat, gear_ratio=GEAR_RATIO):
    """第 k 次循環的起始角錯開量：k/n_repeat 個馬達圈（k = 0..n_repeat−1）。"""
    return [k * motor_rev_rad(gear_ratio) / n_repeat for k in range(n_repeat)]


def stage_config(stage):
    if stage == "QUICK":
        return [ACCEL_LEVELS[-1]], QUICK_N_REPEAT
    if stage == "MID":
        return ACCEL_LEVELS[:2], MID_N_REPEAT
    if stage == "FULL":
        return list(ACCEL_LEVELS), N_REPEAT
    raise ValueError(f"未知的 TEST_STAGE: {stage}")


def safety_check(accels, v_peak, n_repeat, t_hold, qd_max, qdd_max, max_excursion_deg,
                 margin, movej_accel, movej_vel):
    """回傳 (ok, lines, 規劃最遠行程[deg])。最遠行程 = 最大起始錯開 + 單次梯形行程。"""
    lines, ok = [], True

    def chk(name, val, lim, unit):
        nonlocal ok
        lines.append(f"{name} : {val:.4f} {unit}  (限 {lim})")
        if val > lim:
            lines.append("  [FAIL]"); ok = False
        else:
            lines.append("  [OK]")

    chk("峰值速度", v_peak, qd_max, "rad/s")
    chk("最大加速度（加速與減速相同）", max(accels), qdd_max, "rad/s^2")
    chk("movej 回程速度", movej_vel, qd_max, "rad/s")
    chk("movej 回程加速度", movej_accel, qdd_max, "rad/s^2")

    max_off = max(start_offsets_rad(n_repeat))
    excursion = max_off + max(cycle_excursion_rad(a, v_peak, t_hold) for a in accels)
    exc_deg = math.degrees(excursion)
    lines.append(f"單次循環最遠行程（含起始錯開 {math.degrees(max_off):.2f}°） : {exc_deg:.2f}°")
    lines.append(f"x 不確定係數 {margin:.2f} : {exc_deg * margin:.2f}°  (警戒線 {max_excursion_deg}°)")
    if exc_deg * margin > max_excursion_deg:
        lines.append("  [FAIL]"); ok = False
    else:
        lines.append("  [OK]")
    return ok, lines, exc_deg


def check_j0_window(q0_deg, planned_excursion_deg, direction, margin, safe_min_deg, safe_max_deg):
    """同 v2：目前角度 + 規劃行程 × 係數須在 J0 視窗內；未設定一律拒絕。"""
    if safe_min_deg is None or safe_max_deg is None:
        return False, "J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 尚未設定，拒絕執行"
    reach = q0_deg + direction * planned_excursion_deg * margin
    lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
    msg = (f"J0 目前 {q0_deg:+.2f}°，往 {direction:+d} 方向最遠到 {reach:+.2f}°"
           f"（含係數 {margin:.2f}），允許視窗 [{lo:+.1f}°, {hi:+.1f}°]")
    if not (lo <= q0_deg <= hi):
        return False, msg + " -> 起始角度已在視窗外"
    if not (lo <= reach <= hi):
        return False, msg + " -> 行程終點超出視窗"
    return True, msg + " -> OK"


def movej_time(dist_rad, v, a):
    if dist_rad <= 0:
        return 0.0
    t_up = v / a
    d_up = 0.5 * a * t_up ** 2
    if dist_rad <= 2 * d_up:
        return 2.0 * math.sqrt(dist_rad / a)
    return 2.0 * t_up + (dist_rad - 2 * d_up) / v


def estimate_total_time(accels, v_peak, n_repeat, t_hold):
    total = 0.0
    for a in accels:
        d = cycle_excursion_rad(a, v_peak, t_hold)
        t_cycle = 2 * v_peak / a + t_hold + movej_time(d, MOVEJ_RETURN_VEL, MOVEJ_RETURN_ACCEL) + 2 * PAUSE_TIME
        total += n_repeat * t_cycle
    return total


# ============================================================
# 第二部分：URScript
# ============================================================

def build_full_urscript(accels, v_peak, n_repeat, t_hold, direction, joint_index,
                        q_start, movej_accel, movej_vel, pause_time):
    """
    每次循環（逐次展開）：
      movej(起始角 + 錯開) → sleep → speedj(V_PEAK, a, t = 爬升 + 定速) →
      speedj(0, a, t = 減速時間) → stopj(a)（此時速度已近 0，只是確保停止）→ sleep
    最後 movej 回 q_start。
    """
    qd_peak = ["0.0"] * 6
    qd_peak[joint_index] = f"{direction * v_peak:.6f}"
    qd_peak_s = "[" + ", ".join(qd_peak) + "]"
    qd_zero_s = "[" + ", ".join(["0.0"] * 6) + "]"
    offs = start_offsets_rad(n_repeat)

    def q_str(offset):
        q = list(q_start)
        q[joint_index] = q_start[joint_index] + direction * offset
        return "[" + ", ".join(f"{v:.6f}" for v in q) + "]"

    lines = ["def const_accel_v3_program():"]
    for a in accels:
        t_ramp = v_peak / a
        lines.append(f"  # level a = {a}")
        for k, off in enumerate(offs):
            lines.append(f"  movej({q_str(off)}, a={movej_accel}, v={movej_vel})")
            lines.append(f"  sleep({pause_time})")
            lines.append(f"  speedj({qd_peak_s}, a={a}, t={t_ramp + t_hold:.6f})")
            lines.append(f"  speedj({qd_zero_s}, a={a}, t={t_ramp:.6f})")
            lines.append(f"  stopj({a})")
            lines.append(f"  sleep({pause_time})")
    lines.append(f"  movej({q_str(0.0)}, a={movej_accel}, v={movej_vel})")
    lines.append("end")
    return "\n".join(lines) + "\n"


def send_urscript(ip, script_text, port=30002):
    s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
    s.settimeout(5.0)
    s.connect((ip, port))
    s.sendall(script_text.encode("utf-8"))
    s.close()


def send_abort(ip):
    try:
        send_urscript(ip, "def abort_prog():\n  stopj(2.0)\nend\n")
        print("[中止] 已送出 stopj。")
    except Exception as e:
        print(f"[中止] stopj 送出失敗：{e} —— 請立即按下緊急停止。")


def wait_for_motion_complete(logger, planned_s, qd_eps=2e-3,
                             quiet_s=1.0, no_motion_s=10.0, timeout_factor=3.0):
    """同 v2：J0 實際速度連續 quiet_s 秒低於門檻才視為結束。回傳 (status, elapsed, motion_span)。"""
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
# 第三部分：RTDE 記錄（同 v2）
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
        return "\n".join(lines)


# ============================================================
# 第四部分：離線分析
# ============================================================

def _runs(mask):
    idx = np.where(mask)[0]
    return np.split(idx, np.where(np.diff(idx) != 1)[0] + 1) if len(idx) else []


def analyze_const_accel_v3(csv_path, joint_index, kt_out, direction, accel_levels, v_peak,
                           v_lo=0.05, v_hi_frac=0.95, n_bins=6, trim=2, match_tol=0.20,
                           n_boot=500, seed=0):
    """
    以 target_qd 切出每次循環的加速段（斜坡上升）、定速段、減速段（斜坡下降），依加速度歸檔。
    速度分箱（v_lo ～ v_hi_frac·V_PEAK，n_bins 箱）後，對每一檔、每一箱：
        J_bin = [tau_acc − tau_dec] / (|a_acc| + |a_dec|)      （a 用實測速度直線擬合的斜率）
        f_bin = [tau_acc + tau_dec] / 2                         （摩擦曲線）
    結果：
        J（各檔）、J（全部檔位對 a 的過原點回歸與含截距回歸，截距應接近 0）
        B、Tc：由 f(v) 對 v 線性回歸（所有檔位合併，可與定速法交叉驗證）
        tau_hold 與 f(V_PEAK) 的差（一致性檢查）
    所有量以「|tau|、|v|」表示，正反轉相同處理。誤差以「循環」為單位 bootstrap。
    """
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    t = np.array([float(r["timestamp"]) for r in rows])
    d = float(direction)
    qd = d * np.array([float(r[f"actual_qd_{joint_index}"]) for r in rows])
    tq = d * np.array([float(r[f"target_qd_{joint_index}"]) for r in rows])
    tau = d * np.array([float(r[f"actual_current_{joint_index}"]) for r in rows]) * kt_out
    dtq = np.diff(tq, prepend=tq[0])
    v_hi = v_hi_frac * v_peak
    edges = np.linspace(v_lo, v_hi, n_bins + 1)

    acc_runs = [r[trim:] for r in _runs((dtq > 1e-6) & (tq > v_lo) & (tq < v_hi)) if len(r) > trim + 2]
    dec_runs = [r[trim:] for r in _runs((dtq < -1e-6) & (tq > v_lo) & (tq < v_hi)) if len(r) > trim + 2]
    hold_runs = [r for r in _runs(np.abs(tq - v_peak) < 1e-4) if len(r) > 10]

    def classify(run, sign):
        a_cmd = np.polyfit(t[run], tq[run], 1)[0] * sign
        li = int(np.argmin([abs(a_cmd - a) for a in accel_levels]))
        if abs(a_cmd - accel_levels[li]) / accel_levels[li] > match_tol:
            return None, None
        return li, abs(np.polyfit(t[run], qd[run], 1)[0])

    # 依時間把 加速 → 定速 → 減速 配成一次循環
    cycles = []
    for ra in acc_runs:
        li, aa = classify(ra, +1)
        if li is None:
            continue
        t_end = t[ra[-1]]
        rd = next((r for r in dec_runs if t[r[0]] > t_end), None)
        if rd is None:
            continue
        li_d, ad = classify(rd, -1)
        if li_d != li:
            continue
        rh = next((r for r in hold_runs if t_end < t[r[0]] < t[rd[0]]), None)
        tau_hold = float(np.mean(tau[rh[len(rh) // 10:]])) if rh is not None else float("nan")
        b_acc = np.digitize(tq[ra], edges) - 1
        b_dec = np.digitize(tq[rd], edges) - 1

        def bin_mean(x, run, b_idx):
            return np.array([np.mean(x[run][b_idx == b]) if np.any(b_idx == b) else np.nan
                             for b in range(n_bins)])
        cycles.append(dict(level=li, a_acc=aa, a_dec=ad,
                           tau_acc=bin_mean(tau, ra, b_acc), tau_dec=bin_mean(tau, rd, b_dec),
                           v_acc=bin_mean(qd, ra, b_acc), v_dec=bin_mean(qd, rd, b_dec),
                           tau_hold=tau_hold))

    if not cycles:
        print("[警告] 沒有配對成功的循環。")
        return dict(cycles=[])

    def estimate(cs):
        res = {}
        per = []
        for li, a in enumerate(accel_levels):
            c = [x for x in cs if x["level"] == li]
            if not c:
                continue
            per.append(dict(
                TA=np.nanmean([x["tau_acc"] for x in c], axis=0),
                TD=np.nanmean([x["tau_dec"] for x in c], axis=0),
                VA=np.nanmean([x["v_acc"] for x in c], axis=0),
                VD=np.nanmean([x["v_dec"] for x in c], axis=0),
                asum=np.mean([x["a_acc"] + x["a_dec"] for x in c])))
        # 摩擦曲線：f = (tau_acc + tau_dec)/2，速度取兩者實際平均速度的平均
        F = np.concatenate([(p["TA"] + p["TD"]) / 2 for p in per])
        V = np.concatenate([(p["VA"] + p["VD"]) / 2 for p in per])
        ok = ~np.isnan(F) & ~np.isnan(V)
        BT = np.polyfit(V[ok], F[ok], 1)
        # J：同一速度箱內加速、減速的平均速度若略有差異，以 B 修正兩者的摩擦差
        Jl, al = [], []
        for p in per:
            Jb = (p["TA"] - p["TD"] - BT[0] * (p["VA"] - p["VD"])) / p["asum"]
            Jl.append(np.nanmean(Jb)); al.append(p["asum"] / 2)
        res["J_levels"], res["a_levels"] = np.array(Jl), np.array(al)
        h = np.array(Jl) * np.array(al)                       # (tau_acc − tau_dec)/2 各檔平均
        res["J_origin"] = float(np.sum(h * al) / np.sum(np.array(al) ** 2)) if len(al) else float("nan")
        if len(al) >= 2:
            p = np.polyfit(al, h, 1)
            res["J_slope"], res["J_intercept"] = float(p[0]), float(p[1])
        res["B"], res["Tc"] = float(BT[0]), float(BT[1])
        th = [x["tau_hold"] for x in cs if not math.isnan(x["tau_hold"])]
        res["hold_minus_f"] = float(np.mean(th) - (BT[0] * v_peak + BT[1])) if th else float("nan")
        return res

    est = estimate(cycles)
    rng = np.random.default_rng(seed)
    boots = []
    for _ in range(n_boot):
        pick = rng.integers(0, len(cycles), len(cycles))
        try:
            boots.append(estimate([cycles[i] for i in pick]))
        except Exception:
            continue
    def se(key):
        ref = np.shape(est[key])
        vals = [b[key] for b in boots
                if key in b and np.shape(b[key]) == ref and np.all(np.isfinite(b[key]))]
        return np.std(vals, axis=0) if vals else float("nan")
    est["se_J_origin"], est["se_B"], est["se_Tc"] = se("J_origin"), se("B"), se("Tc")
    est["se_J_levels"] = se("J_levels")
    est["cycles"] = cycles

    print(f"配對成功循環數：{len(cycles)}")
    used = [li for li in range(len(accel_levels)) if any(c["level"] == li for c in cycles)]
    se_lv = np.atleast_1d(est["se_J_levels"])
    print(f"{'檔位a':>6}{'循環':>6}{'實測a':>9}{'J':>10}{'±':>8}")
    for k, li in enumerate(used):
        n = sum(1 for x in cycles if x["level"] == li)
        s = se_lv[k] if len(se_lv) == len(used) else float("nan")
        print(f"{accel_levels[li]:6.2f}{n:6d}{est['a_levels'][k]:9.3f}{est['J_levels'][k]:10.4f}{s:8.4f}")
    print(f"J（過原點，各檔合併） = {est['J_origin']:.4f} ± {est['se_J_origin']:.4f}")
    if "J_slope" in est:
        print(f"J（含截距）= {est['J_slope']:.4f}，截距 = {est['J_intercept']:+.4f} N·m（應接近 0）")
    print(f"摩擦曲線：B = {est['B']:.3f} ± {est['se_B']:.3f}，Tc = {est['Tc']:.3f} ± {est['se_Tc']:.3f}")
    print(f"定速段力矩 − f(V_PEAK) = {est['hold_minus_f']:+.4f} N·m（應接近 0）")
    return est


# ============================================================
# 主流程
# ============================================================

def main():
    print("=" * 70)
    print(" 定加速度實驗 v3（梯形速度、起始角錯開）：安全檢查")
    print("=" * 70)
    print(f"方向: {'正轉 (+1)' if DIRECTION > 0 else '反轉 (-1)'}")
    accels, n_repeat = stage_config(TEST_STAGE)
    t_hold = hold_time(V_PEAK, HOLD_MOTOR_REVS)
    print(f"模式: {TEST_STAGE}；加速度檔 {accels} × 各 {n_repeat} 次；峰值 {V_PEAK} rad/s，"
          f"定速 {t_hold:.3f} s（{HOLD_MOTOR_REVS} 個馬達圈）")

    ok, report, planned_deg = safety_check(accels, V_PEAK, n_repeat, t_hold, QD_MAX, QDD_CHECK_LIMIT,
                                           MAX_EXCURSION_DEG, EXCURSION_MARGIN_FACTOR,
                                           MOVEJ_RETURN_ACCEL, MOVEJ_RETURN_VEL)
    print("\n".join(report))
    if not ok:
        print("\n[中止] 安全檢查未通過。"); return
    print("\n[通過安全檢查]")
    if ROBOT_IP is None:
        print("\n[中止] ROBOT_IP 未設定。"); return
    if rtde_receive is None:
        print("\n[中止] 找不到 rtde_receive，僅完成離線設計檢查。"); return

    print("\n執行前姿態確認（J1~J5 診斷原點）")
    if home_pose is None:
        print("[警告] 找不到 ur5_home_pose 模組，請手動確認 J1~J5 為 -90°, 90°, -90°, -90°, 0°")
        if input("已確認，輸入 yes 繼續：").strip().lower() != "yes":
            print("已取消。"); return
    elif not home_pose.ensure_home_pose(ROBOT_IP, mode="CHECK"):
        print("\n[中止] 姿態不符合診斷原點。"); return

    print(f"\n即將對 IP={ROBOT_IP} 送出軌跡，關節 J{JOINT_INDEX}，方向 {DIRECTION:+d}。")
    print("*** 送出後 Ctrl+C 會送出 stopj，但最可靠的仍是緊急停止按鈕 ***")
    input("按 Enter 繼續，或 Ctrl+C 取消...")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_label = "pos" if DIRECTION > 0 else "neg"
    tag = f"{dir_label}_{TEST_STAGE.lower()}_{ts}"
    session_dir = os.path.join(OUTPUT_DIR, f"constaccel_v3_{tag}")
    os.makedirs(session_dir, exist_ok=True)
    csv_path = os.path.join(session_dir, f"constaccel_v3_data_{tag}.csv")
    info_path = os.path.join(session_dir, f"constaccel_v3_info_{tag}.txt")

    logger = RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, JOINT_INDEX)
    logger.start()
    if not logger.started_ok.wait(5.0):
        logger.stop()
        print(f"\n[中止] RTDE 記錄未能啟動：{logger.error}。未送出任何軌跡，手臂未動作。")
        return

    q_start = logger.last_q0
    q_start_full = logger.last_q_full
    w_ok, w_msg = check_j0_window(math.degrees(q_start), planned_deg, DIRECTION,
                                  EXCURSION_MARGIN_FACTOR, J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG)
    print(f"\n[J0 視窗檢查] {w_msg}")
    if not w_ok:
        logger.stop()
        print("\n[中止] J0 視窗檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    script = build_full_urscript(accels, V_PEAK, n_repeat, t_hold, DIRECTION, JOINT_INDEX,
                                 q_start_full, MOVEJ_RETURN_ACCEL, MOVEJ_RETURN_VEL, PAUSE_TIME)
    est_total = estimate_total_time(accels, V_PEAK, n_repeat, t_hold)
    status, elapsed, span = "unknown", 0.0, 0.0
    try:
        print(f"[資訊] 送出 URScript（{len(script.splitlines())} 行，預估 {est_total:.1f} 秒）...")
        send_urscript(ROBOT_IP, script)
        status, elapsed, span = wait_for_motion_complete(logger, est_total)
        if status == "no_motion":
            print(f"\n[警告] 送出後 {elapsed:.1f}s 內未偵測到 J0 運動（腳本可能未被接受）。")
        elif status == "timeout":
            print(f"\n[警告] 超過硬性逾時（{elapsed:.1f}s），主動中止。")
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
        fwd = math.degrees(logger.q0_max - q_start) if DIRECTION > 0 else math.degrees(q_start - logger.q0_min)
        ratio = fwd / planned_deg if planned_deg > 0 else float("nan")
        measured = [
            f"實測行程 / 規劃行程 : {fwd:.3f}° / {planned_deg:.3f}° = {ratio:.3f}",
            f">>> 請將 EXCURSION_MARGIN_FACTOR 回填為 {max(ratio, 1.0) * 1.1:.2f}（實測值 x 1.1）後再跑下一階段",
            f"movej 歸位殘差 : {math.degrees(logger.last_q0 - q_start):+.4f}°（應趨近 0）",
            f"實際運動時長 / 預估 : {span:.2f}s / {est_total:.2f}s",
        ]
        print("\n" + "\n".join(measured))

    with open(info_path, "w", encoding="utf-8") as f:
        f.write(f"方向: {DIRECTION:+d}\n模式: {TEST_STAGE}\n加速度檔位: {accels}\n每檔次數: {n_repeat}\n")
        f.write(f"峰值速度: {V_PEAK}\n定速時間: {t_hold:.4f}s（{HOLD_MOTOR_REVS} 個馬達圈）\n")
        f.write(f"預估總時長: {est_total:.1f}s\n執行狀態: {status}\n")
        f.write(logger.summary() + "\n")
        if measured:
            f.write("\n" + "\n".join(measured) + "\n")
    print(f"\n[完成] 資料存至 {csv_path}")
    print(f"[完成] 摘要存至 {info_path}")
    print("\n下一步：analyze_const_accel_v3(csv_path, JOINT_INDEX, KT_OUT, DIRECTION, ACCEL_LEVELS, V_PEAK)")


if __name__ == "__main__":
    main()
