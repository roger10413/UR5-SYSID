"""
UR5 J0 固定角度窗口定速實驗 —— 鑑別 B, Tc，並觀察摩擦與漣波隨角度的變化
================================================================

★★★ 這支跟原本定速腳本（ur5_const_velocity_ident.py）的差異 ★★★

原本：每個速度檔「定速 6 秒」（以時間為基準）。結果是：
  (1) 每檔掃過的角度長短不同，低速檔只有 1.77 個馬達圈；
  (2) 正轉掃 180°→254°、反轉掃 180°→105°，兩個方向從來沒有量過同一段角度，
      正反轉的差異分不出是方向造成還是位置造成。

這支：以「角度」為基準。所有速度檔、正轉與反轉，都以定速掃過**同一段絕對
角度窗口** [THETA_A, THETA_B]，窗口長度剛好 N_MOTOR_REVS_WINDOW 個馬達圈。

每一段（一個速度檔）的動作：
  1. movej 移到該段起點（窗口起點往回退「加速距離 + 穩定距離」），其他五軸維持不動
  2. speedj 以 SPEEDJ_ACCEL 加速到目標速度，在進入窗口前至少穩定 SETTLE_TIME
  3. 定速掃過整個窗口；**控制器每個週期讀 J0 實際角度，超過終點才停**（角度為基準，
     不受迴圈時序拉伸影響），另設迴圈次數上限作為保險
  4. stopj 減速，停 PAUSE_TIME 後進入下一段
  回程與段間移動都用 movej 絕對定位（同 ur5_const_accel_ident_v2_movej.py，已實測
  200 次往復殘差 −0.0014° / −0.0041°），不會累積角度。

速度順序 SPEED_ORDER：
  "UP"      ：0.02 → 0.20
  "UP_DOWN" ：0.02 → 0.20 → 0.02（共 12 段）。時間效應在往返兩趟互相抵銷，
              兩趟差異也可直接看出有沒有時間效應（暖機、漂移）。

★★★ 使用前必讀 ★★★
1. 【新腳本，尚未經學長審查，也尚未上機】依 TEST_STAGE = QUICK → MID → FULL 階梯執行。
2. 【兩個方向分開執行】DIRECTION=+1 跑完、人工確認後再改 −1。
3. 【URScript 結構與原本不同】迴圈改為「讀 J0 角度判斷是否到終點」，這是新結構，
   需學長確認。迴圈仍遵守「第一行 def 頂格、其餘行至少縮排一格」。
4. 【J0 視窗檢查改為絕對角度】本實驗的軌跡是絕對角度，安全檢查直接檢查
   「所有段的起點、終點、煞車後最遠點 ± 餘裕」是否在 J0 視窗內；未設定視窗一律中止。
5. 【第一次 movej 的距離有上限】目前 J0 與第一段起點相差超過 MAX_INITIAL_MOVE_DEG 就中止，
   避免手臂從意料外的位置大幅轉動。
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
# QUICK: 單一檔 0.05 rad/s、窗口 1 個馬達圈（3.56°），驗證方向、角度停止與 movej 定位
# MID  : 三檔 0.05 / 0.12 / 0.20 rad/s、完整窗口、只升速一趟，驗證多段與較長腳本
# FULL : 六檔、完整窗口、先升後降（12 段）
TEST_STAGE = "QUICK"

SPEED_LEVELS = [0.02, 0.05, 0.08, 0.12, 0.16, 0.20]   # rad/s
SPEED_ORDER = "UP_DOWN"                              # FULL 用；"UP" 或 "UP_DOWN"

GEAR_RATIO = 101
N_MOTOR_REVS_WINDOW = 10        # 窗口長度：10 個馬達圈 = 10 × 360/101 = 35.64°
WINDOW_CENTER_DEG = 180.0       # 窗口中心（絕對角度）
SETTLE_TIME = 0.5               # 進入窗口前，在目標速度至少穩定的時間 [s]
OVERRUN_DEG = 0.5               # 越過窗口終點後再走多少才 stopj [deg]

SPEEDJ_ACCEL = 0.5              # speedj 加速度 [rad/s^2]
STOPJ_DECEL = 1.0               # stopj 減速度 [rad/s^2]
MOVEJ_ACCEL = 0.5               # 段間 movej 加速度 [rad/s^2]
MOVEJ_VEL = 0.15                # 段間 movej 速度 [rad/s]
PAUSE_TIME = 0.5                # 每段 stopj 後、movej 後的靜止時間 [s]（須小於結束判定的 1 s 靜止窗）

QD_MAX = 0.30                   # 速度上限；J0 需低於約 0.32 rad/s，否則每圈 12 次的漣波在 125 Hz 下混疊
QDD_CHECK_LIMIT = 3.0           # 安全檢查用加速度上限
LOOP_TIME_FACTOR = 1.5          # URScript 迴圈次數上限 = 預估時間 × 此係數 / dt
SAFETY_MARGIN_DEG = 2.0         # 絕對角度檢查時，起點與最遠點各再外推的餘裕 [deg]
MAX_INITIAL_MOVE_DEG = 40.0     # 目前 J0 到第一段起點的最大允許距離 [deg]

SAMPLE_HZ = 125.0
DT = 1.0 / SAMPLE_HZ
KT_OUT = 101 * 0.1350

# J0 允許的絕對角度活動視窗。必須由操作者依現場線纜與淨空狀況填寫，
# 未填寫則程式中止 —— 不提供預設值。
J0_SAFE_MIN_DEG = None
J0_SAFE_MAX_DEG = None

OUTPUT_DIR = "."


# ============================================================
# 第一部分：角度窗口與每段規劃
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
    每一段的絕對角度規劃（單位：度）。direction=+1 從 THETA_A 掃到 THETA_B，-1 反之。
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


def estimate_total_time(segs, q0_deg, decel):
    total, pos = 0.0, q0_deg
    for s in segs:
        total += estimate_movej_time(math.radians(abs(s["start_deg"] - pos)), MOVEJ_VEL, MOVEJ_ACCEL)
        total += PAUSE_TIME + s["t_expect"] + s["speed"] / decel + PAUSE_TIME
        pos = s["extreme_deg"]
    return total


def safety_check(segs, qd_max, qdd_max, accel, decel, movej_vel, movej_accel,
                 safe_min_deg, safe_max_deg, margin_deg):
    """
    回傳 (ok, lines, (planned_min_deg, planned_max_deg))。
    絕對角度檢查：所有段的起點與煞車後最遠點，各外推 margin_deg，都必須在 J0 視窗內。
    J0 視窗未設定（None）一律 FAIL（預設拒絕）。
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


# ============================================================
# 第二部分：URScript 產生
# ============================================================

def build_urscript(segs, direction, joint_index, q_others, accel, decel,
                   movej_accel, movej_vel, pause_time, dt):
    """
    q_others：六軸絕對角度 [rad]（取自執行開始時 RTDE 量到的值），只替換 J0 當作各段 movej 目標，
    其他五軸維持原值不動。
    每段：movej → 停 → 迴圈（每週期 speedj 一次，讀 J0 實際角度，越過 stop 角度或達次數上限即離開）
          → stopj → 停。
    """
    lines = ["def const_velocity_angle_program():"]
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
    """同 ur5_const_accel_ident_v2_movej.py：J0 實際速度連續 quiet_s 秒低於門檻才視為結束。"""
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
# 第三部分：RTDE 記錄（同 ur5_const_accel_ident_v2_movej.py）
# ============================================================

class RtdeLogger:
    def __init__(self, robot_ip, sample_hz, csv_path, joint_index):
        self.robot_ip = robot_ip
        self.period = 1.0 / sample_hz
        self.csv_path = csv_path
        self.joint_index = joint_index
        self._stop_flag = threading.Event()
        self._thread = None
        self.rows = 0
        self.loops = 0
        self.gaps = 0
        self.ctrl_times = []
        self.last_q0 = None
        self.last_qd0 = None
        self.last_q_full = None
        self.q0_min = None
        self.q0_max = None
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
        writer = None
        last_ctrl_ts = None
        t0_ctrl = None
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
                q = rtde_r.getActualQ()
                qd = rtde_r.getActualQd()
                i = rtde_r.getActualCurrent()
                tqd = rtde_r.getTargetQd()
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
                self.last_q0 = q[j]
                self.last_qd0 = qd[j]
                self.last_q_full = list(q)
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
        lines = [f"迴圈次數 : {self.loops}, 寫入列數 : {self.rows}"]
        if self.loops > 0:
            dup = self.loops - self.rows
            lines.append(f"重複frame濾掉 : {dup} 次（{100*dup/self.loops:.1f}%）")
        lines.append(f"疑似漏拍 : {self.gaps} 次")
        if len(self.ctrl_times) >= 3:
            d = [b - a for a, b in zip(self.ctrl_times[:-1], self.ctrl_times[1:])]
            lines.append(f"平均間隔 : {statistics.fmean(d)*1000:.4f} ms, "
                         f"std {statistics.pstdev(d)*1000:.4f} ms")
        return "\n".join(lines)


# ============================================================
# 第四部分：離線分析
# ============================================================

def analyze_angle_window(csv_path, joint_index, kt_out, direction, speed_levels,
                         theta_a_deg, theta_b_deg, gear_ratio=GEAR_RATIO):
    """
    每一段：取 target_qd 等於該檔速度、且實際 J0 角度落在窗口 [THETA_A, THETA_B] 內的樣本。
    窗口剛好是整數個馬達圈，漣波在平均中抵銷，不需要再取整數圈。
    以各段平均（實測 qd、tau）回歸 tau = B*qd + Tc*direction；± 為以段為點的回歸標準誤。
    回傳 dict：segments（每段統計）、B、Tc、se_B、se_Tc，以及每段的原始樣本（畫電流對角度用）。
    """
    with open(csv_path, encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    q = np.degrees(np.array([float(r[f"actual_q_{joint_index}"]) for r in rows]))
    qd = np.array([float(r[f"actual_qd_{joint_index}"]) for r in rows])
    tqd = np.array([float(r[f"target_qd_{joint_index}"]) for r in rows])
    tau = np.array([float(r[f"actual_current_{joint_index}"]) for r in rows]) * kt_out
    t = np.array([float(r["timestamp"]) for r in rows])
    d = float(direction)
    in_win = (q >= theta_a_deg) & (q <= theta_b_deg)
    rev_deg = 360.0 / gear_ratio

    segments = []
    for v in speed_levels:
        on = np.abs(tqd - d * v) < 1e-4
        idx = np.where(on & in_win)[0]
        if len(idx) == 0:
            continue
        runs = np.split(idx, np.where(np.diff(idx) != 1)[0] + 1)
        for r in runs:
            if len(r) < 10:
                continue
            cover = (q[r].max() - q[r].min()) / rev_deg
            segments.append(dict(speed=v, idx=r, t0=float(t[r[0]]), n=len(r), revs=float(cover),
                                 qd_mean=float(np.mean(qd[r])), tau_mean=float(np.mean(tau[r]))))
    segments.sort(key=lambda s: s["t0"])
    if len(segments) < 2:
        print("[警告] 有效段數少於 2，無法回歸。")
        return dict(segments=segments)

    x = np.array([s["qd_mean"] for s in segments])
    y = np.array([s["tau_mean"] for s in segments])
    X = np.column_stack([x, np.full(len(x), d)])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    res = y - X @ coef
    dof = max(len(y) - 2, 1)
    se = np.sqrt(np.diag(res @ res / dof * np.linalg.inv(X.T @ X)))

    print(f"{'段':>3}{'速度':>7}{'樣本':>6}{'覆蓋馬達圈':>10}{'平均qd':>10}{'平均tau':>10}")
    for k, s in enumerate(segments):
        print(f"{k+1:3d}{s['speed']:7.2f}{s['n']:6d}{s['revs']:10.2f}{s['qd_mean']:10.4f}{s['tau_mean']:10.3f}")
    print(f"B = {coef[0]:.3f} ± {se[0]:.3f}，Tc = {coef[1]:.3f} ± {se[1]:.3f}（{len(y)} 段）")
    return dict(segments=segments, B=float(coef[0]), Tc=float(coef[1]),
                se_B=float(se[0]), se_Tc=float(se[1]), q_deg=q, tau=tau, t=t)


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
    print(" 固定角度窗口定速實驗：安全檢查")
    print("=" * 70)
    print(f"方向: {'正轉 (+1)' if DIRECTION > 0 else '反轉 (-1)'}")

    levels, order, n_revs = stage_config(TEST_STAGE)
    speeds = speed_sequence(levels, order)
    theta_a, theta_b = window_bounds_deg(WINDOW_CENTER_DEG, n_revs, GEAR_RATIO)
    segs = plan_segments(speeds, DIRECTION, theta_a, theta_b, SPEEDJ_ACCEL, STOPJ_DECEL,
                         SETTLE_TIME, OVERRUN_DEG, LOOP_TIME_FACTOR, DT)
    print(f"模式: {TEST_STAGE}，速度順序: {[f'{v:.2f}' for v in speeds]}")
    print(f"角度窗口: {theta_a:.3f}° 至 {theta_b:.3f}°（{n_revs} 個馬達圈 = {theta_b - theta_a:.3f}°）")

    ok, report, (lo_plan, hi_plan) = safety_check(
        segs, QD_MAX, QDD_CHECK_LIMIT, SPEEDJ_ACCEL, STOPJ_DECEL, MOVEJ_VEL, MOVEJ_ACCEL,
        J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG, SAFETY_MARGIN_DEG)
    print("\n".join(report))
    if not ok:
        print("\n[中止] 安全檢查未通過。")
        return

    print("\n[通過安全檢查]")
    if ROBOT_IP is None:
        print("\n[中止] ROBOT_IP 未設定。")
        return
    if rtde_receive is None:
        print("\n[中止] 找不到 rtde_receive，僅完成離線設計檢查。")
        return

    print("\n" + "=" * 70)
    print(" 執行前姿態確認（診斷原點姿態，J1~J5）")
    print("=" * 70)
    if home_pose is None:
        print("[警告] 找不到 ur5_home_pose 模組，無法自動檢查姿態！")
        print("        請自己手動確認 J1~J5 是否為：-90°, 90°, -90°, -90°, 0°")
        ans = input("已手動確認姿態正確，輸入 yes 繼續，其他任何輸入則取消：").strip().lower()
        if ans != "yes":
            print("已取消。")
            return
    else:
        if not home_pose.ensure_home_pose(ROBOT_IP, mode="CHECK"):
            print("\n[中止] 姿態不符合診斷原點，請先調整姿態後再執行本實驗。")
            return

    print(f"\n即將對 IP={ROBOT_IP} 送出軌跡，關節 J{JOINT_INDEX}，方向 {DIRECTION:+d}，共 {len(segs)} 段。")
    print("*** 送出後 Ctrl+C 會送出 stopj，但最可靠的仍是緊急停止按鈕 ***")
    input("按 Enter 繼續，或 Ctrl+C 取消...")

    timestamp_str = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_label = "pos" if DIRECTION > 0 else "neg"
    file_tag = f"{dir_label}_{TEST_STAGE.lower()}_{timestamp_str}"
    session_dir = os.path.join(OUTPUT_DIR, f"constvel_angle_{file_tag}")
    os.makedirs(session_dir, exist_ok=True)
    csv_path = os.path.join(session_dir, f"constvel_angle_data_{file_tag}.csv")
    info_path = os.path.join(session_dir, f"constvel_angle_info_{file_tag}.txt")

    logger = RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, JOINT_INDEX)
    logger.start()
    if not logger.started_ok.wait(5.0):
        logger.stop()
        print(f"\n[中止] RTDE 記錄未能啟動：{logger.error}")
        print("       未送出任何軌跡，手臂未動作。")
        return

    q0_deg = math.degrees(logger.last_q0)
    pos_ok, pos_msg = check_start_position(q0_deg, segs[0]["start_deg"], MAX_INITIAL_MOVE_DEG,
                                           J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG)
    print(f"\n[起始位置檢查] {pos_msg}")
    if not pos_ok:
        logger.stop()
        print("\n[中止] 起始位置檢查未通過，未送出任何軌跡，手臂未動作。")
        return

    script_text = build_urscript(segs, DIRECTION, JOINT_INDEX, logger.last_q_full,
                                 SPEEDJ_ACCEL, STOPJ_DECEL, MOVEJ_ACCEL, MOVEJ_VEL, PAUSE_TIME, DT)
    est_total = estimate_total_time(segs, q0_deg, STOPJ_DECEL)

    status, elapsed, motion_span = "unknown", 0.0, 0.0
    try:
        print(f"[資訊] 送出 URScript（預估總時長 {est_total:.1f} 秒）...")
        send_urscript(ROBOT_IP, script_text)
        status, elapsed, motion_span = wait_for_motion_complete(logger, est_total)
        if status == "no_motion":
            print(f"\n[警告] 送出後 {elapsed:.1f}s 內未偵測到 J0 運動。")
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
            f"實測 J0 範圍 : {q_lo:.3f}° 至 {q_hi:.3f}°（規劃含餘裕 {lo_plan:.3f}° 至 {hi_plan:.3f}°）",
            f"結束位置     : {math.degrees(logger.last_q0):.3f}°",
            f"實際運動時長 / 預估 : {motion_span:.2f}s / {est_total:.2f}s",
        ]
        if q_lo < lo_plan or q_hi > hi_plan:
            measured_lines.append("[警告] 實測範圍超出規劃範圍，請先查明原因再跑下一階段。")
        print("\n" + "\n".join(measured_lines))

    with open(info_path, "w", encoding="utf-8") as f:
        f.write(f"方向: {DIRECTION:+d}\n模式: {TEST_STAGE}\n")
        f.write(f"速度順序: {speeds}\n")
        f.write(f"角度窗口: {theta_a:.4f} 至 {theta_b:.4f} deg（{n_revs} 個馬達圈）\n")
        f.write(f"預估總時長: {est_total:.1f}s\n執行狀態: {status}\n")
        f.write(logger.summary() + "\n")
        if measured_lines:
            f.write("\n" + "\n".join(measured_lines) + "\n")

    print(f"\n[完成] 資料存至 {csv_path}")
    print(f"[完成] 摘要存至 {info_path}")
    print(f"\n下一步：analyze_angle_window(csv_path, JOINT_INDEX, KT_OUT, DIRECTION, "
          f"SPEED_LEVELS, {theta_a:.4f}, {theta_b:.4f})")


if __name__ == "__main__":
    main()
