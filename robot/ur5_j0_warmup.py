"""
UR5 J0 暖機程式（獨立執行，暖機完成後再執行 ur5_const_velocity_angle_ident.py）
================================================================

動作：J0 以 movej 在 WARMUP_LO_DEG 與 WARMUP_HI_DEG 之間來回轉動，其他五軸不動。
      這段角度涵蓋固定角度窗口實驗兩個方向的規劃範圍（約 152°～208°），
      不會帶手臂去實驗以外的位置。

分批：每批約 BLOCK_MINUTES 分鐘，一次送出一批 URScript；一批跑完、手臂停下後，
      程式計算：
        - J0 溫度（joint_temp_0）
        - 摩擦指標：該批中「以 WARMUP_VEL 定速巡航」樣本的平均 |力矩|，正轉、反轉分開
      再決定要不要送下一批。

結束條件（暖機完成）：
  1. 已暖機至少 MIN_MINUTES 分鐘，且
  2. 最近 TEMP_WINDOW_MIN 分鐘內 J0 溫度變化 < TEMP_STABLE_C，且
  3. 最近 N_STABLE_BLOCKS 批的摩擦指標（兩個方向）最大與最小值相差 < FRICTION_STABLE_FRAC
  未達成但累計超過 MAX_MINUTES 分鐘也結束（並註明「未達穩定」）。
  結束時 movej 回 END_DEG（預設 180°），方便直接接著跑正式實驗。

★★★ 使用前必讀 ★★★
1. 【新腳本，尚未經學長審查，也尚未上機】先用 TEST_STAGE = "QUICK"（只送 1 批、2 次往返）。
2. 【安全檢查預設拒絕】J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG / TEMP_ABORT_C 未設定一律中止。
3. 【溫度上限】任一批結束時 J0 溫度 ≥ TEMP_ABORT_C 就停止，不再送下一批。
4. 【Ctrl+C】會送 stopj；最可靠的仍是緊急停止按鈕。
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

# QUICK：只送 1 批、2 次往返，確認動作與記錄正常
# FULL ：分批暖機直到穩定或達 MAX_MINUTES
TEST_STAGE = "QUICK"

WARMUP_LO_DEG = 155.0         # 來回的兩端（絕對角度）
WARMUP_HI_DEG = 205.0
END_DEG = 180.0               # 暖機結束後停在這裡
WARMUP_VEL = 0.15             # movej 速度 [rad/s]
WARMUP_ACCEL = 0.5            # movej 加速度 [rad/s^2]

BLOCK_MINUTES = 2.0           # 每批長度（約）
MIN_MINUTES = 5.0
MAX_MINUTES = 30.0
TEMP_WINDOW_MIN = 5.0
TEMP_STABLE_C = 0.1
N_STABLE_BLOCKS = 3
FRICTION_STABLE_FRAC = 0.01
CRUISE_TOL = 0.02             # |qd| 與 WARMUP_VEL 相差 2% 以內視為定速巡航

QD_MAX = 0.30
QDD_CHECK_LIMIT = 3.0
SAFETY_MARGIN_DEG = 2.0
MAX_INITIAL_MOVE_DEG = 40.0

SAMPLE_HZ = 125.0
KT_OUT = 101 * 0.1350

# 以下三項未設定則中止（預設拒絕）
J0_SAFE_MIN_DEG = None
J0_SAFE_MAX_DEG = None
TEMP_ABORT_C = None           # J0 溫度上限 [°C]，由操作者依手冊與現場決定

OUTPUT_DIR = "."


# ============================================================
# 規劃與安全檢查
# ============================================================

def movej_time(dist_rad, v, a):
    if dist_rad <= 0:
        return 0.0
    t_up = v / a
    d_up = 0.5 * a * t_up ** 2
    if dist_rad <= 2 * d_up:
        return 2.0 * math.sqrt(dist_rad / a)
    return 2.0 * t_up + (dist_rad - 2 * d_up) / v


def cycles_per_block(block_min, lo_deg, hi_deg, v, a):
    """一次往返（lo→hi→lo）的時間，換算每批幾次往返（至少 1）。"""
    t_cycle = 2 * movej_time(math.radians(hi_deg - lo_deg), v, a)
    return max(1, int(round(block_min * 60.0 / t_cycle))), t_cycle


def safety_check(lo_deg, hi_deg, end_deg, v, a, qd_max, qdd_max,
                 safe_min_deg, safe_max_deg, temp_abort_c, margin_deg):
    lines, ok = [], True

    def chk(name, val, lim, unit):
        nonlocal ok
        lines.append(f"{name} : {val:.4f} {unit}  (限 {lim})")
        if val > lim:
            lines.append("  [FAIL]"); ok = False
        else:
            lines.append("  [OK]")

    chk("movej 速度", v, qd_max, "rad/s")
    chk("movej 加速度", a, qdd_max, "rad/s^2")
    if not lo_deg < hi_deg:
        lines.append("  [FAIL] WARMUP_LO_DEG 必須小於 WARMUP_HI_DEG"); ok = False
    if not lo_deg <= end_deg <= hi_deg:
        lines.append("  [FAIL] END_DEG 必須在來回範圍內"); ok = False
    lo_p, hi_p = lo_deg - margin_deg, hi_deg + margin_deg
    lines.append(f"規劃絕對角度範圍（含餘裕 {margin_deg:.1f}°） : {lo_p:.2f}° 至 {hi_p:.2f}°")
    if safe_min_deg is None or safe_max_deg is None:
        lines.append("  [FAIL] J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 尚未設定，拒絕執行"); ok = False
    else:
        lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
        if lo <= lo_p and hi_p <= hi:
            lines.append(f"  [OK] 在 J0 視窗 [{lo:.1f}°, {hi:.1f}°] 內")
        else:
            lines.append(f"  [FAIL] 超出 J0 視窗 [{lo:.1f}°, {hi:.1f}°]"); ok = False
    if temp_abort_c is None:
        lines.append("  [FAIL] TEMP_ABORT_C 尚未設定，拒絕執行"); ok = False
    else:
        lines.append(f"J0 溫度上限 : {temp_abort_c:.1f} °C  [OK]")
    return ok, lines


def check_start_position(q0_deg, first_target_deg, max_move_deg, safe_min_deg, safe_max_deg):
    if safe_min_deg is None or safe_max_deg is None:
        return False, "J0_SAFE_MIN_DEG / J0_SAFE_MAX_DEG 尚未設定，拒絕執行"
    lo, hi = min(safe_min_deg, safe_max_deg), max(safe_min_deg, safe_max_deg)
    move = abs(first_target_deg - q0_deg)
    msg = (f"J0 目前 {q0_deg:.2f}°，第一個目標 {first_target_deg:.2f}°（需移動 {move:.2f}°，"
           f"上限 {max_move_deg:.1f}°），J0 視窗 [{lo:.1f}°, {hi:.1f}°]")
    if not (lo <= q0_deg <= hi):
        return False, msg + " -> 目前角度在視窗外"
    if move > max_move_deg:
        return False, msg + " -> 第一次 movej 距離過大"
    return True, msg + " -> OK"


# ============================================================
# URScript
# ============================================================

def _q_str(q_others, joint_index, deg):
    q = list(q_others)
    q[joint_index] = math.radians(deg)
    return "[" + ", ".join(f"{v:.6f}" for v in q) + "]"


def build_block_urscript(n_cycles, lo_deg, hi_deg, q_others, joint_index, v, a):
    """一批：先到 lo，再來回 n_cycles 次（lo→hi→lo），只改 J0。"""
    qa, qb = _q_str(q_others, joint_index, lo_deg), _q_str(q_others, joint_index, hi_deg)
    lines = ["def j0_warmup_block():",
             f"  movej({qa}, a={a}, v={v})",
             "  k = 0",
             f"  while k < {n_cycles}:",
             f"    movej({qb}, a={a}, v={v})",
             f"    movej({qa}, a={a}, v={v})",
             "    k = k + 1",
             "  end",
             "end"]
    return "\n".join(lines) + "\n"


def build_end_urscript(end_deg, q_others, joint_index, v, a):
    return ("def j0_warmup_end():\n"
            f"  movej({_q_str(q_others, joint_index, end_deg)}, a={a}, v={v})\n"
            "end\n")


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
    """同其他腳本：J0 實際速度連續 quiet_s 秒低於門檻才視為結束。"""
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
# RTDE 記錄（多一個記憶體緩衝，供每批計算溫度與摩擦指標）
# ============================================================

class RtdeLogger:
    def __init__(self, robot_ip, sample_hz, csv_path, joint_index):
        self.robot_ip = robot_ip
        self.period = 1.0 / sample_hz
        self.csv_path = csv_path
        self.joint_index = joint_index
        self._stop_flag = threading.Event()
        self._thread = None
        self._lock = threading.Lock()
        self.rows = self.loops = self.gaps = 0
        self.ctrl_times = []
        self.last_q0 = self.last_qd0 = self.last_q_full = None
        self.q0_min = self.q0_max = None
        self.buf_t, self.buf_qd, self.buf_cur, self.buf_temp = [], [], [], []
        self.started_ok = threading.Event()
        self.error = None

    def start(self):
        self._thread = threading.Thread(target=self._run, daemon=True)
        self._thread.start()

    def stop(self):
        self._stop_flag.set()
        if self._thread is not None:
            self._thread.join(timeout=10.0)

    def snapshot(self, since_index=0):
        with self._lock:
            return (np.array(self.buf_t[since_index:]), np.array(self.buf_qd[since_index:]),
                    np.array(self.buf_cur[since_index:]), np.array(self.buf_temp[since_index:], dtype=float))

    def n_buffered(self):
        with self._lock:
            return len(self.buf_t)

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
                with self._lock:
                    self.buf_t.append(t_main); self.buf_qd.append(qd[j]); self.buf_cur.append(i[j])
                    self.buf_temp.append(temps[j] if temps[j] is not None else float("nan"))
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
# 每批評估與結束判定（純函式，可離線測試）
# ============================================================

def block_metrics(qd, cur, temp, v, tol, kt_out):
    """回傳 (正轉巡航平均|tau|, 反轉巡航平均|tau|, 批末溫度)；沒有巡航樣本則為 nan。"""
    fwd = np.abs(qd - v) <= tol * v
    rev = np.abs(qd + v) <= tol * v
    f = float(np.mean(np.abs(cur[fwd])) * kt_out) if fwd.any() else float("nan")
    r = float(np.mean(np.abs(cur[rev])) * kt_out) if rev.any() else float("nan")
    tt = temp[~np.isnan(temp)]
    return f, r, (float(tt[-1]) if len(tt) else float("nan"))


def is_warm(history, elapsed_min, min_minutes, temp_window_min, temp_stable_c,
            n_stable_blocks, friction_stable_frac):
    """
    history：每批一筆 dict(t_min=批結束時累計分鐘, temp, fwd, rev)。
    回傳 (是否完成, 說明)。
    """
    if elapsed_min < min_minutes:
        return False, f"未滿最少暖機時間 {min_minutes:.0f} 分鐘"
    recent = [h for h in history if h["t_min"] >= elapsed_min - temp_window_min - 1e-9]
    temps = [h["temp"] for h in recent if not math.isnan(h["temp"])]
    if len(temps) < 2:
        return False, "溫度資料不足"
    d_temp = max(temps) - min(temps)
    if d_temp >= temp_stable_c:
        return False, f"最近 {temp_window_min:.0f} 分鐘溫度變化 {d_temp:.2f} °C ≥ {temp_stable_c} °C"
    if len(history) < n_stable_blocks:
        return False, "批數不足"
    last = history[-n_stable_blocks:]
    for key, name in (("fwd", "正轉"), ("rev", "反轉")):
        vals = [h[key] for h in last]
        if any(math.isnan(x) for x in vals):
            return False, f"{name}摩擦指標缺資料"
        spread = (max(vals) - min(vals)) / np.mean(vals)
        if spread >= friction_stable_frac:
            return False, f"最近 {n_stable_blocks} 批{name}摩擦指標變化 {100*spread:.2f}% ≥ {100*friction_stable_frac:.1f}%"
    return True, f"溫度變化 {d_temp:.2f} °C、摩擦指標兩方向皆穩定"


# ============================================================
# 主流程
# ============================================================

def main():
    print("=" * 70)
    print(" UR5 J0 暖機：安全檢查")
    print("=" * 70)
    n_cyc, t_cycle = cycles_per_block(BLOCK_MINUTES, WARMUP_LO_DEG, WARMUP_HI_DEG, WARMUP_VEL, WARMUP_ACCEL)
    if TEST_STAGE == "QUICK":
        n_cyc, max_blocks = 2, 1
    elif TEST_STAGE == "FULL":
        max_blocks = int(math.ceil(MAX_MINUTES * 60.0 / (n_cyc * t_cycle)))
    else:
        raise ValueError(f"未知的 TEST_STAGE: {TEST_STAGE}")
    print(f"模式: {TEST_STAGE}；來回 {WARMUP_LO_DEG:.1f}° ↔ {WARMUP_HI_DEG:.1f}°，{WARMUP_VEL} rad/s；"
          f"每批 {n_cyc} 次往返（約 {n_cyc * t_cycle:.0f} s），最多 {max_blocks} 批")

    ok, report = safety_check(WARMUP_LO_DEG, WARMUP_HI_DEG, END_DEG, WARMUP_VEL, WARMUP_ACCEL,
                              QD_MAX, QDD_CHECK_LIMIT, J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG,
                              TEMP_ABORT_C, SAFETY_MARGIN_DEG)
    print("\n".join(report))
    if not ok:
        print("\n[中止] 安全檢查未通過。")
        return
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

    print(f"\n即將對 IP={ROBOT_IP} 開始暖機（J0）。Ctrl+C 會送 stopj，最可靠的仍是緊急停止按鈕。")
    input("按 Enter 繼續，或 Ctrl+C 取消...")

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    tag = f"{TEST_STAGE.lower()}_{ts}"
    session_dir = os.path.join(OUTPUT_DIR, f"warmup_{tag}")
    os.makedirs(session_dir, exist_ok=True)
    csv_path = os.path.join(session_dir, f"warmup_data_{tag}.csv")
    info_path = os.path.join(session_dir, f"warmup_info_{tag}.txt")

    logger = RtdeLogger(ROBOT_IP, SAMPLE_HZ, csv_path, JOINT_INDEX)
    logger.start()
    if not logger.started_ok.wait(5.0):
        logger.stop()
        print(f"\n[中止] RTDE 記錄未能啟動：{logger.error}。未送出任何軌跡。")
        return

    q0_deg = math.degrees(logger.last_q0)
    pos_ok, pos_msg = check_start_position(q0_deg, WARMUP_LO_DEG, MAX_INITIAL_MOVE_DEG,
                                           J0_SAFE_MIN_DEG, J0_SAFE_MAX_DEG)
    print(f"\n[起始位置檢查] {pos_msg}")
    if not pos_ok:
        logger.stop()
        print("\n[中止] 起始位置檢查未通過，未送出任何軌跡。")
        return

    q_others = logger.last_q_full
    block_text = build_block_urscript(n_cyc, WARMUP_LO_DEG, WARMUP_HI_DEG, q_others, JOINT_INDEX,
                                      WARMUP_VEL, WARMUP_ACCEL)
    t_block_est = movej_time(math.radians(abs(WARMUP_LO_DEG - q0_deg)), WARMUP_VEL, WARMUP_ACCEL) + n_cyc * t_cycle

    history, result, t0 = [], "未完成", time.time()
    print(f"\n{'批':>3}{'累計[min]':>10}{'溫度[°C]':>10}{'正轉|tau|':>11}{'反轉|tau|':>11}  判定")
    try:
        for b in range(max_blocks):
            i0 = logger.n_buffered()
            send_urscript(ROBOT_IP, block_text)
            status, _, _ = wait_for_motion_complete(logger, t_block_est)
            if status != "done":
                print(f"\n[警告] 第 {b+1} 批狀態 {status}，停止暖機。")
                if status == "timeout":
                    send_abort(ROBOT_IP)
                result = f"中止（{status}）"
                break
            _, qd, cur, temp = logger.snapshot(i0)
            f, r, tmp = block_metrics(qd, cur, temp, WARMUP_VEL, CRUISE_TOL, KT_OUT)
            el = (time.time() - t0) / 60.0
            history.append(dict(t_min=el, temp=tmp, fwd=f, rev=r))
            warm, why = is_warm(history, el, MIN_MINUTES, TEMP_WINDOW_MIN, TEMP_STABLE_C,
                                N_STABLE_BLOCKS, FRICTION_STABLE_FRAC)
            print(f"{b+1:3d}{el:10.1f}{tmp:10.2f}{f:11.3f}{r:11.3f}  {why}")
            if not math.isnan(tmp) and tmp >= TEMP_ABORT_C:
                print(f"\n[停止] J0 溫度 {tmp:.2f} °C ≥ 上限 {TEMP_ABORT_C} °C。")
                result = "停止（溫度達上限）"
                break
            if TEST_STAGE == "FULL" and warm:
                result = "暖機完成"
                break
        else:
            result = "QUICK 完成" if TEST_STAGE == "QUICK" else f"達最長時間 {MAX_MINUTES:.0f} 分鐘，未達穩定"

        if not result.startswith("中止"):
            send_urscript(ROBOT_IP, build_end_urscript(END_DEG, q_others, JOINT_INDEX, WARMUP_VEL, WARMUP_ACCEL))
            wait_for_motion_complete(logger, movej_time(math.radians(abs(WARMUP_HI_DEG - WARMUP_LO_DEG)),
                                                        WARMUP_VEL, WARMUP_ACCEL))
    except KeyboardInterrupt:
        print("\n[中止] 收到 Ctrl+C。")
        send_abort(ROBOT_IP)
        result = "中止（Ctrl+C）"
        raise
    except BaseException:
        send_abort(ROBOT_IP)
        result = "中止（例外）"
        raise
    finally:
        logger.stop()
        with open(info_path, "w", encoding="utf-8") as fi:
            fi.write(f"模式: {TEST_STAGE}\n結果: {result}\n")
            fi.write(f"來回: {WARMUP_LO_DEG} 至 {WARMUP_HI_DEG} deg，{WARMUP_VEL} rad/s，每批 {n_cyc} 次往返\n")
            fi.write("批,累計min,溫度C,正轉|tau|,反轉|tau|\n")
            for k, h in enumerate(history):
                fi.write(f"{k+1},{h['t_min']:.2f},{h['temp']:.2f},{h['fwd']:.4f},{h['rev']:.4f}\n")
            fi.write(logger.summary() + "\n")

    print(f"\n[結果] {result}")
    print(f"[完成] 資料存至 {csv_path}")
    print(f"[完成] 摘要存至 {info_path}")
    print("接著請直接執行 ur5_const_velocity_angle_ident.py（避免間隔過久又變冷）。")


if __name__ == "__main__":
    main()
