# -*- coding: utf-8 -*-
"""
UR5 自動換姿態（取代用教導器手動擺姿勢）
================================================================

用法（在 repo 根目錄）：
    python robot/ur5_goto_pose.py UP        # J1 速度調變（#5、#6）、J4/J5 的 UP
    python robot/ur5_goto_pose.py UP_J3     # J3 定速／速度調變的 UP（#7、#8）：UP 再把 J3 轉到 −130°
    python robot/ur5_goto_pose.py SIX       # 六軸同時鑑別（ur5_six_axis_ident.py）的開始姿態
    python robot/ur5_goto_pose.py DIAG      # 診斷原點（收工或回到 DIAG 系列）
不給參數時用下面的 TARGET。

姿態與各鑑別程式 POSES 裡的 q_deg 完全相同（離線測試會逐一比對），所以到位後鑑別程式的
±0.5° 姿態檢查一定通過。

安全設計
--------
1. 只用一次 movej，慢速（MOVE_VEL 0.3 rad/s、MOVE_ACC 0.3 rad/s²）。
2. 目標角度 = 目前角度 + 「換算到 ±180° 內的差值」：不會因為讀到 −180° 與 180° 而整圈轉。
3. 任一軸要轉超過 MAX_MOVE_DEG 就拒絕（DIAG→UP 是 J2 轉 90°；UP→UP_J3 是 J3 轉 40°）。
4. movej 在關節空間直線插值，程式先沿整條路徑算 J2 之後各座標系（前臂、手腕、法蘭）離底座安裝面的高度，
   最低點低於 MIN_Z 就拒絕。DIAG、UP、UP_J3 之間互換，路徑最低點都是法蘭在 DIAG 時（約 0.43 m）。
5. 動之前印出各軸目前角度、目標、要轉多少與預估時間，輸入 y 才動。
6. 移動中按 Ctrl+C 會送 stopj；等待逾時也會送 stopj。到位判定：六軸都在 ±0.5° 內且速度幾乎為 0。
7. 已經在目標姿態（±0.5°）就不動，直接結束。
"""

import math
import socket
import sys
import time

import numpy as np

try:
    import rtde_receive
except ImportError:  # 離線測試或沒裝 ur_rtde
    rtde_receive = None


# ============================================================
# 參數
# ============================================================
ROBOT_IP = "192.168.50.114"   # 執行前請再次確認
TARGET = "UP"                 # 沒給指令參數時的目標姿態

POSES = {
    "DIAG": (180.0, -90.0, 90.0, -90.0, -90.0, 0.0),    # 診斷原點：前臂水平
    "UP": (180.0, -90.0, 0.0, -90.0, -90.0, 0.0),       # 手臂朝上、工具朝前水平（J1 速度調變、J4/J5 的 UP）
    "UP_J3": (180.0, -90.0, 0.0, -130.0, -90.0, 0.0),   # UP 再把 J3 轉到 −130°（J3 定速／速度調變的 UP）
    "SIX": (180.0, -90.0, 70.0, -70.0, -90.0, 0.0),     # 六軸同時鑑別（ur5_six_axis_ident.py）：前臂往上 20°、工具垂直朝下
}

MOVE_VEL = 0.3          # movej 速度 [rad/s]
MOVE_ACC = 0.3          # movej 加速度 [rad/s^2]
MAX_MOVE_DEG = 100.0    # 任一軸最多轉多少 [deg]
MIN_Z = 0.30            # 路徑上前臂、手腕、法蘭離底座安裝面的最低高度 [m]
POSE_TOL_DEG = 0.5      # 到位容差（同各鑑別程式）
QD_STILL = 2e-3         # 視為靜止的關節速度 [rad/s]
PATH_STEPS = 200        # 檢查路徑高度的取樣點數

DH_D = (0.089159, 0.0, 0.0, 0.10915, 0.09465, 0.0823)
DH_A = (0.0, -0.425, -0.39225, 0.0, 0.0, 0.0)
DH_ALPHA = (math.pi / 2, 0.0, 0.0, math.pi / 2, -math.pi / 2, 0.0)
Z_FRAMES = (2, 3, 4, 5, 6)   # 檢查高度的座標系：J2 之後（肘）、J3、J4、J5 之後、法蘭


# ============================================================
# 計算
# ============================================================
def fk_frames(q_rad):
    """UR5 官方 DH 的各座標系（4×4），frames[k] 為第 k 軸之後的座標系；frames[6] 為法蘭。"""
    T = np.eye(4)
    frames = [T.copy()]
    for th, d, a, al in zip(q_rad, DH_D, DH_A, DH_ALPHA):
        ct, st, ca, sa = math.cos(th), math.sin(th), math.cos(al), math.sin(al)
        T = T @ np.array([[ct, -st * ca, st * sa, a * ct],
                          [st, ct * ca, -ct * sa, a * st],
                          [0.0, sa, ca, d],
                          [0.0, 0.0, 0.0, 1.0]])
        frames.append(T.copy())
    return frames


def wrap180(d):
    """角度差換到 [−180, 180)。"""
    return (d + 180.0) % 360.0 - 180.0


def plan_move(cur_deg, target_name, poses=POSES):
    """
    回傳 (goal_deg, delta_deg)：goal_deg = 目前角度 + 換到 ±180° 內的差值（movej 的實際目標）。
    target_name 不在 poses 或目前角度不是有限值時丟 ValueError。
    """
    if target_name not in poses:
        raise ValueError(f"未知的姿態 {target_name!r}（可用：{', '.join(poses)}）")
    if len(cur_deg) != 6 or not all(math.isfinite(x) for x in cur_deg):
        raise ValueError(f"讀到的六軸角度不正確：{cur_deg}")
    delta = [wrap180(t - c) for c, t in zip(cur_deg, poses[target_name])]
    goal = [c + d for c, d in zip(cur_deg, delta)]
    return goal, delta


def path_min_z(start_deg, goal_deg, steps=PATH_STEPS, frames=Z_FRAMES):
    """movej 關節空間直線路徑上，frames 各座標系原點離底座安裝面的最低高度 [m]。"""
    z_min = float("inf")
    for i in range(steps + 1):
        s = i / steps
        q = [math.radians(a + s * (b - a)) for a, b in zip(start_deg, goal_deg)]
        fr = fk_frames(q)
        z_min = min(z_min, min(fr[k][2, 3] for k in frames))
    return z_min


def estimate_move_time(max_delta_rad, vel=MOVE_VEL, acc=MOVE_ACC):
    """movej 梯形速度的時間估計 [s]。"""
    d = abs(max_delta_rad)
    if d <= vel * vel / acc:
        return 2.0 * math.sqrt(d / acc)
    return d / vel + vel / acc


def safety_check(cur_deg, goal_deg, delta_deg, max_move_deg=MAX_MOVE_DEG, min_z=MIN_Z):
    """回傳 (ok, lines)。任一軸轉太多或路徑太低就拒絕。"""
    lines, ok = [], True
    lines.append(f"{'軸':<4}{'目前[°]':>10}{'目標[°]':>10}{'要轉[°]':>10}")
    for j in range(6):
        lines.append(f"J{j:<3}{cur_deg[j]:>10.2f}{goal_deg[j]:>10.2f}{delta_deg[j]:>10.2f}")
    max_d = max(abs(d) for d in delta_deg)
    if max_d > max_move_deg:
        ok = False
        lines.append(f"[拒絕] 最大轉動 {max_d:.1f}° 超過上限 {max_move_deg}°：請先用教導器移到目標附近")
    z = path_min_z(cur_deg, goal_deg)
    if z < min_z:
        ok = False
        lines.append(f"[拒絕] 路徑上前臂／手腕／法蘭最低 {z:.3f} m，低於 {min_z} m")
    else:
        lines.append(f"路徑上前臂／手腕／法蘭最低 {z:.3f} m（下限 {min_z} m）")
    lines.append(f"最大轉動 {max_d:.1f}°，預估 {estimate_move_time(math.radians(max_d)):.1f} s"
                 f"（{MOVE_VEL} rad/s、{MOVE_ACC} rad/s²，教導器速度滑桿會再等比例放慢）")
    return ok, lines


def at_pose(cur_deg, target_deg, tol_deg=POSE_TOL_DEG):
    return all(abs(wrap180(c - t)) <= tol_deg for c, t in zip(cur_deg, target_deg))


def build_urscript(goal_deg, acc=MOVE_ACC, vel=MOVE_VEL):
    """第一行 def 頂格、其餘行至少縮排一格、最後一行 end 頂格。"""
    q = ", ".join(f"{math.radians(x):.6f}" for x in goal_deg)
    return ("def goto_pose_program():\n"
            f"  movej([{q}], a={acc}, v={vel})\n"
            "  stopj(1.0)\n"
            "end\n")


# ============================================================
# 與手臂通訊
# ============================================================
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


def main(argv=None, robot_ip=ROBOT_IP):
    argv = sys.argv[1:] if argv is None else argv
    target = argv[0].strip().upper() if argv else TARGET
    if target not in POSES:
        print(f"未知的姿態 {target!r}（可用：{', '.join(POSES)}）")
        return 1
    if rtde_receive is None:
        print("找不到 rtde_receive（ur_rtde），無法讀取手臂角度")
        return 1

    try:
        r = rtde_receive.RTDEReceiveInterface(robot_ip)
    except Exception as e:
        print(f"RTDE 連線失敗（{robot_ip}）：{e}")
        return 1
    try:
        cur = [math.degrees(x) for x in r.getActualQ()]
        try:
            goal, delta = plan_move(cur, target)
        except ValueError as e:
            print(e)
            return 1
        print(f"目標姿態 {target}：{list(POSES[target])}°")
        if at_pose(cur, POSES[target]):
            print("[OK] 已經在目標姿態（±0.5°），不用動。")
            return 0
        ok, lines = safety_check(cur, goal, delta)
        print("\n".join(lines))
        if not ok:
            return 1
        print("\n確認：手臂周圍淨空、工具線材有餘裕、緊急停止在手邊。")
        if input("輸入 y 開始移動（其他輸入取消）：").strip().lower() != "y":
            print("已取消，手臂未移動。")
            return 1

        send_urscript(robot_ip, build_urscript(goal))
        t_lim = 3.0 * estimate_move_time(math.radians(max(abs(d) for d in delta))) + 10.0
        t0 = time.time()
        print("移動中…（Ctrl+C 停止）")
        try:
            while True:
                time.sleep(0.2)
                q = [math.degrees(x) for x in r.getActualQ()]
                qd = r.getActualQd()
                if at_pose(q, POSES[target]) and max(abs(x) for x in qd) < QD_STILL:
                    print(f"[完成] 已到 {target}：{[round(x, 2) for x in q]}°")
                    return 0
                if time.time() - t0 > t_lim:
                    print(f"[逾時] {t_lim:.0f} s 內沒有到位（教導器速度滑桿太低或程式被停止？）")
                    send_abort(robot_ip)
                    print(f"目前角度：{[round(x, 2) for x in q]}°")
                    return 1
        except KeyboardInterrupt:
            send_abort(robot_ip)
            return 1
    finally:
        try:
            r.disconnect()
        except Exception:
            pass


if __name__ == "__main__":
    sys.exit(main())
