# 六軸同時鑑別（robot/ur5_six_axis_ident.py）的動畫用軌跡：不連線、不動手臂
#   每段 = movej 到起點（六軸同步梯形）+ 模擬器的 speedj／stopj 段落（sim_six_axis.trajectory，理想追隨）
#   輸出 CSV：t, q0..q5 [deg], qd0..qd5 [rad/s], label —— 給 make_six_axis_animation.m 用
import os, sys, csv, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.path.insert(0, HERE)
import ur5_six_axis_ident as six
import sim_six_axis as sim

DT = six.DT


def movej(q_from, q_to, a=six.MOVEJ_ACCEL, vmax=six.MOVEJ_VEL, dt=DT):
    """六軸同步梯形 movej（最遠的一軸決定時間）。回傳 q 陣列 [rad]、qd 陣列。"""
    d = np.asarray(q_to) - np.asarray(q_from)
    D = float(np.max(np.abs(d)))
    if D < 1e-9:
        return np.empty((0, 6)), np.empty((0, 6))
    vp = min(vmax, math.sqrt(a * D)); ta = vp / a; tc = max(0.0, (D - vp * ta) / vp); T = 2 * ta + tc
    Q, V = [], []
    for k in range(1, int(math.ceil(T / dt)) + 1):
        t = min(k * dt, T)
        if t < ta:
            x, v = 0.5 * a * t * t, a * t
        elif t < ta + tc:
            x, v = 0.5 * a * ta * ta + vp * (t - ta), vp
        else:
            r = T - t; x, v = D - 0.5 * a * r * r, a * r
        Q.append(np.asarray(q_from) + d * x / D); V.append(d * v / D)
    return np.array(Q), np.array(V)


def phase_constvel(s, q, qd):
    if np.max(np.abs(qd)) < 1e-6:
        return "靜止"
    dj = np.array(s["dj"]); ws, we = np.radians(s["win_start_deg"]), np.radians(s["win_end_deg"])
    if np.all(dj * (q - ws) >= 0) and np.all(dj * (q - we) <= 0):
        return "定速【窗口內，分析用】"
    return "加速、穩定" if np.all(dj * (q - ws) < 0) else "越過窗口、stopj 煞車"


def phase_velmod(s, q, qd):
    if np.max(np.abs(qd)) < 1e-6:
        return "靜止"
    d0 = s["dj"][0]; ms, me = math.radians(s["mod_start_deg"][0]), math.radians(s["mod_end_deg"][0])
    if d0 * (q[0] - ms) < 0:
        return "加速到巡航、定速"
    return "速度調變【分析用】" if d0 * (q[0] - me) <= 0 else "stopj 煞車"


def build(mode, stage="FULL"):
    q0 = list(six.POSES["SIX"]["q_deg"])
    if mode == "CONSTVEL":
        segs = six.plan_constvel(six.stage_config(mode, stage), q0)
    else:
        ci, nr = six.stage_config(mode, stage)
        segs = six.plan_velmod(ci, nr, q0)
    t0, q_cur = 0.0, np.radians(q0)
    rows = [(0.0, q_cur.copy(), np.zeros(6), "SIX 姿態")]
    for k, s in enumerate(segs):
        way = "去程" if s["d"] > 0 else "回程"
        if mode == "CONSTVEL":
            tag = f"第 {k + 1}/{len(segs)} 段　v = {s['speed']:.3f} rad/s　{way}"
        else:
            f0, amp = six.TIERS[s["cond"]]
            tag = (f"第 {k + 1}/{len(segs)} 段　條件 {s['cond'] + 1}（基準 {f0:g} Hz、{'大' if amp == 'L' else '小'}振幅）"
                   f"第 {s['repeat'] + 1} 次　{way}")
        Q, V = movej(q_cur, np.radians(s["start_deg"]))
        for qq, vv in zip(Q, V):
            t0 += DT; rows.append((t0, qq, vv, f"{tag}：movej 到起點"))
        T, Qs, QDs, _ = sim.trajectory([s], mode, s["start_deg"])
        ph = phase_constvel if mode == "CONSTVEL" else phase_velmod
        for qq, vv in zip(Qs, QDs):
            t0 += DT; rows.append((t0, qq, vv, f"{tag}：{ph(s, qq, vv)}"))
        q_cur = Qs[-1]
    Q, V = movej(q_cur, np.radians(q0))
    for qq, vv in zip(Q, V):
        t0 += DT; rows.append((t0, qq, vv, "結束：movej 回 SIX"))
    return rows


def save(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["t"] + [f"q{j}" for j in range(6)] + [f"qd{j}" for j in range(6)] + ["label"])
        for t, q, v, lab in rows:
            w.writerow([f"{t:.3f}"] + [f"{x:.4f}" for x in np.degrees(q)] + [f"{x:.4f}" for x in v] + [lab])
    print(path, f"{rows[-1][0]:.1f} s，{len(rows)} 筆")


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    out = sys.argv[1] if len(sys.argv) > 1 else r"D:\cluade_code\ur5\six_axis_anim"
    os.makedirs(out, exist_ok=True)
    for mode in ("CONSTVEL", "VELMOD"):
        save(build(mode), os.path.join(out, f"traj_six_{mode.lower()}_full.csv"))
