# 取放情境鑑別（robot/ur5_pick_place_ident.py）的模擬：不連線、不動手臂
#   照 plan_pick_place 的動作產生 125 Hz 資料：movej 為各軸同步的梯形速度（關節空間直線），
#   movel 為工具垂直線上的梯形速度（沿線 IK 後內插），dwell、停頓為靜止。
#   力矩 = 剛體 Newton-Euler（官方連桿 + 負載 + 轉子）+ 已知摩擦（sim_six_axis.B_TRUE、TC_TRUE），寫成上機相同欄位的 CSV。
import os, sys, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot")); sys.path.insert(0, HERE)
import ur5_pick_place_ident as pp
import sim_six_axis as sim


def _trap(dist, v, a, dt):
    """梯形速度下的位置比例 s(t) ∈ [0, 1] 與速度比例（÷ dist）。"""
    T, t_acc, vp = pp.movej_profile(dist, v, a)
    n = int(math.ceil(T / dt))
    t = np.arange(1, n + 1) * dt
    t = np.minimum(t, T)
    pos = np.where(t < t_acc, 0.5 * a * t ** 2,
                   np.where(t < T - t_acc, 0.5 * a * t_acc ** 2 + vp * (t - t_acc), dist - 0.5 * a * (T - t) ** 2))
    vel = np.where(t < t_acc, a * t, np.where(t < T - t_acc, vp, a * (T - t)))
    vel[-1] = 0.0
    return pos / dist, vel / dist


def trajectory(moves, poses, dt=pp.DT, pause=pp.PAUSE_TIME):
    T, Q, QD = [], [], []
    t, q = 0.0, np.radians(poses["A_up"])

    def hold(dur):
        nonlocal t
        for _ in range(int(round(dur / dt))):
            T.append(t); Q.append(q.copy()); QD.append(np.zeros(6)); t += dt

    hold(pause)
    for m in moves:
        if m["kind"] == "dwell":
            hold(m["t"]); continue
        a, b = np.radians(m["start_deg"]), np.radians(m["target_deg"])
        if m["kind"] == "movej":
            lead = np.max(np.abs(b - a))
            s, sd = _trap(lead, m["v"], m["a"], dt)
            for k in range(len(s)):
                q = a + s[k] * (b - a)
                T.append(t); Q.append(q.copy()); QD.append(sd[k] * (b - a)); t += dt
        else:
            dz = abs(pp.flange_pos(m["start_deg"])[2] - pp.flange_pos(m["target_deg"])[2])
            path = np.radians(np.array(pp.path_samples(m, 41)))           # 沿工具直線的關節角（41 點）
            grid = np.linspace(0.0, 1.0, 41)
            dpath = np.gradient(path, grid, axis=0)
            s, sd = _trap(dz, m["v"], m["a"], dt)
            for k in range(len(s)):
                q = np.array([np.interp(s[k], grid, path[:, j]) for j in range(6)])
                qd = np.array([np.interp(s[k], grid, dpath[:, j]) for j in range(6)]) * sd[k]
                T.append(t); Q.append(q.copy()); QD.append(qd); t += dt
        hold(pause)
    T, Q, QD = np.array(T), np.array(Q), np.array(QD)
    QDD = np.gradient(QD, dt, axis=0)
    return T, Q, QD, QDD


def run(speeds, out_path, noise_a=0.0):
    moves, poses = pp.plan_pick_place(speeds)
    T, Q, QD, QDD = trajectory(moves, poses)
    sim.write_csv(out_path, T, Q, QD, QDD, noise_a=noise_a)
    return moves, poses, out_path


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    out = os.path.join(os.environ.get("TEMP", "."), "pick_place_sim"); os.makedirs(out, exist_ok=True)
    speeds = pp.stage_speeds("FULL")
    _, _, p = run(speeds, os.path.join(out, "sim_pick_place_full.csv"), noise_a=0.01)
    W = pp.transfer_windows()
    rt = pp.analyze_transfer(p, W)
    rf = pp.analyze_full_cycle(p, sim.TOOL_KG, step=2)
    print("\n".join(pp.summary_lines(rt, rf, W)))
    print("\n真值：B " + " ".join(f"{x:.2f}" for x in sim.B_TRUE) + "，Tc " + " ".join(f"{x:.2f}" for x in sim.TC_TRUE)
          + "，轉子 " + " ".join(f"{x:.3f}" for x in sim.ROTOR))
