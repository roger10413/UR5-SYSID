# 六軸同時鑑別（robot/ur5_six_axis_ident.py）的模擬驗證：不連線、不動手臂
#   用 UR5 官方連桿質量／重心 + 末端工具（剛體 Newton-Euler）+ 轉子 + 已知摩擦，照程式規劃的段落產生 125 Hz 資料，
#   寫成和上機相同欄位的 CSV，再用程式內的分析函式算回 B、Tc、J，和真值比。
#   真值：B、Tc 為設定值；J 為 M_jj（含轉子）沿調變路徑的平均。
import os, sys, csv, math, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
import ur5_six_axis_ident as six

G = 9.81
LINK_M = [3.7, 8.393, 2.33, 1.219, 1.219, 0.1879]
LINK_COM = [[0, -0.02561, 0.00193], [0.2125, 0, 0.11336], [0.15, 0, 0.0265],
            [0, -0.0018, 0.01634], [0, 0.0018, 0.01634], [0, 0, -0.001159]]
TOOL_KG, TOOL_COG = 1.30, np.array([0.001, 0.001, 0.074])
ROTOR = np.array([1.915, 1.915, 1.915, 0.20, 0.20, 0.20])        # J3～J5 未知，取控制器 J_model 的量級
B_TRUE = np.array([20.0, 28.0, 30.0, 3.8, 2.8, 4.0])             # 量級同 11/12、11/19 結果
TC_TRUE = np.array([9.0, 7.3, 7.4, 1.8, 1.55, 1.8])
K_TRUE = np.array(six.K_AXIS)


def rnea(q, qd, qdd, pay=TOOL_KG, gravity=True):
    fr = six.fk_frames(q)
    z = [f[:3, 2] for f in fr]; p = [f[:3, 3] for f in fr]
    w, dw = np.zeros(3), np.zeros(3)
    a = np.array([0.0, 0.0, G if gravity else 0.0])
    Fl, Ml = [], []
    for i in range(6):
        w_new = w + z[i] * qd[i]
        dw = dw + z[i] * qdd[i] + np.cross(w, z[i] * qd[i])
        w = w_new
        R = fr[i + 1][:3, :3]
        bodies = [(LINK_M[i], np.array(LINK_COM[i]))] + ([(pay, TOOL_COG)] if i == 5 else [])
        F, M = np.zeros(3), np.zeros(3)
        for m, com in bodies:
            r = R @ com + p[i + 1] - p[i]
            f = m * (a + np.cross(dw, r) + np.cross(w, np.cross(w, r)))
            F += f; M += np.cross(r, f)
        Fl.append(F); Ml.append(M)
        d = p[i + 1] - p[i]
        a = a + np.cross(dw, d) + np.cross(w, np.cross(w, d))
    tau, f, n = np.zeros(6), np.zeros(3), np.zeros(3)
    for i in range(5, -1, -1):
        n = Ml[i] + n + np.cross(p[i + 1] - p[i], f)
        f = Fl[i] + f
        tau[i] = n @ z[i] + ROTOR[i] * qdd[i]
    return tau


def mass_diag(q):
    z6 = np.zeros(6)
    return np.array([rnea(q, z6, e, gravity=False)[j] for j, e in enumerate(np.eye(6))])


def trajectory(segs, mode, q0_deg, dt=six.DT):
    """照 URScript 的邏輯產生目標速度（movej 以靜止跳到起點代替）。回傳 t, q, qd, qdd（rad）。"""
    T, Q, QD = [], [], []
    t, q = 0.0, np.radians(q0_deg)

    def hold(n):
        nonlocal t
        for _ in range(n):
            T.append(t); Q.append(q.copy()); QD.append(np.zeros(6)); t += dt

    for s in segs:
        q = np.radians(s["start_deg"]); hold(int(six.PAUSE_TIME / dt))
        dj = np.array(s["dj"], float)
        cur = np.zeros(6)

        def step(cmd, a_lim):
            nonlocal q, t, cur
            dv = cmd - cur
            m = np.max(np.abs(dv))
            if m > a_lim * dt:                                    # 加速度上限以最快一軸為準、其他軸等比例
                dv = dv * (a_lim * dt / m)
            q = q + (cur + 0.5 * dv) * dt                         # 梯形積分：角度和速度同一時刻（前向積分會讓角度超前半拍 4 ms，J 偏低 B·4 ms）
            cur = cur + dv
            T.append(t); Q.append(q.copy()); QD.append(cur.copy()); t += dt

        if mode == "CONSTVEL":
            stop = np.radians(s["stop_deg"])
            for _ in range(s["n_max"]):
                if np.any(dj * q >= dj * stop):
                    break
                step(np.array(s["qd"]), six.SPEEDJ_ACCEL)
            vmax = s["speed"]
        else:
            guard = np.radians(s["guard_deg"])
            vc = dj * np.array(six.V_CRUISE)
            for _ in range(s["n_pre"]):
                if np.any(dj * q >= dj * guard):
                    break
                step(vc, six.RAMP_ACCEL)
            for n in range(s["n_mod"]):
                if np.any(dj * q >= dj * guard):
                    break
                cmd = np.array([dj[j] * (six.V_CRUISE[j] + s["fa"][j][1] * math.sin(2 * math.pi * s["fa"][j][0] * dt * n))
                                for j in range(6)])
                step(cmd, six.MOD_ACCEL_LIMIT)
            vmax = max(abs(cur))
        n_stop = int(math.ceil(vmax / six.STOPJ_DECEL / dt)) + 1
        for _ in range(n_stop):
            step(np.zeros(6), six.STOPJ_DECEL)
        hold(int(six.PAUSE_TIME / dt))
    T, Q, QD = np.array(T), np.array(Q), np.array(QD)
    QDD = np.gradient(QD, dt, axis=0)
    return T, Q, QD, QDD


def write_csv(path, T, Q, QD, QDD, noise_a=0.01, seed=0):
    rng = np.random.default_rng(seed)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        wr = None
        for k in range(len(T)):
            tau_rb = rnea(Q[k], QD[k], QDD[k])
            tau = tau_rb + B_TRUE * QD[k] + TC_TRUE * np.sign(np.round(QD[k], 9))
            cur = tau / K_TRUE + noise_a * rng.standard_normal(6)
            row = {"timestamp": round(T[k], 6)}
            for j in range(6):
                row[f"actual_q_{j}"] = Q[k, j]; row[f"actual_qd_{j}"] = QD[k, j]
                row[f"actual_current_{j}"] = cur[j]; row[f"target_qd_{j}"] = QD[k, j]
                row[f"target_qdd_{j}"] = QDD[k, j]; row[f"target_moment_{j}"] = tau_rb[j]
                row[f"joint_temp_{j}"] = 30.0
            if wr is None:
                wr = csv.DictWriter(fh, fieldnames=list(row)); wr.writeheader()
            wr.writerow(row)


def run(mode, stage, out_dir, noise_a=0.01):
    q0 = list(six.POSES["SIX"]["q_deg"])
    if mode == "CONSTVEL":
        segs = six.plan_constvel(six.stage_config(mode, stage), q0)
    else:
        ci, nr = six.stage_config(mode, stage)
        segs = six.plan_velmod(ci, nr, q0)
    T, Q, QD, QDD = trajectory(segs, mode, q0)
    path = os.path.join(out_dir, f"sim_six_{mode.lower()}_{stage.lower()}.csv")
    write_csv(path, T, Q, QD, QDD, noise_a)
    if mode == "CONSTVEL":
        return path, six.analyze_constvel(path), None
    moving = np.all(np.abs(QD) > 0.02, axis=1)
    Jtrue = np.mean([mass_diag(q) for q in Q[moving][::25]], axis=0)
    return path, six.analyze_velmod(path, payload_kg=TOOL_KG), Jtrue


if __name__ == "__main__":
    sys.stdout.reconfigure(encoding="utf-8")
    out = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.environ.get("TEMP", "."), "six_axis_sim")
    os.makedirs(out, exist_ok=True)
    rep = {}
    p, rc, _ = run("CONSTVEL", "FULL", out)
    print("===== 定速（FULL）：算回 vs 真值")
    print(f"{'軸':>3}{'B 算回':>10}{'B 真值':>9}{'誤差%':>8}{'Tc 算回':>10}{'Tc 真值':>9}{'誤差%':>8}")
    for j, r in rc.items():
        eb, et = 100 * (r["B"] / B_TRUE[j] - 1), 100 * (r["Tc"] / TC_TRUE[j] - 1)
        print(f"J{j} {r['B']:9.3f}{B_TRUE[j]:9.2f}{eb:8.2f}{r['Tc']:10.3f}{TC_TRUE[j]:9.2f}{et:8.2f}")
        rep[f"J{j}_cv"] = dict(B=r["B"], Tc=r["Tc"], eB=eb, eTc=et)
    p, rv, Jt = run("VELMOD", "FULL", out)
    print("\n===== 速度調變（FULL）：各條件 J（正反平均）vs 真值（M_jj 含轉子，沿路徑平均）")
    for j, r in rv.items():
        vals = []
        for ci in sorted(r["per"]):
            pp = r["per"][ci]
            Jm = np.mean([pp[t]["J"] for t in pp])
            vals.append((r["conds"][ci][0], Jm))
        errs = [100 * (v / Jt[j] - 1) for _, v in vals]
        print(f"J{j} 真值 {Jt[j]:.4f}：" + "  ".join(f"{f:.2f}Hz {v:.4f}" for f, v in vals) +
              f"   誤差 {min(errs):+.1f}%～{max(errs):+.1f}%")
        rep[f"J{j}_vm"] = dict(Jtrue=float(Jt[j]), J=[v for _, v in vals], err=errs)
    json.dump(rep, open(os.path.join(out, "sim_six_axis.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    print("存檔:", out)
