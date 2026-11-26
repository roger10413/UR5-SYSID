"""11/26 J1 共振驗證（RES 條件組，手臂朝上，1.30 與 2.31 kg）。

1. 每個條件的 J（正反轉加權平均）與砝碼的 ΔJ = J(2.31) − J(1.30)。
2. 和共振預測比：若 ΔJ 隨頻率變大是 4.33 Hz 的共振造成，
   ΔJ(f) ≈ ΔJ_model / (1 − (f / f_n)²)，f_n = 4.33 Hz（敲擊測試，手臂朝上）。
   也用資料本身擬合 f_n（ΔJ_0 與 f_n 兩個參數）。
3. 工具端加速度計：扣重力後的 r_app = |acc/qdd| − g/ω²，剛體時 = 0.916 m 且不隨頻率變。
"""
import os, sys, glob, json, math, io, contextlib
import numpy as np, pandas as pd
from scipy.optimize import curve_fit

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "robot"))
import ur5_j1_velocity_modulation_ident as v1  # noqa: E402

D = r"D:\cluade_code\ur5\data_1126\J1"
OUT = r"D:\cluade_code\ur5\analysis_1126"; os.makedirs(OUT, exist_ok=True)
K1 = 12.03
S = K1 / v1.KT_OUT
CONDS = v1.CONDITION_SETS["RES"]["conditions"]
T_MOD = v1.CONDITION_SETS["RES"]["t_mod"]
FN_TAP = 4.33
DT = 0.008


def one(p):
    f = glob.glob(p); assert len(f) == 1, (p, f); return f[0]


def wmean(vals):
    w = np.array([1 / s ** 2 for _, s in vals])
    return float(sum(J * wi for (J, _), wi in zip(vals, w)) / w.sum()), float(1 / math.sqrt(w.sum()))


def analyze():
    res = {}
    for pay in ("p130", "p231"):
        for d, tag in ((+1, "pos"), (-1, "neg")):
            path = one(os.path.join(D, f"j1velmod_up_{pay}_res_{tag}_full_*", "*.csv"))
            with contextlib.redirect_stdout(io.StringIO()):
                r = v1.analyze_j1_velocity_modulation(path, 1, v1.KT_OUT, d, CONDS, v1.V_CRUISE)
            res[f"{pay}_{tag}"] = dict(
                per={str(ci): dict(f=c["freq"], A=c["amp"], J=c["J"] * S, se=c["se_J"] * S, B=c.get("B_local", 0) * S,
                                   gain=c.get("track_gain"), lag=c.get("track_lag_ms"), n=c.get("n_seg"), temp=c.get("temp"))
                     for ci, c in r["conditions"].items()},
                J_model=r.get("J_model_mean"), path=path)
    return res


def toolacc(path):
    """每個條件的 r_app（扣重力）、相位。區塊依序 = 條件 0..6 輪流。"""
    Dd = pd.read_csv(path, usecols=["target_qd_1", "actual_qd_1", "tool_acc_x", "tool_acc_y", "tool_acc_z"])
    qt, qa = Dd.target_qd_1.values, Dd.actual_qd_1.values
    acc = Dd[["tool_acc_x", "tool_acc_y", "tool_acc_z"]].values
    qdd = np.gradient(qa, DT)
    mv = np.abs(qt) > 1e-6
    e = np.flatnonzero(np.diff(mv.astype(int)))
    # 只留調變段（約 5 s）；回起點的 movej（約 2.9 s）不算
    blocks = [(a, b) for a, b in zip(e[::2] + 1, e[1::2] + 1) if (b - a) * DT > T_MOD]
    out = {}
    for k, (a, b) in enumerate(blocks):
        ci = k % len(CONDS)
        f0 = CONDS[ci][0]
        r0 = a + np.flatnonzero(np.abs(qt[a:b]) >= v1.V_CRUISE - 1e-3)[0]   # 第一次到達巡航速度
        i0 = r0 + int(round(0.8 / DT)); i1 = i0 + int(round(T_MOD / DT))
        t = np.arange(i0, i1) * DT
        X = np.column_stack([np.sin(2 * np.pi * f0 * t), np.cos(2 * np.pi * f0 * t), np.ones_like(t), t])
        cq = np.linalg.lstsq(X, qdd[i0:i1], rcond=None)[0][:2]
        ca = np.array([np.linalg.lstsq(X, acc[i0:i1, kk], rcond=None)[0][:2] for kk in range(3)])
        out.setdefault(ci, []).append((ca[:, 0] + 1j * ca[:, 1]) / (cq[0] + 1j * cq[1]))
    rows = {}
    for ci, H in out.items():
        H = np.array(H); Hm = H.mean(axis=0); k = int(np.argmax(np.abs(Hm)))
        mag = np.abs(H[:, k]); w2 = (2 * np.pi * CONDS[ci][0]) ** 2
        rows[ci] = dict(axis="xyz"[k], r_app=float(mag.mean() - 9.81 / w2), se=float(mag.std(ddof=1) / np.sqrt(len(mag))),
                        phase=float(np.degrees(np.angle(Hm[k]))), n=len(H))
    return rows


if __name__ == "__main__":
    res = analyze()
    order = sorted(range(len(CONDS)), key=lambda c: CONDS[c][0])
    jm = {p: res[f"{p}_pos"]["J_model"] for p in ("p130", "p231")}
    print(f"控制器 J_model：1.30 kg {jm['p130']:.3f}、2.31 kg {jm['p231']:.3f}、ΔJ_model {jm['p231'] - jm['p130']:.3f}")
    print(f"{'f':>5}{'A':>7}{'J 1.30':>13}{'J 2.31':>13}{'ΔJ':>14}{'預測(4.33)':>11}{'增益1.30':>9}{'溫度':>7}")
    tab = []
    dJm = jm["p231"] - jm["p130"]
    for ci in order:
        f, A = CONDS[ci]
        Jp = {p: wmean([(res[f"{p}_{t}"]["per"][str(ci)]["J"], res[f"{p}_{t}"]["per"][str(ci)]["se"]) for t in ("pos", "neg")])
              for p in ("p130", "p231")}
        dvals = []
        for t in ("pos", "neg"):
            a, b = res[f"p130_{t}"]["per"][str(ci)], res[f"p231_{t}"]["per"][str(ci)]
            dvals.append((b["J"] - a["J"], math.hypot(a["se"], b["se"])))
        dJ = wmean(dvals)
        pred = dJm / (1 - (f / FN_TAP) ** 2)
        g = np.mean([res[f"p130_{t}"]["per"][str(ci)]["gain"] for t in ("pos", "neg")])
        tp = res["p130_pos"]["per"][str(ci)]["temp"]
        tab.append((f, A, Jp, dJ, pred))
        print(f"{f:5.2f}{A:7.4f}{Jp['p130'][0]:8.3f}±{Jp['p130'][1]:.2f}{Jp['p231'][0]:8.3f}±{Jp['p231'][1]:.2f}"
              f"{dJ[0]:8.3f}±{dJ[1]:.3f}{pred:11.3f}{g:9.3f}{(tp if tp is not None else float('nan')):7.1f}")

    # 擬合 ΔJ(f) = ΔJ0 / (1 − (f/fn)²)
    fs = np.array([r[0] for r in tab]); dj = np.array([r[3][0] for r in tab]); se = np.array([r[3][1] for r in tab])
    model = lambda f, d0, fn: d0 / (1 - (f / fn) ** 2)
    p, cov = curve_fit(model, fs, dj, p0=(0.86, 4.3), sigma=se, absolute_sigma=True, bounds=([0, 3.6], [5, 20]))
    chi_fit = float(np.sum(((dj - model(fs, *p)) / se) ** 2))
    chi_tap = float(np.sum(((dj - model(fs, dJm, FN_TAP)) / se) ** 2))
    chi_const = float(np.sum(((dj - np.average(dj, weights=1 / se ** 2)) / se) ** 2))
    print(f"\n擬合 ΔJ0 = {p[0]:.3f} ± {math.sqrt(cov[0, 0]):.3f}，f_n = {p[1]:.2f} ± {math.sqrt(cov[1, 1]):.2f} Hz，χ² = {chi_fit:.1f}（{len(fs) - 2} 自由度）")
    print(f"固定 ΔJ0 = ΔJ_model、f_n = 4.33 的 χ² = {chi_tap:.1f}（{len(fs)} 自由度）；ΔJ 不隨頻率變的 χ² = {chi_const:.1f}（{len(fs) - 1} 自由度）")

    # J 本身（1.30 kg）也用同一形式擬合
    for pay in ("p130", "p231"):
        Jv = np.array([r[2][pay][0] for r in tab]); Js = np.array([r[2][pay][1] for r in tab])
        pj, cj = curve_fit(model, fs, Jv, p0=(6, 4.3), sigma=Js, absolute_sigma=True, bounds=([0, 3.6], [20, 20]))
        print(f"J {pay}：J0 = {pj[0]:.3f} ± {math.sqrt(cj[0, 0]):.3f}，f_n = {pj[1]:.2f} ± {math.sqrt(cj[1, 1]):.2f} Hz"
              f"（J_model {jm[pay]:.3f}）")

    print("\n工具加速度計 r_app [m]（剛體 0.916，不隨頻率變）")
    ta = {}
    for key in ("p130_pos", "p130_neg", "p231_pos", "p231_neg"):
        ta[key] = toolacc(res[key]["path"])
    print(f"{'f':>5}{'A':>7}" + "".join(f"{k:>16}" for k in ta))
    for ci in order:
        f, A = CONDS[ci]
        print(f"{f:5.2f}{A:7.4f}" + "".join(f"{ta[k][ci]['r_app']:9.3f}±{ta[k][ci]['se']:.3f} " for k in ta))
    with open(os.path.join(OUT, "j1_res_1126.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(res=res, toolacc={k: {str(c): v for c, v in d.items()} for k, d in ta.items()},
                       fit=dict(dJ0=p[0], fn=p[1], chi_fit=chi_fit, chi_tap=chi_tap, chi_const=chi_const)),
                  fh, ensure_ascii=False, indent=1, default=str)
