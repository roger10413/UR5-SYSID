"""10/15 定加速度 v3：三角形（定速段 0）× 2 次與梯形（定速段 4 馬達圈）對照。
含截距回歸 J 與 Δ（減速段比加速段多的摩擦）、換階段跳過時間敏感度、殘差對起步後行程。
摩擦基準：同日定速 3 輪各速度檔 |τ| 對溫度做直線回歸，取該次 v3 的平均溫度（需先執行 constvel_1015.py）。"""
import sys, os, io, contextlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot")); sys.path.insert(0, os.path.join(HERE, "..", "1013")); sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
import ur5_const_accel_ident_v3 as m
import accel_v3_analysis as av            # 10/13 的分析函式（run、fit_intercept、fcv）
import data_1015 as dd

CV = list(np.load(os.path.join(dd.BASE, "analysis_1015", "constvel_1015.npy"), allow_pickle=True))


def cv_levels_at(T, s):
    R = [x for x in CV if x["dir"] == s]; Ts = np.array([x["Tm"] for x in R]); L = np.array([x["levels"] for x in R])
    return np.array([np.polyval(np.polyfit(Ts, L[:, i], 1), T) for i in range(6)])


def residual_vs_travel(F, s, J, edges=np.arange(0, 27, 1.5)):
    """|τ| − 同溫度定速摩擦 − J·a，依加速、定速、減速段對起步後行程分箱（同 10/14 報告圖 5 的定義）。"""
    d = np.genfromtxt(F, delimiter=",", names=True, usecols=(0, 1, 2, 3, 4))
    t = d["timestamp"]; q = np.degrees(d["actual_q_0"]); qd = s * d["actual_qd_0"]; tq = s * d["target_qd_0"]
    tau = s * d["actual_current_0"] * m.KT_OUT
    dtq = np.gradient(tq, t); mvg = tq > 0.05
    ph = {"加速段": mvg & (dtq > 0.1), "定速段": mvg & (np.abs(tq - m.V_PEAK) < 1e-4) & (np.abs(dtq) < 0.05), "減速段": mvg & (dtq < -0.1)}
    for k in ph:
        msk = ph[k].copy()
        for i in np.where(np.diff(msk.astype(int)) == 1)[0] + 1:
            msk[i:i + 5] = False
        ph[k] = msk
    res = tau - av.fcv(qd, s) - J * np.where(ph["定速段"], 0, dtq)
    d1 = np.diff(tq, prepend=tq[0])
    starts = [r[0] - 1 for r in m._runs((d1 > 1e-6) & (tq > 1e-3)) if tq[r].max() > 0.2 and (r[0] < 2 or tq[r[0] - 2] < 1e-3)]
    trav = np.full(len(q), np.nan)
    for i0, i1 in zip(starts, starts[1:] + [len(q)]):
        trav[i0:i1] = s * (q[i0:i1] - q[i0])
    out = {}
    for k, msk in ph.items():
        xs, ys = [], []
        for lo in edges[:-1]:
            sm = msk & (trav >= lo) & (trav < lo + (edges[1] - edges[0]))
            if sm.sum() > 40:
                xs.append(lo + (edges[1] - edges[0]) / 2); ys.append(np.mean(res[sm]))
        out[k] = (np.array(xs), np.array(ys))
    return out


def analyze(runs):
    out = []
    for r in runs:
        if r["kind"] not in ("tri", "trap"):
            continue
        T = np.genfromtxt(r["F"], delimiter=",", names=True, usecols=(5,))["joint_temp_0"]
        av.CV_T[r["dir"]] = cv_levels_at(T.mean(), r["dir"])
        with contextlib.redirect_stdout(io.StringIO()):
            o = av.run(r["label"], r["dir"], r["F"], n_boot=300)
        p = o["res"][0.04]["p"]
        o.update(r=r, Tm=T.mean(), T0=T[0], T1=T[-1], J=p[0], c=p[1], seJ=o["se"][0], sec=o["se"][1],
                 resid=np.array([x["h"] for x in o["L"]]) - np.polyval(p, [x["a"] for x in o["L"]]),
                 sens={k: v["p"] for k, v in o["res"].items()}, travel=residual_vs_travel(r["F"], r["dir"], p[0]))
        del o["res"]
        out.append(o)
    return out


if __name__ == "__main__":
    import warnings; warnings.filterwarnings("ignore")
    out = analyze(dd.run_list())
    for o in out:
        print(f"\n[{o['r']['label']}] {o['r']['ts'][:2]}:{o['r']['ts'][2:4]}，J0 溫度 {o['T0']:.2f}→{o['T1']:.2f}")
        print(f"  含截距：J = {o['J']:.3f} ± {o['seJ']:.3f}，截距 = {o['c']:+.3f} ± {o['sec']:.3f} → Δ = {-2 * o['c']:.2f} N·m；四點殘差 {np.round(o['resid'], 3)}")
        print(f"  各檔 J（標準算法）= {np.round(o['j_levels'], 3)}，a = {np.round(o['a_levels'], 3)}")
        print("  換階段跳過時間：", "  ".join(f"{int(k * 1000)}ms J={v[0]:.3f} Δ={-2 * v[1]:.2f}" for k, v in o["sens"].items()))
    np.save(os.path.join(dd.BASE, "analysis_1015", "accel_v3_1015.npy"), out, allow_pickle=True)
