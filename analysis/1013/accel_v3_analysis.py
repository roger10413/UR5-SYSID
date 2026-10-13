"""10/13 定加速度 v3：含截距回歸 J、只用加速段 J、換階段跳過時間敏感度、殘差對起步後行程。"""
import sys, os, io, contextlib
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "robot"))
sys.stdout.reconfigure(encoding="utf-8")
import ur5_const_accel_ident_v3 as m

BASE = r"D:\cluade_code\ur5"   # 資料資料夾（內含各次 constaccel_v3_* 子資料夾），請改成自己的路徑
RUNS = [("正轉", +1, os.path.join(BASE, "constaccel_v3_pos_full_20261013_101718", "constaccel_v3_data_pos_full_20261013_101718.csv")),
        ("反轉", -1, os.path.join(BASE, "constaccel_v3_neg_full_20261013_102652", "constaccel_v3_data_neg_full_20261013_102652.csv"))]
# 同日定速（第 1 輪，緊接在 v3 之前）各檔 |tau|，升降兩趟平均
CV_V = np.array([0.02, 0.05, 0.08, 0.12, 0.16, 0.20])
CV_T = {+1: np.array([9.388, 10.577, 11.693, 13.026, 14.333, 15.533]),
        -1: np.array([9.623, 10.906, 12.028, 13.456, 14.743, 15.940])}


def fcv(v, s):
    v = np.abs(np.atleast_1d(np.asarray(v, float))); t = CV_T[s]
    out = np.interp(v, CV_V, t); hi = v > 0.2
    out[hi] = t[-1] + (t[-1] - t[-2]) / 0.04 * (v[hi] - 0.2)
    return out


def level_stats(C, B):
    """各檔：h = (tau_acc − tau_dec)/2（以 B 修正速度差）、a、只用加速段 J 需要的量。"""
    out = []
    for li in range(len(m.ACCEL_LEVELS)):
        c = [x for x in C if x["level"] == li]
        TA = np.nanmean([x["tau_acc"] for x in c], 0); TD = np.nanmean([x["tau_dec"] for x in c], 0)
        VA = np.nanmean([x["v_acc"] for x in c], 0); VD = np.nanmean([x["v_dec"] for x in c], 0)
        a = np.mean([x["a_acc"] + x["a_dec"] for x in c]) / 2
        out.append(dict(a=a, h=np.nanmean(TA - TD - B * (VA - VD)) / 2, TA=TA, VA=VA, aa=np.mean([x["a_acc"] for x in c])))
    return out


def fit_intercept(C, B):
    L = level_stats(C, B)
    a = np.array([x["a"] for x in L]); h = np.array([x["h"] for x in L])
    return np.polyfit(a, h, 1), L


def run(name, s, F, settle_list=(0.0, 0.02, 0.04, 0.06, 0.08), n_boot=300):
    res = {}
    for st in settle_list:
        with contextlib.redirect_stdout(io.StringIO()):
            est = m.analyze_const_accel_v3(F, 0, m.KT_OUT, s, m.ACCEL_LEVELS, m.V_PEAK, settle_s=st, n_boot=2)
        p, L = fit_intercept(est["cycles"], est["B"])
        res[st] = dict(est=est, p=p, L=L)
    est = res[0.04]["est"]; C = est["cycles"]; rng = np.random.default_rng(0)
    bs = np.array([fit_intercept([C[i] for i in rng.integers(0, len(C), len(C))], est["B"])[0] for _ in range(n_boot)])
    L = res[0.04]["L"]
    j_acc = [np.nanmean(x["TA"] - fcv(x["VA"], s)) / x["aa"] for x in L]
    return dict(name=name, s=s, res=res, se=bs.std(0), j_levels=est["J_levels"], se_levels=est["se_J_levels"],
                a_levels=est["a_levels"], j_acc=np.array(j_acc), L=L, F=F)


if __name__ == "__main__":
    out = [run(*r) for r in RUNS]
    for o in out:
        p = o["res"][0.04]["p"]
        print(f"{o['name']}：各檔 J = {np.round(o['j_levels'], 3)}（± {np.round(o['se_levels'], 3)}）")
        print(f"  含截距：J = {p[0]:.3f} ± {o['se'][0]:.3f}，截距 = {p[1]:+.3f} ± {o['se'][1]:.3f} N·m → Δ = {-2*p[1]:.2f} N·m")
        print(f"  只用加速段（借同日定速摩擦）：{np.round(o['j_acc'], 3)}")
        print("  換階段跳過時間敏感度：", "  ".join(f"{int(k*1000)}ms: J={v['p'][0]:.3f}, c={v['p'][1]:+.3f}" for k, v in o["res"].items()))
    np.save(os.path.join(BASE, "analysis_1013", "accel_v3.npy"),
            [{k: v for k, v in o.items() if k != "res"} | {"sens": {k: v["p"] for k, v in o["res"].items()}} for o in out],
            allow_pickle=True)
