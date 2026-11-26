"""11/26 六軸定速 4 輪（09:05 冷機、10:10、10:24、10:37）：B、Tc 對溫度。

每輪用 ur5_six_axis_ident.analyze_constvel（正反合併）得到每軸 B、Tc；溫度取窗口內樣本的平均。
另外每輪拆成升速趟、降速趟（six_constvel_drift_1126.py 的方法），每軸得到 8 個點，看 B 對溫度的斜率。
輸出 D:\\cluade_code\\ur5\\analysis_1126\\six_constvel_temp_1126.json 與圖 figs/w1203_six_B_T.png。
"""
import os, sys, glob, json
import numpy as np, pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "robot")); sys.path.insert(0, HERE)
import ur5_six_axis_ident as six  # noqa: E402
from six_constvel_drift_1126 import segment_means, fit  # noqa: E402

D = r"D:\cluade_code\ur5\data_1126\SIX"
OUT = r"D:\cluade_code\ur5\analysis_1126"; FIGS = os.path.join(OUT, "figs"); os.makedirs(FIGS, exist_ok=True)
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 11})
BLUE, ORANGE, INK2 = "#2a78d6", "#d9822b", "#52514e"


def runs():
    return sorted(glob.glob(os.path.join(D, "six_constvel_six_p130_full_*", "six_constvel_*[0-9].csv")))


if __name__ == "__main__":
    win = six.window_bounds(six.POSES["SIX"]["q_deg"])
    res = []
    for path in runs():
        tag = os.path.basename(path)[-10:-4]
        r = six.analyze_constvel(path)
        Dd = pd.read_csv(path)
        per = {}
        for j in range(6):
            df = segment_means(Dd, j, win[j], six.K_AXIS[j])
            df["pass"] = df.groupby(df.v.round(5)).cumcount()
            halves = []
            for p in (0, 1):
                h = df[df["pass"] == p]
                c = fit(h)
                halves.append(dict(B=float(c[1]), Tc=float(c[0]), T=float(h.temp.mean())))
            per[j] = dict(B=r[j]["B"], sd_B=r[j]["sd_B"], Tc=r[j]["Tc"], sd_Tc=r[j]["sd_Tc"], T=float(df.temp.mean()),
                          halves=halves)
        res.append(dict(tag=tag, per=per))
        print(f"== {tag[:2]}:{tag[2:4]}  " + "  ".join(f"J{j} B {per[j]['B']:.2f} Tc {per[j]['Tc']:.3f} T {per[j]['T']:.1f}"
                                                         for j in range(6)))

    print(f"\n{'軸':>3}{'B 冷機':>9}{'B 後三輪':>10}{'變化':>7}{'後三輪 CV':>10}{'Tc 冷機':>9}{'Tc 後三輪':>10}{'T 冷機→後三輪 °C':>20}"
          f"{'dB/dT [%/°C]':>14}")
    summary = {}
    for j in range(6):
        B = np.array([x["per"][j]["B"] for x in res]); Tc = np.array([x["per"][j]["Tc"] for x in res])
        T = np.array([x["per"][j]["T"] for x in res])
        hb = np.array([h["B"] for x in res for h in x["per"][j]["halves"]])
        ht = np.array([h["T"] for x in res for h in x["per"][j]["halves"]])
        slope, icpt = np.polyfit(ht, hb, 1)
        Tref = ht.mean()
        pct = slope / (slope * Tref + icpt) * 100
        cv = B[1:].std(ddof=1) / B[1:].mean() * 100
        summary[j] = dict(B=B.tolist(), Tc=Tc.tolist(), T=T.tolist(), slope=float(slope), pct_per_C=float(pct), cv_last3=float(cv))
        print(f"J{j} {B[0]:9.2f}{B[1:].mean():10.2f}{(B[1:].mean() / B[0] - 1) * 100:+6.0f}%{cv:9.1f}%{Tc[0]:9.3f}{Tc[1:].mean():10.3f}"
              f"{T[0]:12.1f} → {T[1:].mean():.1f}{pct:+12.2f}")

    # 圖：每軸 B 對溫度（左軸）與 Tc 對溫度（右圖）
    fig, axs = plt.subplots(2, 6, figsize=(17, 6.4))
    for j in range(6):
        for row, key, sd, col, lab in ((0, "B", "sd_B", BLUE, "B [N·m·s/rad]"), (1, "Tc", "sd_Tc", ORANGE, "Tc [N·m]")):
            ax = axs[row, j]
            T = [x["per"][j]["T"] for x in res]; y = [x["per"][j][key] for x in res]; e = [x["per"][j][sd] for x in res]
            ax.errorbar(T, y, yerr=e, fmt="o", color=col, ms=7, mec="white", capsize=4, zorder=3)
            for k, x in enumerate(res):
                ax.annotate(f"{x['tag'][:2]}:{x['tag'][2:4]}", (T[k], y[k]), textcoords="offset points", xytext=(6, 4),
                            fontsize=8.5, color=INK2)
            ax.set_title(f"J{j}" if row == 0 else "", loc="left")
            if j == 0:
                ax.set_ylabel(lab)
            if row == 1:
                ax.set_xlabel("關節溫度 [°C]")
            lo, hi = min(y) - max(e), max(y) + max(e)
            pad = max((hi - lo) * 0.25, 0.05 * np.mean(y))
            ax.set_ylim(lo - pad, hi + pad)
            ax.grid(alpha=0.3); [ax.spines[s].set_visible(False) for s in ("top", "right")]
    fig.suptitle("11/26 六軸定速 4 輪：B、Tc 對關節溫度（誤差棒為回歸標準差）", x=0.01, ha="left", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fp = os.path.join(FIGS, "w1203_six_B_T.png"); fig.savefig(fp, dpi=150); plt.close(fig)
    with open(os.path.join(OUT, "six_constvel_temp_1126.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(runs=res, summary=summary), fh, ensure_ascii=False, indent=1)
    print("圖：", fp)
