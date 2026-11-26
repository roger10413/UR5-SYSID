"""12/2 報告的圖（11/26 上機）。
圖 1：六軸定速 4 輪（09:05 冷機、10:10、10:24、10:37）：B、Tc、關節溫度對時間；虛線為 11/19 單軸（J2～J5）。
圖 2：J1 共振驗證：砝碼的 ΔJ 對頻率（含共振擬合），工具端加速度計的等效半徑對頻率。
圖 3：六軸速度調變：J ÷ 控制器 J_model 對頻率；實測電流 vs 控制器模型力矩（target_moment）當電流。
只讀已存的結果（six_constvel_temp_1126.json、j1_res_1126.json）與 11/26 六軸速度調變 CSV，不連網路。
輸出 D:\\cluade_code\\ur5\\analysis_1126\\figs\\w1202_fig*.png。（Microsoft JhengHei 沒有「−」，數值用 ASCII 負號。）"""
import socket
socket.socket = None                                   # 離線分析：不連網路
import os, sys, json, math
from datetime import datetime
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
sys.path.insert(0, os.path.join(ROOT, "robot")); sys.path.insert(0, HERE)
import ur5_six_axis_ident as six  # noqa: E402
from j1_res_1126 import CONDS, wmean  # noqa: E402

plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
OUT = r"D:\cluade_code\ur5\analysis_1126"
FIGS = os.path.join(OUT, "figs"); os.makedirs(FIGS, exist_ok=True)
INK, INK2, MUTED, GRAY = "#0b0b0b", "#52514e", "#8d8b83", "#b9b8b0"
ORANGE, BLUE, GREEN, PURPLE, RED = "#d9822b", "#2a78d6", "#3a9a5b", "#8a5cc2", "#c8463c"
SINGLE_B = {2: 29.85, 3: 3.76, 4: 2.78, 5: 4.02}       # 11/19 單軸定速（診斷姿態）
SINGLE_TC = {2: 7.44, 3: 1.83, 4: 1.55, 5: 1.78}


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", alpha=1, ls=(0, (4, 3))); ax.set_axisbelow(True)


def save(fig, name):
    p = os.path.join(FIGS, name)
    fig.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0.15); plt.close(fig)
    return p


def fig1():
    R = json.load(open(os.path.join(OUT, "six_constvel_temp_1126.json"), encoding="utf-8"))["runs"]
    t = [datetime.strptime("20261126" + r["tag"], "%Y%m%d%H%M%S") for r in R]
    fig, axs = plt.subplots(3, 6, figsize=(18, 8.2), sharex=True)
    rows = (("B", "sd_B", BLUE, "B [N·m·s/rad]", SINGLE_B), ("Tc", "sd_Tc", ORANGE, "Tc [N·m]", SINGLE_TC),
            ("T", None, RED, "關節溫度 [°C]", None))
    for j in range(6):
        for i, (key, sd, col, lab, ref) in enumerate(rows):
            ax = axs[i, j]
            y = [r["per"][str(j)][key] for r in R]
            e = [r["per"][str(j)][sd] for r in R] if sd else None
            ax.errorbar(t, y, yerr=e, fmt="o-", color=col, ms=7, mec="white", capsize=4, lw=1.2, zorder=3)
            ax.plot(t[0], y[0], "o", ms=11, mfc="none", mec=INK2, zorder=4)          # 冷機那輪
            if ref and j in ref:
                ax.axhline(ref[j], color=MUTED, ls="--", lw=1.2)
            if key != "T":
                lo, hi = min(min(y), ref.get(j, min(y)) if ref else min(y)), max(max(y), ref.get(j, max(y)) if ref else max(y))
                pad = max((hi - lo) * 0.3, 0.04 * np.mean(y))
                ax.set_ylim(lo - pad, hi + pad)
            if i == 0:
                ax.set_title(f"J{j}", loc="left", fontsize=13)
            if j == 0:
                ax.set_ylabel(lab)
            style(ax)
            ax.xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
            ax.set_xticks(t); ax.tick_params(axis="x", labelrotation=45, labelsize=9)
    fig.suptitle("11/26 六軸定速 4 輪：B、Tc、關節溫度對時間（09:05 為冷機，圈起來的點；誤差棒為回歸標準差；"
                 "灰色虛線為 11/19 單軸結果，只有 J2～J5）",
                 x=0.01, ha="left", fontsize=14)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    return save(fig, "w1202_fig1_six_friction_time.png")


def fig2():
    A = json.load(open(os.path.join(OUT, "j1_res_1126.json"), encoding="utf-8"))
    res, ta, fit = A["res"], A["toolacc"], A["fit"]
    order = sorted(range(len(CONDS)), key=lambda c: CONDS[c][0])
    fs = np.array([CONDS[c][0] for c in order])
    dj = []
    for ci in order:
        vals = []
        for tg in ("pos", "neg"):
            a, b = res[f"p130_{tg}"]["per"][str(ci)], res[f"p231_{tg}"]["per"][str(ci)]
            vals.append((b["J"] - a["J"], math.hypot(a["se"], b["se"])))
        dj.append(wmean(vals))
    dj = np.array(dj)
    dJm = res["p231_pos"]["J_model"] - res["p130_pos"]["J_model"]
    ff = np.linspace(0.3, 3.8, 200)
    fig, (a1, a2) = plt.subplots(1, 2, figsize=(13, 4.8))
    a1.axhline(dJm, color=MUTED, ls="--", lw=1.2)
    a1.text(2.25, dJm - 0.04, f"控制器 ΔJ_model {dJm:.2f}", color=MUTED, fontsize=10, va="top")
    a1.set_ylim(0.6, 3.8)
    a1.plot(ff, fit["dJ0"] / (1 - (ff / fit["fn"]) ** 2), color=BLUE, lw=1.6,
            label=f"共振擬合：ΔJ0 = {fit['dJ0']:.2f}，f_n = {fit['fn']:.2f} Hz")
    a1.plot(ff, dJm / (1 - (ff / 4.33) ** 2), color=GRAY, lw=1.2, ls=":", label="敲擊測試 4.33 Hz 的預測")
    a1.errorbar(fs, dj[:, 0], yerr=dj[:, 1], fmt="o", color=INK, ms=7, mec="white", capsize=4, zorder=3,
                label="實測 ΔJ = J(2.31 kg) - J(1.30 kg)")
    a1.set_xlabel("調變頻率 [Hz]"); a1.set_ylabel("ΔJ [kg·m²]"); a1.set_title("砝碼的 ΔJ", loc="left")
    a1.legend(frameon=False, fontsize=10, loc="upper left"); style(a1)
    for pay, col, lab in (("p130", BLUE, "1.30 kg"), ("p231", ORANGE, "2.31 kg")):
        y = np.array([np.mean([ta[f"{pay}_{tg}"][str(c)]["r_app"] for tg in ("pos", "neg")]) for c in order])
        e = np.array([math.hypot(ta[f"{pay}_pos"][str(c)]["se"], ta[f"{pay}_neg"][str(c)]["se"]) / 2 for c in order])
        a2.errorbar(fs, y, yerr=e, fmt="o-", color=col, ms=7, mec="white", capsize=4, lw=1.2, label=lab, zorder=3)
    a2.axhline(0.916, color=MUTED, ls="--", lw=1.2); a2.text(0.32, 0.93, "剛體 0.916 m", color=MUTED, fontsize=10)
    a2.set_xlabel("調變頻率 [Hz]"); a2.set_ylabel("等效半徑 [m]")
    a2.set_title("工具端加速度 ÷ 關節角加速度（扣重力）", loc="left")
    a2.legend(frameon=False, fontsize=10, loc="upper left"); style(a2)
    fig.tight_layout()
    return save(fig, "w1202_fig2_j1_resonance.png")


def fig3():
    D = r"D:\cluade_code\ur5\data_1126\SIX\six_velmod_six_p130_full_20261126_092018"
    rm = six.analyze_velmod(os.path.join(D, "six_velmod_six_p130_full_20261126_092018_gcomp.csv"))
    rt = six.analyze_velmod(os.path.join(D, "tm_as_current_gcomp.csv"))
    fig, axs = plt.subplots(1, 6, figsize=(18, 4.2), sharey=True)
    for j, ax in enumerate(axs):
        for r, col, lab, mk in ((rt, GRAY, "控制器模型力矩當電流（驗證方法）", "s"), (rm, BLUE, "實測電流", "o")):
            for A, mfc in ((max(a for _, a in r[j]["conds"]), None), (min(a for _, a in r[j]["conds"]), "white")):
                pts = sorted((r[j]["conds"][ci][0], ci) for ci in r[j]["per"] if r[j]["conds"][ci][1] == A)
                x = [f for f, _ in pts]
                y = [np.mean([r[j]["per"][ci][d]["J"] / r[j]["per"][ci][d]["Jm"] for d in r[j]["per"][ci]]) for _, ci in pts]
                ax.plot(x, y, marker=mk, color=col, mfc=mfc or col, mec=col, ms=7, lw=1.2,
                        label=f"{lab}，{'大' if mfc is None else '小'}振幅" if j == 0 else None)
        ax.axhline(1, color=MUTED, ls="--", lw=1.0); ax.axhline(0, color=GRAY, lw=0.8)
        ax.set_title(f"J{j}", loc="left"); ax.set_xlabel("頻率 [Hz]"); style(ax)
    axs[0].set_ylabel("J ÷ 控制器 J_model")
    fig.legend(loc="lower center", ncol=4, frameon=False, fontsize=10, bbox_to_anchor=(0.5, -0.06))
    fig.suptitle("11/26 六軸速度調變：J ÷ 控制器 J_model（1 = 和模型相同）", x=0.01, ha="left", fontsize=14)
    fig.tight_layout(rect=(0, 0.04, 1, 0.95))
    return save(fig, "w1202_fig3_six_velmod_J.png")


if __name__ == "__main__":
    for fn in (fig1, fig2, fig3):
        print(fn())
