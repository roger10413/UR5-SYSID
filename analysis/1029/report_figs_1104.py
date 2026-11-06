"""11/4 報告第 3 節的圖（改版）：砝碼造成的 J 增加量；校正後的 J 和官方模型、控制器模型比較；比例 c 和頻率。
需先執行 known_mass_1029.py（讀取 OUT/known_mass_1029.json）。只用反轉（EXT45 只做了反轉）。"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
import known_mass_1029 as km

R = json.load(open(os.path.join(km.OUT, "known_mass_1029.json"), encoding="utf-8"))
OUT = os.path.join(km.OUT, "figs"); os.makedirs(OUT, exist_ok=True)
BLUE, GRAY, DARK = "#2a78d6", "#b8b6ad", "#3a3936"
INK, INK2 = "#0b0b0b", "#52514e"

pr = {(p["pose"], p["dir"]): p for p in R["pairs"]}
rev = [pr[("DIAG", -1)], pr[("EXT45", -1)]]
w = np.array([1 / p["sratio_w"] ** 2 for p in rev])
RATIO = float(np.sum(w * [p["ratio_w"] for p in rev]) / np.sum(w))     # 反轉兩組合併的倍數（1.26）
print(f"反轉合併倍數 {RATIO:.4f}")


def label_bars(ax, bars, vals, fmt="{:.2f}", dy=0.02):
    for b, v, d in zip(bars, vals, np.broadcast_to(dy, len(vals))):
        ax.text(b.get_x() + b.get_width() / 2, b.get_height() + d, fmt.format(v), ha="center", va="bottom", fontsize=11, color=INK)


# ---------- 圖 1：砝碼造成的 J 增加量 ----------
fig, ax = plt.subplots(figsize=(9, 5.2))
x = np.arange(2) * 1.3
pred = [p["mr2"] for p in rev]
meas = [p["dJw"] for p in rev]
smeas = [p["sdw"] for p in rev]
W = 0.42
ax.bar(x - W / 2 - 0.02, pred, W, color="#c9c7bf", label="預測 m·r²（只用砝碼的質量與位置）", zorder=2)
ax.bar(x + W / 2 + 0.02, meas, W, color=BLUE, yerr=smeas, capsize=5, ecolor=DARK, error_kw=dict(lw=1.2),
       label="實測（夾砝碼與不夾的 J 差）", zorder=2)
for i, p in enumerate(rev):
    xp, xm = x[i] - W / 2 - 0.02, x[i] + W / 2 + 0.02
    ax.text(xp, pred[i] + 0.015, f"{pred[i]:.2f}", ha="center", va="bottom", fontsize=13, color=INK2)
    ax.text(xm, meas[i] + smeas[i] + 0.015, f"{meas[i]:.2f}", ha="center", va="bottom", fontsize=13, color=BLUE,
            fontweight="bold")
    ytop = meas[i] + smeas[i] + 0.11                       # 預測到實測的括號與倍數
    ax.plot([xp, xp, xm, xm], [pred[i] + 0.07, ytop, ytop, meas[i] + smeas[i] + 0.08], color=INK2, lw=1.1)
    ax.text(x[i], ytop + 0.015, f"× {p['ratio_w']:.2f}", ha="center", va="bottom", fontsize=14, color=INK,
            fontweight="bold")
ax.set_xticks(x, [f"DIAG\n砝碼離 J0 軸 {rev[0]['mr2'] ** 0.5:.2f} m", f"EXT45\n砝碼離 J0 軸 {rev[1]['mr2'] ** 0.5:.2f} m"],
              fontsize=12.5)
ax.tick_params(axis="x", length=0, pad=8)
ax.set_ylabel("J0 增加量 [kg·m²]")
ax.set_ylim(0, 1.08)
ax.set_xlim(-0.75, 2.05)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(alpha=1, axis="y", ls=(0, (4, 3))); ax.set_axisbelow(True)
ax.legend(loc="upper left", fontsize=11, frameon=False, ncol=2, bbox_to_anchor=(0, 1.02), handlelength=1.2)
ax.set_title(f"圖 1　夾上 1 kg 砝碼：砝碼越遠，J 增加越多；兩個姿態都是預測的約 {RATIO:.2f} 倍",
             loc="left", fontsize=13.5, pad=26, color=INK)
fig.text(0.01, 0.005, "反轉；6 個條件依單輪回歸的標準差加權平均，誤差棒為加權平均的標準差；m = 1.000 kg（標稱值，未秤重）",
         ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(os.path.join(OUT, "fig1_v2_known_mass_bars.png"), dpi=150, bbox_inches="tight", pad_inches=0.15)
plt.close(fig)

# ---------- 圖 2：校正後的 J、官方模型、控制器模型 ----------
runs = {(r["pose"], bool(r["w"])): r for r in R["runs"] if r["dir"] < 0}
keys = [("DIAG", False), ("DIAG", True), ("EXT45", False), ("EXT45", True)]
names = ["DIAG\n不夾", "DIAG\n夾砝碼", "EXT45\n不夾", "EXT45\n夾砝碼"]
X = np.array([runs[k]["X"] for k in keys])
Jc = [np.array(runs[k]["J"]) / RATIO for k in keys]
Jm = np.array([j.mean() for j in Jc])
lo = Jm - np.array([j.min() for j in Jc])
hi = np.array([j.max() for j in Jc]) - Jm
CTRL = X + km.J_DRIVE
fig, ax = plt.subplots(figsize=(11, 6.0))
x = np.arange(4)
b1 = ax.bar(x - 0.27, X, 0.26, color=GRAY, label="官方模型：連桿 + 負載（點質量）")
b2 = ax.bar(x, Jm, 0.26, color=BLUE, yerr=[lo, hi], capsize=4, ecolor=INK2, label=f"實測 J ÷ {RATIO:.2f}（砝碼校正後）")
b3 = ax.bar(x + 0.27, CTRL, 0.26, color=DARK, label=f"控制器模型 = 官方模型 + {km.J_DRIVE:.3f}")
label_bars(ax, b1, X)
label_bars(ax, b2, Jm, dy=0.02 + hi)
label_bars(ax, b3, CTRL)
i = 0                                                     # 在 DIAG 不夾標出 1.915 的差
xb = x[i] + 0.27 + 0.17
ax.annotate("", (xb, CTRL[i]), (xb, X[i]), arrowprops=dict(arrowstyle="<->", color=INK, lw=1.3))
ax.text(xb + 0.04, (X[i] + CTRL[i]) / 2, f"{km.J_DRIVE:.3f}\n電流中\n看不到", va="center", fontsize=10.5, color=INK)
ax.set_xticks(x, names)
ax.set_ylabel("J0 慣量 [kg·m²]")
ax.set_ylim(0, 6.6)
ax.set_xlim(-0.6, 3.85)
ax.grid(alpha=1, axis="y"); ax.set_axisbelow(True)
ax.legend(loc="upper left", fontsize=10.5, frameon=False)
ax.text(0.995, -0.13, "反轉；藍色為 6 個條件平均，誤差棒為 6 個條件的最小～最大。控制器模型依第 1 節的組成，以正確負載計算\n"
        "（夾砝碼時教導器負載設 1.40 kg，控制器回報的 J_model 少算 0.9 kg）",
        transform=ax.transAxes, ha="right", va="top", fontsize=9, color=INK2)
fig.suptitle(f"圖 2　校正後的實測 J0 和官方模型一致；控制器模型每次都多 {km.J_DRIVE:.3f}", fontsize=13)
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig2_v2_calibrated_vs_models.png"), dpi=150)
plt.close(fig)

# ---------- 圖 3：比例 c 和頻率 ----------
fig, ax = plt.subplots(figsize=(8.5, 5.4))
for big, mfc, lab, dx in ((True, BLUE, "速度振幅 0.06 rad/s", -0.05), (False, "white", "速度振幅 0.03、0.04 rad/s", 0.05)):
    idx = [k for k, c in enumerate(R["conditions"]) if (c[1] > 0.05) == big]
    idx.sort(key=lambda k: R["conditions"][k][0])
    f = np.array([R["conditions"][k][0] for k in idx])
    c = np.array([R["fit"][k]["c"] for k in idx])
    sc = np.array([R["fit"][k]["sc"] for k in idx])
    ax.errorbar(f + dx, c, sc, fmt="o-" if big else "o--", ms=9, capsize=3, color=BLUE, mfc=mfc, mew=1.6, lw=1.3, label=lab)
    for fi, ci in zip(f, c):
        ax.text(fi + dx + (-0.09 if big else 0.09), ci, f"{ci:.2f}", ha="right" if big else "left", va="center", fontsize=10.5, color=INK2)
ax.set_xticks([1, 2, 3], ["1 Hz", "2 Hz", "3 Hz"])
ax.set_xlim(0.5, 3.5)
ax.set_ylim(1.1, 1.5)
ax.set_xlabel("調變頻率")
ax.set_ylabel("比例 c（實測 ÷ 官方模型）")
ax.grid(alpha=1, axis="y"); ax.set_axisbelow(True)
ax.legend(loc="upper left", fontsize=10.5, frameon=False)
ax.text(0.995, -0.15, "反轉 4 次（2 姿態 × 夾或不夾）合併擬合 J = c·X + J_d；誤差棒為擬合的標準差", transform=ax.transAxes,
        ha="right", va="top", fontsize=9, color=INK2)
fig.suptitle("圖 3　比例 c 隨頻率增加；同頻率、不同振幅幾乎相同", fontsize=13)
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig3_v2_ratio_vs_freq.png"), dpi=150)
plt.close(fig)
print("calibrated", np.round(Jm, 3), "X", np.round(X, 3), "range", np.round(Jm - lo, 2), np.round(Jm + hi, 2))
