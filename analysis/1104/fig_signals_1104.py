"""11/4 報告圖 3：同樣的回歸與 J = c·X + J_d 擬合，換成控制器的 4 種訊號，截距 J_d ÷ c 有沒有 1.915。
需先執行 current_signals_1104.py（讀取 OUT/current_signals_1104.json）。輸出 OUT/figs/w1104_fig3_signals_intercept.png。
（Microsoft JhengHei 沒有「−」，數值用 ASCII 負號。）"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
OUT = os.path.join(r"D:\cluade_code\ur5", "analysis_1104")
FIGS = os.path.join(OUT, "figs"); os.makedirs(FIGS, exist_ok=True)
INK, INK2, DARK = "#0b0b0b", "#52514e", "#3a3936"
ORANGE, BLUE = "#d9822b", "#2a78d6"
J_DRIVE = 1.915

A = json.load(open(os.path.join(OUT, "current_signals_1104.json"), encoding="utf-8"))["A"]
SIG = [("target_moment", "target_moment\n控制器預期的力矩", ORANGE),
       ("target_current", "target_current\n前饋電流", ORANGE),
       ("joint_control_output", "joint_control_output\n最後送出的電流命令", BLUE),
       ("actual_current", "actual_current\n實測電流", BLUE)]

fig, ax = plt.subplots(figsize=(10, 5.6))
x = np.arange(len(SIG))
for i, (name, _, col) in enumerate(SIG):
    d = np.array([r["Jd_over_c"] for r in A[name]])
    m = d.mean()
    ax.bar(x[i], m, 0.58, color=col, alpha=0.9, zorder=2)
    ax.errorbar(x[i], m, yerr=[[m - d.min()], [d.max() - m]], fmt="none", ecolor=DARK, capsize=6, lw=1.2, zorder=3)
    ax.text(x[i], max(m, 0) + (d.max() - m) + 0.06, f"{m:.2f}", ha="center", va="bottom", fontsize=14,
            color=col, fontweight="bold")
    print(f"{name:<22} 平均 {m:+.3f}  範圍 {d.min():+.3f} ～ {d.max():+.3f}")
ax.axhline(J_DRIVE, color=INK2, ls=(0, (5, 3)), lw=1.3, zorder=1)
ax.text(3.45, J_DRIVE + 0.04, "控制器設定的 1.915", ha="right", va="bottom", fontsize=11.5, color=INK2)
ax.axhline(0, color="#8d8b83", lw=1.0, zorder=1)
for xc, t, col in ((0.5, "控制器的前饋（模型算的）", ORANGE), (2.5, "最後的命令與實測（真正送進馬達的）", BLUE)):
    ax.text(xc, 2.72, t, ha="center", va="center", fontsize=12.5, color=col, fontweight="bold")
ax.plot([1.5, 1.5], [-0.3, 2.85], color="#c3c2b7", lw=1.0, ls=":", zorder=0)
ax.set_xticks(x, [s[1] for s in SIG], fontsize=11)
ax.tick_params(axis="x", length=0, pad=8)
ax.set_ylabel("截距 J_d ÷ c [kg·m²]")
ax.set_ylim(-0.3, 2.9)
ax.set_xlim(-0.6, 3.6)
for s in ("top", "right", "left"):
    ax.spines[s].set_visible(False)
ax.tick_params(axis="y", length=0)
ax.grid(alpha=1, axis="y", ls=(0, (4, 3))); ax.set_axisbelow(True)
ax.set_title("圖 3　1.915 在控制器的前饋裡有，在最後送出的電流命令和實測電流裡都沒有", loc="left", fontsize=13.5,
             pad=14, color=INK)
fig.text(0.01, 0.005, "10/29 反轉 4 次（2 姿態 × 夾或不夾），每個條件擬合 J = c·X + J_d；長條為 6 個條件的平均，誤差棒為最小～最大。"
         "控制器訊號的 X 用教導器設定的負載",
         ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(os.path.join(FIGS, "w1104_fig3_signals_intercept.png"), dpi=150, bbox_inches="tight", pad_inches=0.15)
print(os.path.join(FIGS, "w1104_fig3_signals_intercept.png"))
