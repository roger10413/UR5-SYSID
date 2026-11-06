"""11/4 報告 3.3 的圖：J = c·X + J_d 是什麼意思。
橫軸 X = 官方模型（連桿 + 負載 + 砝碼），縱軸 = 量到的 J；4 次實驗畫直線，延長到 X = 0 的高度就是固定項 J_d。
同一個式子用在控制器預期的力矩（有 1.915）與實測電流（沒有）。
需先執行 current_signals_1104.py（讀取 OUT/current_signals_1104.json）。輸出 OUT/figs/w1104_fig3_formula.png。
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
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#8d8b83"
ORANGE, BLUE = "#d9822b", "#2a78d6"
SHOW = 1                                    # 主要畫出的條件：CONDITIONS 順序中的 2 Hz × 0.06
# X：4 次的官方模型（DIAG 不夾、DIAG 夾、EXT45 夾、EXT45 不夾）；控制器只知道教導器設定的負載
X_ACT = np.array([1.134, 1.383, 3.949, 3.317])
X_SET = np.array([1.134, 1.159, 3.380, 3.317])

A = json.load(open(os.path.join(OUT, "current_signals_1104.json"), encoding="utf-8"))["A"]
SIG = [("target_moment", X_SET, ORANGE, "s", "控制器預期的力矩"),
       ("actual_current", X_ACT, BLUE, "o", "實測電流（我們算 J 用的）")]

fig, ax = plt.subplots(figsize=(10, 6.2))
xx = np.linspace(0, 4.25, 50)
ax.axvspan(0, 1.05, color="#f1f0ec", zorder=0)
ax.text(0.52, 6.6, "沒有資料\n直線延長過來", ha="center", va="top", fontsize=10.5, color=INK2)
for name, X, col, mk, lab in SIG:
    J =np.array(A[name][SHOW]["J"])
    p = np.polyfit(X, J, 1)
    ax.plot(xx[xx >= 1.05], np.polyval(p, xx[xx >= 1.05]), color=col, lw=2.2, zorder=3)
    ax.plot(xx[xx <= 1.05], np.polyval(p, xx[xx <= 1.05]), color=col, lw=2.2, ls=(0, (4, 3)), zorder=3)
    ax.scatter(X, J, s=70, marker=mk, color=col, edgecolor="white", linewidth=2, zorder=4)
    ax.scatter([0], [p[1]], s=90, marker=mk, color="white", edgecolor=col, linewidth=2.2, zorder=5)
    jd = [rr["Jd_over_c"] for rr in A[name]]
    print(f"{name:<15} 2 Hz×0.06 斜率 c {p[0]:.3f} 截距 J_d {p[1]:+.3f}；6 條件 J_d÷c {min(jd):+.2f}～{max(jd):+.2f}")
    if name == "target_moment":
        ax.annotate(f"X = 0 時還剩約 {p[1]:.1f}\n→ 固定項有 1.915", (0, p[1]), (0.18, p[1] + 0.95), fontsize=11.5,
                    color=INK, arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
        ax.text(2.3, np.polyval(p, 2.3) + 0.35, lab, ha="center", va="bottom", fontsize=12, color=INK,
                fontweight="bold")
    else:
        ax.annotate("X = 0 時約 0\n→ 沒有固定項", (0, p[1]), (0.18, -0.85), fontsize=11.5, color=INK,
                    arrowprops=dict(arrowstyle="-", color=MUTED, lw=1))
        ax.text(3.0, np.polyval(p, 3.0) - 0.75, lab, ha="center", va="top", fontsize=12, color=INK,
                fontweight="bold")
for x, t in ((1.26, "DIAG\n不夾、夾"), (3.63, "EXT45\n不夾、夾")):
    ax.text(x, -0.95, t, ha="center", va="top", fontsize=10.5, color=INK2)
ax.text(3.1, 1.75, "J = c·X + J_d", ha="center", va="center", fontsize=17, color=INK, fontweight="bold")
ax.text(3.1, 0.85, "c·X：跟著姿態、砝碼變的部分（斜率 c）\nJ_d：不隨姿態改變的固定項（直線在 X = 0 的高度）",
        ha="center", va="center", fontsize=10.5, color=INK2, linespacing=1.5)
ax.axhline(0, color=MUTED, lw=1.0, zorder=1)
ax.set_xlim(-0.12, 4.3); ax.set_ylim(-1.6, 6.9)
ax.set_xlabel("X：官方模型算的慣量（連桿 + 負載 + 砝碼）[kg·m²]")
ax.set_ylabel("量到的 J [kg·m²]")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(alpha=1, ls=(0, (4, 3))); ax.set_axisbelow(True)
ax.set_title("圖 3　同一個式子：控制器的預期力矩算得出 1.915，實測電流算不出來", loc="left", fontsize=13.5, pad=12,
             color=INK)
fig.text(0.01, 0.005, "10/29 反轉 4 次，2 Hz × 0.06（其他 5 個條件結論相同）。轉子慣量不隨姿態改變，若在訊號裡，只會出現在 J_d。"
         "控制器訊號的 X 用教導器設定的負載", ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 1))
fig.savefig(os.path.join(FIGS, "w1104_fig3_formula.png"), dpi=150, bbox_inches="tight", pad_inches=0.15)
print(os.path.join(FIGS, "w1104_fig3_formula.png"))
