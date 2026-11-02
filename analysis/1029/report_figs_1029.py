"""10/29 報告圖：已知質量的比例；實測 J 對官方點質量模型（合併擬合）；與 10/22 比較。
需先執行 known_mass_1029.py（讀取 OUT/known_mass_1029.json）。"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 11,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
import known_mass_1029 as km

R = json.load(open(os.path.join(km.OUT, "known_mass_1029.json"), encoding="utf-8"))
OUT = os.path.join(km.OUT, "figs"); os.makedirs(OUT, exist_ok=True)
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
FREQ_C = {1.0: BLUE, 2.0: ORANGE, 3.0: AQUA}                 # 圖 2：依頻率上色
SMALL = [c[1] < 0.05 for c in R["conditions"]]               # 速度振幅 0.03、0.04 為空心
ORDER = [0, 4, 1, 3, 2, 5]                                   # 依頻率排列：1 Hz×0.06、1 Hz×0.03、2 Hz×0.06、2 Hz×0.03、3 Hz×0.06、3 Hz×0.04
runs = R["runs"]
A_mean = np.mean([r["A"] for r in runs], axis=0)
xt = [f"{R['cond'][i]}\nA = {A_mean[i]:.2f}" for i in ORDER]

# ---------- 圖 1：已知質量的比例 ----------
fig, ax = plt.subplots(figsize=(12.5, 6.0))
xb = np.arange(len(ORDER))
pr = {(p["pose"], p["dir"]): p for p in R["pairs"]}
series = [(("DIAG", -1), BLUE, "o", True, -0.2, "砝碼：DIAG 反轉（r = 0.50 m）"),
          (("EXT45", -1), ORANGE, "s", True, -0.07, "砝碼：EXT45 反轉（r = 0.79 m）"),
          (("DIAG", +1), BLUE, "o", False, 0.06, "砝碼：DIAG 正轉")]
for key, col, mk, filled, dx, lab in series:
    p = pr[key]
    ax.errorbar(xb + dx, [p["ratio"][i] for i in ORDER], [p["sratio"][i] for i in ORDER], fmt=mk, ms=8, capsize=3,
                color=col, mfc=col if filled else "white", mew=1.6, lw=1.1, label=lab)
pose = np.mean([q["ratio"] for q in R["pose_ratio"]], axis=0)
spose = np.mean([q["sratio"] for q in R["pose_ratio"]], axis=0)
ax.errorbar(xb + 0.19, pose[ORDER], spose[ORDER], fmt="D", ms=7, capsize=3, color=INK2, mew=1.4, lw=1.1,
            label="姿態差：(EXT45 - DIAG) ÷ 模型差（10/29 反轉，不夾與夾砝碼平均）")
R22 = json.load(open(km.R1022, encoding="utf-8"))
r22 = [R22["dJ"][i] / R22["dJ_model"] for i in range(4)]
pos22 = {0: 0, 1: 2, 3: 3, 2: 4}                             # 10/22 只有前 4 個條件 → 在 ORDER 中的位置
ax.plot([pos22[i] + 0.3 for i in range(4)], r22, ls="", marker="D", ms=7, color=MUTED, mfc="white", mew=1.4,
        label="姿態差（10/22，10/28 報告的 1.25～1.40 倍）")
ax.axhline(1.0, color=INK, lw=1.2)
ax.text(len(ORDER) - 0.45, 1.0 + 0.03, "實測 = 模型（Kt × 101 = 13.64 N·m/A 正確時）", ha="right", va="bottom", fontsize=9.5, color=INK2)
for x in (1.5, 3.5):
    ax.axvline(x, color="#e1e0d9", lw=1)
ax.set_xticks(xb, xt)
ax.set_xlim(-0.6, len(ORDER) - 0.4)
ax.set_ylim(0.6, 2.0)
ax.set_ylabel("實測 ÷ 模型（Kt = 0.135）")
ax.grid(alpha=1, axis="y")
h, lab = ax.get_legend_handles_labels()
k = [lab.index(s) for s in sorted(lab, key=lambda s: (not s.startswith("砝碼"), "10/22" in s))]   # 砝碼在前，10/22 最後
fig.legend([h[i] for i in k], [lab[i] for i in k], loc="lower center", fontsize=9.5, frameon=False, ncol=3,
           bbox_to_anchor=(0.5, 0.0))
ax.text(0.995, 0.02, "砝碼：ΔJ ÷ m·r²（m = 1.000 kg 標稱值，未秤重）　誤差棒：單輪回歸的標準差", transform=ax.transAxes,
        ha="right", va="bottom", fontsize=9, color=INK2)
fig.suptitle("圖 1　已知質量：砝碼造成的 J 差平均是 m·r² 的 1.23～1.37 倍，與兩姿態差對官方模型的倍數相近", fontsize=12.5)
fig.tight_layout(rect=(0, 0.1, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig1_known_mass_ratio.png"), dpi=150)
plt.close(fig)

# ---------- 圖 2：實測 J 對官方點質量模型 ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.8), gridspec_kw=dict(width_ratios=[1.35, 1]))
rev = [r for r in runs if r["dir"] < 0]
xx = np.linspace(0, 4.3, 50)
ax1.plot(xx, xx + km.J_DRIVE, "--", color=INK, lw=1.6, label="控制器模型：J_model = X + 1.915")
ctrl = [r for r in rev if not r["w"]]
ax1.plot([r["X"] for r in ctrl], [r["J_model"] for r in ctrl], ls="", marker="D", ms=8, color=INK,
         label="控制器 J_model（教導器負載 1.30 kg 的兩次）")
ax1.plot(xx, xx, ":", color=MUTED, lw=1.4, label="J = X（Kt 正確、電流不含驅動端慣量時）")
for i in range(len(R["cond"])):
    f, col = R["fit"][i], FREQ_C[R["conditions"][i][0]]
    ax1.plot(xx, f["c"] * xx + f["Jd"], "-", color=col, lw=1.0, alpha=0.8)
    for r in rev:
        ax1.plot(r["X"], r["J"][i], ls="", marker="s" if r["w"] else "o", ms=7, color=col,
                 mfc="white" if SMALL[i] else col, mew=1.4)
ax1.set_xlim(0, 4.3)
ax1.set_ylim(0, 6.6)
ax1.set_xlabel("官方點質量模型 X = 連桿 + 末端工具 1.30 kg + 砝碼 m·r²  [kg·m²]")
ax1.set_ylabel("實測 J [kg·m²]（Kt = 0.135）")
ax1.grid(alpha=1)
h = [plt.Line2D([], [], ls="-", marker="o", ms=7, color=BLUE, label="1 Hz（c = 1.23）"),
     plt.Line2D([], [], ls="-", marker="o", ms=7, color=ORANGE, label="2 Hz（c = 1.28～1.29）"),
     plt.Line2D([], [], ls="-", marker="o", ms=7, color=AQUA, label="3 Hz（c = 1.40）"),
     plt.Line2D([], [], ls="", marker="o", ms=7, color=INK2, mfc="white", mew=1.4, label="空心：小振幅（0.03、0.04 rad/s）"),
     plt.Line2D([], [], ls="", marker="o", ms=7, color=INK2, label="圓：不夾"),
     plt.Line2D([], [], ls="", marker="s", ms=7, color=INK2, label="方：夾砝碼")]
l1 = ax1.legend(handles=h, loc="lower right", fontsize=9, framealpha=0.95, title="實測（反轉 4 次）與直線 J = c·X + J_d", title_fontsize=9)
ax1.add_artist(l1)
ax1.legend(loc="upper left", fontsize=9, framealpha=0.95)
for x, y, name in ((1.26, 2.45, "DIAG"), (3.63, 6.25, "EXT45")):
    ax1.text(x, y, name, ha="center", fontsize=10, color=INK2)
ax1.set_title("(a) 2 姿態 × 夾或不夾：實測 J 落在通過原點的直線上", fontsize=11.5)

xb = np.arange(len(ORDER))
for k, i in enumerate(ORDER):
    f, col = R["fit"][i], FREQ_C[R["conditions"][i][0]]
    ax2.errorbar(k, f["Jd_c"], f["sJd_c"], fmt="o", ms=8, capsize=3, color=col, mfc="white" if SMALL[i] else col, mew=1.6)
ax2.axhline(km.J_DRIVE, color=INK, ls="--", lw=1.6)
ax2.text(len(ORDER) - 0.5, km.J_DRIVE - 0.06, "控制器模型的驅動端慣量 1.915", ha="right", va="top", fontsize=9.5, color=INK)
ax2.axhline(0, color=INK2, lw=1.0)
ax2.set_xticks(xb, xt, fontsize=9)
ax2.set_xlim(-0.6, len(ORDER) - 0.4)
ax2.set_ylim(-0.4, 2.2)
ax2.set_ylabel("截距 J_d ÷ c [kg·m²]")
ax2.grid(alpha=1, axis="y")
ax2.set_title("(b) 截距 -0.13～+0.05：電流中沒有 1.915 這一項", fontsize=11.5)
ax2.text(0.98, 0.42, "誤差棒：擬合的標準差\n（卡方 ÷ 自由度 > 1 時依殘差放大）", transform=ax2.transAxes, ha="right", fontsize=9, color=INK2)
fig.suptitle("圖 2　實測 J 與官方點質量模型成正比（比例 c），沒有控制器模型中的驅動端慣量", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig2_fit_vs_model.png"), dpi=150)
plt.close(fig)

# ---------- 圖 3：與 10/22 比較 ----------
fig, ax = plt.subplots(figsize=(11, 4.8))
cr = R["cross"]
items = [("DIAG", +1, BLUE, False, "DIAG 正轉（對 10/22 兩次）"), ("DIAG", -1, BLUE, True, "DIAG 反轉"),
         ("EXT45", -1, ORANGE, True, "EXT45 反轉")]
xb = np.arange(4)
for k, (pz, d, col, filled, lab) in enumerate(items):
    sel = [c for c in cr if c["pose"] == pz and c["dir"] == d]
    for j, c in enumerate(sel):
        ax.errorbar(xb + (k - 1) * 0.22 + j * 0.07, c["pct"], c["spct"], fmt="o", ms=8, capsize=3, color=col,
                    mfc=col if filled else "white", mew=1.6, lw=1.1, label=lab if j == 0 else None)
ax.axhspan(-3, 3, color="#e1e0d9", alpha=0.6, lw=0)
ax.axhline(0, color=INK2, lw=1.0)
ax.set_xticks(xb, [f"{R['cond'][i]}\nA = {A_mean[i]:.2f}" for i in range(4)])
ax.set_ylabel("J 的差（10/29 ÷ 10/22 - 1）[%]")
ax.set_ylim(-22, 11)
ax.grid(alpha=1, axis="y")
ax.legend(loc="lower right", fontsize=9.5, framealpha=0.95)
ax.text(0.01, 0.97, "灰帶：±3%　誤差棒：兩天單輪回歸標準差合成　J0 溫度：10/22 33.3 °C，10/29 31.3～32.1 °C",
        transform=ax.transAxes, va="top", fontsize=9, color=INK2)
fig.suptitle("圖 3　隔週重複性（不夾）：2、3 Hz 在 ±6% 內，DIAG 的 1 Hz × 0.06 低了 11～17%", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
fig.savefig(os.path.join(OUT, "fig3_cross_day.png"), dpi=150)
plt.close(fig)
print("圖已存到", OUT)
