"""10/15 報告圖：定速 B、Tc 對溫度；定加速度三角形 vs 梯形（J 回歸、殘差對行程）。
需先執行 constvel_1015.py、accel_v3_1015.py；卡爾曼殘差圖由 kalman_figs_1015.py 產生。"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 11})
import data_1015 as dd

A = os.path.join(dd.BASE, "analysis_1015")
OUT = os.path.join(A, "figs"); os.makedirs(OUT, exist_ok=True)
CB, CR = "#2366b5", "#c8452f"
CV = list(np.load(os.path.join(A, "constvel_1015.npy"), allow_pickle=True))
V3 = list(np.load(os.path.join(A, "accel_v3_1015.npy"), allow_pickle=True))
V6 = np.array([0.02, 0.05, 0.08, 0.12, 0.16, 0.20])

# ---------- 圖 1：定速 B、Tc 對溫度 ----------
fig, axs = plt.subplots(1, 3, figsize=(15, 4.6))
for s, c, lab in ((+1, CB, "正轉"), (-1, CR, "反轉")):
    R = sorted([x for x in CV if x["dir"] == s], key=lambda x: x["n"]); T = np.array([x["Tm"] for x in R])
    for ax, key, se in ((axs[0], "B", "seB"), (axs[1], "Tc", "seTc")):
        ax.errorbar(T, [x[key] for x in R], [x[se] for x in R], fmt="o", color=c, ms=7, capsize=4, label=lab)
        p = np.polyfit(T, [x[key] for x in R], 1); tt = np.linspace(29.5, 32.1, 10)
        ax.plot(tt, np.polyval(p, tt), "-", color=c, lw=1, alpha=0.7)
        for x in R:
            ax.annotate(f"第{x['n']}輪", (x["Tm"], x[key]), textcoords="offset points", xytext=(6, -12 if s > 0 else 6), fontsize=8, color=c)
        if key == "B":
            ax.text(0.03, 0.06 if s > 0 else 0.14, f"{lab}：dB/dT = {p[0]:+.2f} /°C", transform=ax.transAxes, color=c, fontsize=10)
    L = np.array([x["levels"] for x in R]); sl = [np.polyfit(T, L[:, i], 1)[0] for i in range(6)]
    axs[2].plot(V6, sl, "o-", color=c, label=lab)
axs[0].set_ylabel("B [N·m·s/rad]"); axs[1].set_ylabel("Tc [N·m]")
axs[0].set_title("B 對 J0 溫度（± 為單輪回歸的標準差）"); axs[1].set_title("Tc 對 J0 溫度")
axs[2].axhline(0, color="k", lw=0.6); axs[2].set_xlabel("角速度 [rad/s]"); axs[2].set_ylabel("d|τ|/dT [N·m/°C]")
axs[2].set_title("各速度檔力矩隨溫度的變化率")
for ax in axs[:2]:
    ax.set_xlabel("J0 溫度 [°C]（各輪平均）")
for ax in axs:
    ax.grid(alpha=0.2); ax.legend(fontsize=9)
fig.suptitle("10/15 固定角度窗口定速 3 輪：溫度升高約 2 °C，B 下降、Tc 不變", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_constvel_temperature.png"), dpi=140); plt.close(fig)

# ---------- 圖 2：三角形 vs 梯形，(τacc − τdec)/2 對 a ----------
sty = {"tri1": ("o", "-", "三角形 第1次"), "trap": ("s", "--", "梯形"), "tri2": ("^", ":", "三角形 第2次")}
fig, axs = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
for ax, s, c, lab in ((axs[0], +1, CB, "正轉"), (axs[1], -1, CR, "反轉")):
    for o in [x for x in V3 if x["r"]["dir"] == s]:
        key = "trap" if o["r"]["kind"] == "trap" else f"tri{o['r']['n']}"
        mk, ls, nm = sty[key]
        a = np.array([x["a"] for x in o["L"]]); h = np.array([x["h"] for x in o["L"]])
        aa = np.linspace(0, 1.3, 10)
        ax.plot(a, h, mk, color=c, ms=8, mfc=c if key == "trap" else "white")
        ax.plot(aa, o["J"] * aa + o["c"], ls, color=c, lw=1.3,
                label=f"{nm}（{o['Tm']:.1f} °C）：J = {o['J']:.2f}，截距 {o['c']:+.2f} → Δ = {-2 * o['c']:.2f}")
    ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
    ax.set_xlabel("加速度 a [rad/s²]"); ax.set_title(lab); ax.grid(alpha=0.2); ax.legend(fontsize=9, loc="upper left")
axs[0].set_ylabel("(τ_acc - τ_dec) / 2 [N·m]")
fig.suptitle("定加速度 v3：同一天的三角形（定速段 0）與梯形（定速段 4 馬達圈）含截距回歸", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_v3_tri_vs_trap_regression.png"), dpi=140); plt.close(fig)

# ---------- 圖 3：殘差對起步後行程 ----------
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2), sharey=True)
pc = {"加速段": "#3a9e5c", "定速段": "#888888", "減速段": "#b5397a"}
for ax, s, lab in ((axs[0], +1, "正轉"), (axs[1], -1, "反轉")):
    for o in [x for x in V3 if x["r"]["dir"] == s]:
        key = "trap" if o["r"]["kind"] == "trap" else f"tri{o['r']['n']}"
        mk, ls, nm = sty[key]
        for k, (xs, ys) in o["travel"].items():
            if len(xs):
                ax.plot(xs, ys, ls, marker=mk, color=pc[k], ms=4, alpha=0.6 if key == "tri2" else 1, label=f"{nm} {k}")
    ax.axhline(0, color="k", lw=0.6); ax.grid(alpha=0.2); ax.set_title(lab); ax.set_xlabel("每次循環起步後行程 [°]（已扣除起始錯開）")
    ax.legend(fontsize=8, ncol=2)
axs[0].set_ylabel("|τ| - 同溫度定速摩擦 - J·a [N·m]")
fig.suptitle("定加速度 v3 殘差對起步後行程（定速摩擦取 3 輪對溫度回歸在該次溫度的值；J 為各次含截距回歸值）", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_v3_residual_travel.png"), dpi=140); plt.close(fig)
print("完成：", OUT)
