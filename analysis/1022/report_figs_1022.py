"""10/22 報告圖：J 對加速度振幅（DIAG、EXT45）；實測 J 與控制器模型 J_model 的比較、兩姿態差 ΔJ。
需先執行 velmod_1022.py（讀取 OUT/velmod_1022.json）。"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 11})
import velmod_1022 as vv

R = json.load(open(os.path.join(vv.OUT, "velmod_1022.json"), encoding="utf-8"))
OUT = os.path.join(vv.OUT, "figs"); os.makedirs(OUT, exist_ok=True)
FC = {0: "#2366b5", 1: "#e08a1e", 2: "#2a9d5c", 3: "#e08a1e"}          # 依頻率上色：1 Hz 藍、2 Hz 橘、3 Hz 綠
FILLED = {0: True, 1: True, 2: True, 3: False}                          # 速度振幅 0.06 實心、0.03 空心
MK = {"正轉": "o", "反轉": "^", "正轉（重跑）": "s"}
CM = "#c8452f"
POSE_T = {"DIAG": "DIAG（J1/J2 = -90°/90°，法蘭離 J0 軸 0.50 m）", "EXT45": "EXT45（J1/J2 = -45°/45°，0.79 m）"}

# ---------- 圖 1：J 對加速度振幅 ----------
fig, axs = plt.subplots(1, 2, figsize=(14, 5.2))
for ax, pose in zip(axs, ("DIAG", "EXT45")):
    P = R["poses"][pose]
    runs = [r for r in R["runs"] if r["pose"] == pose]
    for k, r in enumerate(runs):
        kind = r["label"].split(" ", 1)[1]
        for i in range(4):
            x = r["A"][i] + (k - (len(runs) - 1) / 2) * 0.012
            ax.errorbar(x, r["J"][i], r["sJ"][i], fmt=MK[kind], ms=7, capsize=3, color=FC[i],
                        mfc=FC[i] if FILLED[i] else "white", mew=1.4, lw=1)
    aa = np.linspace(0.15, 1.3, 100)
    a, b = P["lin"]
    ax.plot(aa, a + b * aa, "-", color="0.35", lw=1.2, label=f"J = a + b·A（A → 0：{a:.2f}）")
    ax.plot(aa, P["J_inf"] + P["k"] / aa, "--", color="0.35", lw=1.2, label=f"J = J∞ + k/A（J∞ = {P['J_inf']:.2f}）")
    ax.set_title(POSE_T[pose], fontsize=11.5)
    ax.set_xlabel("實測加速度振幅 A [rad/s²]")
    ax.set_xlim(0.1, 1.33)
    lo, hi = min(min(r["J"]) for r in runs), max(max(r["J"]) for r in runs)
    ax.set_ylim(lo - 0.25, hi + 0.15)
    ax.grid(alpha=0.3)
    ax.text(0.98, 0.04, f"控制器模型 J_model = {P['J_model']:.2f}（超出圖範圍）", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=9.5, color=CM)
    ax.legend(loc="upper left", fontsize=9.5, framealpha=0.9)
axs[0].set_ylabel("J [kg·m²]（Kt = 0.135）")
h = [plt.Line2D([], [], ls="", marker="o", ms=8, color=FC[0], label="1 Hz × 0.06 rad/s"),
     plt.Line2D([], [], ls="", marker="o", ms=8, color=FC[1], label="2 Hz × 0.06 rad/s"),
     plt.Line2D([], [], ls="", marker="o", ms=8, color=FC[2], label="3 Hz × 0.06 rad/s"),
     plt.Line2D([], [], ls="", marker="o", ms=8, color=FC[3], mfc="white", mew=1.4, label="2 Hz × 0.03 rad/s"),
     plt.Line2D([], [], ls="", marker="o", ms=7, color="0.3", label="正轉"),
     plt.Line2D([], [], ls="", marker="^", ms=7, color="0.3", label="反轉"),
     plt.Line2D([], [], ls="", marker="s", ms=7, color="0.3", label="正轉（EXT45 之後重跑）")]
fig.legend(handles=h, loc="lower center", ncol=7, fontsize=9.5, frameon=False, bbox_to_anchor=(0.5, -0.01))
fig.suptitle("圖 1　速度調變法各條件的 J：隨加速度振幅增加，同振幅不同頻率（1 Hz 實心藍 vs 2 Hz 空心橘）則相同", fontsize=12.5)
fig.tight_layout(rect=(0, 0.06, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig1_J_vs_accel.png"), dpi=150)
plt.close(fig)

# ---------- 圖 2：與控制器模型比較 ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5.4), gridspec_kw=dict(width_ratios=[1, 1.15]))
D, E = R["poses"]["DIAG"], R["poses"]["EXT45"]
xs = np.array([0, 1])
for i in (0, 3, 1, 2):
    y = [D["J"][i], E["J"][i]]
    ax1.plot(xs, y, "-", color=FC[i], lw=1.6, marker="o", ms=7, mfc=FC[i] if FILLED[i] else "white", mew=1.4,
             label=f"實測 {R['cond'][i]}（ΔJ = {y[1] - y[0]:.2f}）")
ax1.plot(xs, [D["J_model"], E["J_model"]], "--", color=CM, lw=2, marker="D", ms=7,
         label=f"控制器模型 J_model（ΔJ = {E['J_model'] - D['J_model']:.2f}）")
ax1.plot(xs - 0.04, [np.mean(D["J_check"]), np.mean(E["J_check"])], ls="", marker="x", ms=9, mew=2, color=CM,
         label="分析流程檢查：以 target_moment 代替電流")
ax1.set_xticks(xs, ["DIAG", "EXT45"])
ax1.set_xlim(-0.3, 1.3)
ax1.set_ylabel("J [kg·m²]（Kt = 0.135）")
ax1.set_title("(a) 實測 J 與控制器模型：實測較低，兩姿態差較大", fontsize=11.5)
ax1.grid(alpha=0.3)
ax1.legend(loc="lower right", fontsize=9, framealpha=0.95)

order = [0, 3, 1, 2]
names = [R["cond"][i] for i in order]
dJ = np.array([R["dJ"][i] for i in order])
lo = dJ - np.array([R["dJ_rng"][i][0] for i in order]); hi = np.array([R["dJ_rng"][i][1] for i in order]) - dJ
xb = np.arange(len(order))
w = 0.38
ax2.bar(xb - w / 2, dJ, w, yerr=[lo, hi], capsize=4, color="#2366b5", label=f"實測（Kt × 101 = {R['KT']:.2f} N·m/A）")
sc = R["K0"] / R["KT"]
ax2.bar(xb + w / 2, dJ * sc, w, yerr=[lo * sc, hi * sc], capsize=4, color="#9cbbe3",
        label=f"實測，改用控制器換算值 {R['K0']:.2f} N·m/A")
ax2.axhline(R["dJ_model"], color=CM, ls="--", lw=2, label=f"控制器模型 ΔJ_model = {R['dJ_model']:.2f}")
for x, v in zip(xb, dJ):
    ax2.text(x - w / 2, v + 0.08, f"{v / R['dJ_model']:.2f} 倍", ha="center", fontsize=9)
for x, v in zip(xb, dJ * sc):
    ax2.text(x + w / 2, v + 0.08, f"{v / R['dJ_model']:.2f} 倍", ha="center", fontsize=9, color="0.35")
ax2.set_xticks(xb, [f"{n}\nA = {D['A'][i]:.2f}" for n, i in zip(names, order)])
ax2.set_ylabel("ΔJ = J(EXT45) - J(DIAG) [kg·m²]")
ax2.set_ylim(0, 4.3)
ax2.set_title("(b) 兩姿態差 ΔJ：為模型的 1.25–1.40 倍（誤差棒為各次組合範圍）", fontsize=11.5)
ax2.grid(alpha=0.3, axis="y")
ax2.legend(loc="upper left", fontsize=9, framealpha=0.95)
fig.suptitle("圖 2　實測 J 與控制器模型 J_model（target_moment 對 target_qdd 回歸）比較", fontsize=12.5)
fig.tight_layout(rect=(0, 0, 1, 0.95))
fig.savefig(os.path.join(OUT, "fig2_model_compare.png"), dpi=150)
plt.close(fig)
print("圖已存到", OUT)
