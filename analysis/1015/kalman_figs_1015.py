"""10/15 卡爾曼濾波殘差分析的圖（需先執行 kalman_residual_1015.py）。
K1 擴增狀態濾波器的模型外力矩 d̂（名義三參數模型 vs 擴充模型）
K2 逐步加入各項後的殘差 RMS
K3 馬達角漣波（部分殘差對馬達角、各次諧波振幅）
K4 起步後行程效應（擬合曲線與各運動類型的部分殘差）
K5 J 依運動類型（只有三參數 vs 完整模型）"""
import sys, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 11, "mathtext.fontset": "dejavusans"})
import data_1015 as dd
import kalman_residual_1015 as kr

FIG = os.path.join(kr.OUT, "figs"); os.makedirs(FIG, exist_ok=True)
CB, CR, CG = "#2366b5", "#c8452f", "#888888"
DHAT = r"$\hat{d}$"
Z = np.load(os.path.join(kr.OUT, "kf_results.npz"))
names = list(Z["names"]); P = dict(zip(names, Z["p"])); h = str(Z["hash"])

runs = dd.run_list(); kn = kr.prepare(runs)
nz_ext = kr.tuned_noise(runs, None, f"ext_{h}", kn)
kf0 = kr.run_kf(runs, kr.NOMINAL, nz_ext, f"nomX_{h}")
kf1 = kr.run_kf(runs, kr.ext_model(P), nz_ext, f"ext_{h}")
cols, y, info, G = kr.build(runs)
X = np.column_stack([cols[n] for n in names]); yhat = X @ Z["p"]


def part(group):
    """部分殘差：實測 − 完整模型 + 該組的擬合值。"""
    idx = [names.index(n) for n in G[group]]
    return y - yhat + X[:, idx] @ Z["p"][idx]


# ---------- K1：濾波器的模型外力矩 ----------
def seg_window(Fe, k_main, pre, post):
    """第 k_main 個主運動段，從起步前 pre 秒到起步後 post 秒。"""
    segs = [g for g in np.unique(Fe["seg"][Fe["main"]]) if g >= 0]
    i0 = np.flatnonzero(Fe["seg"] == segs[k_main])[0]
    return slice(max(i0 - int(pre / dd.DT), 0), min(i0 + int(post / dd.DT), len(Fe["seg"])))


fig, axs = plt.subplots(2, 2, figsize=(15, 7.8), sharex="col", gridspec_kw=dict(height_ratios=[1, 1.7]))
ex = [(next(x for x in runs if x["kind"] == "tri" and x["dir"] > 0 and x["n"] == 1), 45, 1.8, "三角形 v3 第1次 正轉：a = 0.9 rad/s² 的一個循環與 movej 回程"),
      (next(x for x in runs if x["kind"] == "cv" and x["dir"] > 0 and x["n"] == 2), 0, 8.0, "定速第2輪 正轉：0.02 rad/s 段的前 8 s")]
for c, (r, km, post, title) in enumerate(ex):
    i = runs.index(r); D, Fe = r["D"], r["Fe"]; s = Fe["sgn"]
    sl = seg_window(Fe, km, 0.3, post)
    t = D["t"][sl] - D["t"][sl][0]; o0, o1 = kf0[i], kf1[i]; s_ = s[sl]; use = r["keep"][sl]
    axs[0, c].plot(t, D["tq"][sl], color=CG, lw=2.2, label="指令速度")
    axs[0, c].plot(t, r["v"][sl], color=CB, lw=1, label="估計速度（只用編碼器）")
    axs[0, c].set_ylabel(r"$\dot{\theta}$ [rad/s]"); axs[0, c].set_title(title); axs[0, c].legend(fontsize=9); axs[0, c].grid(alpha=0.2)
    Bn = np.where(s_ >= 0, kr.NOMINAL["B"][+1], kr.NOMINAL["B"][-1])
    Tcn = np.where(s_ > 0, kr.NOMINAL["Tc"][+1], np.where(s_ < 0, -kr.NOMINAL["Tc"][-1], 0))
    raw = D["tau"][sl] - kr.NOMINAL["J"] * o0["a"][sl] - Bn * o0["v"][sl] - Tcn
    axs[1, c].plot(t[use], (s_ * raw)[use], ".", color=CG, ms=2.5, alpha=0.6, label="實測力矩 - 名義三參數模型（逐樣本）")
    axs[1, c].plot(t, np.where(s_ != 0, s_ * o0["d"][sl], np.nan), color=CB, lw=1.5, label=f"{DHAT}：名義三參數模型")
    axs[1, c].plot(t, np.where(s_ != 0, s_ * o1["d"][sl], np.nan), color=CR, lw=1.5, label=f"{DHAT}：擴充模型（漣波、溫度、速度非線性、行程、加減速階段）")
    axs[1, c].axhline(0, color="k", lw=0.6); axs[1, c].grid(alpha=0.2); axs[1, c].set_xlabel("時間 [s]")
    axs[1, c].set_ylabel("運動方向的模型外力矩 [N·m]"); axs[1, c].set_ylim(-3.5, 3.5)
axs[1, 0].legend(fontsize=8.5, loc="lower left")
fig.suptitle(f"擴增狀態卡爾曼濾波＋RTS 平滑：模型外力矩 {DHAT}（以運動方向為正；灰點只畫模型適用、有使用電流量測的樣本）", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "figK1_kf_residual_example.png"), dpi=140); plt.close(fig)

# ---------- K2：逐步加入各項 ----------
steps = Z["steps"]; lab = ["三參數\nJ、B±、Tc±", "+ 馬達角漣波\n1～24 次諧波", "+ 溫度\nB、Tc 隨溫度", "+ 速度非線性\nC·v²", "+ 起步後行程", "+ 加減速階段\n常數偏移"]
fig, ax = plt.subplots(figsize=(12, 4.8)); xx = np.arange(len(lab)); w = 0.38
b1 = ax.bar(xx - w / 2, steps[:, 0], w, color=CB, label="殘差 RMS（逐樣本）")
b2 = ax.bar(xx + w / 2, steps[:, 1], w, color=CR, label="慢變殘差 RMS（同段 0.2 s 移動平均後）")
for b in list(b1) + list(b2):
    ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.015, f"{b.get_height():.3f}", ha="center", fontsize=9)
ax.set_xticks(xx); ax.set_xticklabels(lab, fontsize=10); ax.set_ylabel("N·m"); ax.grid(alpha=0.2, axis="y"); ax.legend()
ax.set_title(f"依序加入各項後的殘差（10/15 全部 12 次實驗，模型適用的運動樣本 {len(y)} 個、{len(np.unique(info['seg']))} 段）")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "figK2_nested_rms.png"), dpi=140); plt.close(fig)

# ---------- K3：馬達角漣波 ----------
pr = part("漣波"); ph = np.degrees(np.mod(info["phi"], 2 * np.pi)); edges = np.linspace(0, 360, 181); cen = (edges[1:] + edges[:-1]) / 2
fig, axs = plt.subplots(1, 2, figsize=(15, 5.4), gridspec_kw=dict(width_ratios=[1.6, 1]))
amps = {}
for sg, c, lb in (("+", CB, "正轉"), ("-", CR, "反轉")):
    m = info["s"] > 0 if sg == "+" else info["s"] < 0
    bi = np.clip(np.digitize(ph[m], edges) - 1, 0, 179)
    mean = np.bincount(bi, pr[m], 180) / np.maximum(np.bincount(bi, minlength=180), 1)
    pp = np.radians(cen); curve = sum(P[f"c{k}{sg}"] * np.cos(k * pp) + P[f"s{k}{sg}"] * np.sin(k * pp) for k in range(1, kr.K_RIPPLE + 1))
    axs[0].plot(cen, mean, ".", color=c, ms=4, alpha=0.7, label=f"{lb}：部分殘差（每 2° 平均）")
    axs[0].plot(cen, curve, "-", color=c, lw=1.4, label=f"{lb}：24 次諧波擬合（峰對峰 {np.ptp(curve):.2f} N·m）")
    amps[sg] = np.array([np.hypot(P[f"c{k}{sg}"], P[f"s{k}{sg}"]) for k in range(1, kr.K_RIPPLE + 1)])
axs[0].set_xlabel("馬達角（101 × 關節角，取一圈）[°]"); axs[0].set_ylabel("運動方向的力矩 [N·m]"); axs[0].grid(alpha=0.2)
axs[0].legend(fontsize=8.5, ncol=2, loc="upper center", bbox_to_anchor=(0.5, -0.13))
axs[0].set_title("漣波與馬達角鎖定（每個馬達圈 = 關節 3.564°）")
kk = np.arange(1, kr.K_RIPPLE + 1)
axs[1].bar(kk - 0.2, amps["+"], 0.4, color=CB, label="正轉"); axs[1].bar(kk + 0.2, amps["-"], 0.4, color=CR, label="反轉")
axs[1].set_xlabel("諧波次數（每個馬達圈的週期數）"); axs[1].set_ylabel("振幅 [N·m]"); axs[1].set_xticks(kk[::2])
axs[1].grid(alpha=0.2, axis="y"); axs[1].legend(); axs[1].set_title("各次諧波振幅")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "figK3_ripple_motor_angle.png"), dpi=140); plt.close(fig)

# ---------- K4：起步後行程效應（以行程 10° 為 0） ----------
pr = part("起步後行程"); tr = info["travel"]
gi = [names.index(n) for n in G["起步後行程"]]
ref = (kr.hat_basis(np.array([10.0]), kr.KNOTS) @ Z["p"][gi])[0]
xs = np.linspace(kr.KNOTS[1], 50, 496); gcurve = kr.hat_basis(xs, kr.KNOTS) @ Z["p"][gi] - ref   # 0～0.5° 樣本很少，曲線從 0.5° 畫起
grp = [("三角形 v3（主運動）", (info["kind"] == "tri") & info["main"], "#3a9e5c", "o"),
       ("梯形 v3（主運動）", (info["kind"] == "trap") & info["main"], "#b5397a", "s"),
       ("定速（主運動）", (info["kind"] == "cv") & info["main"], CB, "^"),
       ("movej 回程", ~info["main"], CR, "v")]
fig, axs = plt.subplots(1, 2, figsize=(15, 5), sharey=True)
ed = np.r_[np.arange(0, 4, 0.5), np.arange(4, 10, 1), np.arange(10, 52, 2.5)]; nb = len(ed) - 1; cen = (ed[1:] + ed[:-1]) / 2
for ax, sg, lb in ((axs[0], 1, "正轉"), (axs[1], -1, "反轉")):
    for nm, m, c, mk in grp:
        mm = m & (info["s"] == sg); bi = np.digitize(tr[mm], ed) - 1; ok_ = (bi >= 0) & (bi < nb)
        cnt = np.bincount(bi[ok_], minlength=nb); mean = np.bincount(bi[ok_], pr[mm][ok_], nb) / np.maximum(cnt, 1)
        ok = cnt > 100
        ax.plot(cen[ok], mean[ok] - ref, mk, color=c, ms=5, alpha=0.85, label=nm)
    ax.plot(xs, gcurve, "k-", lw=1.6, label="擬合 g（各類型、正反轉共用）")
    ax.axhline(0, color="k", lw=0.6); ax.grid(alpha=0.2); ax.set_title(lb); ax.set_xlabel("本段起步後行程 [°]")
axs[0].set_ylabel("部分殘差 [N·m]（以行程 10° 為 0）"); axs[1].legend(fontsize=9, loc="lower right")
fig.suptitle("起步後行程效應：摩擦相對於行程 10° 時的差（每箱 > 100 個樣本才畫）", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "figK4_travel_effect.png"), dpi=140); plt.close(fig)

# ---------- K5：J 依運動類型 ----------
js = Z["jsplit"]; jl = ["三角形 v3", "梯形 v3", "定速段（起步、停止）", "movej 回程"]
fig, ax = plt.subplots(figsize=(10, 4.6)); xx = np.arange(4); w = 0.38
ax.bar(xx - w / 2, js[:, 2], w, color=CG, label="只有三參數（漣波等都不扣）")
ax.bar(xx + w / 2, js[:, 0], w, yerr=js[:, 1], capsize=4, color=CB, label="完整模型（± 以運動段重抽樣的標準差）")
for k in range(4):
    ax.text(xx[k] - w / 2, js[k, 2] + 0.03, f"{js[k, 2]:.2f}", ha="center", fontsize=9)
    ax.text(xx[k] + w / 2, js[k, 0] + js[k, 1] + 0.03, f"{js[k, 0]:.2f}", ha="center", fontsize=9)
ax.axhline(P["J"], color=CB, ls="--", lw=1, label=f"完整模型共用一個 J = {P['J']:.2f}")
ax.set_xticks(xx); ax.set_xticklabels(jl); ax.set_ylabel("J [kg·m²]"); ax.grid(alpha=0.2, axis="y"); ax.legend(fontsize=9, loc="lower right")
ax.set_title("各運動類型分開估 J（加速度用指令值）")
fig.tight_layout(); fig.savefig(os.path.join(FIG, "figK5_J_by_motion.png"), dpi=140); plt.close(fig)
print("完成：", FIG)
