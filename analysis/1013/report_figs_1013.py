"""10/14 報告圖：定速 3 輪（回歸、B/Tc、力矩對角度）與定加速度 v3（J 回歸、殘差對行程）。"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot")); sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 11})
import ur5_const_velocity_angle_ident as mv
import accel_v3_analysis as av

BASE = av.BASE                    # 資料資料夾，於 accel_v3_analysis.py 設定
OUT =os.path.join(BASE, "analysis_1013", "figs"); os.makedirs(OUT, exist_ok=True)
CB, CR = "#2366b5", "#c8452f"
REV = 360.0 / 101
rows = list(np.load(os.path.join(BASE, "analysis_1013", "constvel_rounds.npy"), allow_pickle=True))
files = sorted(glob.glob(os.path.join(BASE, "constvel_angle_*_full_20261013_*", "constvel_angle_data_*.csv")), key=lambda p: p[-10:-4])
labels = ["第1輪 正轉", "第1輪 反轉", "第2輪 正轉", "第2輪 反轉", "第3輪 正轉", "第3輪 反轉"]
mk = ["o", "o", "s", "s", "^", "^"]

# ---------- 圖 1：定速回歸（兩象限過原點） ----------
fig, ax = plt.subplots(figsize=(8.5, 6))
for k, x in enumerate(rows):
    s = x["dir"]; c = CB if s > 0 else CR
    v = np.array([0.02, 0.05, 0.08, 0.12, 0.16, 0.20, 0.20, 0.16, 0.12, 0.08, 0.05, 0.02])
    ax.plot(s * v, s * x["levels"], mk[k], color=c, mfc="none" if k < 2 else c, alpha=0.85, ms=6, label=labels[k])
    vv = np.linspace(0, 0.21, 10); ax.plot(s * vv, s * (x["B"] * vv + x["Tc"]), "-", color=c, lw=0.8, alpha=0.6)
ax.axhline(0, color="k", lw=0.6); ax.axvline(0, color="k", lw=0.6)
P = [x for x in rows if x["dir"] > 0]; N = [x for x in rows if x["dir"] < 0]
ax.text(0.03, -6, f"正轉 3 輪\nB = {np.mean([x['B'] for x in P]):.2f}（輪間 std {np.std([x['B'] for x in P], ddof=1):.2f}）\n"
        f"Tc = {np.mean([x['Tc'] for x in P]):.3f}（std {np.std([x['Tc'] for x in P], ddof=1):.3f}）", color=CB, fontsize=10)
ax.text(-0.21, 5, f"反轉 3 輪\nB = {np.mean([x['B'] for x in N]):.2f}（輪間 std {np.std([x['B'] for x in N], ddof=1):.2f}）\n"
        f"Tc = {np.mean([x['Tc'] for x in N]):.3f}（std {np.std([x['Tc'] for x in N], ddof=1):.3f}）", color=CR, fontsize=10)
ax.set_xlabel("J0 角速度 [rad/s]"); ax.set_ylabel("力矩 τ = i·Kt·N [N·m]（窗口內 10 個馬達圈平均）")
ax.set_title("固定角度窗口定速：3 輪 × 正反轉，每輪 12 段（0.02→0.20→0.02）"); ax.grid(alpha=0.2); ax.legend(fontsize=9, ncol=2, loc="lower right")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_constvel_regression.png"), dpi=140); plt.close(fig)

# ---------- 圖 2：各輪 B、Tc ----------
fig, axs = plt.subplots(1, 2, figsize=(12, 4.6))
for ax, key, se, yl in ((axs[0], "B", "seB", "B [N·m·s/rad]"), (axs[1], "Tc", "seTc", "Tc [N·m]")):
    for k, x in enumerate(rows):
        c = CB if x["dir"] > 0 else CR
        ax.errorbar(k + 1, x[key], x[se], fmt="o", color=c, mfc=c, ms=8, capsize=5, lw=1.5)
        ax.plot(k + 1 + 0.18, x["up" + key], "^", color=c, ms=5, alpha=0.5)
        ax.plot(k + 1 + 0.30, x["dn" + key], "v", color=c, ms=5, alpha=0.5)
    ax.set_xticks(range(1, 7)); ax.set_xticklabels(labels, rotation=20); ax.set_ylabel(yl); ax.grid(alpha=0.2)
    ax.set_title(f"{key}：實心 = 12 段回歸（± 以段為點），△ 升速趟、▽ 降速趟")
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_constvel_B_Tc.png"), dpi=140); plt.close(fig)

# ---------- 圖 3：力矩對角度（1 個馬達圈移動平均、扣各段平均） ----------
fig, ax = plt.subplots(figsize=(9, 5))
a_, b_ = mv.window_bounds_deg(180.0, 10, 101)
grid = np.linspace(a_ + REV / 2, b_ - REV / 2, 120)
for k, (F, x) in enumerate(zip(files, rows)):
    s = x["dir"]; d = np.genfromtxt(F, delimiter=",", names=True, usecols=(0, 1, 2, 3, 4))
    q = np.degrees(d["actual_q_0"]); tau = d["actual_current_0"] * mv.KT_OUT
    import io, contextlib
    with contextlib.redirect_stdout(io.StringIO()):
        r = mv.analyze_angle_window(F, 0, mv.KT_OUT, s, mv.SPEED_LEVELS, a_, b_)
    curves = []
    for seg in r["segments"]:
        i = seg["idx"]; qq = q[i]; tt = tau[i]; o = np.argsort(qq); qq, tt = qq[o], tt[o]
        cs = np.concatenate([[0], np.cumsum(tt)]); lo = np.searchsorted(qq, grid - REV / 2); hi = np.searchsorted(qq, grid + REV / 2)
        ma = (cs[hi] - cs[lo]) / np.maximum(hi - lo, 1)
        curves.append(ma - np.mean(tt))
    ax.plot(grid, np.mean(curves, 0), "-", color=CB if s > 0 else CR, lw=1.6, alpha=[1, 1, 0.75, 0.75, 0.5, 0.5][k], label=labels[k])
ax.axhline(0, color="k", lw=0.6)
ax.set_xlabel("J0 角度 [°]"); ax.set_ylabel("τ - 該段平均 [N·m]（帶正負號）")
ax.set_title("力矩隨角度變化（每段先做 1 個馬達圈移動平均消去漣波，再平均 12 段）\n正轉 τ>0、反轉 τ<0，兩方向帶號力矩都隨角度增加")
ax.legend(fontsize=9, ncol=3); ax.grid(alpha=0.2)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_torque_vs_angle.png"), dpi=140); plt.close(fig)

# ---------- 圖 4：v3 J ----------
acc = [av.run(*r) for r in av.RUNS]
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
for o in acc:
    c = CB if o["s"] > 0 else CR; p = o["res"][0.04]["p"]; L = o["L"]
    a = np.array([x["a"] for x in L]); h = np.array([x["h"] for x in L])
    axs[0].plot(a, h, "o", color=c, ms=8, label=f"{o['name']}：J = {p[0]:.2f} ± {o['se'][0]:.2f}，截距 {p[1]:+.2f} N·m")
    aa = np.linspace(0, 1.3, 20); axs[0].plot(aa, np.polyval(p, aa), "-", color=c, lw=1.2)
    axs[0].plot(aa, p[0] * aa, ":", color=c, lw=1)
    axs[1].plot(a, o["j_levels"], "o", color=c, ms=8, label=f"{o['name']}：標準算法（假設加減速摩擦相同）")
    aa = np.linspace(0.25, 1.3, 50); axs[1].plot(aa, p[0] + p[1] / aa, "-", color=c, lw=1.2)
    axs[1].axhline(p[0], color=c, ls=":", lw=1)
axs[0].axhline(0, color="k", lw=0.6); axs[0].set_xlim(0, 1.3)
axs[0].set_xlabel("加速度 a [rad/s²]"); axs[0].set_ylabel("(τ_acc - τ_dec)/2 [N·m]")
axs[0].set_title("各檔平均：斜率 = J，截距 = 加減速段摩擦差的一半（虛線 = 過原點）"); axs[0].legend(fontsize=9); axs[0].grid(alpha=0.2)
axs[1].set_xlabel("加速度 a [rad/s²]"); axs[1].set_ylabel("J [kg·m²]")
axs[1].set_title("標準算法各檔 J；實線 = J + 截距/a，虛線 = 含截距回歸的 J"); axs[1].legend(fontsize=9, loc="lower right"); axs[1].grid(alpha=0.2)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_v3_J.png"), dpi=140); plt.close(fig)

# ---------- 圖 5：殘差對起步後行程 ----------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.6), sharey=True)
for ax, o in zip(axs, acc):
    s = o["s"]; J = o["res"][0.04]["p"][0]
    d = np.genfromtxt(o["F"], delimiter=",", names=True, usecols=(0, 1, 2, 3, 4))
    t = d["timestamp"]; q = np.degrees(d["actual_q_0"]); qd = s * d["actual_qd_0"]; tq = s * d["target_qd_0"]; tau = s * d["actual_current_0"] * av.m.KT_OUT
    dtq = np.gradient(tq, t); mvg = tq > 0.05
    ph = {"加速段": mvg & (dtq > 0.1), "定速段": mvg & (np.abs(tq - 0.25) < 1e-4) & (np.abs(dtq) < 0.05), "減速段": mvg & (dtq < -0.1)}
    for k in ph:
        msk = ph[k].copy()
        for i in np.where(np.diff(msk.astype(int)) == 1)[0] + 1: msk[i:i + 5] = False
        ph[k] = msk
    res = tau - av.fcv(qd, s) - J * np.where(ph["定速段"], 0, dtq)
    trav = np.abs(q - 180.0); edges = np.arange(0, 31, 2)
    for (k, msk), col in zip(ph.items(), ["#3a9e5c", "#888888", "#b5397a"]):
        xs, ys = [], []
        for lo in edges[:-1]:
            sm = msk & (trav >= lo) & (trav < lo + 2)
            if sm.sum() > 40: xs.append(lo + 1); ys.append(np.mean(res[sm]))
        ax.plot(xs, ys, "o-", color=col, label=k, ms=5)
    ax.axhline(0, color="k", lw=0.6); ax.grid(alpha=0.2); ax.set_xlabel("起步後行程 |q - 180°| [°]")
    ax.set_title(f"{o['name']}（J = {J:.2f}）"); ax.legend(fontsize=9)
axs[0].set_ylabel("|τ| - 同日定速摩擦 - J·a [N·m]")
fig.suptitle("定加速度 v3 殘差：同一行程下，減速段仍比定速段高約 0.5 N·m", fontsize=12)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_v3_residual_travel.png"), dpi=140); plt.close(fig)
print("figs ->", OUT)
