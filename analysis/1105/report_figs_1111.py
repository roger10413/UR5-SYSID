"""11/11 報告的圖（11/5 上機，手臂朝上 UP）。
圖 1：J1 定速，正轉與反轉（取絕對值）的窗口平均力矩對速度：重疊表示窗口內重力抵消。
圖 2：J1 速度調變，各條件兩個方向的 J，加控制器 J_model。
另存（報告不放，自己記錄）：J0 手臂朝上的 J，和上週直線 J = c·X + J_d（10/29 反轉 4 次，actual_current）
在 X_UP 的預測並排，加控制器 J_model。
分析用 robot/ 內的分析函式（不執行任何 main，不連網路）；結果存 OUT/results_1105.json（已存在時直接讀取）。
輸出 OUT/figs/w1111_fig*.png、w1111_j0_up_prediction.png。（Microsoft JhengHei 沒有「−」，數值用 ASCII 負號。）"""
import socket
socket.socket = None                                   # 離線分析：不連網路
import os, sys, io, json, glob, math, contextlib
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.path.insert(0, os.path.join(HERE, "..", "1104"))
import ur5_velocity_modulation_ident as v0
import ur5_j1_velocity_modulation_ident as v1
import ur5_j1_const_velocity_ident as c1
import arm_model_1104 as am

sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
DATA = r"D:\cluade_code\ur5\data_1105"
OUT = r"D:\cluade_code\ur5\analysis_1105"
FIGS = os.path.join(OUT, "figs"); os.makedirs(FIGS, exist_ok=True)
INK, INK2, MUTED, GRAY = "#0b0b0b", "#52514e", "#8d8b83", "#b9b8b0"
ORANGE, BLUE = "#d9822b", "#2a78d6"                    # 反轉、正轉（圖 1 的實測為正轉）
PAY = 1.30


def csv_of(pattern):
    return glob.glob(os.path.join(DATA, pattern, "*.csv"))[0]


def quiet(fn, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()):
        return fn(*a, **k)


def compute():
    R = {}
    r = quiet(v0.analyze_velocity_modulation, csv_of("velmod_up_p130_pos_full_*"), 0, v0.KT_OUT, +1,
              v0.CONDITIONS, v0.V_CRUISE)
    R["j0_up"] = dict(J=[r["conditions"][i]["J"] for i in range(6)], se=[r["conditions"][i]["se_J"] for i in range(6)],
                      temp=r["conditions"][0]["temp"], J_model=r["J_model_mean"],
                      J_mean=r["J_mean"], se_J_mean=r["se_J_mean"])
    for d, tag in ((+1, "pos"), (-1, "neg")):
        r = quiet(v1.analyze_j1_velocity_modulation, csv_of(f"j1velmod_up_p130_{tag}_full_*"), 1, v1.KT_OUT, d,
                  v1.CONDITIONS, v1.V_CRUISE)
        cs = [r["conditions"][i] for i in range(6)]
        R[f"j1_velmod_{tag}"] = dict(J=[c["J"] for c in cs], se=[c["se_J"] for c in cs],
                                     J_spread=[c["J_spread"] for c in cs], G_s=[c["G_s"] for c in cs],
                                     G_c=[c["G_c"] for c in cs], gain=[c["track_gain"] for c in cs],
                                     temp=cs[0]["temp"], J_model=r["J_model_mean"])
    a, b = c1.window_bounds_deg(c1.POSES["UP"]["center_deg"], c1.N_MOTOR_REVS_WINDOW, c1.GEAR_RATIO)
    res = {}
    for d, tag in ((+1, "pos"), (-1, "neg")):
        res[tag] = quiet(c1.analyze_angle_window, csv_of(f"j1constvel_up_{tag}_p130_full_*"), 1, c1.KT_OUT, d,
                         c1.SPEED_LEVELS, a, b)
        R[f"j1_constvel_{tag}"] = dict(qd=[s["qd_mean"] for s in res[tag]["segments"]],
                                       tau=[s["tau_mean"] for s in res[tag]["segments"]],
                                       B=res[tag]["B"], C=res[tag]["C"], sd_B=res[tag]["sd_B"], sd_C=res[tag]["sd_C"])
    R["j1_constvel"] = {k: float(x) for k, x in quiet(c1.combine_directions, res["pos"], res["neg"]).items()}
    return R


path = os.path.join(OUT, "results_1105.json")
if os.path.exists(path):
    R = json.load(open(path, encoding="utf-8"))
else:
    R = compute()
    json.dump(R, open(path, "w", encoding="utf-8"), ensure_ascii=False, indent=1)


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", alpha=1, ls=(0, (4, 3))); ax.set_axisbelow(True)


# ---------------- J0 手臂朝上，預測 vs 實測（報告不放，自己記錄） ----------------
X_UP = float(am.mass_matrix(np.radians(am.POSES["UP"]), PAY, rotor=np.zeros(6))[0, 0])
A = json.load(open(r"D:\cluade_code\ur5\analysis_1104\current_signals_1104.json", encoding="utf-8"))["A"]
ac = A["actual_current"]
pred = np.array([r["c"] * (X_UP + r["Jd_over_c"]) for r in ac])
pred_rotor = np.array([r["c"] * (X_UP + 1.915) for r in ac])
meas, se = np.array(R["j0_up"]["J"]), np.array(R["j0_up"]["se"])
jm0 = R["j0_up"]["J_model"]
labels0 = [f"條件 {i + 1}\n{f:g} Hz × {amp:g}" for i, (f, amp) in enumerate(v0.CONDITIONS)]

fig, ax = plt.subplots(figsize=(10, 5.9))
x, w = np.arange(6), 0.36
ax.bar(x - w / 2 - 0.01, pred, w, color=GRAY, edgecolor="white", linewidth=2, zorder=3)
ax.bar(x + w / 2 + 0.01, meas, w, color=BLUE, edgecolor="white", linewidth=2, zorder=3)
ax.errorbar(x + w / 2 + 0.01, meas, yerr=se, fmt="none", ecolor=INK, elinewidth=1.2, capsize=3, zorder=4)
for i in range(6):
    ax.text(x[i] - w / 2 - 0.01, pred[i] + 0.04, f"{pred[i]:.2f}", ha="center", va="bottom", fontsize=10, color=INK2)
    ax.text(x[i] + w / 2 + 0.01, meas[i] + se[i] + 0.04, f"{meas[i]:.2f}", ha="center", va="bottom", fontsize=10,
            color=INK, fontweight="bold")
ax.axhline(jm0, color=INK, lw=1.8, ls=(0, (6, 4)), zorder=2)
ax.text(5.55, jm0 + 0.06, f"控制器 J_model = {jm0:.2f}（連桿與負載 {X_UP:.2f} + 1.915）", ha="right", va="bottom",
        fontsize=11.5, color=INK)
ax.text(-0.45, 1.55, "灰：上週的直線預測（電流不含 1.915）\n藍：11/5 實測", ha="left", va="center", fontsize=11.5,
        color=INK, linespacing=1.6)
ax.set_xticks(x); ax.set_xticklabels(labels0, fontsize=10.5)
ax.set_xlim(-0.6, 5.6); ax.set_ylim(0, 2.45)
ax.set_ylabel("J0 [kg·m²]")
style(ax)
ax.set_title("手臂朝上的 J0：實測和上週預測一致，遠低於控制器的 2.10", loc="left", fontsize=13.5, pad=12, color=INK)
fig.text(0.01, 0.005, f"預測：上週 4 次（DIAG、EXT45 × 夾或不夾，10/29 反轉）的直線 J = c·X + J_d 在 X = {X_UP:.2f} 的值。"
         f"實測：11/5 正轉，每個條件 7 段，誤差棒為單輪回歸的標準差。\n電流若含 1.915，J 會在 "
         f"{pred_rotor.min():.1f}～{pred_rotor.max():.1f}（= c × ({X_UP:.2f} + 1.915)）。",
         ha="left", va="bottom", fontsize=9, color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.06, 1, 1))
f0 = os.path.join(FIGS, "w1111_j0_up_prediction.png")
fig.savefig(f0, dpi=150, bbox_inches="tight", pad_inches=0.15); plt.close(fig)

# ---------------- 圖 1：J1 定速，正反轉重疊 ----------------
cv = R["j1_constvel"]
fig, ax = plt.subplots(figsize=(10, 5.9))
for tag, col, mk, lab in (("pos", BLUE, "o", "正轉"), ("neg", ORANGE, "s", "反轉（取絕對值）")):
    d = R[f"j1_constvel_{tag}"]
    ax.scatter(np.abs(d["qd"]), np.abs(d["tau"]), s=70, marker=mk, color=col, edgecolor="white", linewidth=1.8,
               zorder=4, label=lab)
vv = np.linspace(0, 0.21, 20)
ax.plot(vv, cv["Tc"] + cv["B"] * vv, color=INK, lw=1.6, zorder=3)
ax.text(0.205, 9.35, f"零速截距：正轉 {cv['Tc'] + cv['H']:.2f}、反轉 {cv['Tc'] - cv['H']:.2f}\n"
        f"→ 窗口內重力平均 {cv['H']:+.2f} ± {cv['sd_H']:.2f} N·m", ha="right", va="bottom", fontsize=11.5,
        color=INK, linespacing=1.6)
ax.text(0.004, 12.75, f"B = {cv['B']:.1f} ± {cv['sd_B']:.1f} N·m·s/rad\nTc = {cv['Tc']:.2f} ± {cv['sd_Tc']:.2f} N·m",
        ha="left", va="top", fontsize=12.5, color=INK, fontweight="bold", linespacing=1.6)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.02), frameon=False, fontsize=11.5)
ax.set_xlim(0, 0.215); ax.set_ylim(9.0, 13.6)
ax.set_xlabel("J1 速度 [rad/s]"); ax.set_ylabel("窗口平均力矩（電流 × 13.64）[N·m]")
style(ax); ax.grid(axis="x", alpha=1, ls=(0, (4, 3)))
ax.set_title("圖 1　J1 定速：正、反轉的力矩重疊 → 重力在窗口內抵消，量到的就是摩擦", loc="left", fontsize=13.5,
             pad=12, color=INK)
fig.text(0.01, 0.005, "手臂朝上，窗口 -102.7°～-81.3°（6 個馬達圈，中心 -92° 為重力零點）；"
         "每個方向 6 檔速度先升後降共 12 段，每點為一段的平均。±：單輪回歸的標準差合成。",
         ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 1))
f1 = os.path.join(FIGS, "w1111_fig1_j1_constvel.png")
fig.savefig(f1, dpi=150, bbox_inches="tight", pad_inches=0.15); plt.close(fig)

# ---------------- 圖 2：J1 速度調變 ----------------
order = sorted(range(6), key=lambda i: (v1.CONDITIONS[i][0], v1.CONDITIONS[i][1]))
jp, jn = R["j1_velmod_pos"], R["j1_velmod_neg"]
jm1 = 0.5 * (jp["J_model"] + jn["J_model"])
X1 = float(am.mass_matrix(np.radians(am.POSES["UP"]), PAY, rotor=np.zeros(6))[1, 1])
fig, ax = plt.subplots(figsize=(10, 5.9))
x = np.arange(6)
for k, (d, col, lab, sh) in enumerate(((jp, BLUE, "正轉", -1), (jn, ORANGE, "反轉", +1))):
    J, s = np.array(d["J"])[order], np.array(d["se"])[order]
    xx = x + sh * (w / 2 + 0.01)
    ax.bar(xx, J, w, color=col, edgecolor="white", linewidth=2, zorder=3, label=lab)
    ax.errorbar(xx, J, yerr=s, fmt="none", ecolor=INK, elinewidth=1.2, capsize=3, zorder=4)
    for i in range(6):
        ax.text(xx[i], J[i] - 0.3, f"{J[i]:.2f}", ha="center", va="top", fontsize=10.5, color="white",
                fontweight="bold", zorder=6)
ax.axhline(jm1, color=INK, lw=1.8, ls=(0, (6, 4)), zorder=5)
ax.text(-0.45, jm1 + 0.2, f"控制器 J_model = {jm1:.2f}", ha="left", va="bottom", fontsize=11.5, color=INK)
ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), frameon=False, fontsize=11.5, ncol=2)
ax.set_xticks(x); ax.set_xticklabels([f"{v1.CONDITIONS[i][0]:g} Hz × {v1.CONDITIONS[i][1]:g}" for i in order],
                                     fontsize=10.5)
ax.set_xlim(-0.6, 5.6); ax.set_ylim(0, 9.4)
ax.set_ylabel("J1 [kg·m²]")
style(ax)
ax.set_title("圖 2　J1 速度調變：兩個方向的 J 一致", loc="left", fontsize=13.5, pad=12, color=INK)
fig.text(0.01, 0.005, "手臂朝上，調變段中心在重力零點 -92°，巡航 0.08 rad/s，每個條件 7 段。電流 × 13.64；"
         "分析含重力項。誤差棒為單輪回歸的標準差。", ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 1))
f2 = os.path.join(FIGS, "w1111_fig2_j1_velmod.png")
fig.savefig(f2, dpi=150, bbox_inches="tight", pad_inches=0.15); plt.close(fig)

# ---------------- 數字（講解頁、講稿用） ----------------
print("X_UP J0 =", round(X_UP, 4), " J0 J_model UP =", round(jm0, 4), " J0 溫度", round(R["j0_up"]["temp"], 2))
for i in range(6):
    print(f"J0 UP 條件 {i + 1}: 預測 {pred[i]:.3f}  實測 {meas[i]:.3f} ± {se[i]:.3f}  差 {meas[i] - pred[i]:+.3f}"
          f"  若含 1.915 {pred_rotor[i]:.2f}")
print("J0 UP 加權平均", round(R["j0_up"]["J_mean"], 4), "±", round(R["j0_up"]["se_J_mean"], 4))
print("J1 X_UP（點質量）", round(X1, 4), " 含連桿自身慣量",
      round(float(am.mass_matrix(np.radians(am.POSES["UP"]), PAY, self_inertia=True, rotor=np.zeros(6))[1, 1]), 4),
      " J1 J_model", round(jm1, 4), " 溫度", round(jp["temp"], 2), round(jn["temp"], 2))
q0 = math.radians(-92.0)
for tag, d, sgn in (("pos", jp, +1), ("neg", jn, -1)):
    sl = [sgn * (gs * math.cos(q0) - gc * math.sin(q0)) * math.pi / 180 for gs, gc in zip(d["G_s"], d["G_c"])]
    print(f"J1 速度調變 {tag}: J {[round(j, 3) for j in d['J']]}  ± {[round(s, 3) for s in d['se']]}"
          f"  追隨增益 {[round(g, 3) for g in d['gain']]}  重力斜率 {min(sl):+.3f}～{max(sl):+.3f} N·m/°")
print("J1 定速", {k: round(x, 4) for k, x in cv.items()})
for f in (f1, f2, f0):
    print(f)
