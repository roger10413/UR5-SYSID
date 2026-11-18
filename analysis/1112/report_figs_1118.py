"""11/18 報告的圖（11/12 上機，手臂朝上 UP）。
圖 1：J1 速度調變，J 對頻率：11/5 與 11/12 重複（1.30 kg），加砝碼（2.31 kg）。
圖 2：砝碼造成的 ΔJ 對頻率：J1（加單一彈簧模型）、J2，和控制器 ΔJ_model 比。
圖 3：換算比例（N·m/A）的各種估計，和文獻值 13.64（減速比 101 × 轉矩常數 0.135）比。
圖 4：敲擊測試的餘振頻譜（11/12 已畫好的 w1112_tap_spectrum.png，複製改名）。
圖 5：J2 定速（B、Tc）與 J2 速度調變各條件的 J。
只讀 11/12 分析已存的結果（analysis_1112/*.json，由同資料夾的 j1_1112*.py、j2_*.py、tap_*.py 產生），不連網路。
輸出 OUT/figs/w1118_fig*.png，並複製到 docs/report_1118/。（Microsoft JhengHei 沒有「−」，數值用 ASCII 負號。）"""
import socket
socket.socket = None                                   # 離線分析：不連網路
import os, re, sys, json, math, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r"D:\cluade_code\ur5\analysis_1112"
FIGS = os.path.join(OUT, "figs"); os.makedirs(FIGS, exist_ok=True)
DOCS = os.path.join(HERE, "..", "..", "docs", "report_1118"); os.makedirs(DOCS, exist_ok=True)
INK, INK2, MUTED, GRAY = "#0b0b0b", "#52514e", "#8d8b83", "#b9b8b0"
ORANGE, BLUE, GREEN = "#d9822b", "#2a78d6", "#3a9a5b"
KT = 13.64                                             # 文獻值 101 × 0.135（Raviola 2021，無出處）；analysis_1112/*.json 都用它換算
KT_CB = 12.63                                          # 文獻值 101 × 0.125（Kovincic 2019、Porcelli 2020：CB 系列文件值）
M_CUP, G = 1.008, 9.81

J1 = json.load(open(os.path.join(OUT, "j1_1112.json"), encoding="utf-8"))
J1B = json.load(open(os.path.join(OUT, "j1_1112_b.json"), encoding="utf-8"))
J2 = json.load(open(os.path.join(OUT, "j2_1112.json"), encoding="utf-8"))
J2C = json.load(open(os.path.join(OUT, "j2_cup_1112.json"), encoding="utf-8"))

# 11/18 起改用各軸砝碼校正的換算比例：J、ΔJ、B、Tc 都和換算比例成正比，JSON（× 13.64）的值乘 K / 13.64。
# 控制器的 J_model、ΔJ_model 是控制器自己的單位（N·m），不換算。
mgh1 = M_CUP * G * J1B["h"]                            # 砝碼離 J1 軸的高度 h（重心由控制器 ΔJ_model 反推，d ≈ 56 mm）
dS1 = float(np.mean(J1B["dS"]))
K1 = KT * mgh1 / abs(dS1)                              # J1 砝碼重力斜率（與 KT 無關：dS1 也是 × KT 算的）
K2 = float(J2C["Kcal"])                                # J2 砝碼重力斜率
S1, S2 = K1 / KT, K2 / KT


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", alpha=1, ls=(0, (4, 3))); ax.set_axisbelow(True)


def save(fig, name):
    p = os.path.join(FIGS, name)
    fig.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0.15); plt.close(fig)
    shutil.copy(p, os.path.join(DOCS, name))
    return p


def cond_mean(run, ci):
    """同一條件兩個方向的加權平均 J。"""
    v = [J1[f"velmod_{run}_{t}"]["per"][str(ci)] for t in ("pos", "neg")]
    w = np.array([1 / x["se"] ** 2 for x in v])
    J, se = float(sum(x["J"] * wi for x, wi in zip(v, w)) / w.sum()), float(1 / math.sqrt(w.sum()))
    return S1 * J, S1 * se, v[0]["f"], v[0]["A"]


# ---------------- 圖 1：J1 速度調變，J 對頻率（兩種振幅各畫一張，含誤差棒） ----------------
runs = (("1105_p130", MUTED, "11/5　1.30 kg", "^"), ("1112_p130", BLUE, "11/12　1.30 kg", "o"),
        ("1112_p231", ORANGE, "11/12　2.31 kg（加砝碼 1.008 kg）", "s"))
AMP1 = ((True, "大振幅 0.04 rad/s"), (False, "小振幅 0.02 rad/s（3 Hz 為 0.027）"))
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.8), sharey=True)
pts_all = {run: [cond_mean(run, ci) for ci in range(6)] for run, *_ in runs}
fmean = {run: [np.mean([J for J, s, f, A in pts if abs(f - ff) < 0.2]) for ff in (1, 2, 3)]
         for run, pts in pts_all.items()}
jm130 = np.mean([J1["velmod_1112_p130_pos"]["per"][str(c)]["Jm"] for c in range(6)])
jm231 = np.mean([J1["velmod_1112_p231_pos"]["per"][str(c)]["Jm"] for c in range(6)])
for ax, (big, ttl) in zip(axs, AMP1):
    y3 = {}
    for run, col, lab, mk in runs:
        pp = sorted((f, J, s) for J, s, f, A in pts_all[run] if (A >= 0.04) == big)
        ax.errorbar([p[0] for p in pp], [p[1] for p in pp], yerr=[p[2] for p in pp], fmt="-" + mk, color=col, lw=1.8,
                    ms=8 if mk != "^" else 9, mfc=col, mec="white", mew=1.2, ecolor=col, elinewidth=1.3, capsize=4,
                    capthick=1.3, zorder=4, label=lab)
        y3[run] = pp[-1][1]
    for jm, col in ((jm130, BLUE), (jm231, ORANGE)):
        ax.axhline(jm, color=col, lw=1.2, ls=(0, (6, 4)), zorder=2)
    g0, g1 = y3["1105_p130"], y3["1112_p130"]
    # 箭頭放在 3 Hz 點的右邊，數字標在箭頭右側，不壓到任何線
    ax.annotate("", xy=(3.1, g1), xytext=(3.1, g0), arrowprops=dict(arrowstyle="<->", color=INK2, lw=1.1))
    ax.text(3.14, 0.5 * (g0 + g1), f"兩週\n{g1 - g0:+.2f}", ha="left", va="center", fontsize=10.5, color=INK,
            linespacing=1.3, zorder=6)
    ax.set_title(ttl, loc="left", fontsize=12.5, color=INK)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["1 Hz", "2 Hz", "3 Hz"])
    ax.set_xlim(0.7, 3.45); ax.set_ylim(4.5, 10.2); ax.set_xlabel("調變頻率")
    style(ax)
axs[0].set_ylabel(f"J1 [kg·m²]（電流 × {K1:.2f}，砝碼校正）")
axs[0].legend(loc="upper left", frameon=False, fontsize=10.5)
for jm, col, txt in ((jm130, BLUE, "1.30 kg"), (jm231, ORANGE, "2.31 kg")):
    axs[1].text(3.5, jm, f"控制器 J_model {jm:.2f}\n（{txt}）", ha="left", va="center", fontsize=10, color=col,
                linespacing=1.3)
fig.suptitle("圖 1　J1 速度調變：J 隨頻率上升（兩週都是），加 1 kg 後上升得更快", x=0.01, ha="left", fontsize=13.5,
             color=INK)
fig.text(0.01, 0.005, "手臂朝上，調變段中心在重力零點 -92°，巡航 0.08 rad/s。每點為一個條件正反轉的加權平均；"
         f"左右兩張是兩種振幅，同一條件只出現一次；誤差棒為單輪回歸的標準差。電流 × {K1:.2f}（J1 砝碼校正）。\n"
         "虛線為控制器的 J_model；箭頭為 1.30 kg 在 3 Hz 的兩週差。11/12 比 11/5 低約 0.8 °C（29.9 vs 30.7 °C）。", ha="left", va="bottom", fontsize=9,
         color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.06, 0.93, 0.95))
f1 = save(fig, "w1118_fig1_j1_velmod_freq.png")

# ---------------- 圖 2：砝碼的 ΔJ 對頻率（列：J1、J2；欄：大、小振幅） ----------------
# J1 單一彈簧模型（j1_1112_fit.py 的「整個慣量在彈簧後」）：J_meas = s·J/(1 − J·ω²/k)，砝碼加在 J 上
S_FIT, J_FIT, K_FIT, DJ_FIT = 1.095, 5.5467, 6897.9, 0.905   # j1_1112_fit.py（× 13.64）的擬合值
S_FIT1 = S_FIT * S1                                    # 換成 J1 砝碼校正值後的 s
def dj_spring(f):
    w2 = (2 * math.pi * f) ** 2
    app = lambda J: S_FIT1 * J / (1 - J * w2 / K_FIT)
    return app(J_FIT + DJ_FIT) - app(J_FIT)

dj1 = J1B["dJ"]; jmd1 = J1B["Jm_d"]
j1 = {True: [], False: []}
for lab, (m, s, jm) in dj1.items():
    f = float(lab.split()[0]); A = float(lab.split("×")[1])
    j1[A >= 0.04].append((f, S1 * m, S1 * s))
j2 = {True: [], False: []}
for r in J2C["rows"]:
    w = [1 / s ** 2 for s in r["sdJ"]]
    j2[r["A"] >= 0.05].append((r["f"], S2 * sum(d * wi for d, wi in zip(r["dJ"], w)) / sum(w), S2 / math.sqrt(sum(w))))
jmd2 = float(np.mean([r["dJm"] for r in J2C["rows"]]))
fig, axs = plt.subplots(2, 2, figsize=(12.5, 8.8), sharex=True)
rows = (("J1（肩）", j1, BLUE, jmd1, (0, 2.8), {True: "大振幅 0.04 rad/s", False: "小振幅 0.02 rad/s（3 Hz 為 0.027）"}),
        ("J2（肘）", j2, GREEN, jmd2, (0, 0.8), {True: "大振幅 0.06 rad/s", False: "小振幅 0.03 rad/s（3 Hz 為 0.04）"}))
ff = np.linspace(0.6, 3.3, 60)
for i, (jn, dd, col, jmd, yl, ttl) in enumerate(rows):
    for k, big in enumerate((True, False)):
        ax = axs[i, k]
        pp = sorted(dd[big])
        ax.errorbar([p[0] for p in pp], [p[1] for p in pp], yerr=[p[2] for p in pp], fmt="-o", color=col, lw=1.8, ms=8,
                    mfc=col, mec="white", mew=1.2, ecolor=col, elinewidth=1.3, capsize=4, capthick=1.3, zorder=4)
        dy = 0.035 * (yl[1] - yl[0])
        for n, (f, m, s) in enumerate(pp):              # 數值標在點的右邊，避開誤差棒和往右的折線
            if n == len(pp) - 1:
                y, va = m, "center"
            elif pp[n + 1][1] >= m:                    # 折線往上 → 標在線下方
                y, va = m - dy, "top"
            else:                                       # 折線往下 → 標在線上方
                y, va = m + dy, "bottom"
            ax.text(f + 0.07, y, f"{m:.2f}", ha="left", va=va, fontsize=10.5, color=col,
                    zorder=6, bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
        ax.axhline(jmd, color=INK, lw=1.4, ls=(0, (6, 4)), zorder=2)
        ax.text(3.3, jmd - 0.02 * (yl[1] - yl[0]), f"控制器 ΔJ_model {jmd:.2f}", ha="right", va="top", fontsize=10,
                color=INK)
        if i == 0:
            ax.plot(ff, [dj_spring(x) for x in ff], color=INK2, lw=1.3, ls=(0, (2, 2)), zorder=3)
            if k == 0:
                ax.text(1.95, dj_spring(1.95) + 0.3, "單一彈簧模型\n（固有頻率約 5.6 Hz）", ha="right", va="bottom",
                        fontsize=9.5, color=INK2, linespacing=1.3)
        ax.set_title(f"{jn}　{ttl[big]}", loc="left", fontsize=12, color=INK)
        ax.set_xlim(0.7, 3.4); ax.set_ylim(*yl)
        ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["1 Hz", "2 Hz", "3 Hz"])
        if k == 0:
            ax.set_ylabel("ΔJ [kg·m²]")
        style(ax)
fig.suptitle("圖 2　同一個 1 kg 砝碼：頻率越高，量到的 ΔJ 越大（J1 只有 1 Hz 接近控制器的值）", x=0.01, ha="left",
             fontsize=13.5, color=INK)
fig.text(0.01, 0.005, "ΔJ = J(2.31 kg) - J(1.30 kg)，手臂朝上，同一天、同一方向相減後兩個方向加權平均。左右兩欄是兩種振幅，"
         f"誤差棒為標準差。電流 × 砝碼校正值（J1 {K1:.2f}、J2 {K2:.2f}）。\n控制器 ΔJ_model：target_moment 對 target_qdd 回歸，2.31 kg 減 1.30 kg。"
         f"J1 的模型：J_量 = s·J / (1 - J·ω²/k)，12 個條件擬合（s = {S_FIT1:.3f}、J = 5.55、k = 6898 N·m/rad；χ² = 42，自由度 9）。",
         ha="left", va="bottom", fontsize=9, color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.05, 1, 0.96))
f2 = save(fig, "w1118_fig2_cup_dJ.png")

# ---------------- 圖 3：換算比例的各種估計 ----------------
c1d = J1["constvel_diag"]; c2 = J2["constvel"]
est = [  # (標籤, 值, 顏色, 說明)
    ("J1 砝碼重力斜率", K1, BLUE, f"m·g·h {mgh1:.2f} ÷ 量測 Δ斜率 {abs(dS1):.2f}"),
    ("J2 砝碼重力斜率", J2C["Kcal"], GREEN, f"m·g·h {J2C['mgh']:.2f} ÷ 量測 Δ斜率 5.43"),
    ("J1 彈簧模型擬合", KT / S_FIT, BLUE, "圖 2 左的模型參數 s"),
    ("J1 定速 DIAG 斜率", c1d["K_slope"], GRAY, "力矩對角度斜率，假設控制器重力模型正確"),
    ("J2 定速 UP 斜率", c2["K_slope"], GRAY, "同上"),
    ("J1 定速 DIAG 平均", c1d["K_mean"], GRAY, "窗口平均重力，假設控制器重力模型正確"),
]
fig, ax = plt.subplots(figsize=(10, 5.4))
y = np.arange(len(est))[::-1]
for yi, (lab, v, col, note) in zip(y, est):
    ax.barh(yi, v, 0.6, color=col, edgecolor="white", zorder=3)
    ax.text(v - 0.12, yi, f"{v:.2f}", va="center", ha="right", fontsize=11.5, color="white", fontweight="bold",
            zorder=4)                                   # 數值放在長條內，避開兩條文獻值虛線
    ax.text(0.15, yi, note, va="center", ha="left", fontsize=9.5, color="white", zorder=4)
ax.axvline(KT, color=INK, lw=1.8, ls=(0, (6, 4)), zorder=5)
ax.text(KT + 0.08, y[0] + 0.5, f"文獻值 {KT}\n（Raviola 2021）", ha="left", va="bottom", fontsize=10, color=INK,
        linespacing=1.2)
ax.axvline(KT_CB, color=INK2, lw=1.4, ls=(0, (2, 2)), zorder=5)
ax.text(KT_CB - 0.08, y[0] + 0.5, f"文獻值 {KT_CB}\n（Kovincic 2019）", ha="right", va="bottom", fontsize=10,
        color=INK2, linespacing=1.2)
ax.set_yticks(y); ax.set_yticklabels([e[0] for e in est], fontsize=11)
ax.set_xlim(0, 16.2); ax.set_ylim(-0.6, len(est) + 0.5)
ax.set_xlabel("換算比例（力矩 ÷ 電流）[N·m/A]")
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(axis="x", alpha=1, ls=(0, (4, 3))); ax.set_axisbelow(True)
ax.set_title("圖 3　換算比例：用砝碼的重力斜率量到 12.0～12.1，比文獻值 13.64 低約 11%、比 12.63 低約 4%", loc="left",
             fontsize=13.5, pad=12, color=INK)
fig.text(0.01, 0.005, "彩色：只用已知質量（砝碼 1.008 kg）與幾何，不需要控制器的模型；灰色：要假設控制器的重力模型正確。"
         "本報告其他圖改用砝碼校正值（J1 %.2f、J2 %.2f）。"
         "\n控制器自己的重力斜率與 m·g·h 相比：J1 1.00、J2 1.03（控制器知道砝碼重量）。"
         "兩個文獻值：13.64 = 101 × 0.135 沒寫出處；12.63 = 101 × 0.125 標註來自手臂文件。" % (K1, K2), ha="left", va="bottom",
         fontsize=9, color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.06, 1, 1))
f3 = save(fig, "w1118_fig3_conversion_ratio.png")

# ---------------- 圖 4：敲擊測試（11/12 已畫好） ----------------
src = os.path.join(OUT, "w1112_tap_spectrum.png")
shutil.copy(src, os.path.join(FIGS, "w1118_fig4_tap_spectrum.png"))
shutil.copy(src, os.path.join(DOCS, "w1118_fig4_tap_spectrum.png"))

# ---------------- 圖 5：J2 定速與速度調變 ----------------
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4), gridspec_kw=dict(width_ratios=[1, 1.15]))
ax = axs[0]
# 定速：各段平均力矩取自 j2_1112_log.txt（j2_1112.json 只存回歸結果）；每檔速度畫升、降兩段的平均
seg = {"pos": [], "neg": []}; cur = None
for line in open(os.path.join(OUT, "j2_1112_log.txt"), encoding="utf-8"):
    if line.startswith("--- pos"): cur = "pos"
    elif line.startswith("--- neg"): cur = "neg"
    elif line.startswith("---") or line.startswith("====="): cur = None if "J2 速度調變" in line or "合併" in line else cur
    elif cur and re.match(r"^\s+\d+\s+0\.\d+\s+\d+", line):
        t = line.split(); seg[cur].append((abs(float(t[4])), abs(float(t[5]))))
vv = np.linspace(0, 0.21, 20)
# 扣掉窗口內的平均重力 H 再畫：正轉 |τ| = Tc + B·v + H，反轉 |τ| = Tc + B·v − H，
# 所以 |τ| − sgn·H 兩個方向都剩下摩擦 Tc + B·v，應該疊在同一條線上
for d, col, sgn, lab, kw in (("pos", BLUE, +1, "正轉", dict(marker="o", s=55, color=BLUE, zorder=5)),
                             ("neg", ORANGE, -1, "反轉（取絕對值）",
                              dict(marker="s", s=150, facecolors="none", edgecolors=ORANGE, linewidths=2, zorder=4))):
    a = np.array(seg[d]); assert len(a) == 12, (d, len(a))
    v = np.round(a[:, 0], 2); vs = np.unique(v)
    ax.scatter(vs, [S2 * (a[v == x, 1].mean() - sgn * c2["H"]) for x in vs], label=lab, **kw)
B2, sB2, Tc2, sTc2, H2 = (S2 * c2[k] for k in ("B", "sd_B", "Tc", "sd_Tc", "H"))
ax.plot(vv, Tc2 + B2 * vv, color=INK2, lw=1.4, zorder=3, label="摩擦 Tc + B·v")
ax.text(0.205, 6.6, f"B = {B2:.1f} ± {sB2:.1f} N·m·s/rad\nTc = {Tc2:.2f} ± {sTc2:.2f} N·m\n"
        f"已扣窗口內重力平均 {H2:+.2f} N·m", ha="right", va="bottom", fontsize=11, color=INK, linespacing=1.5)
ax.legend(loc="upper left", frameon=False, fontsize=10.5)
ax.set_xlim(0, 0.215); ax.set_ylim(6, 14.5)
ax.set_xlabel("J2 速度 [rad/s]"); ax.set_ylabel("扣掉重力後的力矩大小 [N·m]")
ax.set_title("定速（-16.4°～+5.0°，6 檔速度），扣掉重力", loc="left", fontsize=12, color=INK)
style(ax)
ax = axs[1]
vm = J2["velmod"]
order = sorted(range(6), key=lambda i: (vm["pos"][str(i)]["f"], vm["pos"][str(i)]["A"]))
x, w = np.arange(6), 0.36
for d, col, lab, sh in (("pos", BLUE, "正轉", -1), ("neg", ORANGE, "反轉", +1)):
    J = S2 * np.array([vm[d][str(i)]["J"] for i in order]); s = S2 * np.array([vm[d][str(i)]["se"] for i in order])
    Jall = J if d == "pos" else np.concatenate([Jall, J])
    xx = x + sh * (w / 2 + 0.01)
    ax.bar(xx, J, w, color=col, edgecolor="white", linewidth=2, zorder=3, label=lab)
    ax.errorbar(xx, J, yerr=s, fmt="none", ecolor=INK, elinewidth=1.2, capsize=3, zorder=4)
    for i in range(6):
        ax.text(xx[i], J[i] - 0.08, f"{J[i]:.2f}", ha="center", va="top", fontsize=9, color="white",
                fontweight="bold", zorder=6, rotation=90)
jm2 = float(np.mean([J2C["R"]["p130_pos"][str(i)]["Jm"] for i in range(6)]))
ax.axhline(jm2, color=INK, lw=1.8, ls=(0, (6, 4)), zorder=5)
ax.text(5.55, jm2 + 0.05, f"控制器 J_model = {jm2:.2f}", ha="right", va="bottom", fontsize=11, color=INK)
ax.axhline(1.01, color=MUTED, lw=1.3, ls=(0, (2, 2)), zorder=5)
ax.annotate("官方連桿 + 末端工具 1.01\n（不含轉子）", xy=(5.72, 1.02), xytext=(5.85, 2.05), ha="right", va="bottom",
            fontsize=9.5, color=INK2, linespacing=1.3, zorder=7,
            arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.0, shrinkA=2, shrinkB=1,
                            relpos=(0.97, 0)))         # 箭頭從文字右下角往下，走最右邊的空白，不壓長條
ax.legend(loc="upper left", frameon=False, fontsize=10.5, ncol=2)
ax.set_xticks(x)
ax.set_xticklabels([f"{vm['pos'][str(i)]['f']:g} Hz\n× {vm['pos'][str(i)]['A']:g}" for i in order], fontsize=10)
ax.set_xlim(-0.6, 5.95); ax.set_ylim(0, 3.6)
ax.set_ylabel(f"J2 [kg·m²]（電流 × {K2:.2f}，砝碼校正）")
ax.set_title("速度調變（調變段中心 -5.4°）", loc="left", fontsize=12, color=INK)
style(ax)
fig.suptitle(f"圖 5　J2 第一次上機：扣掉重力後定速正反轉重合；速度調變的 J 也隨頻率上升，約為控制器的 "
             f"{Jall.min() / jm2:.2f}～{Jall.max() / jm2:.2f} 倍", x=0.01, ha="left", fontsize=13.5, color=INK)
fig.text(0.01, 0.005, f"11/12 手臂朝上（J2 = 0°，前臂朝上），末端工具 1.30 kg。電流 × {K2:.2f}（J2 砝碼校正）；分析含重力項。"
         f"定速每點為同一檔速度升、降兩段的平均，已扣窗口內平均重力（正轉加 {abs(H2):.2f}、反轉減 {abs(H2):.2f}）；"
         "右圖誤差棒為單輪回歸的標準差。", ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 0.94))
f5 = save(fig, "w1118_fig5_j2.png")

# ---------------- 數字（講稿用） ----------------
print("J1 各頻率平均 J：", {k: [round(v, 3) for v in fm] for k, fm in fmean.items()})
print(f"J1 J_model 1.30 {jm130:.3f}、2.31 {jm231:.3f}")
print("J1 ΔJ：", {k: [round(S1 * v[0], 3), round(S1 * v[1], 3)] for k, v in dj1.items()}, " ΔJ_model", round(jmd1, 4))
print("J1 彈簧模型 ΔJ：", [round(dj_spring(f), 3) for f in (1, 2, 3)])
print(f"J1 換算比例：m g h {mgh1:.3f}，Δ斜率 {dS1:.3f} → {K1:.2f}；J2 {J2C['Kcal']:.2f} ± {J2C['sdK']:.2f}；"
      f"彈簧擬合 {KT / S_FIT:.2f}；換算 S1 {S1:.4f}、S2 {S2:.4f}")
print(f"J2 定速（× {K2:.2f}）B {B2:.3f} ± {sB2:.3f}、Tc {Tc2:.3f} ± {sTc2:.3f}、H {H2:+.3f}；J2 J_model {jm2:.3f}；"
      f"速度調變 J {Jall.min():.3f}～{Jall.max():.3f}")
print(f"J1 定速 DIAG（× {K1:.2f}）B {S1 * c1d['B']:.3f} ± {S1 * c1d['sd_B']:.3f}、Tc {S1 * c1d['Tc']:.3f} ± "
      f"{S1 * c1d['sd_Tc']:.3f}、H {S1 * c1d['H']:+.3f}")
for p in (f1, f2, f3, f5):
    print(p)
