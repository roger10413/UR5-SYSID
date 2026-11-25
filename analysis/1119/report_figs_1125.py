"""11/25 報告的圖（11/19 上機）。
圖 2：J1 速度調變，J 對頻率：11/5、11/12、11/19 三次（1.30 kg），11/12、11/19 加砝碼（2.31 kg）。
圖 3：J1 砝碼的 ΔJ 對頻率：11/12 與 11/19，和控制器 ΔJ_model 比。
圖 1：J2～J5 定速 DIAG：扣掉窗口內平均重力後的力矩對速度，Tc + B·v。
圖 4：手腕 J3～J5 速度調變的 J（負值），和控制器 J_model 比。
圖 5：敲擊測試的餘振頻譜，UP（11/12）與 DIAG（11/19）比。
只讀已存的結果（analysis_1112/j1_1112*.json、analysis_1119/analyze_1119.json）與敲擊 CSV，不連網路。
輸出 OUT/figs/w1125_fig*.png，並複製到 docs/report_1125/。（Microsoft JhengHei 沒有「−」，數值用 ASCII 負號。）"""
import socket
socket.socket = None                                   # 離線分析：不連網路
import os, sys, json, math, glob, shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"F:\cluade_code\UR5-SYSID\robot")
import ur5_tap_test_record as tp

plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12,
                     "axes.edgecolor": "#c3c2b7", "axes.labelcolor": "#0b0b0b", "xtick.color": "#52514e",
                     "ytick.color": "#52514e", "grid.color": "#e1e0d9", "grid.linewidth": 0.8})
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = r"D:\cluade_code\ur5\analysis_1119"
FIGS = os.path.join(OUT, "figs"); os.makedirs(FIGS, exist_ok=True)
DOCS = os.path.join(HERE, "..", "..", "docs", "report_1125"); os.makedirs(DOCS, exist_ok=True)
INK, INK2, MUTED, GRAY = "#0b0b0b", "#52514e", "#8d8b83", "#b9b8b0"
ORANGE, BLUE, GREEN, PURPLE, RED = "#d9822b", "#2a78d6", "#3a9a5b", "#8a5cc2", "#c8463c"
KT = 13.64                                             # 文獻值 101 × 0.135；analysis_1112/*.json 用它換算
K1, K2, KW = 12.03, 12.10, 9.31                        # J1、J2 砝碼校正；手腕暫用 101 × 0.0922（Kovincic 2019）
S1 = K1 / KT

A12 = json.load(open(r"D:\cluade_code\ur5\analysis_1112\j1_1112.json", encoding="utf-8"))
R = json.load(open(os.path.join(OUT, "analyze_1119.json"), encoding="utf-8"))   # 已換成各軸換算比例


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(axis="y", alpha=1, ls=(0, (4, 3))); ax.set_axisbelow(True)


def save(fig, name):
    p = os.path.join(FIGS, name)
    fig.savefig(p, dpi=150, bbox_inches="tight", pad_inches=0.15); plt.close(fig)
    shutil.copy(p, os.path.join(DOCS, name))
    return p


def wmean(vals):
    """[(J, se), ...] 的加權平均與標準差。"""
    w = np.array([1 / s ** 2 for _, s in vals])
    return float(sum(J * wi for (J, _), wi in zip(vals, w)) / w.sum()), float(1 / math.sqrt(w.sum()))


def j1_points(run):
    """run = '1105_p130' 等（11/5、11/12 用 13.64 存的 JSON）或 '1119_p130'。回傳 [(f, A, J, se)]。"""
    out = []
    for ci in range(6):
        if run.startswith("1119"):
            v = [R[f"J1_velmod_up_{run[5:]}_{t}"]["per"][str(ci)] for t in ("pos", "neg")]
            J, se = wmean([(x["J"], x["se"]) for x in v])
        else:
            v = [A12[f"velmod_{run}_{t}"]["per"][str(ci)] for t in ("pos", "neg")]
            J, se = wmean([(S1 * x["J"], S1 * x["se"]) for x in v])
        out.append((v[0]["f"], v[0]["A"], J, se))
    return out


def label_right(ax, pp, col, yl, fmt="{:.2f}"):
    """數值標在點的右邊：折線往上就標在線下方，往下就標在線上方。"""
    dy = 0.035 * (yl[1] - yl[0])
    for n, (f, m, s) in enumerate(pp):
        if n == 0:                                      # 第一點標在左邊空白，避開另一條線
            ax.text(f - 0.07, m, fmt.format(m), ha="right", va="center", fontsize=10.5, color=col, zorder=6,
                    bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))
            continue
        if n == len(pp) - 1:
            y, va = m, "center"
        elif pp[n + 1][1] >= m:
            y, va = m - dy, "top"
        else:
            y, va = m + dy, "bottom"
        ax.text(f + 0.07, y, fmt.format(m), ha="left", va=va, fontsize=10.5, color=col, zorder=6,
                bbox=dict(boxstyle="square,pad=0.1", fc="white", ec="none"))


AMP_TTL = ((True, "大振幅 0.04 rad/s"), (False, "小振幅 0.02 rad/s（3 Hz 為 0.027）"))

# ---------------- 圖 1：J1 三次上機 ----------------
runs = (("1105_p130", GRAY, "11/5　1.30 kg", "^"), ("1112_p130", BLUE, "11/12　1.30 kg", "o"),
        ("1119_p130", GREEN, "11/19　1.30 kg", "D"), ("1112_p231", ORANGE, "11/12　2.31 kg（加砝碼）", "s"),
        ("1119_p231", RED, "11/19　2.31 kg（加砝碼）", "P"))
pts = {r: j1_points(r) for r, *_ in runs}
jm130 = R["J1_velmod_up_p130_pos"]["J_model_mean"]; jm231 = R["J1_velmod_up_p231_pos"]["J_model_mean"]
fig, axs = plt.subplots(1, 2, figsize=(12.5, 6.0), sharey=True)
for ax, (big, ttl) in zip(axs, AMP_TTL):
    for run, col, lab, mk in runs:
        pp = sorted((f, J, s) for f, A, J, s in pts[run] if (A >= 0.04) == big)
        ax.errorbar([p[0] for p in pp], [p[1] for p in pp], yerr=[p[2] for p in pp], fmt="-" + mk, color=col, lw=1.8,
                    ms=8, mfc=col, mec="white", mew=1.0, ecolor=col, elinewidth=1.3, capsize=4, capthick=1.3,
                    zorder=4, label=lab)
    for jm, col in ((jm130, BLUE), (jm231, ORANGE)):
        ax.axhline(jm, color=col, lw=1.2, ls=(0, (6, 4)), zorder=2)
    ax.set_title(ttl, loc="left", fontsize=12.5, color=INK)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["1 Hz", "2 Hz", "3 Hz"])
    ax.set_xlim(0.7, 3.3); ax.set_ylim(4.8, 10.2); ax.set_xlabel("調變頻率")
    style(ax)
axs[0].set_ylabel(f"J1 [kg·m²]（電流 × {K1}，砝碼校正）")
axs[0].legend(loc="upper left", frameon=False, fontsize=10.5)
for jm, col, txt in ((jm130, BLUE, "1.30 kg"), (jm231, ORANGE, "2.31 kg")):
    axs[1].text(3.35, jm, f"控制器 J_model {jm:.2f}\n（{txt}）", ha="left", va="center", fontsize=10, color=col,
                linespacing=1.3)
fig.suptitle("圖 2　J1 速度調變：11/12 與 11/19 幾乎重合，11/5 偏低", x=0.01, ha="left", fontsize=13.5, color=INK)
fig.text(0.01, 0.005, "手臂朝上，調變段中心 -92°（重力零點），巡航 0.08 rad/s，STD 條件組。每點為一個條件正反轉的加權平均，"
         f"誤差棒為單輪回歸的標準差（大多比點小）。電流 × {K1}（J1 砝碼校正）。\n"
         "虛線為控制器的 J_model。J1 溫度：11/5 約 30.7 °C、11/12 約 29.9 °C、11/19 約 33.0 °C。", ha="left", va="bottom",
         fontsize=9, color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.06, 0.93, 0.95))
f1 = save(fig, "w1125_fig2_j1_velmod_3runs.png")

# ---------------- 圖 2：J1 砝碼 ΔJ，11/12 與 11/19 ----------------
def dj_points(date):
    out = []
    for ci in range(6):
        vals, f, A = [], None, None
        for t in ("pos", "neg"):
            if date == "1119":
                a = R[f"J1_velmod_up_p130_{t}"]["per"][str(ci)]; b = R[f"J1_velmod_up_p231_{t}"]["per"][str(ci)]
                vals.append((b["J"] - a["J"], math.hypot(a["se"], b["se"])))
            else:
                a = A12[f"velmod_1112_p130_{t}"]["per"][str(ci)]; b = A12[f"velmod_1112_p231_{t}"]["per"][str(ci)]
                vals.append((S1 * (b["J"] - a["J"]), S1 * math.hypot(a["se"], b["se"])))
            f, A = a["f"], a["A"]
        m, s = wmean(vals)
        out.append((f, A, m, s))
    return out

djp = {"1112": dj_points("1112"), "1119": dj_points("1119")}
djm = jm231 - jm130
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
yl = (0, 2.6)
for ax, (big, ttl) in zip(axs, AMP_TTL):
    for date, col, mk, lab in (("1112", ORANGE, "s", "11/12"), ("1119", RED, "P", "11/19")):
        pp = sorted((f, m, s) for f, A, m, s in djp[date] if (A >= 0.04) == big)
        ax.errorbar([p[0] for p in pp], [p[1] for p in pp], yerr=[p[2] for p in pp], fmt="-" + mk, color=col, lw=1.8,
                    ms=8, mfc=col, mec="white", mew=1.0, ecolor=col, elinewidth=1.3, capsize=4, capthick=1.3,
                    zorder=4, label=lab)
        if date == "1119":
            label_right(ax, pp, col, yl)
    ax.axhline(djm, color=INK, lw=1.4, ls=(0, (6, 4)), zorder=2)
    ax.text(3.3, djm - 0.05, f"控制器 ΔJ_model {djm:.2f}", ha="right", va="top", fontsize=10, color=INK)
    ax.set_title(ttl, loc="left", fontsize=12.5, color=INK)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["1 Hz", "2 Hz", "3 Hz"])
    ax.set_xlim(0.55, 3.4); ax.set_ylim(*yl); ax.set_xlabel("調變頻率")
    style(ax)
axs[0].set_ylabel("ΔJ [kg·m²]")
axs[0].legend(loc="upper left", frameon=False, fontsize=10.5)
fig.suptitle("圖 3　J1 加砝碼的 ΔJ：兩週一致，1 Hz 接近控制器的值，頻率越高越偏大", x=0.01, ha="left", fontsize=13.5,
             color=INK)
fig.text(0.01, 0.005, "ΔJ = J(2.31 kg) - J(1.30 kg)，同一天、同一方向相減後兩個方向加權平均；數值標的是 11/19。"
         f"誤差棒為標準差。電流 × {K1}。控制器 ΔJ_model：target_moment 對 target_qdd 回歸，2.31 kg 減 1.30 kg。",
         ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
f2 = save(fig, "w1125_fig3_j1_dJ.png")

# ---------------- 圖 1：J2～J5 定速 DIAG ----------------
fig, axs = plt.subplots(1, 4, figsize=(15, 4.8))
vv = np.linspace(0, 0.21, 20)
for ax, j, yl in zip(axs, (2, 3, 4, 5), ((6, 15), (1.4, 3.0), (1.2, 2.4), (1.4, 2.8))):
    c = R[f"J{j}_constvel_diag"]
    for d, col, sgn, lab, kw in (("pos", BLUE, +1, "正轉", dict(marker="o", s=40, color=BLUE, zorder=5)),
                                 ("neg", ORANGE, -1, "反轉（絕對值）",
                                  dict(marker="s", s=110, facecolors="none", edgecolors=ORANGE, linewidths=1.8, zorder=4))):
        seg = [s for s in c["seg"] if s[0] == d]
        # 段內力矩是原始正負號；乘方向後是「運動方向的力矩」= Tc + B·v + sgn·H
        v = np.round(np.abs([s[1] for s in seg]), 3); tau = sgn * np.array([s[3] for s in seg])
        vs = np.unique(v)
        ax.scatter(vs, [tau[v == x].mean() - sgn * c["H"] for x in vs], label=lab, **kw)
        # 各方向自己的直線（扣重力後截距都是 Tc，差別只在斜率）：C_pos = Tc + H、C_neg = Tc - H（C 已乘方向）
        b_d, c_d = c[f"B_{d}"], c[f"C_{d}"] - sgn * c["H"]
        ax.plot(vv, c_d + b_d * vv, color=col, lw=1.3, ls="--", zorder=2,
                label=("正轉擬合" if d == "pos" else "反轉擬合"))
    ax.plot(vv, c["Tc"] + c["B"] * vv, color=INK2, lw=1.2, zorder=3, alpha=0.7, label="兩方向合併 Tc + B·v")
    dB = 100 * (c["B_pos"] / c["B_neg"] - 1)
    ax.text(0.01, yl[1] - 0.03 * (yl[1] - yl[0]),
            f"B 正轉 {c['B_pos']:.2f}、反轉 {c['B_neg']:.2f}（差 {dB:+.0f}%）\n"
            f"合併 B = {c['B']:.2f}，Tc = {c['Tc']:.2f}\n已扣重力 {c['H']:+.2f}",
            ha="left", va="top", fontsize=10.5, color=INK, linespacing=1.4)
    ax.set_xlim(0, 0.215); ax.set_ylim(*yl)
    ax.set_xlabel(f"J{j} 速度 [rad/s]")
    ax.set_title(f"J{j}（× {c['K_used']}）", loc="left", fontsize=12.5, color=INK)
    style(ax)
axs[0].set_ylabel("扣掉重力後的力矩大小 [N·m]")
axs[3].legend(loc="lower right", frameon=False, fontsize=10)
fig.suptitle("圖 1　J2～J5 定速（診斷姿態）：沒有重力的 J4、J5 正反重合；有重力的 J2、J3 高速時分開，重力幫忙的正轉 B 較大",
             x=0.01, ha="left", fontsize=13.5, color=INK)
fig.text(0.01, 0.005, "11/19 診斷姿態，末端工具 1.30 kg，9 檔速度（0.005～0.20 rad/s）先升後降，每點為同一檔升、降兩段的平均。"
         f"B [N·m·s/rad]、Tc [N·m]。電流 × 各軸換算比例：J2 {K2}（砝碼校正），J3～J5 {KW}（暫定，101 × 0.0922）。\n"
         "已扣重力 = 兩方向截距差的一半，所以扣完後兩方向的截距必然相同，要比的是兩條虛線的斜率。"
         "J2、J3 的正轉是重力幫忙轉的方向。",
         ha="left", va="bottom", fontsize=9, color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.075, 1, 0.94))
f3 = save(fig, "w1125_fig1_friction_j2_j5.png")

# ---------------- 圖 4：手腕 J3～J5 速度調變 ----------------
AMPW = ((True, "大振幅 0.09 rad/s"), (False, "小振幅 0.045 rad/s（3 Hz 為 0.06）"))
fig, axs = plt.subplots(1, 2, figsize=(12.5, 5.4), sharey=True)
yl = (-0.16, 0.30)
for ax, (big, ttl) in zip(axs, AMPW):
    for j, col, mk in ((3, BLUE, "o"), (4, GREEN, "s"), (5, PURPLE, "D")):
        P, N = R[f"J{j}_velmod_diag_p130_pos"]["per"], R[f"J{j}_velmod_diag_p130_neg"]["per"]
        pp = []
        for ci in P:
            if (P[ci]["A"] >= 0.09) == big:
                m, s = wmean([(P[ci]["J"], P[ci]["se"]), (N[ci]["J"], N[ci]["se"])])
                pp.append((P[ci]["f"], m, s))
        pp.sort()
        ax.errorbar([p[0] for p in pp], [p[1] for p in pp], yerr=[p[2] for p in pp], fmt="-" + mk, color=col, lw=1.8,
                    ms=8, mfc=col, mec="white", mew=1.0, ecolor=col, elinewidth=1.3, capsize=4, capthick=1.3,
                    zorder=4, label=f"J{j} 量測")
        jm = R[f"J{j}_velmod_diag_p130_pos"]["J_model_mean"]
        ax.axhline(jm, color=col, lw=1.1, ls=(0, (6, 4)), zorder=2)
        if not big:
            ax.text(3.35, jm, f"J{j} 控制器 {jm:.2f}", ha="left", va="center", fontsize=10, color=col)
    ax.axhline(0, color=INK, lw=1.2, zorder=3)
    ax.set_title(ttl, loc="left", fontsize=12.5, color=INK)
    ax.set_xticks([1, 2, 3]); ax.set_xticklabels(["1 Hz", "2 Hz", "3 Hz"])
    ax.set_xlim(0.7, 3.3); ax.set_ylim(*yl); ax.set_xlabel("調變頻率")
    style(ax)
axs[0].set_ylabel(f"J [kg·m²]（電流 × {KW}，暫定）")
axs[0].legend(loc="center left", frameon=False, fontsize=10.5)
fig.suptitle("圖 4　手腕 J3～J5 速度調變：量到的 J 是負的，頻率越高越負（控制器約 0.21～0.27）", x=0.01, ha="left",
             fontsize=13.5, color=INK)
fig.text(0.01, 0.005, "11/19 診斷姿態，巡航 0.12 rad/s，每段調變 3 s；分析含重力項。每點為正反轉加權平均，誤差棒為標準差。"
         "虛線為控制器的 J_model（target_moment 的單位）。\n拿掉漣波項、v² 項或改寫重力項，J 都差不到 0.015。",
         ha="left", va="bottom", fontsize=9, color=INK2, linespacing=1.5)
fig.tight_layout(rect=(0, 0.06, 0.92, 0.95))
f4 = save(fig, "w1125_fig4_wrist_J.png")

# ---------------- 圖 5：敲擊頻譜 UP vs DIAG ----------------
DT, NF = 0.008, 8192
fr = np.fft.rfftfreq(NF, DT)


def spec(paths):
    """工具加速度三軸功率相加（不受工具方向影響），各下敲擊平均。"""
    P = np.zeros(len(fr)); n = 0
    for F in paths:
        rows, col = tp._load_csv(F); t = col("timestamp")
        acc = np.column_stack([col(f"tool_acc_{a}") for a in "xyz"])
        taps = tp.find_taps(t, np.linalg.norm(acc - np.median(acc, axis=0), axis=1), DT, floor=0.2)
        for a, k in enumerate(taps):
            b = k + 5; e = min(len(t), b + 312, taps[a + 1] - 5 if a + 1 < len(taps) else len(t))
            if e - b < 100:
                continue
            for i in range(3):
                y = acc[b:e, i]; x = np.arange(len(y)); y = y - np.polyval(np.polyfit(x, y, 1), x)
                P += np.abs(np.fft.rfft(y * np.hanning(len(y)), NF)) ** 2
            n += 1
    return P / n, n


g = lambda base, s: glob.glob(os.path.join(base, f"*{s}", "*.csv"))[0]
U, Dg = r"D:\cluade_code\ur5\data_1112\TAP", r"D:\cluade_code\ur5\data_1119\TAP"
cases = (("前後敲", [g(U, "111246"), g(U, "111353")], [g(Dg, "085527")]),
         ("左右敲", [g(U, "111624")], [g(Dg, "085743")]))
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8))
peaks = {}
for ax, (ttl, up, dg) in zip(axs, cases):
    Pu, nu = spec(up); Pd, nd = spec(dg)
    m = (fr >= 1) & (fr <= 20)
    for P, col, lab in ((Pu, BLUE, f"手臂朝上 UP（11/12，{nu} 下）"), (Pd, RED, f"診斷姿態 DIAG（11/19，{nd} 下）")):
        db = 10 * np.log10(P[m] / P[m].max())
        ax.plot(fr[m], db, color=col, lw=2, label=lab)
        lo = (fr[m] >= 3) & (fr[m] <= 6)
        fpk = fr[m][lo][np.argmax(db[lo])]
        peaks[(ttl, lab[:4])] = round(float(fpk), 2)
        ax.annotate(f"{fpk:.1f} Hz", xy=(fpk, db[lo].max()), xytext=(fpk + (1.6 if col == RED else -0.4), 2.0),
                    ha="left" if col == RED else "right", fontsize=10.5, color=col,
                    arrowprops=dict(arrowstyle="-", color=col, lw=0.9))
    ax.set_title(ttl, loc="left", fontsize=12.5, color=INK)
    ax.set_xlabel("頻率 [Hz]"); ax.set_xlim(1, 20); ax.set_ylim(-45, 6)
    style(ax)
axs[0].set_ylabel("功率 [dB]（相對最大峰）")
axs[1].legend(loc="lower left", frameon=False, fontsize=10.5)
fig.suptitle("圖 5　敲擊後的餘振：換到診斷姿態，最低的峰只從 4.3 Hz 升到 4.4～4.7 Hz", x=0.01, ha="left", fontsize=13.5,
             color=INK)
fig.text(0.01, 0.005, "末端工具 1.30 kg。工具加速度計三軸功率相加（兩個姿態工具方向不同）；每一下敲擊後最長 2.5 s、漢寧窗、各下平均。",
         ha="left", va="bottom", fontsize=9, color=INK2)
fig.tight_layout(rect=(0, 0.05, 1, 0.94))
f5 = save(fig, "w1125_fig5_tap_up_diag.png")

# ---------------- 數字 ----------------
for r, _, lab, _ in runs:
    print(lab, [(f, A, round(J, 3), round(s, 3)) for f, A, J, s in pts[r]])
print("ΔJ_model", round(djm, 3))
for d in djp:
    print("ΔJ", d, [(f, A, round(m, 3), round(s, 3)) for f, A, m, s in djp[d]])
print("tap peaks", peaks)
for p in (f1, f2, f3, f4, f5):
    print(p)
