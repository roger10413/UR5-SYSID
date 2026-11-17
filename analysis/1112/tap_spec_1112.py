# 離線：敲擊後餘振頻譜，有無砝碼比較；不連線、不動手臂
import sys, glob, socket
socket.socket = None
import numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
sys.path.insert(0, r"D:\cluade_code\UR5-SYSID\robot")
import ur5_tap_test_record as tp
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12})
DT = 0.008; NF = 8192; fr = np.fft.rfftfreq(NF, DT)
def spec(stamps, ax_i):
    P = np.zeros(len(fr)); n = 0
    for s in stamps:
        F = glob.glob(rf"D:\cluade_code\ur5\data_1112\TAP\*{s}\*.csv")[0]
        rows, col = tp._load_csv(F); t = col("timestamp")
        acc = np.column_stack([col(f"tool_acc_{a}") for a in "xyz"])
        taps = tp.find_taps(t, np.linalg.norm(acc - np.median(acc, axis=0), axis=1), DT, floor=0.2)
        for a, k in enumerate(taps):
            b = k + 5; e = min(len(t), b + 312, taps[a + 1] - 5 if a + 1 < len(taps) else len(t))
            if e - b < 100: continue
            y = acc[b:e, ax_i]; x = np.arange(len(y)); y = y - np.polyval(np.polyfit(x, y, 1), x)
            P += np.abs(np.fft.rfft(y * np.hanning(len(y)), NF)) ** 2; n += 1
    return P / n
C0, C1 = "#1f5fbf", "#d62728"
fig, axs = plt.subplots(1, 2, figsize=(12.5, 4.8), sharey=False)
cases = [("前後方向（工具加速度計 z）：前後敲", ["111246", "111353"], ["111944", "112031"], 2),
         ("左右方向（工具加速度計 x）：左右敲", ["111624"], ["111832"], 0)]
for ax, (title, a, b, i) in zip(axs, cases):
    Pa, Pb = spec(a, i), spec(b, i); m = (fr >= 1) & (fr <= 20); ref = Pa[m].max()
    ax.plot(fr[m], 10*np.log10(Pa[m] / ref), color=C0, lw=2, label="只有末端工具（1.30 kg）")
    ax.plot(fr[m], 10*np.log10(Pb[m] / ref), color=C1, lw=2, label="加砝碼（2.31 kg）")
    ax.set_title(title, loc="left", fontsize=12.5); ax.set_xlabel("頻率 [Hz]"); ax.set_xlim(1, 20); ax.set_ylim(-50, 3)
    ax.grid(True, which="major", color="#e5e5e5"); ax.spines[["top", "right"]].set_visible(False)
axs[0].set_ylabel("功率 [dB]（相對最大峰）")
axs[0].annotate("4.3 Hz：加砝碼幾乎不動\n（降不到 3%）", xy=(4.3, 0), xytext=(7.2, -6), fontsize=11, arrowprops=dict(arrowstyle="->", color="#444"))
axs[1].annotate("4.3 → 4.2 Hz", xy=(4.25, 0), xytext=(6.0, -4), fontsize=11, arrowprops=dict(arrowstyle="->", color="#444"))
axs[1].annotate("10.7 → 9.7 Hz\n（降 9%）", xy=(9.8, -10.5), xytext=(12.8, -8), fontsize=11, arrowprops=dict(arrowstyle="->", color="#444"))
axs[1].legend(loc="lower left", frameon=False, fontsize=11)
fig.suptitle("敲擊後的餘振：4.3 Hz 的振動加 1 kg 幾乎不變，10.7 Hz 的振動明顯變低", x=0.01, ha="left", fontsize=14, fontweight="bold")
fig.text(0.01, 0.01, "11/12 手臂朝上（UP）；每一下敲擊後最長 2.5 s、漢寧窗、各下平均；用拳頭敲。", fontsize=10, color="#555")
fig.tight_layout(rect=(0, 0.04, 1, 0.93)); fig.savefig("w1112_tap_spectrum.png", dpi=150); print("ok")
