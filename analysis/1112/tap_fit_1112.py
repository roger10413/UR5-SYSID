# 離線：11/12 敲擊測試，每一下餘振用阻尼正弦擬合頻率與阻尼比；不連線、不動手臂
import sys, glob, os, socket
socket.socket = None
sys.stdout.reconfigure(encoding="utf-8")
import numpy as np
from scipy.optimize import least_squares
sys.path.insert(0, r"D:\cluade_code\UR5-SYSID\robot")
import ur5_tap_test_record as tp
DT = 0.008
NOTE = {"110802": "無", "111155": "亂敲", "111246": "FB", "111353": "FB 輕敲", "111624": "LR",
        "111832": "LR 砝碼", "111944": "FB 砝碼", "112031": "FB 輕敲砝碼（漏一下）"}
def fit(y, f0):
    t = np.arange(len(y)) * DT
    y = y - np.polyval(np.polyfit(t, y, 1), t)
    A0 = np.max(np.abs(y[:20]))
    def r(p):
        A, s, f, ph, c, d = p
        return A * np.exp(-s * t) * np.cos(2 * np.pi * f * t + ph) + c + d * t - y
    best = None
    for ph in (0, np.pi / 2, np.pi, 3 * np.pi / 2):
        o = least_squares(r, [A0, 1.0, f0, ph, 0, 0], bounds=([0, 0, f0 - 1.0, -10, -np.inf, -np.inf], [np.inf, 30, f0 + 1.0, 10, np.inf, np.inf]))
        if best is None or o.cost < best.cost:
            best = o
    A, s, f, ph, c, d = best.x
    r2 = 1 - 2 * best.cost / np.sum((y - y.mean()) ** 2)
    return f, s / np.sqrt(s ** 2 + (2 * np.pi * f) ** 2), r2
out = {}
for F in sorted(glob.glob(r"D:\cluade_code\ur5\data_1112\TAP\*\*.csv")):
    stamp = F.split("_")[-1][:6]
    rows, col = tp._load_csv(F)
    t = col("timestamp"); acc = np.column_stack([col(f"tool_acc_{a}") for a in "xyz"])
    if stamp == "110802":
        print(f"[{NOTE[stamp]}] tool_acc 中位數 x,y,z = {np.round(np.median(acc, axis=0), 3)}")
        continue
    trig = np.linalg.norm(acc - np.median(acc, axis=0), axis=1)
    taps = tp.find_taps(t, trig, DT, floor=0.2)
    print(f"\n[{NOTE[stamp]}] {os.path.basename(os.path.dirname(F))}  敲擊 {len(taps)} 下 @ {[round(t[k],2) for k in taps]}，峰值 {[round(trig[k:k+10].max(),2) for k in taps]} m/s²")
    sig = {"acc_x": acc[:, 0], "acc_z": acc[:, 2], "q1": col("actual_q_1"), "q0": col("actual_q_0"), "q2": col("actual_q_2")}
    for name, y in sig.items():
        fs, zs, r2s = [], [], []
        for a, k in enumerate(taps):
            s = k + 5; e = min(len(t), s + 312, taps[a + 1] - 5 if a + 1 < len(taps) else len(t))
            if e - s < 100: continue
            seg = y[s:e]
            fr = np.fft.rfftfreq(4096, DT); P = np.abs(np.fft.rfft(seg - seg.mean(), 4096))
            m = (fr > 3.0) & (fr < 6.0); f0 = fr[m][np.argmax(P[m])]
            f, z, r2 = fit(seg, f0)
            fs.append(f); zs.append(z); r2s.append(r2)
        if fs:
            print(f"   {name:6s} f = " + ", ".join(f"{v:.2f}" for v in fs) + f"  → 平均 {np.mean(fs):.2f} ± {np.std(fs, ddof=1) if len(fs)>1 else 0:.2f} Hz；ζ 平均 {np.mean(zs):.3f}；R² " + ", ".join(f"{v:.2f}" for v in r2s))
            out.setdefault(stamp, {})[name] = dict(f=fs, z=zs, r2=r2s)
    # 低頻（0.5～3 Hz）有沒有峰：整段扣掉敲擊後的頻譜
    fr = np.fft.rfftfreq(4096, DT)
    for name in ("acc_x", "acc_z", "q1"):
        y = sig[name]; P = np.zeros(len(fr)); n = 0
        for a, k in enumerate(taps):
            s = k + 5; e = min(len(t), s + 312, taps[a + 1] - 5 if a + 1 < len(taps) else len(t))
            if e - s < 100: continue
            seg = y[s:e]; x = np.arange(len(seg)); seg = seg - np.polyval(np.polyfit(x, seg, 1), x)
            P += np.abs(np.fft.rfft(seg * np.hanning(len(seg)), 4096)) ** 2; n += 1
        P /= max(n, 1); m = (fr >= 0.8) & (fr <= 45)
        idx = np.flatnonzero(m); pk = [k for k in idx[1:-1] if P[k] > P[k-1] and P[k] >= P[k+1]]
        pk.sort(key=lambda k: -P[k])
        print(f"   {name:6s} 漢寧窗頻譜峰（0.8–45 Hz）: " + ", ".join(f"{fr[k]:.2f} ({P[k]/P[pk[0]]:.2f})" for k in pk[:4]))
import json; json.dump(out, open(r"D:\cluade_code\ur5\analysis_1112\tap_fit_1112.json", "w"), indent=1)
