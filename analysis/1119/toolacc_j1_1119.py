"""11/19 J1 速度調變：工具端加速度計 vs 編碼器算的加速度。剛體、底座不動時，
工具加速度在切線方向 = r·qdd（r 固定），比值 K = acc/qdd 不隨頻率變；底座晃時 K 會隨頻率變。
加速度計含重力：角度振盪 δ = −qdd/ω² 讓重力在切線方向的分量變 g·δ，所以 K = r + g/ω²（剛體）。
扣掉 g/ω² 後的 r_app 若隨頻率上升，表示工具實際動得比編碼器算的多（平台或連桿在晃）。"""
import sys, glob
import numpy as np, pandas as pd
sys.stdout.reconfigure(encoding="utf-8")
DT = 0.008
COND = [(1, .04), (2, .04), (3, .04), (2, .02), (1, .02), (3, .027)]
FILES = {f"{p}_{d}": glob.glob(rf"D:\cluade_code\ur5\data_1119\J1\j1velmod_up_{p}_{d}_full_*\j1velmod_data_*.csv")[0]
         for p in ("p130", "p231") for d in ("pos", "neg")}


def mod_segments(qt, _):
    """target_qd 非零的區塊依序 = 條件 1..6 輪流；到達巡航速度後 0.8 s 開始調變 3 s。"""
    mv = np.abs(qt) > 1e-6
    edges = np.flatnonzero(np.diff(mv.astype(int)))
    blocks = [(a, b) for a, b in zip(edges[::2] + 1, edges[1::2] + 1) if b - a > 300]
    segs = []
    for k, (a, b) in enumerate(blocks):
        r = a + np.flatnonzero(np.abs(np.abs(qt[a:b]) - 0.08) < 1e-4)[0]
        i0 = r + int(round(0.8 / DT)); i1 = i0 + int(round(3.0 / DT))
        segs.append((i0, i1, k % 6))
    return segs

for key, fn in FILES.items():
    D = pd.read_csv(fn)
    qt, qa = D.target_qd_1.values, D.actual_qd_1.values
    acc = D[["tool_acc_x", "tool_acc_y", "tool_acc_z"]].values
    qdd = np.gradient(qa, DT)
    res = {}
    for i0, i1, ci in mod_segments(qt, None):
        f0 = COND[ci][0]
        t = np.arange(i0, i1) * DT
        X = np.column_stack([np.sin(2 * np.pi * f0 * t), np.cos(2 * np.pi * f0 * t), np.ones_like(t), t])
        cq = np.linalg.lstsq(X, qdd[i0:i1], rcond=None)[0][:2]
        ca = np.array([np.linalg.lstsq(X, acc[i0:i1, k], rcond=None)[0][:2] for k in range(3)])  # 3×2
        zq = cq[0] + 1j * cq[1]
        za = ca[:, 0] + 1j * ca[:, 1]
        res.setdefault(ci, []).append((za / zq))           # 各軸 acc/qdd 的複數比 [m]
    print(f"== {key}")
    for ci in sorted(res):
        H = np.array(res[ci])                               # n×3
        Hm = H.mean(axis=0)
        k = np.argmax(np.abs(Hm))
        mag = np.abs(H[:, k]); ph = np.degrees(np.angle(H[:, k] * np.conj(Hm[k]) / np.abs(Hm[k]))) + np.degrees(np.angle(Hm[k]))
        print(f"  {COND[ci][0]} Hz A={COND[ci][1]:.3f}  n={len(H)}  主軸 {'xyz'[k]}  |acc/qdd| = {mag.mean():.3f} ± {mag.std(ddof=1)/np.sqrt(len(mag)):.3f} m"
              f"  扣重力 r_app = {mag.mean() - 9.81 / (2 * np.pi * COND[ci][0]) ** 2:.3f} m  相位 {np.degrees(np.angle(Hm[k])):+.1f}°")


