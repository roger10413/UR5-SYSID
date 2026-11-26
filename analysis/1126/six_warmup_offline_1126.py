"""11/26 六軸暖機程式（robot/ur5_six_axis_warmup.py）的離線驗證。

用今天 4 輪六軸定速（09:05 冷機、10:10、10:24、10:37），每輪切成升速趟、降速趟兩半，
每半只取暖機程式用的速度（0.05～0.20），當成一批「暖機批」，用暖機程式的 block_result 算六軸 B。
1. 只用 0.05～0.20 算的 B 和全部 9 檔算的 B 差多少（暖機批的 B 是否代表鑑別結果）。
2. 熱機後（後三輪 6 個半輪）相鄰半輪之間 B 的差異 = 批與批之間的雜訊，用來定 B_STABLE_FRAC。
3. 把 8 個半輪依時間餵給 is_converged，看會在哪一批判定暖機完成。
"""
import os, sys, glob
from datetime import datetime
import numpy as np, pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "robot"))
import ur5_six_axis_ident as six  # noqa: E402
import ur5_six_axis_warmup as wu  # noqa: E402

D = r"D:\cluade_code\ur5\data_1126\SIX"
OUT = r"D:\cluade_code\ur5\analysis_1126\warmup_blocks"; os.makedirs(OUT, exist_ok=True)


def split_halves(path):
    """以 J0 的 0.20 rad/s 段（升速趟去、回，降速趟去、回）切開：第 2 段結束與第 3 段開始的中點。"""
    Dd = pd.read_csv(path)
    on = np.abs(np.abs(Dd.target_qd_0.values) - 0.20) < 1e-4
    e = np.flatnonzero(np.diff(on.astype(int))) + 1
    starts, ends = e[::2], e[1::2]
    assert len(starts) == 4, len(starts)
    cut = (ends[1] + starts[2]) // 2
    tag = os.path.basename(path)[-10:-4]
    t0 = datetime.strptime(os.path.basename(path)[-19:-4], "%Y%m%d_%H%M%S")
    out = []
    for k, part in enumerate((Dd.iloc[:cut], Dd.iloc[cut:])):
        fp = os.path.join(OUT, f"{tag}_{'up' if k == 0 else 'dn'}.csv")
        part.to_csv(fp, index=False)
        mid = t0.timestamp() + float(part.timestamp.iloc[len(part) // 2] - Dd.timestamp.iloc[0])
        out.append((fp, f"{tag[:2]}:{tag[2:4]} {'升' if k == 0 else '降'}", mid))
    return out


if __name__ == "__main__":
    runs = sorted(glob.glob(os.path.join(D, "six_constvel_six_p130_full_*", "six_constvel_*[0-9].csv")))
    full = [six.analyze_constvel(p) for p in runs]
    sub = [wu.block_result(p) for p in runs]
    print("1. 整輪：只用 0.05～0.20 的 B ÷ 全部 9 檔的 B − 1")
    for p, f, s in zip(runs, full, sub):
        print(f"  {os.path.basename(p)[-10:-6]}  " + "  ".join(f"J{j} {(s['B'][j] / f[j]['B'] - 1) * 100:+5.1f}%" for j in range(6)))

    halves = [h for p in runs for h in split_halves(p)]
    res = []
    for fp, lab, mid in halves:
        r = wu.block_result(fp)
        res.append((lab, mid, r))
    print("\n2. 每半輪（暖機批）的 B")
    print(f"{'':>10}" + "".join(f"{'J' + str(j):>8}" for j in range(6)))
    for lab, _, r in res:
        print(f"{lab:>10}" + "".join(f"{b:8.3f}" for b in r["B"]))
    warm = np.array([r["B"] for _, _, r in res[2:]])          # 後三輪的 6 個半輪
    d = np.abs(np.diff(warm, axis=0)) / warm.mean(axis=0)
    spread3 = np.array([(warm[k:k + 3].max(0) - warm[k:k + 3].min(0)) / warm[k:k + 3].mean(0) for k in range(len(warm) - 2)])
    print("\n   熱機後相鄰半輪差 |ΔB|/B：中位數 " + "  ".join(f"J{j} {np.median(d[:, j]) * 100:.1f}%" for j in range(6)))
    print("   熱機後連續 3 批的 spread：最大  " + "  ".join(f"J{j} {spread3[:, j].max() * 100:.1f}%" for j in range(6)))
    cold = np.array([r["B"] for _, _, r in res[:2]])
    print("   冷機第 1 輪兩半 → 熱機平均：  " + "  ".join(
        f"J{j} {(cold[0, j] / warm[:, j].mean() - 1) * 100:+.0f}%/{(cold[1, j] / warm[:, j].mean() - 1) * 100:+.0f}%"
        for j in range(6)))

    print("\n3. 整輪（升＋降，等同暖機程式的一批：先升後降）只取 0.05～0.20 的 B")
    Bs = np.array([s["B"] for s in sub])
    for p, b in zip(runs, Bs):
        print(f"  {os.path.basename(p)[-10:-6]}  " + "".join(f"{x:8.3f}" for x in b))
    w = Bs[1:]
    print("   熱機後三批 spread  " + "  ".join(f"J{j} {(w[:, j].max() - w[:, j].min()) / w[:, j].mean() * 100:.1f}%" for j in range(6)))
    print("   冷機 → 熱機平均    " + "  ".join(f"J{j} {(Bs[0, j] / w[:, j].mean() - 1) * 100:+.0f}%" for j in range(6)))

    print(f"\n4. 整輪依時間餵給 is_converged（判定軸 {wu.CONVERGE_JOINTS}、門檻 {wu.B_STABLE_FRAC:.0%}，最少時間設 0）")
    hist = []
    for p, s in zip(runs, sub):
        hist.append(dict(B=s["B"]))
        ok, why, sp = wu.is_converged(hist, elapsed_min=99, min_minutes=0)
        print(f"  {os.path.basename(p)[-10:-6]}  {'完成' if ok else '    '}  {why}")
