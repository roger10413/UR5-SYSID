# 用力矩漣波的相位估「電流相對角度的延遲」（不受摩擦影響；不連線、不動手臂）
#   漣波只跟馬達角 φ = N·q 有關。電流若比角度晚 Δ，量到的是 r(φ − N·v·Δ)，
#   第 k 次諧波的相位會隨速度 v 線性移動，斜率 = k·N·Δ。正反轉、9 檔速度一起回歸。
#   資料：11/19 J2～J5 定速 DIAG；11/12 J1 定速 DIAG。
import sys, os, glob, json, math, csv
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")

N = 101
KMAX = 12
SRC = {1: r"D:\cluade_code\ur5\data_1112", 2: r"D:\cluade_code\ur5\data_1119\J2", 3: r"D:\cluade_code\ur5\data_1119\J3",
       4: r"D:\cluade_code\ur5\data_1119\J4", 5: r"D:\cluade_code\ur5\data_1119\J5"}
OUT = r"D:\cluade_code\ur5\analysis_1119"


def load(path, j):
    with open(path, encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
    g = lambda k: np.array([float(r[k]) for r in rr])
    return g("timestamp"), g(f"actual_q_{j}"), g(f"target_qd_{j}"), g(f"actual_current_{j}")


def segments(t, tqd, min_s=1.5, skip_s=0.3):
    """目標速度固定、|v| > 0.004 的連續區段（去掉頭尾 skip_s）。"""
    out, i, n = [], 0, len(t)
    while i < n:
        if abs(tqd[i]) < 0.004:
            i += 1; continue
        k = i
        while k + 1 < n and abs(tqd[k + 1] - tqd[i]) < 1e-6:
            k += 1
        a, b = t[i] + skip_s, t[k] - skip_s
        ii = np.where((t >= a) & (t <= b))[0]
        if len(ii) and t[ii[-1]] - t[ii[0]] >= min_s:
            out.append((float(tqd[i]), ii))
        i = k + 1
    return out


def fit_seg(q, cur, ii, v):
    """cur = 三次多項式(q) + Σ_k (a_k cos kφ + b_k sin kφ)；只留 k 次諧波在 125 Hz 下不混疊的。"""
    qq = q[ii] - q[ii].mean(); phi = N * q[ii]
    f_mot = abs(v) * N / (2 * math.pi)                      # 馬達每秒轉幾圈
    ks = [k for k in range(1, KMAX + 1) if k * f_mot < 0.4 * 125]
    cols = [np.ones_like(qq), qq, qq ** 2, qq ** 3]
    for k in ks:
        cols += [np.cos(k * phi), np.sin(k * phi)]
    c, *_ = np.linalg.lstsq(np.column_stack(cols), cur[ii], rcond=None)
    e = cur[ii] - np.column_stack(cols) @ c
    return {k: complex(c[4 + 2 * n_], c[5 + 2 * n_]) for n_, k in enumerate(ks)}, float(np.std(e)), len(ii)


res = {}
for j in (1, 2, 3, 4, 5):
    rows = []
    for tag in ("pos", "neg"):
        path = glob.glob(os.path.join(SRC[j], f"j{j}constvel_diag_{tag}_p130_full_*", "*.csv"))[0]
        t, q, tqd, cur = load(path, j)
        for v, ii in segments(t, tqd):
            h, sd, n = fit_seg(q, cur, ii, v)
            rows.append((v, h, sd, n))
    # 每個 k：相位 θ_k(v) = θ0 + k·N·Δ·v；用相量 h_k 的相位（展開）對 v 加權回歸
    print(f"\n===== J{j}：{len(rows)} 段，速度 {sorted({round(abs(r[0]), 3) for r in rows})}")
    print("  k   平均振幅[A]  可用段數   Δ [ms]（± 標準誤差）   殘差相位 RMS [°]")
    resk = {}
    for k in range(1, KMAX + 1):
        pts = [(v, h[k], sd / math.sqrt(n)) for v, h, sd, n in rows if k in h]
        if len(pts) < 6:
            continue
        amp = np.array([abs(p[1]) for p in pts]); amp_m = float(np.median(amp))
        noise = np.array([p[2] for p in pts]) * math.sqrt(2)
        ok = amp > 4 * noise                                  # 振幅明顯大於雜訊才用
        if ok.sum() < 6:
            print(f"{k:>3}   {amp_m:>9.4f}     {ok.sum():>3}      （振幅太小，略過）")
            continue
        v = np.array([p[0] for p in pts])[ok]; z = np.array([p[1] for p in pts])[ok]
        w = (amp[ok] / noise[ok]) ** 2
        # 先以 v≈0 附近的平均相位為基準，相位相對基準展開到 ±π，再線性回歸
        ref = np.angle(np.sum(z[np.argsort(np.abs(v))[:4]]))
        o = np.argsort(v)                                     # 依速度排序後展開（相鄰速度檔的相位差 < π）
        ph = np.empty_like(v); ph[o] = np.unwrap(np.angle(z[o] * np.exp(-1j * ref)))
        ph -= 2 * math.pi * round(float(np.median(ph[np.argsort(np.abs(v))[:4]])) / (2 * math.pi))
        A = np.column_stack([np.ones_like(v), v])
        W = np.diag(w)
        cov = np.linalg.inv(A.T @ W @ A)
        p = cov @ A.T @ W @ ph
        r = ph - A @ p
        s2 = float(r @ W @ r / max(len(v) - 2, 1))
        se = math.sqrt(cov[1, 1] * s2)
        # 量到的相位 = θ0 − k·N·v·Δ（電流晚 Δ 時，量到 r(φ − NvΔ)；cos(kφ − θ) 形式下 θ 增加 kNvΔ）
        # 本程式相量 a − ... 取 complex(a, b) → 相位 = atan2(b, a) = θ，故 Δ = 斜率 / (k·N)
        D = p[1] / (k * N) * 1000; seD = se / (k * N) * 1000
        resk[k] = dict(amp=amp_m, n=int(ok.sum()), delay_ms=float(D), se_ms=float(seD), rms_deg=float(np.degrees(np.sqrt(np.mean(r ** 2)))))
        print(f"{k:>3}   {amp_m:>9.4f}     {ok.sum():>3}      {D:>7.1f} ± {seD:<6.1f}          {np.degrees(np.sqrt(np.mean(r**2))):>6.1f}")
    if resk:
        ww = np.array([1 / v["se_ms"] ** 2 for v in resk.values()]); dd = np.array([v["delay_ms"] for v in resk.values()])
        good = [k for k, v in resk.items() if v["se_ms"] < 5]
        if good:
            wg = np.array([1 / resk[k]["se_ms"] ** 2 for k in good]); dg = np.array([resk[k]["delay_ms"] for k in good])
            print(f"  → 標準誤差 < 5 ms 的諧波（k = {good}）加權平均 Δ = {np.sum(wg*dg)/np.sum(wg):.1f} ± {1/math.sqrt(np.sum(wg)):.1f} ms")
    res[f"J{j}"] = resk

json.dump(res, open(os.path.join(OUT, "ripple_delay_1119.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
