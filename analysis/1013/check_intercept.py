"""10/13 定加速度 v3：含截距回歸 J 的離線檢查（不需上機）。

1. 電流與速度的時間對齊：換階段（加速→定速、定速→減速）前後 ±0.2 s 做 20 次循環的集成平均，
   以 tau = J·qdd(t+k) + B·v + c 擬合，掃描錯位 k（樣本），殘差最小的 k 即為時間差；
   同時得到「換階段瞬間」的 J（只看加速度的突變，與常數摩擦差 Δ 無關）。
2. J 對時間錯位的敏感度：把電流相對速度平移 −2～+2 個樣本後重做含截距回歸。
3. 共同速度範圍：v_lo 提高到 0.10 rad/s，使 4 檔在跳過 0～80 ms 時用到相同速度範圍。
4. Δ 對速度、加速度：Δ_bin = 2(J·a − h_bin)，檢查是否為常數。
5. 拿掉一檔的 J（leave-one-level-out）。
6. 換階段步階：前 [−0.2, −0.064] s、後 [+0.064, +0.2] s 的平均力矩差（不微分），各檔 J_step·a 對 a 做含截距回歸，
   與斜坡分箱的回歸是不同段的資料（加速→定速在 a = 1.2 時前窗口會碰到起步暫態，不用）。
7. 速度範圍組合：v_lo ∈ {0.05, 0.10}、v_hi ∈ {0.95, 0.82}·V_PEAK、跳過 20/40/60 ms。
"""
import sys, os, io, csv, contextlib, tempfile
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot")); sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
import ur5_const_accel_ident_v3 as m
import accel_v3_analysis as av

FS = 125.0; DT = 1 / FS; KT = m.KT_OUT


def load(F, s):
    d = np.genfromtxt(F, delimiter=",", names=True)
    return d["timestamp"], s * d["actual_qd_0"], s * d["target_qd_0"], s * d["actual_current_0"] * KT


def transitions(t, qd, tq, tau):
    """回傳 [(level, 'acc_hold' | 'hold_dec', idx)]：idx = target 開始不變 / 開始下降的樣本。"""
    dtq = np.diff(tq, prepend=tq[0]); out = []
    acc = [r for r in m._runs((dtq > 1e-6) & (tq > 1e-3)) if len(r) >= 10 and tq[r].max() > 0.2]
    dec = [r for r in m._runs((dtq < -1e-6) & (tq > 1e-3)) if len(r) >= 10 and tq[r].max() > 0.2]
    for ra in acc:
        a = np.polyfit(t[ra], tq[ra], 1)[0]; li = int(np.argmin(np.abs(np.array(m.ACCEL_LEVELS) - a)))
        rd = next((r for r in dec if r[0] > ra[-1]), None)
        if rd is None:
            continue
        out.append((li, "acc_hold", ra[-1] + 1)); out.append((li, "hold_dec", rd[0] - 1))
    return out


def step_fit(t, qd, tau, idxs, B_fix, W=25, kmax=4):
    """集成平均後，對錯位 k 擬合 tau = J·qdd(k) + B·v + c；回傳每個 k 的 (J, rms)。"""
    seg = lambda x: np.mean([x[i - W - kmax:i + W + kmax + 1] for i in idxs], 0)
    T, V = seg(tau), seg(qd)
    A = np.gradient(V, DT)
    res = {}
    core = slice(kmax, kmax + 2 * W + 1)
    for k in range(-kmax, kmax + 1):
        a_k = A[kmax + k: kmax + k + 2 * W + 1]           # k > 0：電流對應較晚的加速度（電流超前）
        y = T[core] - B_fix * V[core]
        X = np.column_stack([a_k, np.ones_like(a_k)])
        p, *_ = np.linalg.lstsq(X, y, rcond=None)
        res[k] = (p[0], np.sqrt(np.mean((y - X @ p) ** 2)))
    return res


def step_J(qd, tau, idxs, B_fix=30.0, W0=8, W1=25):
    """換階段前後窗口的平均力矩差 / 加速度差（加速度以速度直線擬合，不微分）。"""
    b = np.arange(-W1, -W0 + 1); f = np.arange(W0, W1 + 1); off = W1 + 1
    T = np.mean([tau[i - off:i + off + 1] for i in idxs], 0); V = np.mean([qd[i - off:i + off + 1] for i in idxs], 0)
    ab = np.polyfit(b * DT, V[off + b], 1)[0]; af = np.polyfit(f * DT, V[off + f], 1)[0]
    dT = T[off + f].mean() - T[off + b].mean(); dV = V[off + f].mean() - V[off + b].mean()
    return (dT - B_fix * dV) / (af - ab), abs(af - ab)


def shifted_csv(F, k, tmpdir):
    """電流欄位相對其他欄位平移 k 個樣本（k > 0：電流往後移，即原本電流超前）。"""
    with open(F, encoding="utf-8") as f:
        rows = list(csv.reader(f))
    head, body = rows[0], rows[1:]
    ci = head.index("actual_current_0"); cur = [r[ci] for r in body]
    n = len(body)
    out = os.path.join(tmpdir, f"shift_{k:+d}.csv")
    with open(out, "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f); w.writerow(head)
        for i in range(max(0, k), min(n, n + k)):
            r = list(body[i]); r[ci] = cur[i - k]; w.writerow(r)
    return out


def intercept_J(F, s, **kw):
    with contextlib.redirect_stdout(io.StringIO()):
        est = m.analyze_const_accel_v3(F, 0, KT, s, m.ACCEL_LEVELS, m.V_PEAK, n_boot=2, **kw)
    p, L = av.fit_intercept(est["cycles"], est["B"])
    return p, L, est


if __name__ == "__main__":
    tmp = tempfile.mkdtemp()
    for name, s, F in av.RUNS:
        print(f"\n==================== {name} ====================")
        t, qd, tq, tau = load(F, s)
        tr = transitions(t, qd, tq, tau)

        # ---- 1. 換階段瞬間 J 與時間錯位 ----
        print("[1] 換階段前後 ±0.2 s 集成平均：tau = J·qdd(t+k) + B·v + c（B 固定 30），各 k 的 J 與殘差 RMS")
        for kind in ("acc_hold", "hold_dec"):
            for li, a in enumerate(m.ACCEL_LEVELS):
                idxs = [i for (l, k_, i) in tr if l == li and k_ == kind]
                r = step_fit(t, qd, tau, idxs, B_fix=30.0)
                kbest = min(r, key=lambda k: r[k][1])
                row = "  ".join(f"k={k:+d}: J={r[k][0]:.2f} rms={r[k][1]:.3f}" for k in (-2, -1, 0, 1, 2))
                print(f"  {kind:8s} a={a:.1f} (n={len(idxs)})  最佳 k={kbest:+d}, J={r[kbest][0]:.2f} | {row}")

        # ---- 2. 時間錯位敏感度 ----
        print("[2] 電流平移 k 個樣本後的含截距回歸 J（k > 0：電流往後移）")
        for k in (-2, -1, 0, 1, 2):
            Fk = F if k == 0 else shifted_csv(F, k, tmp)
            p, _, _ = intercept_J(Fk, s)
            print(f"  k={k:+d}（{k * 8:+d} ms）：J = {p[0]:.3f}，截距 = {p[1]:+.3f}")

        # ---- 3. 共同速度範圍 ----
        print("[3] 跳過時間敏感度：原設定 v_lo = 0.05 vs 共同範圍 v_lo = 0.10")
        for v_lo in (0.05, 0.10):
            row = []
            for st in (0.0, 0.02, 0.04, 0.06, 0.08):
                p, _, _ = intercept_J(F, s, v_lo=v_lo, settle_s=st)
                row.append(f"{int(st * 1000)}ms: {p[0]:.3f}")
            print(f"  v_lo={v_lo:.2f}：" + "  ".join(row))

        # ---- 4. Δ 對速度、加速度 ----
        p, L, est = intercept_J(F, s)
        J = p[0]; edges = np.linspace(0.05, 0.95 * m.V_PEAK, 7)
        print(f"[4] Δ_bin = 2(J·a − h_bin) [N·m]，J = {J:.3f}；速度箱中心 " + " ".join(f"{(edges[i] + edges[i + 1]) / 2:.3f}" for i in range(6)))
        for li, a in enumerate(m.ACCEL_LEVELS):
            c = [x for x in est["cycles"] if x["level"] == li]
            TA = np.nanmean([x["tau_acc"] for x in c], 0); TD = np.nanmean([x["tau_dec"] for x in c], 0)
            VA = np.nanmean([x["v_acc"] for x in c], 0); VD = np.nanmean([x["v_dec"] for x in c], 0)
            aa = L[li]["a"]; h = (TA - TD - est["B"] * (VA - VD)) / 2
            print(f"  a={a:.1f}: " + " ".join(f"{v:6.2f}" for v in 2 * (J * aa - h)))

        # ---- 5. 拿掉一檔 ----
        a = np.array([x["a"] for x in L]); h = np.array([x["h"] for x in L])
        print("[5] 拿掉一檔的含截距 J：" + "  ".join(
            f"去掉 {m.ACCEL_LEVELS[i]:.1f}: {np.polyfit(np.delete(a, i), np.delete(h, i), 1)[0]:.3f}" for i in range(4)))
        res = h - np.polyval(p, a)
        print(f"    四檔殘差 [N·m]：{np.round(res, 3)}")

        # ---- 6. 換階段步階 ----
        print("[6] 換階段步階的 J_step 與其含截距回歸")
        for kind in ("acc_hold", "hold_dec"):
            aa, hh = [], []
            for li, a in enumerate(m.ACCEL_LEVELS):
                if kind == "acc_hold" and a > 1.0:
                    continue
                Js, da = step_J(qd, tau, [i for (l, k_, i) in tr if l == li and k_ == kind])
                aa.append(da); hh.append(Js * da)
            q = np.polyfit(aa, hh, 1)
            print(f"  {kind}: J_step = {np.round(np.array(hh) / np.array(aa), 2)} → 含截距 J = {q[0]:.2f}，截距 {q[1]:+.2f} N·m")

        # ---- 7. 速度範圍組合 ----
        print("[7] 速度範圍 × 跳過時間的含截距 J")
        for v_lo in (0.05, 0.10):
            for vh in (0.95, 0.82):
                row = [f"{int(st * 1000)}ms {intercept_J(F, s, v_lo=v_lo, v_hi_frac=vh, settle_s=st)[0][0]:.3f}" for st in (0.02, 0.04, 0.06)]
                print(f"  v {v_lo:.2f}～{vh * m.V_PEAK:.3f}：" + "  ".join(row))
