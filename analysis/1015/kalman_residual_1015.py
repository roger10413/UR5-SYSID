"""10/15 卡爾曼濾波殘差分析：找出三參數模型（J、B、Tc）漏掉的力矩，並量化各成分。

流程
 1. 只用編碼器的卡爾曼濾波＋RTS 平滑（kalman_1015.kin_rts）估計每個樣本的速度 v̂；
    雜訊參數以新息概似最大化（三角形 v3 正轉前 8000 樣本）
 2. 回歸分解：以實測力矩對 [指令加速度、v̂、…] 回歸，依序加入候選項，看殘差下降多少
      三參數 → 馬達角漣波 → 溫度 → 速度非線性 → 起步後行程 → 加減速階段
    加速度用指令值（估計的加速度含雜訊，當回歸變數會低估 J）。
    以運動段為單位 bootstrap 估計標準差（殘差在時間上相關，逐樣本的回歸誤差會低估）。
 3. 檢查：J 是否隨運動類型改變、絕對角度與行程能否分開、行程與時間哪個較能解釋
 4. 擴增狀態卡爾曼濾波（kalman_1015.kf_rts，狀態 [q, v, a, d]）：名義三參數模型與擴充模型
    （把鑑別出的各項當成已知力矩）各跑一次，比較模型外力矩 d̂ 的大小與概似
所有數值以運動方向為正（τ、v、a 乘上指令速度的正負號）。
"""
import sys, os, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
import data_1015 as dd
import kalman_1015 as kf

OUT = os.path.join(dd.BASE, "analysis_1015"); os.makedirs(OUT, exist_ok=True)
CACHE = os.path.join(OUT, "kf_cache"); os.makedirs(CACHE, exist_ok=True)
T_REF = 31.0                                        # 溫度項的參考溫度 [°C]
NOMINAL = dict(J=1.75, B={+1: 33.80, -1: 33.79}, Tc={+1: 8.92, -1: 9.16})   # 10/15 定速第 2 輪（31.1 °C）；J 取 10/13 含截距回歸
K_RIPPLE = 24                                       # 漣波：馬達角 1～24 次諧波，正反轉各一組
KNOTS = np.array([0, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10, 13, 16, 20, 25, 30, 40, 50.0])   # 起步後行程分段線性基底 [°]，g(0) = 0
V_MIN_STEADY, V_MIN_RAMP = 0.015, 0.05
ORDER = ["三參數", "漣波", "溫度", "速度非線性", "起步後行程", "加減速階段"]
SMOOTH_N = 25                                       # 「慢變殘差」：同一運動段內 25 個樣本（0.2 s）移動平均
TUNE = slice(0, 8000)                               # 調雜訊用的樣本（三角形 v3 正轉前 64 s）
NOISE_KEYS = ("sig_q", "sig_tau", "sig_j", "sig_d")


def tuning_run(runs):
    return next(x for x in runs if x["kind"] == "tri" and x["dir"] > 0)


# ---------------- 運動學（只用編碼器） ----------------
def kin_noise(runs):
    fn = os.path.join(CACHE, "noise_kin.npz")
    if os.path.exists(fn):
        z = np.load(fn); return {k: float(z[k]) for k in ("sig_q", "sig_j")}
    best = kf.tune_kin(tuning_run(runs)["D"]["q"][TUNE])
    np.savez(fn, **best)
    return best


def run_kin(runs, kn):
    for r in runs:
        fn = os.path.join(CACHE, f"kin_{r['kind']}_{r['ts']}.npz")
        if not os.path.exists(fn):
            xs = kf.kin_rts(r["D"]["q"], kn["sig_q"], kn["sig_j"])["xs"]
            np.savez(fn, v=xs[:, 1], a=xs[:, 2])
        z = np.load(fn); r["v"] = z["v"]; r["a"] = z["a"]
        r["keep"] = keep_mask(r["Fe"], r["Fe"]["sgn"] * r["v"])


# ---------------- 擴增狀態濾波器 ----------------
def kf_inputs(r, model, known=None):
    D, Fe = r["D"], r["Fe"]; sgn = Fe["sgn"]
    hv = np.where(sgn > 0, model["B"][+1], model["B"][-1]) * np.ones(len(sgn))
    if "kB" in model:
        hv = hv + model["kB"] * (D["T"] - T_REF)
    off = np.where(sgn > 0, model["Tc"][+1], np.where(sgn < 0, -model["Tc"][-1], 0.0))
    if known is not None:
        off = off + known
    return D["q"], D["tau"] - off, r["keep"], hv


def tuned_noise(runs, model, tag, kn, known_fn=None):
    """sig_q、sig_j 沿用運動學的值，只調 sig_tau、sig_d。"""
    fn = os.path.join(CACHE, f"noise_{tag}.npz")
    if os.path.exists(fn):
        z = np.load(fn); return {k: float(z[k]) for k in NOISE_KEYS}
    r = tuning_run(runs)
    q, zt, use, hv = kf_inputs(r, model, known_fn(r) if known_fn else None)
    init = dict(sig_q=kn["sig_q"], sig_tau=0.7, sig_j=kn["sig_j"], sig_d=1.0)
    best, _ = kf.tune_noise(q[TUNE], zt[TUNE], use[TUNE], hv[TUNE], model["J"], init)
    np.savez(fn, **best)
    return best


def run_kf(runs, model, noise, tag, known_fn=None):
    res = []
    for r in runs:
        fn = os.path.join(CACHE, f"{tag}_{r['kind']}_{r['ts']}.npz")
        if not os.path.exists(fn):
            q, zt, use, hv = kf_inputs(r, model, known_fn(r) if known_fn else None)
            o = kf.kf_rts(q, zt, use, hv, model["J"], noise)
            np.savez(fn, v=o["xs"][:, 1], a=o["xs"][:, 2], d=o["xs"][:, 3], nu=o["nu_tau"], S=o["S_tau"], ll=o["loglik"])
        res.append(dict(np.load(fn)))
    return res


# ---------------- 回歸 ----------------
def hat_basis(x, knots):
    x = np.clip(x, knots[0], knots[-1]); Bm = np.zeros((len(x), len(knots)))
    for j in range(len(knots)):
        if j > 0:
            m = (x >= knots[j - 1]) & (x <= knots[j]); Bm[m, j] = (x[m] - knots[j - 1]) / (knots[j] - knots[j - 1])
        if j < len(knots) - 1:
            m = (x >= knots[j]) & (x < knots[j + 1]); Bm[m, j] = (knots[j + 1] - x[m]) / (knots[j + 1] - knots[j])
    return Bm[:, 1:]                                 # 去掉 g(0)：以起步瞬間為 0


def keep_mask(Fe, vm):
    """模型適用的樣本：運動中、換階段 40 ms 後、速度不太低（定速段 ≥ 0.015、加減速段 ≥ 0.05 rad/s）。"""
    s, ph = Fe["sgn"], Fe["phase"]
    return (s != 0) & (Fe["since"] >= dd.SETTLE_N) & (vm > 0) & \
           (((ph == 2) & (vm >= V_MIN_STEADY)) | (((ph == 1) | (ph == 3)) & (vm >= V_MIN_RAMP)))


def build(runs, j_split=False):
    """回傳 cols（欄名 → 陣列）、y（運動方向的力矩）、info（各樣本的分組資訊）、G（各組欄名）。"""
    keys = ("seg", "run", "kind", "main", "phase", "travel", "s", "v", "a", "T", "phi", "tsince", "q")
    cols, Y, info = {}, [], {k: [] for k in keys}
    for ri, r in enumerate(runs):
        D, Fe = r["D"], r["Fe"]; s = Fe["sgn"]
        idx = np.flatnonzero(r["keep"])
        vals = dict(seg=ri * 10000 + Fe["seg"][idx], run=np.full(len(idx), ri), kind=np.full(len(idx), r["kind"]),
                    main=Fe["main"][idx], phase=Fe["phase"][idx], travel=Fe["travel"][idx], s=s[idx], v=s[idx] * r["v"][idx],
                    a=Fe["am_cmd"][idx],        # 加速度用指令值（無量測雜訊，避免回歸變數含誤差而低估 J）
                    T=D["T"][idx], phi=Fe["phi"][idx], tsince=Fe["tsince"][idx], q=D["q"][idx])
        for k in keys:
            info[k].append(vals[k])
        Y.append(s[idx] * D["tau"][idx])
    info = {k: np.concatenate(v) for k, v in info.items()}
    y = np.concatenate(Y)
    s, v, a, T, phi = info["s"], info["v"], info["a"], info["T"], info["phi"]
    pos, neg = (s > 0).astype(float), (s < 0).astype(float)
    G = {}
    if j_split:
        for nm, m in (("J 三角形", (info["kind"] == "tri") & info["main"]), ("J 梯形", (info["kind"] == "trap") & info["main"]),
                      ("J 定速段", (info["kind"] == "cv") & info["main"]), ("J movej", ~info["main"])):
            cols[nm] = a * m
        G["三參數"] = list(cols)
    else:
        cols["J"] = a; G["三參數"] = ["J"]
    cols.update({"B+": v * pos, "B-": v * neg, "Tc+": pos, "Tc-": neg}); G["三參數"] += ["B+", "B-", "Tc+", "Tc-"]
    G["漣波"] = []
    for sg, m in (("+", pos), ("-", neg)):
        for k in range(1, K_RIPPLE + 1):
            cols[f"c{k}{sg}"] = np.cos(k * phi) * m; cols[f"s{k}{sg}"] = np.sin(k * phi) * m
            G["漣波"] += [f"c{k}{sg}", f"s{k}{sg}"]
    cols["kB"] = (T - T_REF) * v; cols["kTc"] = T - T_REF; G["溫度"] = ["kB", "kTc"]
    cols["C+"] = v ** 2 * pos; cols["C-"] = v ** 2 * neg; G["速度非線性"] = ["C+", "C-"]
    Hb = hat_basis(info["travel"], KNOTS); G["起步後行程"] = []
    for j in range(Hb.shape[1]):
        cols[f"g{KNOTS[j + 1]:g}"] = Hb[:, j]; G["起步後行程"].append(f"g{KNOTS[j + 1]:g}")
    cols["δ加速"] = (info["phase"] == 1).astype(float); cols["δ減速"] = (info["phase"] == 3).astype(float)
    G["加減速階段"] = ["δ加速", "δ減速"]
    cols["角度"] = s * np.degrees(info["q"] - np.pi)          # 關節座標中與方向無關、隨角度線性變化的力矩
    return cols, y, info, G


def lsq(cols, names, y):
    X = np.column_stack([cols[n] for n in names])
    p, *_ = np.linalg.lstsq(X, y, rcond=None)
    return p, y - X @ p


def seg_bootstrap(cols, names, y, seg, n_boot=200, seed=0):
    """以運動段為單位重抽樣：先算各段的 XᵀX、Xᵀy，再加權相加。"""
    X = np.column_stack([cols[n] for n in names])
    u, inv = np.unique(seg, return_inverse=True); p = X.shape[1]
    XtX = np.zeros((len(u), p, p)); Xty = np.zeros((len(u), p))
    order = np.argsort(inv, kind="stable"); bounds = np.searchsorted(inv[order], np.arange(len(u) + 1))
    for g in range(len(u)):
        ii = order[bounds[g]:bounds[g + 1]]; Xg = X[ii]
        XtX[g] = Xg.T @ Xg; Xty[g] = Xg.T @ y[ii]
    rng = np.random.default_rng(seed); out = []
    for _ in range(n_boot):
        w = np.bincount(rng.integers(0, len(u), len(u)), minlength=len(u)).astype(float)
        out.append(np.linalg.solve(np.einsum("g,gij->ij", w, XtX) + 1e-9 * np.eye(p), w @ Xty))
    return np.array(out)


def slow_rms(e, seg, n=SMOOTH_N):
    """同一運動段內做 n 點移動平均後的 RMS（去掉白雜訊後的系統性殘差）。"""
    order = np.argsort(seg, kind="stable"); cut = np.flatnonzero(np.diff(seg[order])) + 1; out = []
    for x in np.split(e[order], cut):
        if len(x) >= n:
            out.append(np.convolve(x, np.ones(n) / n, "valid"))
    return np.sqrt(np.mean(np.concatenate(out) ** 2))


def known_extended(r, P):
    """完整模型中 J、B、Tc、kB 以外的各項 → 關節座標的已知力矩（給擴充模型的濾波器用）。"""
    D, Fe = r["D"], r["Fe"]; s = Fe["sgn"]; phi = Fe["phi"]; vm = s * r["v"]
    c = np.zeros(len(s))
    for sg, m in ((+1, s > 0), (-1, s < 0)):
        tag = "+" if sg > 0 else "-"
        rip = sum(P[f"c{k}{tag}"] * np.cos(k * phi) + P[f"s{k}{tag}"] * np.sin(k * phi) for k in range(1, K_RIPPLE + 1))
        c[m] += (rip + P[f"C{tag}"] * vm ** 2)[m]
    c += P["kTc"] * (D["T"] - T_REF)
    c += hat_basis(np.nan_to_num(Fe["travel"]), KNOTS) @ np.array([P[f"g{k:g}"] for k in KNOTS[1:]])
    c += P["δ加速"] * (Fe["phase"] == 1) + P["δ減速"] * (Fe["phase"] == 3)
    return np.where(s != 0, s * c, 0.0)


def prepare(runs):
    """讀檔、算運動特徵、跑運動學濾波器；回傳運動學雜訊參數。"""
    for r in runs:
        r["D"] = dd.load(r["F"]); r["Fe"] = dd.motion_features(r["D"], r["dir"])
    kn = kin_noise(runs)
    run_kin(runs, kn)
    return kn


def ext_model(P):
    return dict(J=P["J"], B={+1: P["B+"], -1: P["B-"]}, Tc={+1: P["Tc+"], -1: P["Tc-"]}, kB=P["kB"])


def param_tag(p):
    return hashlib.md5(np.round(p, 6).tobytes()).hexdigest()[:8]


def kf_summary(runs, kres):
    """使用電流量測的樣本上：d̂ RMS、慢變 d̂ RMS、電流通道新息 RMS、每個樣本的電流通道對數概似。"""
    keep = [r["keep"] for r in runs]
    d = np.concatenate([r["Fe"]["sgn"][k] * o["d"][k] for r, o, k in zip(runs, kres, keep)])
    sg = np.concatenate([ri * 10000 + r["Fe"]["seg"][k] for ri, (r, k) in enumerate(zip(runs, keep))])
    nu = np.concatenate([o["nu"][k] for o, k in zip(kres, keep)]); S = np.concatenate([o["S"][k] for o, k in zip(kres, keep)])
    return np.sqrt(np.mean(d ** 2)), slow_rms(d, sg), np.sqrt(np.mean(nu ** 2)), -0.5 * np.mean(np.log(2 * np.pi * S) + nu ** 2 / S)


if __name__ == "__main__":
    runs = dd.run_list(); kn = prepare(runs)

    # ---- 0. 運動學濾波器設定 ----
    r0 = tuning_run(runs); D0, F0 = r0["D"], r0["Fe"]
    acmd = np.gradient(D0["tq"], dd.DT); k0 = (F0["sgn"] != 0) & (F0["since"] >= dd.SETTLE_N); st = (F0["phase"] == 2) & (F0["since"] > 25)
    print(f"[運動學濾波器] 概似最大化：sig_q = {kn['sig_q']:.3g} rad（編碼器 LSB 1.19e-5），sig_j = {kn['sig_j']:.3g} rad/s³/√Hz")
    kin_chk = []
    for nm, (sq, sj) in (("sig_q 固定為量化誤差 3.4e-6、sig_j = 30", (3.4e-6, 30.0)), ("概似最大化", (kn["sig_q"], kn["sig_j"]))):
        xs = kf.kin_rts(D0["q"], sq, sj)["xs"]
        kin_chk.append((np.sqrt(np.mean((xs[k0, 2] - acmd[k0]) ** 2)), np.std(xs[st, 1] - D0["tq"][st])))
        print(f"  {nm}：â − 指令加速度 RMS {kin_chk[-1][0]:.3f} rad/s²；定速時 v̂ − 指令速度 標準差 {kin_chk[-1][1]:.4f} rad/s"
              f"（控制器的 actual_qd：{np.std(D0['qd'][st] - D0['tq'][st]):.4f}）")

    cols, y, info, G = build(runs)
    print(f"\n回歸樣本數 {len(y)}，運動段 {len(np.unique(info['seg']))}")

    # ---- 1. 逐步加入各項 ----
    print("\n[逐步加入] 殘差 RMS [N·m]；慢變殘差 = 同段 0.2 s 移動平均後的 RMS")
    names, steps = [], []
    for g in ORDER:
        names += G[g]; p, e = lsq(cols, names, y)
        steps.append((np.sqrt(np.mean(e ** 2)), slow_rms(e, info["seg"])))
        print(f"  + {g:6s}: RMS {steps[-1][0]:.4f}，慢變殘差 {steps[-1][1]:.4f}")
    full = list(names); p_full, e_full = lsq(cols, full, y); P = dict(zip(full, p_full))

    # ---- 2. 只拿掉一組 ----
    print("\n[只拿掉一組] 殘差 RMS / 慢變殘差 增加")
    base, base_s = steps[-1]; drop = []
    for g in ORDER[1:]:
        _, e = lsq(cols, [n_ for n_ in full if n_ not in G[g]], y)
        drop.append((np.sqrt(np.mean(e ** 2)) - base, slow_rms(e, info["seg"]) - base_s))
        print(f"  − {g:6s}: +{drop[-1][0]:.4f} / +{drop[-1][1]:.4f}")

    # ---- 3. 係數與 bootstrap ----
    bf = seg_bootstrap(cols, full, y, info["seg"], n_boot=200); sd = dict(zip(full, bf.std(0)))
    print("\n[完整模型係數 ± 以運動段重抽樣的標準差]")
    for n_ in ["J", "B+", "B-", "Tc+", "Tc-", "kB", "kTc", "C+", "C-", "δ加速", "δ減速"]:
        print(f"  {n_:6s} = {P[n_]:+.4f} ± {sd[n_]:.4f}")
    print("  起步後行程 g [N·m]：", "  ".join(f"{k[1:]}°:{P[k]:+.2f}±{sd[k]:.2f}" for k in G["起步後行程"]))
    ph_ = np.linspace(0, 2 * np.pi, 721); rip = {}
    for sg in ("+", "-"):
        rip[sg] = sum(P[f"c{k}{sg}"] * np.cos(k * ph_) + P[f"s{k}{sg}"] * np.sin(k * ph_) for k in range(1, K_RIPPLE + 1))
        amp = {k: np.hypot(P[f"c{k}{sg}"], P[f"s{k}{sg}"]) for k in range(1, K_RIPPLE + 1)}
        top = sorted(amp, key=amp.get, reverse=True)[:5]
        print(f"  漣波（{'正轉' if sg == '+' else '反轉'}）：峰對峰 {np.ptp(rip[sg]):.2f} N·m，RMS {rip[sg].std():.2f}；最大諧波 "
              + "、".join(f"{k} 次 {amp[k]:.2f}" for k in top))
    print(f"  正轉與反轉漣波（運動方向）的相關係數 {np.corrcoef(rip['+'], rip['-'])[0, 1]:+.2f}（−1 表示關節座標中兩方向相同）")
    for sg, lab in (("+", "正轉"), ("-", "反轉")):
        print(f"  速度非線性（{lab}）：局部斜率 B + 2Cv 在 0.02 rad/s = {P['B' + sg] + 2 * P['C' + sg] * 0.02:.1f}，"
              f"0.20 rad/s = {P['B' + sg] + 2 * P['C' + sg] * 0.2:.1f}")

    # ---- 4. J 是否隨運動類型改變 ----
    cj, yj, ij, Gj = build(runs, j_split=True)
    nj = Gj["三參數"] + [n_ for g in ORDER[1:] for n_ in Gj[g]]
    pj, _ = lsq(cj, nj, yj); bj = seg_bootstrap(cj, nj, yj, ij["seg"], n_boot=100)
    p3, _ = lsq(cj, Gj["三參數"], yj)
    jsplit = [(n_, pj[i], bj[:, i].std(), p3[i]) for i, n_ in enumerate(nj[:4])]
    print("\n[J 依運動類型]（完整模型 ± 重抽樣標準差；括號內為只有三參數時）")
    for n_, v_, s_, v3 in jsplit:
        print(f"  {n_:8s} = {v_:.3f} ± {s_:.3f}（{v3:.3f}）")

    # ---- 5. 角度 vs 行程、時間 vs 行程 ----
    off = info["s"] * np.degrees(info["q"] - np.pi) - info["travel"]
    order = np.argsort(info["seg"], kind="stable"); cut = np.flatnonzero(np.diff(info["seg"][order])) + 1
    within = max(x.std() for x in np.split(off[order], cut))
    pa, ea = lsq(cols, full + ["角度"], y); Pa = dict(zip(full + ["角度"], pa))
    print(f"\n[角度 vs 行程] 每段內 s·(q−π) − 行程 的標準差最大 {within:.1e}°（= 常數）；"
          f"加入角度項：RMS {base:.4f} → {np.sqrt(np.mean(ea ** 2)):.4f}，角度係數 {pa[-1]:+.4f} N·m/°")
    print("  加入角度項後的 g：", "  ".join(f"{k[1:]}°:{Pa[k]:+.2f}" for k in G["起步後行程"]))
    tk = np.array([0, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0, 40.0])
    Ht = hat_basis(info["tsince"], tk); ct = dict(cols); tn = []
    for j in range(Ht.shape[1]):
        ct[f"t{j}"] = Ht[:, j]; tn.append(f"t{j}")
    e_t = lsq(ct, [n_ for n_ in full if n_ not in G["起步後行程"]] + tn, y)[1]
    tvt = (np.sqrt(np.mean(e_t ** 2)), slow_rms(e_t, info["seg"]))
    print(f"[時間 vs 行程] 起步後行程 RMS {base:.4f}（慢變 {base_s:.4f}）；改成起步後時間 RMS {tvt[0]:.4f}（慢變 {tvt[1]:.4f}）")

    # ---- 6. 擴增狀態濾波器：名義模型 vs 擴充模型 ----
    h = param_tag(p_full); kfn = lambda r: known_extended(r, P); mext = ext_model(P)
    nz_nom = tuned_noise(runs, NOMINAL, "nom", kn)
    nz_ext = tuned_noise(runs, mext, f"ext_{h}", kn, kfn)
    print("\n[擴增狀態濾波器] sig_q、sig_j 沿用運動學的值，只調 sig_tau、sig_d（概似最大化）")
    print(f"  名義三參數模型：sig_tau {nz_nom['sig_tau']:.3f} N·m，sig_d {nz_nom['sig_d']:.3f} N·m/√s")
    print(f"  擴充模型　　　：sig_tau {nz_ext['sig_tau']:.3f} N·m，sig_d {nz_ext['sig_d']:.3f} N·m/√s")
    print("  兩個模型都用擴充模型的雜訊參數跑全部 12 次實驗（只算使用電流量測的樣本）：")
    kfcmp = []
    for nm, kres in (("三參數名義模型", run_kf(runs, NOMINAL, nz_ext, f"nomX_{h}")), ("擴充模型", run_kf(runs, mext, nz_ext, f"ext_{h}", kfn))):
        kfcmp.append(kf_summary(runs, kres))
        print(f"  {nm}：d̂ RMS {kfcmp[-1][0]:.3f}，慢變 {kfcmp[-1][1]:.3f}；電流通道新息 RMS {kfcmp[-1][2]:.3f}，"
              f"每個樣本的電流通道對數概似 {kfcmp[-1][3]:.3f}")

    np.savez(os.path.join(OUT, "kf_results.npz"), names=np.array(full), p=p_full, sd=np.array([sd[n_] for n_ in full]),
             steps=np.array(steps), drop=np.array(drop), jsplit=np.array([(v_, s__, v3) for _, v_, s__, v3 in jsplit]), hash=h,
             kin=np.array([kn["sig_q"], kn["sig_j"]]), kin_chk=np.array(kin_chk),
             noise_nom=np.array([nz_nom[k] for k in NOISE_KEYS]), noise_ext=np.array([nz_ext[k] for k in NOISE_KEYS]),
             kfcmp=np.array(kfcmp), angle=np.array([np.sqrt(np.mean(ea ** 2)), pa[-1]]), time_vs_travel=np.array(tvt))
