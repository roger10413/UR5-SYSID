"""1.915 為什麼在電流中看不到：用 10/29 的速度調變資料做兩個離線檢查（不用上機）。

A. 換訊號：同樣的回歸（沿用 robot/ur5_velocity_modulation_ident.py 的完整模型：J·a + B·v + C·v² + 漣波 + 各段常數 + 各段行程斜率），
   力矩改用 4 種訊號 × 13.64：actual_current（實測電流）、target_current、joint_control_output、target_moment（不乘 13.64）。
   反轉 4 次（2 姿態 × 夾或不夾）每個條件擬合 J = c·X + J_d（X = 官方點質量模型，取自 known_mass_1029.json）。
   J_d ÷ c 約 1.915 → 這個訊號含驅動端慣量；約 0 → 不含。target_moment 是控制器的模型，應約 1.915（檢查方法本身）。
B. 分開兩種加速度：actual_current 的回歸多加一欄目標加速度 target_qdd：τ = J_a·a（實際）+ J_t·target_qdd + …
   推測：若回報電流前先用目標加速度扣掉了轉子項，J_t ≈ −1.915·c、J_a ≈ (X + 1.915)·c；若電流單純不含轉子，J_t ≈ 0。
   實際與目標加速度只差追隨增益 1.00～1.09 與延遲，J_t 的標準差會比 J 大很多，看標準差再判斷。
只做離線分析，執行時封鎖網路。結果存到 OUT/current_signals_1104.log 與 .json。
"""
import os, sys, io, csv, json, socket
import numpy as np

socket.socket = None                                   # 封鎖網路：只做離線分析
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.path.insert(0, os.path.join(HERE, "..", "1029"))
sys.stdout.reconfigure(encoding="utf-8")
import ur5_velocity_modulation_ident as vm
import known_mass_1029 as km

OUT = os.path.join(os.path.dirname(km.BASE), "analysis_1104"); os.makedirs(OUT, exist_ok=True)
# (訊號, 換算, 是否為控制器的目標值)：目標值的 X 用教導器設定的負載（控制器只知道設定值；夾砝碼時設 1.40、實際 2.30 kg）
SIGNALS = [("actual_current", km.KT, False), ("target_current", km.KT, True), ("joint_control_output", km.KT, False),
           ("target_moment", 1.0, True)]
K_RIPPLE = 24
LOG = io.StringIO()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def load(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) if r[k] not in ("", "None") else np.nan for r in rows]) for k in rows[0]}


def fit_conditions(D, s, regressors):
    """
    每個條件：y = Σ regressors·係數 + B·v + C·v² + 漣波 + 各段常數 + 各段行程斜率（同 analyze_velocity_modulation 的「完整」寫法）。
    regressors：{名稱: 函式(段) -> 欄}；y 欄由 regressors["_y"] 給。回傳 {條件: {名稱: (係數, 標準差)}}。
    """
    t, q = D["timestamp"], D["actual_q_0"]
    tq = s * D["target_qd_0"]
    segs, _ = vm.find_modulation_segments(t, tq, vm.V_CRUISE, vm.CONDITIONS, vm.DT)
    for sg in segs:
        ii = sg["idx"]
        sg["v"], sg["a"], *_ = vm._fit_kinematics(t[ii], s * q[ii], sg["freq"])
    names_x = [n for n in regressors if n != "_y"]
    out = {}
    for ci in range(len(vm.CONDITIONS)):
        cs = [sg for sg in segs if sg["cond"] == ci]
        if len(cs) < 2:
            continue
        ng = len(cs)
        Xs, ys = [], []
        for g, sg in enumerate(cs):
            ii = sg["idx"]
            phi = vm.GEAR_RATIO * q[ii]
            cols = [regressors[n](sg) for n in names_x] + [sg["v"], sg["v"] ** 2]
            for k in range(1, K_RIPPLE + 1):
                cols += [np.cos(k * phi), np.sin(k * phi)]
            one, trv, zero = np.ones(len(ii)), s * (q[ii] - q[ii[0]]), np.zeros(len(ii))
            cols += [one if h == g else zero for h in range(ng)] + [trv if h == g else zero for h in range(ng)]
            Xs.append(np.column_stack(cols))
            ys.append(regressors["_y"](sg))
        XtX = sum(X.T @ X for X in Xs)
        Xty = sum(X.T @ y for X, y in zip(Xs, ys))
        p = np.linalg.solve(XtX, Xty)
        es = [y - X @ p for X, y in zip(Xs, ys)]
        cov = vm._hac_cov(Xs, es, XtX, int(round(0.25 / vm.DT)))
        out[ci] = {n: (float(p[k]), float(np.sqrt(cov[k, k]))) for k, n in enumerate(names_x)}
    return out


def wls_line(X, y, s):
    """y = c·X + J_d 的加權最小平方；回傳 (c, J_d, 共變異)（卡方 ÷ 自由度 > 1 時放大）。"""
    A = np.column_stack([X, np.ones(len(X))])
    p, cov, chi2, dof = km.wls(A, np.asarray(y), np.asarray(s))
    return p[0], p[1], cov


def ratio_se(c, d, cov):
    """J_d ÷ c 的標準差（一階誤差傳播）。"""
    g = np.array([-d / c ** 2, 1 / c])
    return float(np.sqrt(g @ cov @ g))


def main():
    R = json.load(open(os.path.join(km.OUT, "known_mass_1029.json"), encoding="utf-8"))
    runs = {r["ts"]: r for r in R["runs"]}
    res = {}
    for r in km.velmod_runs():
        D = load(r["F"])
        s = float(r["dir"])
        reg = {}
        for name, k, _ in SIGNALS:
            col = s * D[f"{name}_0"] * k
            reg[name] = fit_conditions(D, s, {"J": lambda sg: sg["a"], "_y": lambda sg, c=col: c[sg["idx"]]})
        cur = s * D["actual_current_0"] * km.KT
        tqdd = s * D["target_qdd_0"]
        reg["split"] = fit_conditions(D, s, {"J_a": lambda sg: sg["a"], "J_t": lambda sg: tqdd[sg["idx"]],
                                             "_y": lambda sg: cur[sg["idx"]]})
        rr = runs[r["ts"]]
        X_ctrl = rr["X_link"] + rr["X_tool"] * rr["pay"] / km.TOOL_KG          # 控制器看到的官方模型（負載 = 教導器設定）
        res[r["ts"]] = dict(pose=r["pose"], w=r["w"], dir=r["dir"], X=rr["X"], X_ctrl=X_ctrl, label=r["label"], reg=reg)
        a, b = D["actual_current_0"], D["joint_control_output_0"]
        ok = np.isfinite(a) & np.isfinite(b)
        g, o = np.polyfit(b[ok], a[ok], 1)
        res[r["ts"]]["jco"] = dict(gain=float(g), offset=float(o), resid_rms=float(np.std(a[ok] - (g * b[ok] + o))),
                                   std_actual=float(np.std(a[ok])), corr=float(np.corrcoef(a[ok], b[ok])[0, 1]))
        chk = max(abs(reg["actual_current"][ci]["J"][0] - runs[r["ts"]]["J"][ci]) for ci in reg["actual_current"])
        say(f"{r['label']}：X = {runs[r['ts']]['X']:.3f}；actual_current 的 J 和 known_mass_1029 最大差 {chk:.1e}")

    # ---------------- A. 換訊號 ----------------
    rev = [v for v in res.values() if v["dir"] < 0]
    say("\n" + "=" * 90)
    say("A. 各訊號在反轉 4 次的 J = c·X + J_d（每個條件分開擬合；X = 官方點質量模型；控制器目標值用教導器設定的負載）")
    say("   J_d ÷ c 約 1.915 → 含驅動端慣量；約 0 → 不含")
    say("=" * 90)
    out_a = {}
    for name, _, ctrl in SIGNALS:
        X = np.array([v["X_ctrl" if ctrl else "X"] for v in rev])
        say(f"\n[{name}]" + ("（× 13.64）" if name != "target_moment" else "（N·m，不換算）")
            + ("；X 用教導器設定的負載" if ctrl else "；X 用實際負載"))
        say(f"  {'條件':<10}" + "".join(f"{v['pose'] + (' 夾' if v['w'] else ' 不夾'):>11}" for v in rev)
            + f"{'c':>9}{'J_d ÷ c':>16}")
        say(f"  {'X':<10}" + "".join(f"{x:11.3f}" for x in X))
        out_a[name] = []
        for ci, cn in enumerate(km.COND):
            J = np.array([v["reg"][name][ci]["J"][0] for v in rev])
            sJ = np.array([v["reg"][name][ci]["J"][1] for v in rev])
            c, d, cov = wls_line(X, J, np.maximum(sJ, 1e-6))
            say(f"  {cn:<10}" + "".join(f"{j:11.3f}" for j in J) + f"{c:9.3f}{d / c:+9.3f} ± {ratio_se(c, d, cov):.3f}")
            out_a[name].append(dict(J=J.tolist(), c=c, Jd_over_c=d / c, se=ratio_se(c, d, cov)))

    say("\n  actual_current 對 joint_control_output（整段資料）：actual = 斜率 × jco + 截距")
    for v in res.values():
        j = v["jco"]
        say(f"    {v['label']:<16} 斜率 {j['gain']:.4f}  截距 {j['offset']:+.4f} A  相關係數 {j['corr']:.5f}  "
            f"殘差 RMS {j['resid_rms']:.4f} A（actual 本身的標準差 {j['std_actual']:.4f} A）")

    # ---------------- B. 分開兩種加速度 ----------------
    say("\n" + "=" * 90)
    say("B. actual_current：τ = J_a·a（實際加速度）+ J_t·target_qdd + …（各次、各條件）")
    say("   推測：電流先用目標加速度扣掉轉子項 → J_t ≈ −1.915 × c（約 −2.4）；電流單純不含轉子 → J_t ≈ 0")
    say("=" * 90)
    out_b = []
    for v in res.values():
        say(f"  {v['label']}（X = {v['X']:.3f}）")
        for ci, cn in enumerate(km.COND):
            if ci not in v["reg"]["split"]:
                continue
            (ja, sja), (jt, sjt) = v["reg"]["split"][ci]["J_a"], v["reg"]["split"][ci]["J_t"]
            j1 = v["reg"]["actual_current"][ci]["J"][0]
            say(f"    {cn:<10} J_a {ja:+7.3f} ± {sja:.3f}   J_t {jt:+7.3f} ± {sjt:.3f}   J_a + J_t {ja + jt:6.3f}（單一 J {j1:.3f}）")
            out_b.append(dict(label=v["label"], cond=cn, J_a=ja, se_J_a=sja, J_t=jt, se_J_t=sjt))
    # 合併：各次各條件的 J_t 依標準差加權平均
    jt = np.array([b["J_t"] for b in out_b]); sjt = np.array([b["se_J_t"] for b in out_b])
    w = 1 / sjt ** 2
    m = float(np.sum(w * jt) / np.sum(w)); sm = float(np.sqrt(1 / np.sum(w)))
    chi2 = float(np.sum(w * (jt - m) ** 2))
    say(f"\n  J_t 全部 {len(jt)} 個加權平均 = {m:+.3f} ± {sm:.3f}（卡方 {chi2:.1f}，自由度 {len(jt) - 1}）")

    with open(os.path.join(OUT, "current_signals_1104.log"), "w", encoding="utf-8") as fh:
        fh.write(LOG.getvalue())
    with open(os.path.join(OUT, "current_signals_1104.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(A=out_a, B=out_b, J_t_mean=m, se_J_t_mean=sm), fh, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
