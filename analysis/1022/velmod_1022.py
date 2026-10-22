"""10/22 速度調變實驗分析（DIAG、EXT45 兩個姿態，末端負載 1.3 kg）。

1. 各次、各條件的 J：沿用 robot/ur5_velocity_modulation_ident.py 的 analyze_velocity_modulation。
2. 分析流程檢查：把電流換成控制器的 target_moment（不含摩擦、沒有量測雜訊），用同一套分析重算；
   流程沒有偏差時，各條件的 J 都應等於 J_model。
3. 姿態合併：各條件平均 J、輪間標準差；J 對 1/A 回歸；同加速度振幅、不同頻率的 J 差（1 Hz vs 2 Hz 小振幅），
   同頻率、不同振幅的 J 差（2 Hz vs 2 Hz 小振幅）；兩姿態的 J 差 ΔJ 與控制器模型的比較。
4. 控制器常數：J0 的 target_current 對 target_moment、target_qd、sign(target_qd) 回歸（力矩換電流的比值、摩擦前饋）；
   J1 在兩個姿態的 target_moment、target_current、實測電流。
5. 同日其他資料：定速（速度滑桿未回 100%）實際速度比例、暖機各批結果。

結果存到 OUT/velmod_1022.json，文字輸出存到 OUT/velmod_1022.log。只做離線分析，執行時封鎖網路。
"""
import os, sys, glob, csv, io, json, math, socket, contextlib
import numpy as np

socket.socket = None                                   # 封鎖網路：只做離線分析
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.stdout.reconfigure(encoding="utf-8")
import ur5_velocity_modulation_ident as vm

BASE = r"D:\cluade_code\ur5\data_1022"   # 10/22 資料（各次 velmod_*、constvel_angle_*、warmup_* 子資料夾），請改成自己的路徑
OUT = os.path.join(os.path.dirname(BASE), "analysis_1022")
KT = vm.KT_OUT                           # 本研究的輸出端力矩常數 0.135 × 101 [N·m/A]
COND = ["1 Hz", "2 Hz", "3 Hz", "2 Hz 小振幅"]
NC = len(vm.CONDITIONS)
LOG = io.StringIO()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def run_list():
    """依時間排序的速度調變 FULL：dict(pose, dir, ts, F, label)。"""
    out = []
    for d in glob.glob(os.path.join(BASE, "velmod_*_full_*")):
        _, pose, dr, _, _, ts = os.path.basename(d).split("_")
        out.append(dict(pose=pose.upper(), dir=+1 if dr == "pos" else -1, ts=ts,
                        F=glob.glob(os.path.join(d, "velmod_data_*.csv"))[0]))
    out.sort(key=lambda r: r["ts"])
    seen = {}
    for r in out:
        k = (r["pose"], r["dir"])
        seen[k] = seen.get(k, 0) + 1
        r["label"] = f"{r['pose']} {'正轉' if r['dir'] > 0 else '反轉'}" + ("（重跑）" if seen[k] > 1 else "")
    return out


def analyze(F, direction, kt=KT):
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = vm.analyze_velocity_modulation(F, 0, kt, direction, vm.CONDITIONS, vm.V_CRUISE)
    return r, buf.getvalue()


def moment_as_current(F, tmp):
    """複製 CSV，actual_current_0 換成 target_moment_0 / KT（分析後的力矩即為控制器模型力矩）。"""
    with open(F, encoding="utf-8") as f:
        rows = list(csv.reader(f))
    ic, im = rows[0].index("actual_current_0"), rows[0].index("target_moment_0")
    for r in rows[1:]:
        r[ic] = repr(float(r[im]) / KT)
    with open(tmp, "w", encoding="utf-8", newline="") as f:
        csv.writer(f).writerows(rows)


def controller_constants(F):
    """J0：target_current = target_moment / K + b·target_qd + c·sign(target_qd) + c0（運動中樣本）。J1、J2：平均值。"""
    d = np.genfromtxt(F, delimiter=",", names=True)
    tq, tm, tc = d["target_qd_0"], d["target_moment_0"], d["target_current_0"]
    m = np.abs(tq) > 0.05
    X = np.column_stack([tm[m], tq[m], np.sign(tq[m]), np.ones(m.sum())])
    p, *_ = np.linalg.lstsq(X, tc[m], rcond=None)
    K = 1.0 / p[0]
    out = dict(K0=float(K), B_ff=float(p[1] * K), Tc_ff=float(p[2] * K), c0=float(p[3]),
               resid=float(np.std(tc[m] - X @ p)), T0=float(np.mean(d["joint_temp_0"])))
    for j in (1, 2):
        for key, col in (("tm", "target_moment"), ("tc", "target_current"), ("ac", "actual_current")):
            out[f"{key}{j}"] = float(np.mean(d[f"{col}_{j}"]))
    return out


def wls(X, y, s):
    w = 1.0 / np.asarray(s) ** 2
    cov = np.linalg.inv((X.T * w) @ X)
    p = cov @ ((X.T * w) @ y)
    return p, np.sqrt(np.diag(cov)), float(np.sum(w * (y - X @ p) ** 2))


def main():
    os.makedirs(OUT, exist_ok=True)
    runs = run_list()
    tmp = os.path.join(OUT, "_moment_as_current.csv")
    for r in runs:
        res, txt = analyze(r["F"], r["dir"])
        c = res["conditions"]
        r["J"] = [c[i]["J"] for i in range(NC)]
        r["sJ"] = [c[i]["se_J"] for i in range(NC)]
        r["A"] = [c[i]["a_amp"] for i in range(NC)]
        r["B_local"] = [c[i]["B_local"] for i in range(NC)]
        r["sB_local"] = [c[i]["se_B_local"] for i in range(NC)]
        r["gain"] = [c[i]["track_gain"] for i in range(NC)]
        r["lag_ms"] = [c[i]["track_lag_ms"] for i in range(NC)]
        r["J_spread"] = [c[i]["J_spread"] for i in range(NC)]
        r["n_seg"] = [c[i]["n_seg"] for i in range(NC)]
        r["J_model"] = res["J_model_mean"]
        r["J_inf"], r["sJ_inf"], r["k"], r["sk"] = res["J_inf"], res["se_J_inf"], res["k_invA"], res["se_k_invA"]
        r["delay_ms"], r["sdelay_ms"] = res["delay_ms"], res["se_delay_ms"]
        moment_as_current(r["F"], tmp)
        res_m, _ = analyze(tmp, r["dir"])
        r["J_check"] = [res_m["conditions"][i]["J"] for i in range(NC)]
        r["ctrl"] = controller_constants(r["F"])
        with open(os.path.join(OUT, f"velmod_{r['pose'].lower()}_{r['ts']}_analysis.txt"), "w", encoding="utf-8") as f:
            f.write(txt)
    os.remove(tmp)

    say("=" * 100)
    say("各次、各條件（J ± 單輪回歸的標準差；Kt = 0.135）")
    for r in runs:
        say(f"{r['label']:<14} T0 = {r['ctrl']['T0']:.2f} °C  J_model = {r['J_model']:.3f}")
        for i in range(NC):
            say(f"   {COND[i]:<9} A = {r['A'][i]:.3f}  J = {r['J'][i]:.3f} ± {r['sJ'][i]:.3f}"
                f"（模型寫法最大差 {r['J_spread'][i]:.3f}）  B_local = {r['B_local'][i]:.2f} ± {r['sB_local'][i]:.2f}"
                f"  追隨增益 {r['gain'][i]:.3f}  延遲 {r['lag_ms'][i]:+.1f} ms  流程檢查 J = {r['J_check'][i]:.3f}")
        say(f"   J∞ = {r['J_inf']:.3f} ± {r['sJ_inf']:.3f}，k = {r['k']:+.3f} ± {r['sk']:.3f}"
            f"（減速段多的摩擦 {-math.pi * r['k'] / 2:.2f} N·m），延遲 {r['delay_ms']:+.1f} ± {r['sdelay_ms']:.1f} ms")
        say(f"   同加速度振幅不同頻率（1 Hz − 2 Hz 小振幅）：{r['J'][0] - r['J'][3]:+.3f} ± {math.hypot(r['sJ'][0], r['sJ'][3]):.3f}；"
            f"同頻率振幅加倍（2 Hz − 2 Hz 小振幅）：{r['J'][1] - r['J'][3]:+.3f} ± {math.hypot(r['sJ'][1], r['sJ'][3]):.3f}")

    say("=" * 100)
    say("姿態合併")
    P = {}
    for pose in ("DIAG", "EXT45"):
        rr = [r for r in runs if r["pose"] == pose]
        J = np.array([r["J"] for r in rr]); sJ = np.array([r["sJ"] for r in rr]); A = np.array([r["A"] for r in rr])
        X = np.column_stack([np.ones(J.size), 1.0 / A.ravel()])
        p, s, chi2 = wls(X, J.ravel(), sJ.ravel())
        XL = np.column_stack([np.ones(J.size), A.ravel()])
        pL, sL, chi2L = wls(XL, J.ravel(), sJ.ravel())
        P[pose] = dict(n=len(rr), J=J.mean(0).tolist(), J_sd=J.std(0, ddof=1).tolist(), sJ=sJ.mean(0).tolist(),
                       A=A.mean(0).tolist(), J_inf=float(p[0]), sJ_inf=float(s[0]), k=float(p[1]), sk=float(s[1]),
                       chi2=chi2, dof=J.size - 2, lin=[float(pL[0]), float(pL[1])], slin=[float(sL[0]), float(sL[1])],
                       chi2_lin=chi2L, J_model=float(np.mean([r["J_model"] for r in rr])),
                       J_check=np.mean([r["J_check"] for r in rr], 0).tolist(),
                       J_all=J.tolist(), A_all=A.tolist(), labels=[r["label"] for r in rr])
        q = P[pose]
        say(f"{pose}（{len(rr)} 次）：J_model = {q['J_model']:.3f}")
        for i in range(NC):
            say(f"   {COND[i]:<9} A = {q['A'][i]:.3f}  J = {q['J'][i]:.3f}（輪間標準差 {q['J_sd'][i]:.3f}，"
                f"單輪回歸的標準差 {q['sJ'][i]:.3f}）  J ÷ J_model = {q['J'][i] / q['J_model']:.3f}  流程檢查 {q['J_check'][i]:.3f}")
        say(f"   J = J∞ + k/A：J∞ = {q['J_inf']:.3f} ± {q['sJ_inf']:.3f}，k = {q['k']:+.3f} ± {q['sk']:.3f}，"
            f"卡方 {q['chi2']:.0f}（自由度 {q['dof']}）；J = a + b·A：a = {q['lin'][0]:.3f}，b = {q['lin'][1]:+.3f}，卡方 {q['chi2_lin']:.0f}")
    D, E = P["DIAG"], P["EXT45"]
    dJm = E["J_model"] - D["J_model"]
    dJ = [E["J"][i] - D["J"][i] for i in range(NC)]
    pairs = [(e, d) for e in E["J_all"] for d in D["J_all"]]
    dJ_rng = [(min(e[i] - d[i] for e, d in pairs), max(e[i] - d[i] for e, d in pairs)) for i in range(NC)]
    dJinf = E["J_inf"] - D["J_inf"]
    say(f"兩姿態差 ΔJ（EXT45 − DIAG）；控制器模型 ΔJ_model = {dJm:.3f}")
    for i in range(NC):
        say(f"   {COND[i]:<9} ΔJ = {dJ[i]:.3f}（各次組合 {dJ_rng[i][0]:.3f}～{dJ_rng[i][1]:.3f}）  ÷ ΔJ_model = {dJ[i] / dJm:.3f}")
    say(f"   J∞        ΔJ = {dJinf:.3f}  ÷ ΔJ_model = {dJinf / dJm:.3f}")
    cs = [dJ[i] / dJm for i in range(NC)] + [dJinf / dJm]
    say("   直線 J = c·J_model + d（兩姿態）：" + "，".join(
        f"{n} c = {c:.3f}、d = {(D['J'][i] if i < NC else D['J_inf']) - c * D['J_model']:+.3f}"
        for i, (n, c) in enumerate(zip(COND + ['J∞'], cs))))

    say("=" * 100)
    say("控制器常數")
    K0 = float(np.mean([r["ctrl"]["K0"] for r in runs]))
    for r in runs:
        c = r["ctrl"]
        say(f"{r['label']:<14} J0：target_current = target_moment / {c['K0']:.3f} + 摩擦前饋"
            f"（黏滯 {c['B_ff']:.2f} N·m·s/rad、庫倫 {c['Tc_ff']:.2f} N·m，以 {c['K0']:.2f} 換算），殘差 {c['resid']:.4f} A"
            f" | J1：target_moment {c['tm1']:+.2f}，target_current {c['tc1']:+.3f}，實測 {c['ac1']:+.3f} A"
            f" | J2：{c['tm2']:+.2f}、{c['tc2']:+.3f}、{c['ac2']:+.3f}")
    say(f"J0 力矩換電流比值 K0 = {K0:.3f} N·m/A（本研究 KT = {KT:.3f}，比值 {KT / K0:.3f}）")
    cD = [r["ctrl"] for r in runs if r["pose"] == "DIAG"]; cE = [r["ctrl"] for r in runs if r["pose"] == "EXT45"]
    tmD, tmE = np.mean([c["tm1"] for c in cD]), np.mean([c["tm1"] for c in cE])
    tcD, tcE = np.mean([c["tc1"] for c in cD]), np.mean([c["tc1"] for c in cE])
    K1 = (tmE - tmD) / (tcE - tcD)
    acE = np.mean([c["ac1"] for c in cE])
    K1_act = [(tmE - tmD) / (acE - c["ac1"]) for c in cD]
    say(f"J1 控制器比值 K1 = Δtarget_moment / Δtarget_current = {K1:.3f}；J2 比值 {np.mean([c['tm2'] / c['tc2'] for c in cD]):.3f}")
    ac_txt = ", ".join(f"{c['ac1']:+.3f}" for c in cD)
    k_txt = ", ".join(f"{k:.2f}" for k in K1_act)
    say(f"J1 實測電流：DIAG 各次 {ac_txt}，EXT45 平均 {acE:+.3f} A；"
        f"Δtarget_moment ÷ Δ實測電流 = {k_txt} N·m/A（靜摩擦使靜止電流隨進入方向改變）")
    say(f"以控制器 K0 換算時 ΔJ ÷ ΔJ_model：" + "，".join(f"{n} {c * K0 / KT:.3f}" for n, c in zip(COND + ['J∞'], cs)))
    say(f"兩姿態差要等於模型時的 K（= KT × ΔJ_model ÷ ΔJ）：" + "，".join(f"{n} {KT / c:.2f}" for n, c in zip(COND + ['J∞'], cs)))

    say("=" * 100)
    say("同日其他資料")
    other = {}
    for d in sorted(glob.glob(os.path.join(BASE, "constvel_angle_*_full_*"))):
        F = glob.glob(os.path.join(d, "*_data_*.csv"))[0]
        x = np.genfromtxt(F, delimiter=",", names=True, usecols=(0, 1, 2, 3, 4, 5))
        ratio = np.max(np.abs(x["target_qd_0"])) / 0.20
        other[os.path.basename(d)] = dict(speed_ratio=float(ratio), dur=float(x["timestamp"][-1] - x["timestamp"][0]),
                                         T=[float(x["joint_temp_0"][0]), float(x["joint_temp_0"][-1])])
        say(f"{os.path.basename(d)}：最高目標速度 = 規劃 0.20 rad/s 的 {ratio:.2f} 倍，歷時 {other[os.path.basename(d)]['dur']:.0f} s，"
            f"J0 {x['joint_temp_0'][0]:.2f} → {x['joint_temp_0'][-1]:.2f} °C")
    for d in glob.glob(os.path.join(BASE, "warmup_full_*")):
        info = open(glob.glob(os.path.join(d, "warmup_info_*.txt"))[0], encoding="utf-8").read()
        blocks = [ln.split(",") for ln in info.splitlines() if ln[:1].isdigit() and ln.count(",") == 4]
        other["warmup"] = [[float(v) for v in b] for b in blocks]
        say("暖機各批（批, 累計 min, 溫度, 正轉|τ|, 反轉|τ|）：" + "；".join(",".join(b) for b in blocks))

    res = dict(KT=KT, K0=K0, K1=float(K1), K1_act=[float(k) for k in K1_act], dJ_model=dJm, dJ=dJ, dJ_rng=dJ_rng,
               dJ_inf=dJinf, poses=P, other=other, cond=COND,
               runs=[{k: v for k, v in r.items() if k != "F"} for r in runs])
    with open(os.path.join(OUT, "velmod_1022.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1)
    with open(os.path.join(OUT, "velmod_1022.log"), "w", encoding="utf-8") as f:
        f.write(LOG.getvalue())


if __name__ == "__main__":
    main()
