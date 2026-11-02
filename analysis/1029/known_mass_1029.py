"""10/29 已知質量實驗分析（夾爪橫夾 1 kg 圓柱砝碼；DIAG、EXT45 兩個姿態；末端工具 FT300、夾爪、3D 相機 1.30 kg）。

1. 速度調變各次、各條件的 J：沿用 robot/ur5_velocity_modulation_ident.py 的 analyze_velocity_modulation（Kt = 0.135）。
2. 官方點質量模型 X：UR 官方 DH 頁面的連桿質量與重心（點質量，不含連桿自身慣量）+ 末端工具 1.30 kg 點質量
   （重心 (1, 1, 74) mm）+ 砝碼 m·r²。姿態用每次資料的實測平均角度。
3. 控制器模型：J_model − 官方連桿點質量 − 教導器負載點質量 = 剩下的常數（Gnad 2026 表 3 的驅動端慣量為 1.915）。
   這次夾砝碼時教導器負載設為 1.40 kg（實際 2.30 kg），ΔJ_model 應為 0.10·r²。
4. 夾砝碼 − 不夾（同姿態、同方向）：ΔJ ÷ m·r²，以及換算比例 K = KT × m·r² ÷ ΔJ。
   另列只用姿態變化的比例 (J_EXT45 − J_DIAG) ÷ (X_EXT45 − X_DIAG)。
5. 合併擬合（反轉 4 次 = 2 姿態 × 夾或不夾，每個條件分開）：J = c·X + J_d。
   c 為實測對官方點質量模型的比例（K = KT ÷ c）；J_d ÷ c 為不隨姿態、負載改變的慣量
   （控制器的驅動端慣量若在電流中，應約為 1.915）。另以 J = c·(X + 1.915) 擬合比較卡方，
   並以 X 加上連桿自身的轉動慣量（CB3 未公開，借用 Clochiatti 2024 的 UR5e 值）重算作對照。
   係數標準差在卡方 ÷ 自由度 > 1 時乘上 √(卡方 ÷ 自由度)。
6. J1 靜止電流：夾砝碼 − 不夾 → K（兩方向靠近的中心電流）；並列各靠近方向單獨算的 K 與靜摩擦帶寬。
7. 與 10/22 比較（不夾，前 4 個條件；需先執行 ../1022/velmod_1022.py）。
8. 同日定速 B、Tc（固定角度窗口），並列 10/13、10/15 各輪對照。

砝碼質量先用標稱 1.000 kg（未秤重；可用第一個參數改，例如 python known_mass_1029.py 1.002）。
砝碼重心假設在工具軸線上：工具垂直朝下時 r = 法蘭中心離 J0 軸的水平距離，與離法蘭多遠無關（另列工具軸的實際傾角）。
結果存到 OUT/known_mass_1029.json，文字輸出存到 OUT/known_mass_1029.log。只做離線分析，執行時封鎖網路。
"""
import os, sys, glob, io, json, math, socket, contextlib
import numpy as np

socket.socket = None                                   # 封鎖網路：只做離線分析
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.stdout.reconfigure(encoding="utf-8")
import ur5_velocity_modulation_ident as vm
import ur5_j1_static_current as js
import ur5_const_velocity_angle_ident as cv

BASE = r"D:\cluade_code\ur5\data_1029"   # 10/29 資料（各次 velmod_*、j1static_*、constvel_angle_* 子資料夾），請改成自己的路徑
OUT = os.path.join(os.path.dirname(BASE), "analysis_1029")
R1022 = os.path.join(os.path.dirname(BASE), "analysis_1022", "velmod_1022.json")
CV_OLD = os.path.dirname(BASE)            # 10/13、10/15 定速資料（constvel_angle_*_full_2026101[35]_*）
KT = vm.KT_OUT                           # 本研究的輸出端力矩常數 0.135 × 101 [N·m/A]
M_W = float(sys.argv[1]) if len(sys.argv) > 1 else 1.000   # 砝碼質量 [kg]
TOOL_KG, TOOL_COG = 1.30, np.array([0.001, 0.001, 0.074])  # 末端工具（教導器負載設定）質量與重心（法蘭座標）
J_DRIVE = 1.915                          # Gnad 2026 表 3：UR5e J1～J3 驅動端慣量（設定檔轉子慣量 × 101²）
COND = ["1 Hz×0.06", "2 Hz×0.06", "3 Hz×0.06", "2 Hz×0.03", "1 Hz×0.03", "3 Hz×0.04"]
NC = len(vm.CONDITIONS)
# UR 官方 DH 頁面（UR5，CB3）：各連桿質量 [kg] 與重心（各連桿 DH 座標系）[m]
LINK_M = [3.7, 8.393, 2.33, 1.219, 1.219, 0.1879]
LINK_COM = [[0, -0.02561, 0.00193], [0.2125, 0, 0.11336], [0.15, 0, 0.0265],
            [0, -0.0018, 0.01634], [0, 0.0018, 0.01634], [0, 0, -0.001159]]
# 連桿自身的轉動慣量（重心處，各連桿 DIAG 座標系）[kg·m²]：CB3 官方未公開，借用 Clochiatti 2024 的 UR5e 值，只作對照
LINK_IC = [np.diag([0.0067, 0.0064, 0.0067]), np.diag([0.0149, 0.3564, 0.3553]),
           np.array([[0.0025, 0, 0.0034], [0, 0.0551, 0], [0.0034, 0, 0.0546]]),
           np.diag([0.0012, 0.0012, 0.0009]), np.diag([0.0012, 0.0012, 0.009]), np.diag([0.001, 0.001, 0.001])]
LOG = io.StringIO()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def point_mass_model(q, pay_kg):
    """J0 轉軸的點質量慣量：官方連桿、教導器負載（pay_kg，重心 TOOL_COG）；並回傳法蘭離 J0 軸的水平距離與工具軸傾角。"""
    fr = js.fk_frames(q)
    j_link = 0.0
    for i in range(6):
        c = fr[i + 1][:3, :3] @ np.array(LINK_COM[i]) + fr[i + 1][:3, 3]
        j_link += LINK_M[i] * (c[0] ** 2 + c[1] ** 2)
    cp = fr[6][:3, :3] @ TOOL_COG + fr[6][:3, 3]
    r_fl = float(np.hypot(*fr[6][:2, 3]))
    tilt = math.degrees(math.acos(min(1.0, -fr[6][2, 2])))          # 工具 z 軸與鉛直向下的夾角
    return j_link, pay_kg * (cp[0] ** 2 + cp[1] ** 2), r_fl, tilt


def link_self_inertia(q):
    """各連桿自身轉動慣量（LINK_IC，借用 UR5e 值）對 J0 轉軸的總和。"""
    fr = js.fk_frames(q)
    return float(sum((fr[i + 1][:3, :3] @ LINK_IC[i] @ fr[i + 1][:3, :3].T)[2, 2] for i in range(6)))


def wls(X, y, s):
    """加權最小平方；回傳係數、共變異（卡方 ÷ 自由度 > 1 時放大）、卡方、自由度。"""
    w = 1.0 / np.asarray(s) ** 2
    cov = np.linalg.inv((X.T * w) @ X)
    p = cov @ ((X.T * w) @ y)
    chi2 = float(np.sum(w * (y - X @ p) ** 2))
    dof = len(y) - X.shape[1]
    if dof > 0 and chi2 / dof > 1:
        cov = cov * chi2 / dof
    return p, cov, chi2, dof


def velmod_runs():
    """依時間排序的速度調變 FULL：資料夾 velmod_{pose}_p{kg×100}_{pos|neg}_full_{date}_{time}。"""
    runs = []
    for d in glob.glob(os.path.join(BASE, "velmod_*_full_*")):
        _, pose, pay, dr, _, day, tm = os.path.basename(d).split("_")
        pay = int(pay[1:]) / 100
        runs.append(dict(pose=pose.upper(), pay=pay, dir=+1 if dr == "pos" else -1, ts=f"{day}_{tm}",
                         w=pay > TOOL_KG + 0.005, F=glob.glob(os.path.join(d, "velmod_data_*.csv"))[0]))
    runs.sort(key=lambda r: r["ts"])
    for r in runs:
        r["label"] = f"{r['pose']} {'夾砝碼' if r['w'] else '不夾'} {'正轉' if r['dir'] > 0 else '反轉'}"
    return runs


def main():
    os.makedirs(OUT, exist_ok=True)
    say(f"砝碼質量 m = {M_W:.3f} kg；KT = {KT:.3f} N·m/A；驅動端慣量對照值 {J_DRIVE}")

    # ---------------- 1～3. 速度調變各次 ----------------
    say("=" * 100)
    say("速度調變各次（J ± 單輪回歸的標準差；Kt = 0.135）")
    runs = velmod_runs()
    for r in runs:
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            res = vm.analyze_velocity_modulation(r["F"], 0, KT, r["dir"], vm.CONDITIONS, vm.V_CRUISE)
        with open(os.path.join(OUT, f"velmod_{r['pose'].lower()}_p{round(r['pay'] * 100)}_{r['ts']}_analysis.txt"),
                  "w", encoding="utf-8") as f:
            f.write(buf.getvalue())
        c = res["conditions"]
        r["J"] = [c[i]["J"] for i in range(NC)]
        r["sJ"] = [c[i]["se_J"] for i in range(NC)]
        r["A"] = [c[i]["a_amp"] for i in range(NC)]
        r["J_model"] = res["J_model_mean"]
        r["J_inf"], r["k"] = res["J_inf"], res["k_invA"]
        d = np.genfromtxt(r["F"], delimiter=",", names=True)
        r["T"] = [float(d["joint_temp_0"][0]), float(d["joint_temp_0"][-1])]
        q = [float(np.mean(d[f"actual_q_{k}"])) for k in range(6)]
        r["q_deg"] = [math.degrees(v) for v in q]
        r["X_link"], r["X_tool"], r["r_fl"], r["tilt"] = point_mass_model(q, TOOL_KG)
        r["X_w"] = M_W * r["r_fl"] ** 2 if r["w"] else 0.0
        r["X"] = r["X_link"] + r["X_tool"] + r["X_w"]
        r["X_self"] = link_self_inertia(q)
        j_link_c, j_pay_c, _, _ = point_mass_model(q, r["pay"])          # 控制器看到的負載 = 教導器設定
        r["drive"] = r["J_model"] - j_link_c - j_pay_c
        del d
        say(f"{r['label']:<16} {r['ts']}  J0 {r['T'][0]:.2f}→{r['T'][1]:.2f} °C  教導器負載 {r['pay']:.2f} kg"
            f"  法蘭離 J0 軸 {r['r_fl']:.4f} m  工具軸傾角 {r['tilt']:.2f}°")
        say(f"   官方點質量：連桿 {r['X_link']:.4f} + 末端工具 {r['X_tool']:.4f} + 砝碼 {r['X_w']:.4f} = X {r['X']:.4f}；"
            f"控制器 J_model {r['J_model']:.4f} − 連桿 − 教導器負載 {j_pay_c:.4f} = {r['drive']:.4f}")
        for i in range(NC):
            say(f"   {COND[i]:<10} A = {r['A'][i]:.3f}  J = {r['J'][i]:.4f} ± {r['sJ'][i]:.4f}  J ÷ X = {r['J'][i] / r['X']:.3f}")

    # ---------------- 4. 夾砝碼 − 不夾 ----------------
    say("=" * 100)
    say(f"夾砝碼 − 不夾（同姿態、同方向）；m·r²：m = {M_W:.3f} kg，r = 法蘭離 J0 軸水平距離")
    pairs = []
    for pose in ("DIAG", "EXT45"):
        for dr in (+1, -1):
            w = [r for r in runs if r["pose"] == pose and r["dir"] == dr and r["w"]]
            o = [r for r in runs if r["pose"] == pose and r["dir"] == dr and not r["w"]]
            if not (w and o):
                continue
            w, o = w[0], o[0]
            mr2 = M_W * ((w["r_fl"] + o["r_fl"]) / 2) ** 2
            dJ = np.array(w["J"]) - np.array(o["J"])
            sdJ = np.hypot(w["sJ"], o["sJ"])
            wt = 1 / sdJ ** 2
            dJw, sdw = float(np.sum(wt * dJ) / np.sum(wt)), float(1 / math.sqrt(np.sum(wt)))
            dJm = w["J_model"] - o["J_model"]
            say(f"{pose} {'正轉' if dr > 0 else '反轉'}：m·r² = {mr2:.4f}；控制器 ΔJ_model = {dJm:.4f}"
                f"（教導器負載差 {w['pay'] - o['pay']:.2f} kg × r² = {(w['pay'] - o['pay']) * mr2 / M_W:.4f}）")
            for i in range(NC):
                say(f"   {COND[i]:<10} ΔJ = {dJ[i]:.4f} ± {sdJ[i]:.4f}  ÷ m·r² = {dJ[i] / mr2:.3f} ± {sdJ[i] / mr2:.3f}"
                    f"  → K = {KT * mr2 / dJ[i]:.2f} N·m/A")
            say(f"   加權平均 ΔJ = {dJw:.4f} ± {sdw:.4f}  ÷ m·r² = {dJw / mr2:.3f}  → K = {KT * mr2 / dJw:.2f} N·m/A")
            pairs.append(dict(pose=pose, dir=dr, mr2=mr2, dJ=dJ.tolist(), sdJ=sdJ.tolist(),
                              ratio=(dJ / mr2).tolist(), sratio=(sdJ / mr2).tolist(),
                              dJw=dJw, sdw=sdw, ratio_w=dJw / mr2, sratio_w=sdw / mr2, dJ_model=dJm))
    say("只用姿態變化（反轉，同負載）：(J_EXT45 − J_DIAG) ÷ (X_EXT45 − X_DIAG)")
    pose_ratio = []
    for wflag in (False, True):
        D = [r for r in runs if r["pose"] == "DIAG" and r["dir"] < 0 and r["w"] == wflag]
        E = [r for r in runs if r["pose"] == "EXT45" and r["dir"] < 0 and r["w"] == wflag]
        if not (D and E):
            continue
        D, E = D[0], E[0]
        dX = E["X"] - D["X"]
        rr = [(E["J"][i] - D["J"][i]) / dX for i in range(NC)]
        sr = [math.hypot(E["sJ"][i], D["sJ"][i]) / dX for i in range(NC)]
        pose_ratio.append(dict(w=wflag, dX=dX, ratio=rr, sratio=sr))
        say(f"   {'夾砝碼' if wflag else '不夾'}：ΔX = {dX:.4f}；" + "，".join(f"{COND[i]} {rr[i]:.3f}" for i in range(NC)))

    # ---------------- 5. 合併擬合 ----------------
    say("=" * 100)
    say("合併擬合（反轉 4 次）：J = c·X + J_d；對照 J = c·(X + 1.915)")
    rev = [r for r in runs if r["dir"] < 0]
    say("   使用：" + "、".join(f"{r['label']}（X = {r['X']:.3f}）" for r in rev))
    Xv = np.array([r["X"] for r in rev])
    fit = []
    for i in range(NC):
        y = np.array([r["J"][i] for r in rev]); s = np.array([r["sJ"][i] for r in rev])
        p, cov, chi2, dof = wls(np.column_stack([Xv, np.ones_like(Xv)]), y, s)
        c, Jd = p
        g = np.array([-Jd / c ** 2, 1 / c])                      # J_d ÷ c 的梯度
        s_ratio = float(math.sqrt(g @ cov @ g))
        p1, cov1, chi2_1, dof1 = wls((Xv + J_DRIVE)[:, None], y, s)
        Xs = Xv + np.array([r["X_self"] for r in rev])                # 加上連桿自身慣量（UR5e 值）
        p2, cov2, chi2_2, _ = wls(np.column_stack([Xs, np.ones_like(Xs)]), y, s)
        resid = y - (c * Xv + Jd)
        fit.append(dict(c=float(c), sc=float(math.sqrt(cov[0, 0])), Jd=float(Jd), sJd=float(math.sqrt(cov[1, 1])),
                        Jd_c=float(Jd / c), sJd_c=s_ratio, chi2=chi2, dof=dof, rms=float(np.sqrt(np.mean(resid ** 2))),
                        c_drive=float(p1[0]), chi2_drive=chi2_1, dof_drive=dof1, K=float(KT / c), sK=float(KT * math.sqrt(cov[0, 0]) / c ** 2),
                        c_self=float(p2[0]), sc_self=float(math.sqrt(cov2[0, 0])), Jd_c_self=float(p2[1] / p2[0]), chi2_self=chi2_2))
        f = fit[-1]
        say(f"   {COND[i]:<10} c = {c:.3f} ± {f['sc']:.3f}（K = {f['K']:.2f} ± {f['sK']:.2f}）  J_d = {Jd:+.3f} ± {f['sJd']:.3f}"
            f"  J_d ÷ c = {Jd / c:+.3f} ± {s_ratio:.3f}  殘差 RMS {f['rms']:.4f}  卡方 {chi2:.1f}（自由度 {dof}）"
            f" | 含 1.915：c = {p1[0]:.3f}，卡方 {chi2_1:.0f}（自由度 {dof1}）"
            f" | 加連桿自身慣量：c = {p2[0]:.3f} ± {f['sc_self']:.3f}，J_d ÷ c = {p2[1] / p2[0]:+.3f}，卡方 {chi2_2:.1f}")
    pw = [p for p in pairs if p["dir"] < 0]
    if len(pw) == 2:
        say("   只用砝碼（反轉兩組合併）÷ m·r²：" + "，".join(
            f"{COND[i]} {np.average([p['ratio'][i] for p in pw], weights=[p['sratio'][i] ** -2 for p in pw]):.3f}"
            for i in range(NC)))

    # ---------------- 6. J1 靜止電流 ----------------
    say("=" * 100)
    say("J1 靜止電流（FULL；夾砝碼 − 不夾）")
    J1 = {}
    for d in sorted(glob.glob(os.path.join(BASE, "j1static_*_full_*"))):
        _, pose, pay, _, _, _ = os.path.basename(d).split("_")
        F = glob.glob(os.path.join(d, "j1static_data_*.csv"))[0]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            res = js.analyze_static_current(F)
        say(f"--- {os.path.basename(d)}\n" + buf.getvalue().rstrip())
        for sgn, key in ((+1, "pos"), (-1, "neg")):
            ii = [h["i"] for h in res["holds"] if h["dir"] == sgn]
            res[f"sd_{key}"] = float(np.std(ii, ddof=1) / math.sqrt(len(ii)))
        J1[(pose.upper(), int(pay[1:]) > round(TOOL_KG * 100) + 0.5)] = res
    j1 = {}
    for pose in ("DIAG", "EXT45"):
        if (pose, True) not in J1 or (pose, False) not in J1:
            continue
        a, b = J1[(pose, True)], J1[(pose, False)]
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            cm = js.compare_known_mass(a, b, M_W)
        out = {k: (float(v) if isinstance(v, (float, np.floating)) else v) for k, v in cm.items()}
        for key in ("pos", "neg"):
            di = a[f"i_{key}"] - b[f"i_{key}"]
            out[f"K_{key}"] = abs(cm["d_tau"]) / abs(di)
            out[f"sK_{key}"] = out[f"K_{key}"] * math.hypot(a[f"sd_{key}"], b[f"sd_{key}"]) / abs(di)
        out.update(half_without=abs(b["i_half"]), half_with=abs(a["i_half"]), T=[b["temp"], a["temp"]])
        j1[pose] = out
        say(f"### {pose}\n" + buf.getvalue().rstrip())
        say(f"   只用正向靠近 K = {out['K_pos']:.2f} ± {out['sK_pos']:.2f}，只用負向靠近 K = {out['K_neg']:.2f} ± {out['sK_neg']:.2f}；"
            f"靜摩擦帶寬 ±{out['half_without']:.3f} A（不夾）→ ±{out['half_with']:.3f} A（夾砝碼）")

    # ---------------- 7. 與 10/22 比較 ----------------
    say("=" * 100)
    cross = []
    if os.path.exists(R1022):
        R = json.load(open(R1022, encoding="utf-8"))
        say("與 10/22 比較（不夾，前 4 個條件；10/22 J0 " + f"{np.mean([x['ctrl']['T0'] for x in R['runs']]):.2f} °C）")
        for x in R["runs"]:
            t = [r for r in runs if r["pose"] == x["pose"] and r["dir"] == x["dir"] and not r["w"]]
            if not t:
                continue
            t = t[0]
            pct = [(t["J"][i] / x["J"][i] - 1) * 100 for i in range(4)]
            spct = [100 * t["J"][i] / x["J"][i] * math.hypot(t["sJ"][i] / t["J"][i], x["sJ"][i] / x["J"][i]) for i in range(4)]
            cross.append(dict(label=x["label"], pose=x["pose"], dir=x["dir"], J22=x["J"][:4], J29=t["J"][:4], pct=pct, spct=spct))
            say(f"   {x['label']:<14} 10/22 " + " ".join(f"{v:.3f}" for v in x["J"][:4]) + "  10/29 "
                + " ".join(f"{v:.3f}" for v in t["J"][:4]) + "  差 % " + " ".join(f"{v:+.1f}±{s:.1f}" for v, s in zip(pct, spct)))
    else:
        say(f"找不到 {R1022}，略過與 10/22 比較（先執行 ../1022/velmod_1022.py）")

    # ---------------- 8. 定速 B、Tc ----------------
    say("=" * 100)
    say("定速（固定角度窗口）B、Tc；10/22 的定速因速度滑桿未回 100% 不列入")
    lo, hi = cv.window_bounds_deg(180.0, 10, 101)
    files = sorted(glob.glob(os.path.join(BASE, "constvel_angle_*_full_*", "constvel_angle_data_*.csv"))
                   + glob.glob(os.path.join(CV_OLD, "constvel_angle_*_full_2026101[35]_*", "constvel_angle_data_*.csv")),
                   key=lambda p: os.path.basename(p)[-19:-4])
    cvr = []
    for F in files:
        s = +1 if "_pos_" in os.path.basename(F) else -1
        with contextlib.redirect_stdout(io.StringIO()):
            r = cv.analyze_angle_window(F, 0, KT, s, cv.SPEED_LEVELS, lo, hi)
        T = np.genfromtxt(F, delimiter=",", names=True, usecols=("joint_temp_0",))["joint_temp_0"]
        cvr.append(dict(name=os.path.basename(F)[-19:-4], dir=s, B=float(r["B"]), sB=float(r["se_B"]),
                        Tc=float(r["Tc"]), sTc=float(r["se_Tc"]), T=float(np.mean(T))))
        x = cvr[-1]
        say(f"   {x['name']}  {'正轉' if s > 0 else '反轉'}  J0 {x['T']:.2f} °C  B = {x['B']:.2f} ± {x['sB']:.2f}  Tc = {x['Tc']:.3f} ± {x['sTc']:.3f}")
    for s in (+1, -1):
        old = [x for x in cvr if x["dir"] == s and "20261029" not in x["name"]]
        new = [x for x in cvr if x["dir"] == s and "20261029" in x["name"]]
        if len(old) >= 3 and new:
            pB = np.polyfit([x["T"] for x in old], [x["B"] for x in old], 1)
            for x in new:
                say(f"   {'正轉' if s > 0 else '反轉'} 10/29：B = {x['B']:.2f}；10/13、10/15 的 B 對溫度回歸（斜率 {pB[0]:+.2f} /°C）"
                    f"在 {x['T']:.2f} °C 的預測 {np.polyval(pB, x['T']):.2f}")

    res = dict(M_W=M_W, KT=KT, J_DRIVE=J_DRIVE, cond=COND, conditions=[list(c) for c in vm.CONDITIONS],
               runs=[{k: v for k, v in r.items() if k != "F"} for r in runs], pairs=pairs, pose_ratio=pose_ratio,
               fit=fit, fit_runs=[r["label"] for r in rev], j1=j1, cross=cross, constvel=cvr)
    with open(os.path.join(OUT, "known_mass_1029.json"), "w", encoding="utf-8") as f:
        json.dump(res, f, ensure_ascii=False, indent=1, default=float)
    with open(os.path.join(OUT, "known_mass_1029.log"), "w", encoding="utf-8") as f:
        f.write(LOG.getvalue())


if __name__ == "__main__":
    main()
