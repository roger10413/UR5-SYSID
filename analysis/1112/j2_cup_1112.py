# 離線：J2 速度調變 UP，夾紙杯砝碼（2.31 kg）vs 只有工具（1.30 kg）；ΔJ 與重力斜率 → J2 換算比例（不連線、不動手臂）
import sys, os, glob, json, math, socket
socket.socket = None
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\cluade_code\UR5-SYSID\robot")
import ur5_j2_velocity_modulation_ident as v2
import io, contextlib
D = r"D:\cluade_code\ur5\data_1112\J2"; K = v2.KT_OUT; G = 9.81; M = 1.008
def one(p):
    f = glob.glob(p); assert len(f) == 1, (p, f); return f[0]
def tm_fit_full(cs, tm, tqdd, tq, q):
    Xs, ys = [], []
    for g, sg in enumerate(cs):
        ii = sg["idx"]; Dm = np.zeros((len(ii), len(cs))); Dm[:, g] = 1.0
        Xs.append(np.column_stack([tqdd[ii], tq[ii], np.sin(q[ii]), np.cos(q[ii]), Dm])); ys.append(tm[ii])
    p, *_ = np.linalg.lstsq(np.vstack(Xs), np.concatenate(ys), rcond=None); return p
R = {}
for pay in ("p130", "p231"):
    for d, tag in ((+1, "pos"), (-1, "neg")):
        path = one(os.path.join(D, f"j2velmod_up_{pay}_{tag}_full_*", "*.csv"))
        with contextlib.redirect_stdout(io.StringIO()):
            r = v2.analyze_j2_velocity_modulation(path, 2, K, d, v2.CONDITIONS, v2.V_CRUISE)
        Dd = np.genfromtxt(path, delimiter=",", names=True)
        q = Dd["actual_q_2"]; tm = Dd["target_moment_2"]; tqdd = d * Dd["target_qdd_2"]; tq = d * Dd["target_qd_2"]
        per = {}
        for ci, c in r["conditions"].items():
            cs = [sg for sg in r["segments"] if sg["cond"] == ci]
            p = tm_fit_full(cs, d * tm, tqdd, tq, q)
            idx = np.concatenate([sg["idx"] for sg in cs])
            per[ci] = dict(f=c["freq"], A=c["amp"], J=c["J"], se=c["se_J"], G_s=d * c["G_s"], G_c=d * c["G_c"],
                           Jm=float(p[0]), Gs_m=d * float(p[2]), Gc_m=d * float(p[3]), q0=float(np.mean(q[idx])),
                           temp=c["temp"], gain=c["track_gain"])
        R[f"{pay}_{tag}"] = per
        print(f"{pay} {tag}: J2 溫度 {Dd['joint_temp_2'][0]:.2f} → {Dd['joint_temp_2'][-1]:.2f}")
slope = lambda p, ks, kc, q0: p[ks] * math.cos(q0) - p[kc] * math.sin(q0)
q0 = float(np.mean([R["p130_pos"][ci]["q0"] for ci in R["p130_pos"]]))
print(f"\n調變段平均角度 {math.degrees(q0):.2f}°")
print("\n條件         J 1.30（+/−）       J 2.31（+/−）       ΔJ（+/−）           ΔJ_model    Δ斜率量測（+/−）  Δ斜率控制器")
rows = []
for ci in sorted(R["p130_pos"]):
    a = R["p130_pos"][ci]; lab = f"{a['f']:.0f} Hz A{a['A']:.2f}"
    dJ = [R[f"p231_{t}"][ci]["J"] - R[f"p130_{t}"][ci]["J"] for t in ("pos", "neg")]
    sdJ = [math.hypot(R[f"p231_{t}"][ci]["se"], R[f"p130_{t}"][ci]["se"]) for t in ("pos", "neg")]
    dJm = np.mean([R[f"p231_{t}"][ci]["Jm"] - R[f"p130_{t}"][ci]["Jm"] for t in ("pos", "neg")])
    dS = [slope(R[f"p231_{t}"][ci], "G_s", "G_c", q0) - slope(R[f"p130_{t}"][ci], "G_s", "G_c", q0) for t in ("pos", "neg")]
    dSm = np.mean([slope(R[f"p231_{t}"][ci], "Gs_m", "Gc_m", q0) - slope(R[f"p130_{t}"][ci], "Gs_m", "Gc_m", q0) for t in ("pos", "neg")])
    rows.append(dict(ci=ci, lab=lab, f=a["f"], A=a["A"], dJ=dJ, sdJ=sdJ, dJm=dJm, dS=dS, dSm=dSm,
                     J130=[R[f"p130_{t}"][ci]["J"] for t in ("pos", "neg")], J231=[R[f"p231_{t}"][ci]["J"] for t in ("pos", "neg")],
                     Jm130=R["p130_pos"][ci]["Jm"]))
    print(f"{lab:12s} {rows[-1]['J130'][0]:.3f}/{rows[-1]['J130'][1]:.3f}     {rows[-1]['J231'][0]:.3f}/{rows[-1]['J231'][1]:.3f}     "
          f"{dJ[0]:+.3f}±{sdJ[0]:.3f}/{dJ[1]:+.3f}±{sdJ[1]:.3f}   {dJm:+.4f}    {dS[0]:+.2f}/{dS[1]:+.2f}        {dSm:+.3f}")
dS_all = np.array([r["dS"] for r in rows]); dSm_all = np.array([r["dSm"] for r in rows])
mgh = M * G * 0.4869
print(f"\n量測 Δ斜率平均 {dS_all.mean():+.3f}（各條件、方向標準差 {dS_all.std(ddof=1):.3f}，n = {dS_all.size}）；控制器 Δ斜率 {dSm_all.mean():+.3f}；紙杯 m g h = {mgh:.3f} N·m/rad")
Kcal = K * mgh / abs(dS_all.mean()); sdK = Kcal * dS_all.std(ddof=1) / math.sqrt(dS_all.size) / abs(dS_all.mean())
print(f"→ J2 換算比例（已知質量，重力斜率）K = {K:.2f} × {mgh:.3f} ÷ {abs(dS_all.mean()):.3f} = {Kcal:.2f}（± {sdK:.2f}，各條件平均的標準差，推測各條件不完全獨立）")
print(f"   控制器自己：Δ斜率 {abs(dSm_all.mean()):.3f} vs m g h {mgh:.3f}（比例 {abs(dSm_all.mean())/mgh:.3f}）")
for dd in (0.056, 0.10, 0.15, 0.20):
    x = 0.0823 + dd; print(f"   紙杯中心離法蘭 {dd*1000:.0f} mm：ρ = {math.hypot(x, 0.4869):.4f}，m ρ² = {M*(x*x+0.4869**2):.4f}")
json.dump(dict(R={k: {str(c): v for c, v in val.items()} for k, val in R.items()}, rows=rows, q0=q0, Kcal=Kcal, sdK=sdK, mgh=mgh),
          open("j2_cup_1112.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1, default=float)
