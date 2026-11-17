# 11/12 J1 第二步：彙整 j1_1112.json —— 定速 DIAG vs UP、速度調變重複性、紙杯砝碼校正（不連線）
import sys, json, math, socket
socket.socket = None
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\cluade_code\UR5-SYSID\robot")
import ur5_j1_velocity_modulation_ident as v1
R = json.load(open(r"D:\cluade_code\ur5\analysis_1112\j1_1112.json", encoding="utf-8"))
K0 = v1.KT_OUT; G = 9.81; M = 1.008
ORDER = [4, 0, 3, 1, 5, 2]                       # 1Hz×0.02, 1Hz×0.04, 2Hz×0.02, 2Hz×0.04, 3Hz×0.027, 3Hz×0.04
LAB = {i: f"{f:.0f} Hz×{a}" for i, (f, a) in enumerate(v1.CONDITIONS)}
out = {}

# ---------- 幾何：紙杯在 UP 姿態（q1 = −92°）離 J1 軸的距離 ----------
q = np.radians([180, -92, 0, -90, -90, 0]); fr = v1.fk_frames(q)
o1, z1, F = fr[1][:3, 3], fr[1][:3, 2], fr[6]
def geo(d_m):
    v = F[:3, :3] @ np.array([0.001, 0.001, d_m]) + F[:3, 3] - o1
    v = v - (v @ z1) * z1
    return float(np.linalg.norm(v)), float(v[2]), float(v[0])     # ρ, 高度 h, 水平 x（世界 x）

# ---------- 1. 定速 DIAG vs 11/5 UP ----------
c = R["constvel_diag"]
print("=== 1. J1 定速：DIAG（11/12）vs UP（11/5：B 19.775 ± 0.263，Tc 9.164 ± 0.032，H̄ +0.004，30.58 °C）")
print(f"DIAG: B = {c['B']:.3f} ± {c['sd_B']:.3f}（正 {c['B_pos']:.3f}、反 {c['B_neg']:.3f}）  Tc = {c['Tc']:.3f} ± {c['sd_Tc']:.3f}  "
      f"H̄ = {c['H']:+.3f}  C₊ {c['C_pos']:.3f}  C₋ {c['C_neg']:.3f}  溫度 {c['temps']}")
print(f"      target_moment 窗口平均 {c['tm_win']:+.3f} → K_mean {c['K_mean']:.2f}；斜率 tau {c['slope_tau']:+.4f}、tm {c['slope_tm']:+.4f} → K_slope {c['K_slope']:.2f}")

# ---------- 2. 速度調變重複：11/5 vs 11/12（p130） ----------
print("\n=== 2. J1 速度調變 UP 重複（×13.64）")
print(f"{'條件':<11}" + "".join(f"{k:>16}" for k in ("11/5 正", "11/12 正", "11/5 反", "11/12 反")) + f"{'差(正)/σ':>10}{'差(反)/σ':>10}")
for ci in ORDER:
    row, dd = [], []
    for tag in ("pos", "neg"):
        a, b = R[f"velmod_1105_p130_{tag}"]["per"][str(ci)], R[f"velmod_1112_p130_{tag}"]["per"][str(ci)]
        row += [f"{a['J']:.3f}±{a['se']:.3f}", f"{b['J']:.3f}±{b['se']:.3f}"]
        dd.append((b["J"] - a["J"]) / math.hypot(a["se"], b["se"]))
    print(f"{LAB[ci]:<11}" + "".join(f"{x:>16}" for x in (row[0], row[1], row[2], row[3])) + f"{dd[0]:>10.2f}{dd[1]:>10.2f}")
for nm in ("1105_p130", "1112_p130", "1112_p231"):
    for tag in ("pos", "neg"):
        r = R[f"velmod_{nm}_{tag}"]
        print(f"  {nm} {tag}: 溫度 {r['temp'][0]:.2f}→{r['temp'][1]:.2f}  追隨增益 " +
              " ".join(f"{r['per'][str(ci)]['gain']:.3f}" for ci in ORDER) + "  J_model " +
              " ".join(f"{r['per'][str(ci)]['Jm']:.3f}" for ci in ORDER))

# ---------- 3. 紙杯砝碼：ΔJ ----------
print("\n=== 3. 紙杯砝碼 1.008 kg：ΔJ = J(2.31) − J(1.30)，同一天同方向（×13.64）")
dJ = {}
for tag in ("pos", "neg"):
    for ci in ORDER:
        a, b = R[f"velmod_1112_p130_{tag}"]["per"][str(ci)], R[f"velmod_1112_p231_{tag}"]["per"][str(ci)]
        dJ[(tag, ci)] = (b["J"] - a["J"], math.hypot(a["se"], b["se"]), b["Jm"] - a["Jm"])
print(f"{'條件':<11}{'ΔJ 正':>16}{'ΔJ 反':>16}{'平均':>16}{'控制器 ΔJ_model':>16}")
dJm = []
for ci in ORDER:
    p, n = dJ[("pos", ci)], dJ[("neg", ci)]
    w = np.array([1 / p[1] ** 2, 1 / n[1] ** 2]); m = (p[0] * w[0] + n[0] * w[1]) / w.sum(); s = 1 / math.sqrt(w.sum())
    dJm.append((ci, m, s, 0.5 * (p[2] + n[2])))
    print(f"{LAB[ci]:<11}{p[0]:>9.3f}±{p[1]:.3f}{n[0]:>9.3f}±{n[1]:.3f}{m:>9.3f}±{s:.3f}{0.5*(p[2]+n[2]):>16.4f}")
out["dJ"] = {LAB[ci]: [m, s, jm] for ci, m, s, jm in dJm}

# 控制器的 ΔJ_model → 推回教導器上填的重心（點質量：2.31 kg 在 z_e 與 1.30 kg 在 74 mm 之差）
Jm_d = float(np.mean([x[3] for x in dJm]))
def dJ_ctrl(z_e):
    return 2.31 * geo(z_e)[0] ** 2 - 1.30 * geo(0.074)[0] ** 2
zz = np.linspace(0.0, 0.4, 4001); jj = np.array([dJ_ctrl(z) for z in zz])
z_e = float(zz[np.argmin(np.abs(jj - Jm_d))])
d_e = (2.308 * z_e - 1.30 * 0.074) / 1.008
print(f"\n控制器 ΔJ_model 平均 {Jm_d:.4f} → 教導器重心 z ≈ {z_e*1000:.0f} mm → 紙杯中心離法蘭 d ≈ {d_e*1000:.0f} mm（由公式反推）")

# ---------- 4. 重力：斜率與位準 ----------
print("\n=== 4. 重力（絕對座標，×13.64）：q = −92° 的斜率 dτ/dq [N·m/rad]、整段平均位準 [N·m]")
def slope(p, key_s, key_c, qdeg):
    qq = math.radians(qdeg); return p[key_s] * math.cos(qq) - p[key_c] * math.sin(qq)
for tag in ("pos", "neg"):
    print(f"  {tag}: " + "  ".join(
        f"{LAB[ci]} {slope(R[f'velmod_1112_p231_{tag}']['per'][str(ci)], 'G_s', 'G_c', -92) - slope(R[f'velmod_1112_p130_{tag}']['per'][str(ci)], 'G_s', 'G_c', -92):+.2f}"
        for ci in ORDER))
dS = np.array([[slope(R[f'velmod_1112_p231_{t}']['per'][str(ci)], 'G_s', 'G_c', -92) - slope(R[f'velmod_1112_p130_{t}']['per'][str(ci)], 'G_s', 'G_c', -92)
                for ci in ORDER] for t in ("pos", "neg")])
dSm = np.array([[slope(R[f'velmod_1112_p231_{t}']['per'][str(ci)], 'Gs_m', 'Gc_m', -92) - slope(R[f'velmod_1112_p130_{t}']['per'][str(ci)], 'Gs_m', 'Gc_m', -92)
                 for ci in ORDER] for t in ("pos", "neg")])
print(f"  量測 Δ斜率：平均 {dS.mean():+.3f}，各條件標準差 {dS.std(ddof=1):.3f}；控制器 Δ斜率 {dSm.mean():+.3f}")
lv = {nm: 0.5 * (R[f"velmod_{nm}_pos"]["lvl_cur"] + R[f"velmod_{nm}_neg"]["lvl_cur"]) for nm in ("1112_p130", "1112_p231")}
lt = {nm: 0.5 * (R[f"velmod_{nm}_pos"]["lvl_tm"] + R[f"velmod_{nm}_neg"]["lvl_tm"]) for nm in ("1112_p130", "1112_p231")}
fr_ = {nm: 0.5 * (R[f"velmod_{nm}_pos"]["lvl_cur"] - R[f"velmod_{nm}_neg"]["lvl_cur"]) for nm in ("1112_p130", "1112_p231")}
qm = {nm: 0.5 * (R[f"velmod_{nm}_pos"]["q_mean"] + R[f"velmod_{nm}_neg"]["q_mean"]) for nm in ("1112_p130", "1112_p231")}
print(f"  位準（正反平均 = 重力）：電流 1.30 {lv['1112_p130']:+.3f} → 2.31 {lv['1112_p231']:+.3f}（Δ {lv['1112_p231']-lv['1112_p130']:+.3f}）；"
      f"target_moment {lt['1112_p130']:+.3f} → {lt['1112_p231']:+.3f}（Δ {lt['1112_p231']-lt['1112_p130']:+.3f}）；平均角度 {qm}")
print(f"  摩擦（正反差 / 2）：1.30 {fr_['1112_p130']:.3f}，2.31 {fr_['1112_p231']:.3f}")
for dm in (0.10, 0.15, 0.20, 0.25):
    rho, h, x = geo(dm)
    print(f"  幾何 d = {dm*1000:.0f} mm：ρ = {rho:.4f}，m ρ² = {M*rho**2:.4f}，m g h = {M*G*h:.3f} N·m/rad，m g x = {M*G*x:.3f} N·m")
rho_e, h_e, x_e = geo(d_e)
out.update(Jm_d=Jm_d, z_e=z_e, d_e=d_e, rho=rho_e, h=h_e, x=x_e, dS=dS.tolist(), dSm=dSm.tolist(), lv=lv, lt=lt, fr=fr_)
json.dump(out, open(r"D:\cluade_code\ur5\analysis_1112\j1_1112_b.json", "w", encoding="utf-8"), ensure_ascii=False, indent=1)
