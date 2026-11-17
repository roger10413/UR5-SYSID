# 11/12 J1：用紙杯砝碼的 ΔJ 檢驗「彈性」模型（推測用，不連線）
import sys, json, math, socket
socket.socket = None
import numpy as np
from scipy.optimize import least_squares
sys.stdout.reconfigure(encoding="utf-8")
R = json.load(open(r"D:\cluade_code\ur5\analysis_1112\j1_1112.json", encoding="utf-8"))
DJT = 0.905                                     # 紙杯 1.008 kg × ρ²（d ≈ 178 mm，ρ ≈ 0.948 m）
def cond_mean(nm, ci):
    v = [R[f"velmod_{nm}_{t}"]["per"][str(ci)] for t in ("pos", "neg")]
    w = np.array([1 / x["se"] ** 2 for x in v]); J = sum(x["J"] * wi for x, wi in zip(v, w)) / w.sum()
    return J, 1 / math.sqrt(w.sum()), v[0]["f"]
rows = []
for nm, dj in (("1112_p130", 0.0), ("1112_p231", DJT)):
    for ci in range(6):
        J, s, f = cond_mean(nm, ci); rows.append((nm, ci, f, dj, J, s))
f = np.array([r[2] for r in rows]); dj = np.array([r[3] for r in rows]); Jm = np.array([r[4] for r in rows]); sJ = np.array([r[5] for r in rows])
w2 = (2 * math.pi * f) ** 2
def m1(p):                                      # 整個慣量都在彈簧後：s·J/(1 − J ω²/K)
    s, J, K = p; JJ = J + dj; return s * JJ / (1 - JJ * w2 / K)
def m2(p):                                      # 只有 Jb 在彈簧後，砝碼加在 Jb
    s, Ja, Jb, K = p; JJ = Jb + dj; return s * (Ja + JJ / (1 - JJ * w2 / K))
for name, fun, p0 in (("整個慣量在彈簧後", m1, [1.1, 5.5, 3000.0]), ("部分慣量 Jb 在彈簧後", m2, [1.1, 4.0, 1.5, 800.0])):
    r = least_squares(lambda p: (fun(p) - Jm) / sJ, p0)
    J_ = np.linalg.inv(r.jac.T @ r.jac); se = np.sqrt(np.diag(J_)); chi2 = float(np.sum(r.fun ** 2)); dof = len(Jm) - len(p0)
    print(f"\n[{name}] 參數 {np.round(r.x, 4)} ± {np.round(se, 4)}；χ² = {chi2:.1f}（自由度 {dof}）")
    if fun is m1:
        s, J, K = r.x; print(f"  換算比例 K = 13.64 / s = {13.64 / s:.2f}；固有頻率 1.30 kg {math.sqrt(K / J) / 2 / math.pi:.2f} Hz，2.31 kg {math.sqrt(K / (J + DJT)) / 2 / math.pi:.2f} Hz")
    else:
        s, Ja, Jb, K = r.x; print(f"  換算比例 K = 13.64 / s = {13.64 / s:.2f}；固有頻率 1.30 kg {math.sqrt(K / Jb) / 2 / math.pi:.2f} Hz，2.31 kg {math.sqrt(K / (Jb + DJT)) / 2 / math.pi:.2f} Hz")
    for (nm, ci, ff, d_, J, s_), mv in zip(rows, fun(r.x)):
        print(f"  {nm} {ff:.0f} Hz 條件{ci+1}: 量測 {J:.3f} ± {s_:.3f}  模型 {mv:.3f}  差/σ {(J - mv) / s_:+.1f}")
# 1/J 對 ω²（單一彈簧時為直線，兩種負載斜率相同）
print("\n1/J 對 ω²（各頻率兩種振幅、兩方向平均）：")
for nm in ("1105_p130", "1112_p130", "1112_p231"):
    out = []
    for fr_ in (1.0, 2.0, 3.0):
        vals = [cond_mean(nm, ci)[0] for ci in range(6) if abs(R[f"velmod_{nm}_pos"]["per"][str(ci)]["f"] - fr_) < 0.2]
        out.append((fr_, float(np.mean(vals))))
    sl = [(1 / out[k + 1][1] - 1 / out[k][1]) / ((2 * math.pi) ** 2 * (out[k + 1][0] ** 2 - out[k][0] ** 2)) for k in range(2)]
    print(f"  {nm}: J = " + "、".join(f"{a:.0f} Hz {b:.3f}" for a, b in out) + f"；斜率 1→2 Hz {sl[0]:.2e}、2→3 Hz {sl[1]:.2e}")
