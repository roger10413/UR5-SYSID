# 11/19 J1 換算比例重算：砝碼（1.008 kg）造成的重力斜率 vs m·g·h（做法同 analysis/1112/j1_1112_b.py 第 4 節）
#   換算比例 = 13.64 × m·g·h / 量到的 Δ斜率（電流 × 13.64）。11/12 得到 12.03。
import sys, os, glob, json, math, socket
socket.socket = None                                   # 封鎖網路
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"F:\cluade_code\UR5-SYSID\robot")
import ur5_j1_velocity_modulation_ident as v1

D19 = r"D:\cluade_code\ur5\data_1119\J1"
A12 = r"D:\cluade_code\ur5\analysis_1112"
OUT = r"D:\cluade_code\ur5\analysis_1119"
K0, G, M = v1.KT_OUT, 9.81, 1.008
ORDER = [4, 0, 3, 1, 5, 2]
LAB = {i: f"{f:.0f} Hz×{a}" for i, (f, a) in enumerate(v1.CONDITIONS)}
d_e = json.load(open(os.path.join(A12, "j1_1112_b.json"), encoding="utf-8"))["d_e"]   # 砝碼中心離法蘭（11/12 由控制器 ΔJ_model 反推）


def one(pattern):
    f = glob.glob(pattern)
    assert len(f) == 1, (pattern, f)
    return f[0]


def h_of(q1_deg):
    """砝碼在 J1 軸上方的高度（UP 姿態，J1 角度 q1）。"""
    fr = v1.fk_frames(np.radians([180, q1_deg, 0, -90, -90, 0]))
    o1, z1, F = fr[1][:3, 3], fr[1][:3, 2], fr[6]
    v = F[:3, :3] @ np.array([0.001, 0.001, d_e]) + F[:3, 3] - o1
    v = v - (v @ z1) * z1
    return float(v[2])


def tm_fit_full(cs, tm, tqdd, tq, q):
    Xs, ys = [], []
    for g, sg in enumerate(cs):
        ii = sg["idx"]
        Dm = np.zeros((len(ii), len(cs))); Dm[:, g] = 1.0
        Xs.append(np.column_stack([tqdd[ii], tq[ii], np.sin(q[ii]), np.cos(q[ii]), Dm])); ys.append(tm[ii])
    p, *_ = np.linalg.lstsq(np.vstack(Xs), np.concatenate(ys), rcond=None)
    return p


def run(pay_tag, tag, d):
    path = one(os.path.join(D19, f"j1velmod_up_{pay_tag}_{tag}_full_*", "*.csv"))
    r = v1.analyze_j1_velocity_modulation(path, 1, K0, d, v1.CONDITIONS, v1.V_CRUISE)
    Dd = np.genfromtxt(path, delimiter=",", names=True)
    q = Dd["actual_q_1"]; tm = Dd["target_moment_1"]
    tqdd = d * Dd["target_qdd_1"]; tq = d * Dd["target_qd_1"]
    per = {}
    for ci, c in r["conditions"].items():
        cs = [sg for sg in r["segments"] if sg["cond"] == ci]
        p = tm_fit_full(cs, d * tm, tqdd, tq, q)
        idx = np.concatenate([sg["idx"] for sg in cs])
        per[ci] = dict(G_s=d * c["G_s"], G_c=d * c["G_c"], Gs_m=d * float(p[2]), Gc_m=d * float(p[3]),
                       q_mean=float(np.degrees(np.mean(q[idx]))))
    return per


def slope(p, ks, kc, qdeg):
    qq = math.radians(qdeg)
    return p[ks] * math.cos(qq) - p[kc] * math.sin(qq)


import io, contextlib
with contextlib.redirect_stdout(io.StringIO()):
    R = {(pt, t): run(pt, t, d) for pt in ("p130", "p231") for t, d in (("pos", +1), ("neg", -1))}

qref = float(np.mean([R[k][ci]["q_mean"] for k in R for ci in ORDER]))
mgh = M * G * h_of(qref)
dS = np.array([[slope(R[("p231", t)][ci], "G_s", "G_c", qref) - slope(R[("p130", t)][ci], "G_s", "G_c", qref) for ci in ORDER]
               for t in ("pos", "neg")])
dSm = np.array([[slope(R[("p231", t)][ci], "Gs_m", "Gc_m", qref) - slope(R[("p130", t)][ci], "Gs_m", "Gc_m", qref) for ci in ORDER]
                for t in ("pos", "neg")])
K = K0 * mgh / np.abs(dS)                    # 斜率的正負號只是座標定義，取大小

print(f"11/19 J1 平均角度 {qref:.2f}°，砝碼中心離法蘭 {d_e*1000:.0f} mm，m·g·h = {mgh:.3f} N·m/rad")
print(f"{'條件':<11}{'Δ斜率 正':>10}{'Δ斜率 反':>10}{'換算比例 正':>12}{'換算比例 反':>12}")
for j, ci in enumerate(ORDER):
    print(f"{LAB[ci]:<11}{dS[0, j]:>10.2f}{dS[1, j]:>10.2f}{K[0, j]:>12.2f}{K[1, j]:>12.2f}")
Km = K0 * mgh / abs(dS.mean())
se = K0 * mgh * dS.std(ddof=1) / math.sqrt(dS.size) / dS.mean() ** 2
print(f"量測 Δ斜率 平均 {dS.mean():.3f}（12 個條件標準差 {dS.std(ddof=1):.3f}）→ 換算比例 {Km:.2f} ± {se:.2f}（平均的標準誤差）")
print(f"控制器 Δ斜率 {dSm.mean():.3f}，與 m·g·h 比值 {abs(dSm.mean())/mgh:.3f}")

# 11/12 同法（同一角度、同一 m·g·h 定義）對照
R12 = json.load(open(os.path.join(A12, "j1_1112.json"), encoding="utf-8"))
q12 = float(np.mean([R12[f"velmod_1112_{n}_{t}"]["q_mean"] for n in ("p130", "p231") for t in ("pos", "neg")]))
mgh12 = M * G * h_of(q12)
dS12 = np.array([[slope(R12[f"velmod_1112_p231_{t}"]["per"][str(ci)], "G_s", "G_c", q12) -
                  slope(R12[f"velmod_1112_p130_{t}"]["per"][str(ci)], "G_s", "G_c", q12) for ci in ORDER] for t in ("pos", "neg")])
print(f"\n11/12 同法：平均角度 {q12:.2f}°，m·g·h = {mgh12:.3f}，Δ斜率 {dS12.mean():.3f}（標準差 {dS12.std(ddof=1):.3f}）→ 換算比例 {K0*mgh12/abs(dS12.mean()):.2f}")

json.dump(dict(qref=qref, d_e=d_e, mgh=mgh, dS=dS.tolist(), dSm=dSm.tolist(), K=Km, se_K=se,
               q12=q12, mgh12=mgh12, dS12=dS12.tolist(), K12=K0 * mgh12 / abs(dS12.mean())),
          open(os.path.join(OUT, "j1_gravity_1119.json"), "w", encoding="utf-8"), ensure_ascii=False, indent=1)
