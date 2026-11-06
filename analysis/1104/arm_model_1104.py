"""UR5（CB3）6 軸剛體模型：用 UR 官方 DH 頁面的連桿質量與重心 + 教導器負載，算任意姿態的力矩、慣量矩陣與重力矩。

模型（Newton-Euler 遞迴，世界座標）：
  τ = M(q)·q̈ + C(q, q̇)·q̇ + g(q) + 轉子項（ROTOR[j]·q̈_j）
  - 連桿：官方點質量（LINK_M、LINK_COM）；可選擇加上連桿自身的轉動慣量（CB3 未公開，借用 Clochiatti 2024 的 UR5e 值，只作對照）。
  - 負載：點質量，重心在法蘭座標 TOOL_COG。
  - 轉子 + 減速機：J0～J2 用設定檔的 1.915（Gnad 2026 表 3）；J3～J5 未知，設 0（只有 J0 動作或靜止時不影響比較）。
這是名目模型（和控制器同一組資料），不是量測；沒有摩擦、沒有關節彈性、沒有電流換算比例 K。

用 10/29 資料對照：
1. 自我檢查：M[0,0] 與 known_mass_1029 的點質量 J0 慣量相同；DIAG 的 J1 重力保持力矩約 −23.3 N·m；M 對稱。
2. 控制器的 target_moment（6 軸）：J1 靜止 4 次（重力）、速度調變 6 次（J0 動作時 6 軸）。負載用教導器設定值。
   若控制器用的就是這組資料，差應接近 0。
3. 實測電流 × K 和模型比：J1、J2 靜止時的重力（負載用實際值：不夾 1.30、夾砝碼 2.30 kg）；K 用 13.64 與砝碼校正的 10.8。
4. 常用姿態的表：J0～J2 的慣量對角項與 J1、J2 重力保持力矩。
只做離線分析，執行時封鎖網路。結果存到 OUT/arm_model_1104.log 與 .json。
"""
import os, sys, glob, io, csv, json, socket, contextlib
import numpy as np

socket.socket = None                                   # 封鎖網路：只做離線分析
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.path.insert(0, os.path.join(HERE, "..", "1029"))
sys.stdout.reconfigure(encoding="utf-8")
import ur5_j1_static_current as js
import known_mass_1029 as km

BASE = km.BASE                                         # 10/29 資料
OUT = os.path.join(os.path.dirname(BASE), "analysis_1104"); os.makedirs(OUT, exist_ok=True)
G = 9.81
ROTOR = np.array([km.J_DRIVE] * 3 + [0.0] * 3)         # 轉子 + 減速機換算到關節側 [kg·m²]；J3～J5 未知
K_NOM, K_CAL = km.KT, 10.8                             # 電流換算比例：文獻 0.135 × 101、砝碼校正（11/4 報告 3.1）
W_KG = km.M_W                                          # 砝碼質量（標稱 1.000 kg）
POSES = {"DIAG": (180, -90, 90, -90, -90, 0), "EXT45": (180, -45, 45, -90, -90, 0), "UP": (180, -90, 0, -90, -90, 0)}
LOG = io.StringIO()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def bodies(pay_kg, self_inertia=False):
    """每根連桿上的剛體：(連桿編號, 質量, 重心（連桿座標）, 重心處慣量（連桿座標）)；負載掛在第 6 根連桿上。"""
    out = [(i, km.LINK_M[i], np.array(km.LINK_COM[i]), km.LINK_IC[i] if self_inertia else np.zeros((3, 3)))
           for i in range(6)]
    out.append((5, pay_kg, km.TOOL_COG, np.zeros((3, 3))))
    return out


def rnea(q, qd, qdd, pay_kg, self_inertia=False, rotor=ROTOR, gravity=True):
    """Newton-Euler 逆動力學：回傳 6 軸力矩 [N·m]（正負號依 DH 轉軸方向，和 target_moment 相同）。"""
    fr = js.fk_frames(q)
    z = [f[:3, 2] for f in fr]
    p = [f[:3, 3] for f in fr]
    bs = bodies(pay_kg, self_inertia)
    w, dw = np.zeros(3), np.zeros(3)
    a = np.array([0.0, 0.0, G if gravity else 0.0])    # 基座往上加速 g，等效於重力
    Fl, Ml = [], []                                    # 每根連桿的合力、對其轉軸點 p[i] 的合力矩
    for i in range(6):
        w_new = w + z[i] * qd[i]
        dw = dw + z[i] * qdd[i] + np.cross(w, z[i] * qd[i])
        w = w_new
        F, M = np.zeros(3), np.zeros(3)
        R = fr[i + 1][:3, :3]
        for li, m, com, Ic in bs:
            if li != i:
                continue
            r = R @ com + p[i + 1] - p[i]
            ac = a + np.cross(dw, r) + np.cross(w, np.cross(w, r))
            Iw = R @ Ic @ R.T
            f = m * ac
            F += f
            M += np.cross(r, f) + Iw @ dw + np.cross(w, Iw @ w)
        Fl.append(F); Ml.append(M)
        d = p[i + 1] - p[i]
        a = a + np.cross(dw, d) + np.cross(w, np.cross(w, d))
    tau = np.zeros(6)
    f, n = np.zeros(3), np.zeros(3)
    for i in range(5, -1, -1):
        n = Ml[i] + n + np.cross(p[i + 1] - p[i], f)
        f = Fl[i] + f
        tau[i] = n @ z[i] + rotor[i] * qdd[i]
    return tau


def mass_matrix(q, pay_kg, self_inertia=False, rotor=ROTOR):
    z6 = np.zeros(6)
    return np.column_stack([rnea(q, z6, e, pay_kg, self_inertia, rotor, gravity=False) for e in np.eye(6)])


def holding(q, pay_kg):
    """撐住重力的力矩 g(q)（靜止時控制器要給的力矩）。"""
    return rnea(q, np.zeros(6), np.zeros(6), pay_kg)


def read_csv(path, step=1):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))[::step]
    keys = rows[0].keys()
    return {k: np.array([float(r[k]) if r[k] not in ("", "None") else np.nan for r in rows]) for k in keys}


def cols(D, name):
    return np.column_stack([D[f"{name}_{j}"] for j in range(6)])


def run_tag(path):
    """資料夾名稱 {種類}_{姿態}_p{教導器負載×100}_..."""
    parts = os.path.basename(os.path.dirname(path)).split("_")
    pose, pay = parts[1].upper(), int(parts[2][1:]) / 100
    return pose, pay, pay > km.TOOL_KG + 0.005


def fmt(v, f="{:+7.3f}"):
    return " ".join(f.format(x) for x in v)


def main():
    out = {}
    # ---------------- 1. 自我檢查 ----------------
    say("=" * 70); say("1. 自我檢查（負載 1.30 kg）"); say("=" * 70)
    for name in ("DIAG", "EXT45"):
        q = np.radians(POSES[name])
        Mq = mass_matrix(q, km.TOOL_KG, rotor=np.zeros(6))
        jl, jt, _, _ = km.point_mass_model(q, km.TOOL_KG)
        say(f"{name}: M[0,0] = {Mq[0, 0]:.4f}，known_mass 點質量 {jl + jt:.4f}；"
            f"M 不對稱最大 {np.max(np.abs(Mq - Mq.T)):.1e}；J1 重力保持 {holding(q, km.TOOL_KG)[1]:+.3f} N·m")
        Ms = mass_matrix(q, km.TOOL_KG, self_inertia=True, rotor=np.zeros(6))
        say(f"   加連桿自身慣量後 M[0,0] = {Ms[0, 0]:.4f}（差 {Ms[0, 0] - Mq[0, 0]:.4f}；known_mass {km.link_self_inertia(q):.4f}）")

    # ---------------- 2. 控制器 target_moment ----------------
    say("\n" + "=" * 70); say("2. 模型 vs 控制器 target_moment（負載 = 教導器設定值）"); say("=" * 70)
    say("2a. J1 靜止（每 25 筆取 1 筆；用 target_q）：各軸 target_moment 平均、模型平均、差的 RMS [N·m]")
    out["static_tm"] = {}
    for F in sorted(glob.glob(os.path.join(BASE, "j1static_*_full_*", "j1static_data_*.csv"))):
        pose, pay, w = run_tag(F)
        D = read_csv(F, step=25)
        Q, TM = cols(D, "target_q"), cols(D, "target_moment")
        Hm = np.array([holding(q, pay) for q in Q])
        rms = np.sqrt(np.mean((TM - Hm) ** 2, axis=0))
        say(f"  {pose:5s} 設定 {pay:.2f} kg  控制器 {fmt(TM.mean(0))}")
        say(f"  {'':5s} {'':10s}  模型   {fmt(Hm.mean(0))}   差 RMS {fmt(rms, '{:.3f}')}")
        out["static_tm"][f"{pose}_{pay:.2f}"] = dict(tm=TM.mean(0).tolist(), model=Hm.mean(0).tolist(), rms=rms.tolist())

    say("\n2b. 速度調變（J0 動作；每 5 筆取 1 筆；用 actual_q、target_qd、target_qdd）：差的 RMS 與 target_moment 本身的 RMS [N·m]")
    out["velmod_tm"] = {}
    for F in sorted(glob.glob(os.path.join(BASE, "velmod_*_full_*", "velmod_data_*.csv"))):
        pose, pay, w = run_tag(F)
        dr = "正轉" if "_pos_" in F else "反轉"
        D = read_csv(F, step=5)
        Q, QD, QDD, TM = cols(D, "actual_q"), cols(D, "target_qd"), cols(D, "target_qdd"), cols(D, "target_moment")
        Tm = np.array([rnea(*x, pay) for x in zip(Q, QD, QDD)])
        Tn = np.array([rnea(q, np.zeros(6), np.zeros(6), pay) for q in Q])     # 只有重力（不含動態項）
        e, en = TM - Tm, TM - Tn
        rms = lambda x: np.sqrt(np.mean(x ** 2, axis=0))
        say(f"  {pose:5s} {dr} 設定 {pay:.2f} kg  差 RMS {fmt(rms(e), '{:.3f}')}   只用重力時 {fmt(rms(en), '{:.3f}')}")
        say(f"  {'':5s} {'':4s} {'':10s}  J0 動態部分 RMS（target_moment − 重力）{rms(en)[0]:.3f}，J0 最大 |target_moment| {np.max(np.abs(TM[:, 0])):.2f}")
        out["velmod_tm"][f"{pose}_{dr}_{pay:.2f}"] = dict(rms=rms(e).tolist(), rms_gravity_only=rms(en).tolist())

    # ---------------- 3. 實測電流 vs 模型（重力） ----------------
    say("\n" + "=" * 70); say("3. 實測電流 × K vs 模型重力（負載 = 實際值：不夾 1.30、夾砝碼 2.30 kg）"); say("=" * 70)
    say("   J1：兩方向靠近的中心電流（消去靜摩擦）；J2：整段平均（J2 沒動，含未知的靜摩擦）")
    out["static_cur"] = {}
    for F in sorted(glob.glob(os.path.join(BASE, "j1static_*_full_*", "j1static_data_*.csv"))):
        pose, pay, w = run_tag(F)
        real = km.TOOL_KG + (W_KG if w else 0.0)
        with contextlib.redirect_stdout(io.StringIO()):          # 不印 analyze_static_current 的逐項輸出
            res = js.analyze_static_current(F)
        D = read_csv(F, step=25)
        H = holding(res["q_mean"], real)
        i2 = float(np.mean(D["actual_current_2"]))
        tm = holding(res["q_mean"], pay)
        say(f"  {pose:5s} {'夾砝碼' if w else '不夾  '}（實際 {real:.2f} kg）")
        for j, i_m in ((1, res["i_mid"]), (2, i2)):
            say(f"     J{j}: 模型 {H[j]:+7.2f}  控制器 {tm[j]:+7.2f}（設定負載）  電流 {i_m:+.4f} A × 13.64 = {i_m * K_NOM:+7.2f}"
                f"（÷ 模型 {i_m * K_NOM / H[j]:.3f}）  × 10.8 = {i_m * K_CAL:+7.2f}  → 電流換算 K = 模型 ÷ 電流 = {H[j] / i_m:.2f}")
        out["static_cur"][f"{pose}_{'w' if w else 'nw'}"] = dict(model=H.tolist(), i1=res["i_mid"], i2=i2,
                                                                  K1=H[1] / res["i_mid"], K2=H[2] / i2)

    # ---------------- 4. J0 實測 vs 模型 ----------------
    say("\n" + "=" * 70); say("4. J0 慣量：模型 M[0,0] vs 實測（反轉 6 條件平均；÷ 1.26 為砝碼校正）"); say("=" * 70)
    R = json.load(open(os.path.join(km.OUT, "known_mass_1029.json"), encoding="utf-8"))
    out["j0"] = []
    for r in R["runs"]:
        if r["dir"] > 0:
            continue
        real = km.TOOL_KG
        q = np.radians(r["q_deg"])
        M0 = mass_matrix(q, real, rotor=np.zeros(6))[0, 0] + W_KG * r["r_fl"] ** 2 * r["w"]
        Jm = float(np.mean(r["J"]))
        say(f"  {r['pose']:5s} {'夾砝碼' if r['w'] else '不夾  '}  模型 {M0:.3f}（+1.915 = {M0 + km.J_DRIVE:.3f}）  "
            f"實測 {Jm:.3f}  ÷ 1.26 = {Jm / 1.26:.3f}")
        out["j0"].append(dict(pose=r["pose"], w=r["w"], model=M0, J=Jm))

    # ---------------- 5. 常用姿態 ----------------
    say("\n" + "=" * 70); say("5. 常用姿態（負載 1.30 kg；不含轉子 1.915、不含連桿自身慣量）"); say("=" * 70)
    say("  姿態    M00     M11     M22     M01     M12   | J1 重力保持  J2 重力保持 [N·m]")
    out["poses"] = {}
    for name, qd in POSES.items():
        q = np.radians(qd)
        Mq = mass_matrix(q, km.TOOL_KG, rotor=np.zeros(6))
        H = holding(q, km.TOOL_KG)
        say(f"  {name:5s} {Mq[0, 0]:6.3f}  {Mq[1, 1]:6.3f}  {Mq[2, 2]:6.3f}  {Mq[0, 1]:+6.3f}  {Mq[1, 2]:+6.3f}  | {H[1]:+8.2f}     {H[2]:+8.2f}")
        out["poses"][name] = dict(M=Mq.tolist(), holding=H.tolist())

    with open(os.path.join(OUT, "arm_model_1104.log"), "w", encoding="utf-8") as fh:
        fh.write(LOG.getvalue())
    with open(os.path.join(OUT, "arm_model_1104.json"), "w", encoding="utf-8") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
