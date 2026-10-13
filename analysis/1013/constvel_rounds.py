"""10/13 固定角度窗口定速：3 輪（正反轉各 3 次）比較 B、Tc、溫度、升降趟、力矩對角度斜率。"""
import sys, glob, os
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "robot"))
sys.stdout.reconfigure(encoding="utf-8")
import io, contextlib
import ur5_const_velocity_angle_ident as m

BASE = r"D:\cluade_code\ur5"   # 資料資料夾（內含各次 constvel_angle_* 子資料夾），請改成自己的路徑
a, b = m.window_bounds_deg(180.0, 10, 101)
files = sorted(glob.glob(os.path.join(BASE, "constvel_angle_*_full_20261013_*", "constvel_angle_data_*.csv")),
               key=lambda p: p[-10:-4])
rows = []
for F in files:
    s = +1 if "_pos_" in os.path.basename(F) else -1
    with contextlib.redirect_stdout(io.StringIO()):
        r = m.analyze_angle_window(F, 0, m.KT_OUT, s, m.SPEED_LEVELS, a, b)
    S = r["segments"]
    v = np.abs([x["qd_mean"] for x in S]); t = np.abs([x["tau_mean"] for x in S])
    d = np.genfromtxt(F, delimiter=",", names=True, usecols=(0, 1, 2, 3, 4, 5))
    q = np.degrees(d["actual_q_0"]); tau = d["actual_current_0"] * m.KT_OUT
    up = np.polyfit(v[:6], t[:6], 1); dn = np.polyfit(v[6:], t[6:], 1)
    lo = np.polyfit(v[v < 0.13], t[v < 0.13], 1); hi = np.polyfit(v[v > 0.1], t[v > 0.1], 1)
    sl = np.mean([np.polyfit(q[x["idx"]], np.abs(tau[x["idx"]]), 1)[0] for x in S])
    T = d["joint_temp_0"]
    rows.append(dict(name=os.path.basename(F)[-19:-4], dir=s, B=r["B"], seB=r["se_B"], Tc=r["Tc"], seTc=r["se_Tc"],
                     T0=T[0], T1=T[-1], upB=up[0], upTc=up[1], dnB=dn[0], dnTc=dn[1], loB=lo[0], hiB=hi[0],
                     slope=sl, levels=t, nseg=len(S)))
print(f"{'run':16s}{'dir':>4}{'B':>8}{'±':>6}{'Tc':>7}{'±':>6}{'T0':>7}{'T1':>7}{'upB':>7}{'dnB':>7}{'upTc':>7}{'dnTc':>7}{'B<=.12':>8}{'B>=.12':>8}{'d|t|/dq':>9}")
for x in rows:
    print(f"{x['name']:16s}{x['dir']:+4d}{x['B']:8.2f}{x['seB']:6.2f}{x['Tc']:7.3f}{x['seTc']:6.3f}{x['T0']:7.2f}{x['T1']:7.2f}"
          f"{x['upB']:7.2f}{x['dnB']:7.2f}{x['upTc']:7.3f}{x['dnTc']:7.3f}{x['loB']:8.2f}{x['hiB']:8.2f}{x['slope']:+9.4f}")
print("\n各段 |tau| [N*m]（12 段，速度 0.02→0.20→0.02）")
for x in rows: print(f"{x['name']:16s}", np.round(x["levels"], 3))
P = [x for x in rows if x["dir"] > 0]; N = [x for x in rows if x["dir"] < 0]
print("\n反轉 − 正轉（同一輪）|tau| 平均到 6 檔：")
for p, n in zip(P, N):
    dp = (p["levels"][:6] + p["levels"][6:][::-1]) / 2; dn_ = (n["levels"][:6] + n["levels"][6:][::-1]) / 2
    print(f"  {p['name'][-6:]} vs {n['name'][-6:]}:", np.round(dn_ - dp, 3))
for lab, G in (("正轉", P), ("反轉", N)):
    B = np.array([x["B"] for x in G]); Tc = np.array([x["Tc"] for x in G])
    print(f"{lab} 3 輪：B = {B.mean():.2f}（輪間 std {B.std(ddof=1):.2f}），Tc = {Tc.mean():.3f}（std {Tc.std(ddof=1):.3f}）")
os.makedirs(os.path.join(BASE, "analysis_1013"), exist_ok=True)
np.save(os.path.join(BASE, "analysis_1013", "constvel_rounds.npy"), rows, allow_pickle=True)
