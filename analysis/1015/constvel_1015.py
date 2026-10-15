"""10/15 固定角度窗口定速 3 輪：B、Tc 對溫度，正反轉差，各速度檔力矩隨溫度的變化。"""
import sys, os, io, contextlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot")); sys.path.insert(0, HERE)
sys.stdout.reconfigure(encoding="utf-8")
import ur5_const_velocity_angle_ident as m
import data_1015 as dd

V6 = np.array([0.02, 0.05, 0.08, 0.12, 0.16, 0.20])


def analyze(runs):
    a, b = m.window_bounds_deg(180.0, 10, 101)
    rows = []
    for r in runs:
        if r["kind"] != "cv":
            continue
        with contextlib.redirect_stdout(io.StringIO()):
            x = m.analyze_angle_window(r["F"], 0, m.KT_OUT, r["dir"], m.SPEED_LEVELS, a, b)
        S = x["segments"]; v = np.abs([s["qd_mean"] for s in S]); t = np.abs([s["tau_mean"] for s in S])
        T = np.genfromtxt(r["F"], delimiter=",", names=True, usecols=(5,))["joint_temp_0"]
        up = np.polyfit(v[:6], t[:6], 1); dn = np.polyfit(v[6:], t[6:], 1)
        rows.append(dict(label=r["label"], ts=r["ts"], dir=r["dir"], n=r["n"], B=x["B"], seB=x["se_B"], Tc=x["Tc"], seTc=x["se_Tc"],
                         T0=T[0], T1=T[-1], Tm=T.mean(), upB=up[0], dnB=dn[0], upTc=up[1], dnTc=dn[1],
                         levels=(t[:6] + t[6:][::-1]) / 2))
    return rows


if __name__ == "__main__":
    rows = analyze(dd.run_list())
    print(f"{'實驗':14s}{'時間':>8}{'溫度':>14}{'B':>8}{'±':>6}{'Tc':>8}{'±':>7}{'升速B':>8}{'降速B':>8}")
    for x in rows:
        print(f"{x['label']:14s}{x['ts'][:2]}:{x['ts'][2:4]:>5}{x['T0']:7.2f}→{x['T1']:5.2f}{x['B']:8.2f}{x['seB']:6.2f}{x['Tc']:8.3f}{x['seTc']:7.3f}{x['upB']:8.2f}{x['dnB']:8.2f}")
    print("\n各速度檔 |τ| [N·m]（升降兩趟平均）：速度", V6)
    for x in rows:
        print(f"  {x['label']:14s}", np.round(x["levels"], 3))
    print("\n對溫度的線性回歸（3 輪）：")
    for s, lab in ((+1, "正轉"), (-1, "反轉")):
        R = [x for x in rows if x["dir"] == s]; T = np.array([x["Tm"] for x in R])
        pB = np.polyfit(T, [x["B"] for x in R], 1); pT = np.polyfit(T, [x["Tc"] for x in R], 1)
        L = np.array([x["levels"] for x in R]); pl = [np.polyfit(T, L[:, i], 1)[0] for i in range(6)]
        print(f"  {lab}：dB/dT = {pB[0]:+.2f} N·m·s/rad/°C（31 °C 時 B = {np.polyval(pB, 31):.2f}），dTc/dT = {pT[0]:+.3f} N·m/°C")
        print(f"        各速度檔 d|τ|/dT [N·m/°C]：{np.round(pl, 3)}")
    P = sorted([x for x in rows if x["dir"] > 0], key=lambda x: x["n"]); N = sorted([x for x in rows if x["dir"] < 0], key=lambda x: x["n"])
    print("\n反轉 − 正轉（同一輪）各速度檔 [N·m]：")
    for p, n_ in zip(P, N):
        print(f"  第 {p['n']} 輪：", np.round(n_["levels"] - p["levels"], 3))
    os.makedirs(os.path.join(dd.BASE, "analysis_1015"), exist_ok=True)
    np.save(os.path.join(dd.BASE, "analysis_1015", "constvel_1015.npy"), rows, allow_pickle=True)
