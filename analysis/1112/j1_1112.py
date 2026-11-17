# 離線分析：11/12 J1 定速 DIAG、速度調變 UP 重複、速度調變 UP 夾紙杯砝碼（不連線、不動手臂）
import sys, os, glob, json, math, socket, io, contextlib
socket.socket = None                                   # 封鎖網路
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\cluade_code\UR5-SYSID\robot")
import ur5_j1_velocity_modulation_ident as v1
import ur5_j1_const_velocity_ident as c1

D12 = r"D:\cluade_code\ur5\data_1112"
D05 = r"D:\cluade_code\ur5\data_1105"
OUT = r"D:\cluade_code\ur5\analysis_1112"
K = v1.KT_OUT
res = {}

def one(pattern):
    f = glob.glob(pattern)
    assert len(f) == 1, (pattern, f)
    return f[0]

# ---------- 1. J1 定速 DIAG ----------
lo, hi = c1.window_bounds_deg(c1.POSES["DIAG"]["center_deg"], 6, c1.GEAR_RATIO)
print(f"===== J1 定速 DIAG，窗口 {lo:.4f} ~ {hi:.4f}")
cv = {}
for d, tag in ((+1, "pos"), (-1, "neg")):
    print(f"--- {tag}")
    cv[tag] = c1.analyze_angle_window(one(os.path.join(D12, f"j1constvel_diag_{tag}_p130_full_*", "*.csv")),
                                      1, K, d, c1.SPEED_LEVELS, lo, hi)
comb = c1.combine_directions(cv["pos"], cv["neg"])
temps = {}
for tag in ("pos", "neg"):
    import csv
    with open(one(os.path.join(D12, f"j1constvel_diag_{tag}_p130_full_*", "*.csv")), encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
    temps[tag] = (float(rr[0]["joint_temp_1"]), float(rr[-1]["joint_temp_1"]))
print("J1 溫度（開始, 結束）:", temps)
res["constvel_diag"] = dict(B=comb["B"], sd_B=comb["sd_B"], B_pos=comb["B_pos"], B_neg=comb["B_neg"],
                            Tc=comb["Tc"], sd_Tc=comb["sd_Tc"], H=comb["H"], C_pos=cv["pos"]["C"], C_neg=cv["neg"]["C"],
                            tm_win=comb["tm_win"], K_mean=comb["K_mean"], K_slope=comb["K_slope"],
                            slope_tau=comb["slope_tau"], slope_tm=comb["slope_tm"], temps=temps,
                            seg=[(s["speed"], s["qd_mean"], s["tau_mean"], s["tm_mean"]) for t in ("pos", "neg") for s in cv[t]["segments"]])

# ---------- 2. J1 速度調變 ----------
def tm_fit_full(cs, tm, tqdd, tq, q):
    """同 v1._fit_target_moment，但回傳全部係數（J、B、G_s、G_c、各段常數）。"""
    Xs, ys = [], []
    for g, sg in enumerate(cs):
        ii = sg["idx"]
        Dm = np.zeros((len(ii), len(cs))); Dm[:, g] = 1.0
        Xs.append(np.column_stack([tqdd[ii], tq[ii], np.sin(q[ii]), np.cos(q[ii]), Dm])); ys.append(tm[ii])
    X, y = np.vstack(Xs), np.concatenate(ys)
    p, *_ = np.linalg.lstsq(X, y, rcond=None)
    return p

runs = {
    "1105_p130_pos": (D05, "j1velmod_up_p130_pos_full_*", +1, 1.30),
    "1105_p130_neg": (D05, "j1velmod_up_p130_neg_full_*", -1, 1.30),
    "1112_p130_pos": (D12, "j1velmod_up_p130_pos_full_*", +1, 1.30),
    "1112_p130_neg": (D12, "j1velmod_up_p130_neg_full_*", -1, 1.30),
    "1112_p231_pos": (D12, "j1velmod_up_p231_pos_full_*", +1, 2.31),
    "1112_p231_neg": (D12, "j1velmod_up_p231_neg_full_*", -1, 2.31),
}
for name, (base, pat, d, pay) in runs.items():
    path = one(os.path.join(base, pat, "*.csv"))
    print(f"\n===== {name}  {os.path.basename(path)}")
    r = v1.analyze_j1_velocity_modulation(path, 1, K, d, v1.CONDITIONS, v1.V_CRUISE)
    # 原始欄位（絕對座標），算整體平均位準與控制器模型的完整係數
    Dd = np.genfromtxt(path, delimiter=",", names=True)
    q = Dd["actual_q_1"]; cur = Dd["actual_current_1"] * K; tm = Dd["target_moment_1"]
    tqdd = d * Dd["target_qdd_1"]; tq = d * Dd["target_qd_1"]
    segs = r["segments"]
    allidx = np.concatenate([sg["idx"] for sg in segs])
    per = {}
    for ci in sorted(r["conditions"]):
        c = r["conditions"][ci]
        cs = [sg for sg in segs if sg["cond"] == ci]
        p = tm_fit_full(cs, d * tm, tqdd, tq, q)
        idx = np.concatenate([sg["idx"] for sg in cs])
        per[ci] = dict(f=c["freq"], A=c["amp"], J=c["J"], se=c["se_J"], G_s=d * c["G_s"], G_c=d * c["G_c"],
                       a_amp=c["a_amp"], gain=c["track_gain"], lag=c["track_lag_ms"], temp=c["temp"], B=c["B_local"],
                       Jm=float(p[0]), Gs_m=d * float(p[2]), Gc_m=d * float(p[3]),
                       lvl_cur=float(np.mean(cur[idx])), lvl_tm=float(np.mean(tm[idx])), q_mean=float(np.degrees(np.mean(q[idx]))))
    res["velmod_" + name] = dict(d=d, pay=pay, per=per, lvl_cur=float(np.mean(cur[allidx])), lvl_tm=float(np.mean(tm[allidx])),
                                 q_mean=float(np.degrees(np.mean(q[allidx]))),
                                 temp=(float(Dd["joint_temp_1"][0]), float(Dd["joint_temp_1"][-1])))

with open(os.path.join(OUT, "j1_1112.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1, default=float)
print("\n存檔:", os.path.join(OUT, "j1_1112.json"))
