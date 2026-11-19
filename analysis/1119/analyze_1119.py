# 離線分析：11/19 上機全部資料（不連線、不動手臂）
#   J2 定速 DIAG、J3 定速／速度調變 DIAG、J4／J5 定速／速度調變 DIAG、J1 速度調變 UP（STD，1.30 與 2.31 kg）、敲擊測試 DIAG
# 用程式內的 KT_OUT（13.64，文獻值）算，最後再依各軸換算比例縮放（B、Tc、J 都和換算比例成正比）。
import sys, os, glob, json, csv, socket, io, contextlib
socket.socket = None                                   # 封鎖網路
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
ROBOT = r"F:\cluade_code\UR5-SYSID\robot"
sys.path.insert(0, ROBOT)
import importlib

D = r"D:\cluade_code\ur5\data_1119"
OUT = r"D:\cluade_code\ur5\analysis_1119"
K_AXIS = {1: 12.03, 2: 12.10, 3: 9.31, 4: 9.31, 5: 9.31}   # 各軸換算比例（J3～J5 暫定）


def one(pattern):
    f = glob.glob(pattern)
    assert len(f) == 1, (pattern, f)
    return f[0]


def temps(path, j):
    with open(path, encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
    return float(rr[0][f"joint_temp_{j}"]), float(rr[-1][f"joint_temp_{j}"])


def scal(x, s):
    return None if x is None else float(x) * s


res, log = {}, io.StringIO()


def run(label, fn):
    """執行分析並把印出的內容同時存進 log。"""
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        r = fn()
    log.write(f"\n\n########## {label}\n" + buf.getvalue())
    return r


# ---------- 定速（J2～J5，DIAG） ----------
for j in (2, 3, 4, 5):
    m = importlib.import_module(f"ur5_j{j}_const_velocity_ident")
    K = m.KT_OUT
    lo, hi = m.window_bounds_deg(m.POSES["DIAG"]["center_deg"], 6, m.GEAR_RATIO)
    cv, tp = {}, {}
    for d, tag in ((+1, "pos"), (-1, "neg")):
        path = one(os.path.join(D, f"J{j}", f"j{j}constvel_diag_{tag}_p130_full_*", "*.csv"))
        cv[tag] = run(f"J{j} 定速 DIAG {tag}", lambda: m.analyze_angle_window(path, j, K, d, m.SPEED_LEVELS, lo, hi))
        tp[tag] = temps(path, j)
    comb = run(f"J{j} 定速 DIAG 合併", lambda: m.combine_directions(cv["pos"], cv["neg"]))
    s = K_AXIS[j] / K
    res[f"J{j}_constvel_diag"] = dict(
        window=(lo, hi), K_used=K_AXIS[j],
        B=scal(comb["B"], s), sd_B=scal(comb["sd_B"], s), Tc=scal(comb["Tc"], s), sd_Tc=scal(comb["sd_Tc"], s),
        H=scal(comb["H"], s), B_pos=scal(comb["B_pos"], s), B_neg=scal(comb["B_neg"], s),
        C_pos=scal(cv["pos"]["C"], s), C_neg=scal(cv["neg"]["C"], s),
        tm_win=comb.get("tm_win"), K_mean=comb.get("K_mean"), K_slope=comb.get("K_slope"), temps=tp,
        seg=[(t, sg["speed"], sg["qd_mean"], sg["tau_mean"] * s, sg.get("tm_mean")) for t in ("pos", "neg") for sg in cv[t]["segments"]])

# ---------- 速度調變（J3～J5 DIAG、J1 UP） ----------
vm_runs = [(3, "diag_p130", 3), (4, "diag_p130", 4), (5, "diag_p130", 5),
           (1, "up_p130", 1), (1, "up_p231", 1)]
for j, tag0, _ in vm_runs:
    m = importlib.import_module(f"ur5_j{j}_velocity_modulation_ident")
    K = m.KT_OUT
    fn = getattr(m, f"analyze_j{j}_velocity_modulation")
    conds = m.CONDITIONS          # 今天 J1 也是 STD（info：條件組 STD）
    s = K_AXIS[j] / K
    for d, tag in ((+1, "pos"), (-1, "neg")):
        path = one(os.path.join(D, f"J{j}", f"j{j}velmod_{tag0}_{tag}_full_*", "*.csv"))
        r = run(f"J{j} 速度調變 {tag0} {tag}", lambda: fn(path, j, K, d, conds, m.V_CRUISE))
        per = {}
        for ci, c in r["conditions"].items():
            per[str(ci)] = dict(f=c["freq"], A=c["amp"], J=scal(c["J"], s), se=scal(c["se_J"], s),
                                B=scal(c.get("B_local"), s), gain=c.get("track_gain"), lag=c.get("track_lag_ms"),
                                Jm=c.get("J_model"), temp=c.get("temp"))
        res[f"J{j}_velmod_{tag0}_{tag}"] = dict(K_used=K_AXIS[j], per=per, temps=temps(path, j),
                                                 J_model_mean=r.get("J_model_mean"))

# ---------- 敲擊測試 DIAG ----------
import ur5_tap_test_record as tt
for lab in ("none", "fb", "lr"):
    path = one(os.path.join(D, "TAP", f"taptest_diag_p130_{lab}_*", "*.csv"))
    r = run(f"敲擊 DIAG {lab}", lambda: tt.analyze_tap_test(path))
    log.write("\n".join(tt.tap_analysis_lines(r)) + "\n")
    res[f"tap_diag_{lab}"] = {k: v for k, v in r.items() if isinstance(v, (int, float, str, list, dict))}

with open(os.path.join(OUT, "analyze_1119.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
with open(os.path.join(OUT, "analyze_1119_log.txt"), "w", encoding="utf-8") as fh:
    fh.write(log.getvalue())
print("存檔:", OUT)
