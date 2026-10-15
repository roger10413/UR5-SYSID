"""10/15 實驗資料：檔案清單、讀檔、運動段特徵（方向、起步後行程、加減速階段、馬達角）。"""
import os, glob
import numpy as np

BASE = r"D:\cluade_code\ur5"   # 資料資料夾（內含各次 constvel_angle_*、constaccel_v3_*、constaccel_v3trap_* 子資料夾），請改成自己的路徑
KT_OUT = 101 * 0.1350           # 輸出端力矩常數 [N·m/A]（Kt × 減速比）
GEAR = 101
DT = 1 / 125.0
A_PHASE = 0.05                  # 指令加速度門檻 [rad/s²]：大於為加速段、小於負值為減速段
SETTLE_N = 5                    # 換階段後不用的樣本數（40 ms）

KIND_NAME = {"cv": "定速", "tri": "三角形 v3", "trap": "梯形 v3"}


def run_list(day="20261015"):
    """依時間排序的 FULL 實驗：dict(kind = cv / tri / trap, dir = ±1, ts = 時間戳, F = CSV 路徑, label)。"""
    pats = [("cv", f"constvel_angle_*_full_{day}_*"), ("tri", f"constaccel_v3_*_full_{day}_*"),
            ("trap", f"constaccel_v3trap_*_full_{day}_*")]
    out = []
    for kind, p in pats:
        for d in glob.glob(os.path.join(BASE, p)):
            F = glob.glob(os.path.join(d, "*_data_*.csv"))[0]
            out.append(dict(kind=kind, dir=+1 if "_pos_" in os.path.basename(d) else -1, ts=os.path.basename(d)[-6:], F=F))
    out.sort(key=lambda r: r["ts"])
    cnt = {}
    for r in out:
        cnt[(r["kind"], r["dir"])] = cnt.get((r["kind"], r["dir"]), 0) + 1
        r["n"] = cnt[(r["kind"], r["dir"])]
        nth = {"cv": f"第{r['n']}輪", "tri": f" 第{r['n']}次", "trap": ""}[r["kind"]]
        r["label"] = f"{KIND_NAME[r['kind']]}{nth} {'正轉' if r['dir'] > 0 else '反轉'}"
    return out


def load(F):
    d = np.genfromtxt(F, delimiter=",", names=True, usecols=(0, 1, 2, 3, 4, 5))
    return dict(t=d["timestamp"], q=d["actual_q_0"], qd=d["actual_qd_0"], tau=d["actual_current_0"] * KT_OUT,
                tq=d["target_qd_0"], T=d["joint_temp_0"])


def motion_features(D, run_dir):
    """逐樣本特徵。sgn：指令速度方向（靜止為 0）；travel：本段起步後行程 [°]；tsince：起步後時間 [s]；
    phase：0 靜止、1 加速、2 定速、3 減速（依指令加速度）；since：距上次換階段的樣本數；
    main：是否為本次實驗的主運動（方向與實驗方向相同；反向的是段間 movej 回程）。"""
    tq, q = D["tq"], D["q"]; n = len(tq)
    mv = np.abs(tq) > 1e-6
    sgn = np.where(mv, np.sign(tq), 0.0)
    travel = np.full(n, np.nan); tsince = np.full(n, np.nan); seg = np.full(n, -1)
    e = np.flatnonzero(np.diff(np.r_[0, mv.astype(int), 0]))
    for k, (a, b) in enumerate(zip(e[0::2], e[1::2])):
        a0 = max(a - 1, 0)
        seg[a:b] = k
        travel[a:b] = np.degrees(np.sign(tq[a]) * (q[a:b] - q[a0]))
        tsince[a:b] = (np.arange(a, b) - a0) * DT
    am = sgn * np.gradient(tq, DT)
    phase = np.where(~mv, 0, np.where(am > A_PHASE, 1, np.where(am < -A_PHASE, 3, 2)))
    since = np.zeros(n, int)
    for i in range(1, n):
        since[i] = since[i - 1] + 1 if phase[i] == phase[i - 1] else 0
    return dict(sgn=sgn, travel=travel, tsince=tsince, seg=seg, phase=phase, since=since,
                main=(sgn == run_dir), phi=GEAR * q, am_cmd=am)
