"""11/26 六軸定速：升速趟 vs 降速趟的 B、Tc（看暖機造成的漂移）。

9 檔速度先升後降，每檔去、回各一段。每軸每段取窗口內平均力矩（電流 × K），
每個方向分「升速趟」「降速趟」各自回歸 τ = d·Tc + B·v + g（g 由正反合併消掉）。
"""
import os, sys
import numpy as np, pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "robot"))
import ur5_six_axis_ident as six  # noqa: E402

SRC = r"D:\cluade_code\ur5\data_1126\SIX\six_constvel_six_p130_full_20261126_090545\six_constvel_six_p130_full_20261126_090545.csv"


def segment_means(D, j, win, k):
    qt, q = D[f"target_qd_{j}"].values, np.degrees(D[f"actual_q_{j}"].values)
    tau = D[f"actual_current_{j}"].values * k
    temp = D[f"joint_temp_{j}"].values
    lo, hi = win
    moving = np.abs(qt) > 1e-6
    edges = np.flatnonzero(np.diff(moving.astype(int))) + 1
    starts, ends = edges[::2], edges[1::2]
    rows = []
    for a, b in zip(starts, ends):
        v = np.median(qt[a:b])
        m = np.zeros(len(qt), bool); m[a:b] = True
        m &= (q >= lo) & (q <= hi) & (np.abs(qt - v) < 1e-6)
        if m.sum() < 5:
            continue
        rows.append(dict(t=a, v=v, tau=tau[m].mean(), temp=temp[m].mean()))
    return pd.DataFrame(rows)


def fit(df):
    """τ = Tc·sign(v) + B·v + g"""
    X = np.column_stack([np.sign(df.v), df.v, np.ones(len(df))])
    c, *_ = np.linalg.lstsq(X, df.tau, rcond=None)
    return c  # Tc, B, g


if __name__ == "__main__":
    D = pd.read_csv(SRC)
    win = six.window_bounds(six.POSES["SIX"]["q_deg"])
    print(f"{'軸':>3}{'B 升':>8}{'B 降':>8}{'差':>7}{'Tc 升':>8}{'Tc 降':>8}{'溫度 升→降 °C':>16}")
    for j in range(6):
        df = segment_means(D, j, win[j], six.K_AXIS[j])
        # 每個 (速度, 方向) 第一次出現 = 升速趟，第二次 = 降速趟
        df["key"] = df.v.round(5)
        df["pass"] = df.groupby("key").cumcount()
        up, dn = df[df["pass"] == 0], df[df["pass"] == 1]
        cu, cd = fit(up), fit(dn)
        print(f"J{j} {cu[1]:8.2f}{cd[1]:8.2f}{(cd[1] / cu[1] - 1) * 100:+6.0f}%{cu[0]:8.3f}{cd[0]:8.3f}"
              f"{up.temp.mean():9.1f} → {dn.temp.mean():.1f}   段數 {len(up)}/{len(dn)}")
