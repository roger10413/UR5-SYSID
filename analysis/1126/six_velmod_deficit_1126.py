"""11/26 六軸速度調變：實測 J 比控制器模型（target_moment 算出的 J）少多少，和局部摩擦斜率 B_local 的關係。

假設摩擦落後速度 τ 秒：F(v(t−τ)) ≈ F(v) − B_local·τ·a，會把 J 低估 B_local·τ。
若各軸的 ΔJ / B_local 差不多，表示同一種機制。
"""
import os, sys
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "robot"))
import ur5_six_axis_ident as six  # noqa: E402

D = r"D:\cluade_code\ur5\data_1126\SIX\six_velmod_six_p130_full_20261126_092018"
MEAS = os.path.join(D, "six_velmod_six_p130_full_20261126_092018_gcomp.csv")
TM = os.path.join(D, "tm_as_current_gcomp.csv")      # six_axis_check_1126.py 產生

if __name__ == "__main__":
    rm = six.analyze_velmod(MEAS)
    rt = six.analyze_velmod(TM)
    print(f"{'軸':>3}{'f':>6}{'A':>7}{'J實測':>9}{'J模型力矩':>10}{'ΔJ':>8}{'B_local':>9}{'ΔJ/B [ms]':>11}")
    for j in range(6):
        taus = []
        for ci in sorted(rm[j]["per"]):
            f, A = rm[j]["conds"][ci]
            pm, pt = rm[j]["per"][ci], rt[j]["per"][ci]
            Jm = np.mean([pm[d]["J"] for d in pm]); Jt = np.mean([pt[d]["J"] for d in pt])
            B = np.mean([pm[d]["B"] for d in pm])
            dJ = Jt - Jm
            taus.append(dJ / B * 1e3)
            print(f"J{j} {f:6.2f}{A:7.3f}{Jm:9.3f}{Jt:10.3f}{dJ:8.3f}{B:9.2f}{taus[-1]:11.0f}")
        print(f"   J{j} ΔJ/B 中位數 {np.median(taus):.0f} ms")
