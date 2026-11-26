"""11/26 六軸速度調變：用控制器自己的模型力矩（target_moment）跑同一套分析，檢查方法本身。

target_moment 是 UR 控制器依自己的動力學模型算出的力矩，含慣量項、不含摩擦。
把它當成「電流 × K」丟進同一個分析流程：
  - 若算出的 J ≈ J_model（名目模型），代表六軸同時動時分析流程沒問題，實測 J 偏低／為負是真的；
  - 若連 target_moment 都算不出 J_model，問題在分析流程（耦合、修邊、重力補償）。
"""
import os, sys, csv
import numpy as np

sys.stdout.reconfigure(encoding="utf-8")
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(ROOT, "robot"))
import ur5_six_axis_ident as six  # noqa: E402

SRC = r"D:\cluade_code\ur5\data_1126\SIX\six_velmod_six_p130_full_20261126_092018\six_velmod_six_p130_full_20261126_092018.csv"
OUT = os.path.join(os.path.dirname(SRC), "tm_as_current.csv")


def make_tm_csv(src=SRC, out=OUT, k_axis=six.K_AXIS):
    """actual_current_j ← target_moment_j / K_j（分析會再乘回 K_j）"""
    with open(src, encoding="utf-8") as fi, open(out, "w", newline="", encoding="utf-8") as fo:
        rd = csv.DictReader(fi)
        wr = csv.DictWriter(fo, fieldnames=rd.fieldnames)
        wr.writeheader()
        for r in rd:
            for j in range(6):
                r[f"actual_current_{j}"] = float(r[f"target_moment_{j}"]) / k_axis[j]
            wr.writerow(r)
    return out


if __name__ == "__main__":
    if not os.path.exists(OUT):
        make_tm_csv()
    res = six.analyze_velmod(OUT, payload_kg=1.30)
    print("\n".join(six.summary_lines("VELMOD", res)))
