# 11/19 結果摘要（讀 analyze_1119.json，另與 11/12 的 J1、J2 結果並列）
import json, sys, os, math
sys.stdout.reconfigure(encoding="utf-8")
H = r"D:\cluade_code\ur5\analysis_1119"
R = json.load(open(os.path.join(H, "analyze_1119.json"), encoding="utf-8"))

print("===== 定速 DIAG（換算比例見 K）")
for j in (2, 3, 4, 5):
    c = R[f"J{j}_constvel_diag"]
    print(f"J{j}: K={c['K_used']}  B={c['B']:.3f}±{c['sd_B']:.3f}  Tc={c['Tc']:.3f}±{c['sd_Tc']:.3f}  "
          f"H={c['H']:.3f}  B+={c['B_pos']:.3f} B-={c['B_neg']:.3f}  C+={c['C_pos']:.3f} C-={c['C_neg']:.3f}  "
          f"tm_win={c['tm_win']}  K_mean={c['K_mean']}  K_slope={c['K_slope']}  temp={c['temps']}")

print("\n===== 速度調變（各條件 J ± se，方向 pos / neg）")
for key in ("J3_velmod_diag_p130", "J4_velmod_diag_p130", "J5_velmod_diag_p130", "J1_velmod_up_p130", "J1_velmod_up_p231"):
    P, N = R[key + "_pos"], R[key + "_neg"]
    print(f"-- {key}  K={P['K_used']}  temp pos {P['temps']} neg {N['temps']}  J_model_mean {P['J_model_mean']} / {N['J_model_mean']}")
    for ci in sorted(P["per"], key=int):
        p, n = P["per"][ci], N["per"].get(ci, {})
        print(f"   {p['f']:.1f} Hz A={p['A']:.3f}: J+ {p['J']:.4f}±{p['se']:.4f}  J- {n.get('J', float('nan')):.4f}±{n.get('se', float('nan')):.4f}"
              f"  Jm {p.get('Jm')} / {n.get('Jm')}  gain {p.get('gain')}")

print("\n===== J1 ΔJ（2.31 − 1.30，同條件同方向平均）")
for ci in sorted(R["J1_velmod_up_p130_pos"]["per"], key=int):
    a = [R[f"J1_velmod_up_p130_{t}"]["per"][ci]["J"] for t in ("pos", "neg")]
    b = [R[f"J1_velmod_up_p231_{t}"]["per"][ci]["J"] for t in ("pos", "neg")]
    f = R["J1_velmod_up_p130_pos"]["per"][ci]["f"]; A = R["J1_velmod_up_p130_pos"]["per"][ci]["A"]
    print(f"   {f:.1f} Hz A={A:.3f}: J(1.30)={sum(a)/2:.4f}  J(2.31)={sum(b)/2:.4f}  ΔJ={(sum(b)-sum(a))/2:.4f}")

print("\n===== 敲擊 DIAG")
for lab in ("none", "fb", "lr"):
    t = R[f"tap_diag_{lab}"]
    print(lab, {k: v for k, v in t.items() if not isinstance(v, (list, dict))})
