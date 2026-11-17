# 離線分析：11/12 J2 定速 UP、速度調變 UP（不連線、不動手臂）
import sys, os, glob, json, csv, socket
socket.socket = None
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"D:\cluade_code\UR5-SYSID\robot")
import ur5_j2_const_velocity_ident as c2
import ur5_j2_velocity_modulation_ident as v2
D = r"D:\cluade_code\ur5\data_1112\J2"
one = lambda p: (lambda f: (f[0] if len(f) == 1 else (_ for _ in ()).throw(AssertionError(f))))(glob.glob(p))
K = c2.KT_OUT
lo, hi = c2.window_bounds_deg(c2.POSES["UP"]["center_deg"], 6, c2.GEAR_RATIO)
print(f"===== J2 定速 UP 窗口 {lo:.3f} ~ {hi:.3f}，K = {K:.3f}")
cv = {}
for d, tag in ((+1, "pos"), (-1, "neg")):
    print(f"--- {tag}")
    cv[tag] = c2.analyze_angle_window(one(os.path.join(D, f"j2constvel_up_{tag}_p130_full_*", "*.csv")), 2, K, d, c2.SPEED_LEVELS, lo, hi)
print("\n--- 合併")
comb = c2.combine_directions(cv["pos"], cv["neg"])
print({k: v for k, v in comb.items() if not isinstance(v, (list, dict))})
for tag in ("pos", "neg"):
    with open(one(os.path.join(D, f"j2constvel_up_{tag}_p130_full_*", "*.csv")), encoding="utf-8") as fh:
        rr = list(csv.DictReader(fh))
    print(tag, "J2 溫度", rr[0]["joint_temp_2"], "→", rr[-1]["joint_temp_2"], "；J1", rr[0]["joint_temp_1"], "→", rr[-1]["joint_temp_1"])
vm = {}
for d, tag in ((+1, "pos"), (-1, "neg")):
    path = one(os.path.join(D, f"j2velmod_up_p130_{tag}_full_*", "*.csv"))
    print(f"\n===== J2 速度調變 {tag}")
    r = v2.analyze_j2_velocity_modulation(path, 2, K, d, v2.CONDITIONS, v2.V_CRUISE)
    vm[tag] = {str(ci): dict(f=c["freq"], A=c["amp"], J=c["J"], se=c["se_J"], B=c.get("B_local"), gain=c.get("track_gain"),
                             temp=c.get("temp")) for ci, c in r["conditions"].items()}
    vm[tag]["J_model_mean"] = r.get("J_model_mean")
    vm[tag]["model"] = {str(k): v for k, v in r.get("model", {}).items()} if isinstance(r.get("model"), dict) else None
json.dump(dict(constvel={k: v for k, v in comb.items() if not isinstance(v, (list, dict))}, velmod=vm),
          open("j2_1112.json", "w", encoding="utf-8"), indent=1, ensure_ascii=False, default=float)
