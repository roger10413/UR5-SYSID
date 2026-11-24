# 11/19 手腕 J 為負：離線檢查「電流比角度晚」的假設（不連線、不動手臂）
#   1. 調變頻率上的相位差 → actual_current 相對 target_current、target_moment 晚多少 ms（每段擬合正弦）。
#   2. 把電流往前移 Δ（用 t + Δ 的電流，線性內插）再跑原本的速度調變分析，看 J 怎麼變。
#   資料每 8 ms 一筆（125 Hz）。J1（UP，1.30 與 2.31 kg）一起算，當對照。
import sys, os, glob, json, csv, io, contextlib, socket, types, importlib, math
socket.socket = None
import numpy as np
sys.stdout.reconfigure(encoding="utf-8")
sys.path.insert(0, r"F:\cluade_code\UR5-SYSID\robot")

D = r"D:\cluade_code\ur5\data_1119"
OUT = r"D:\cluade_code\ur5\analysis_1119"
K_AXIS = {1: 12.03, 3: 9.31, 4: 9.31, 5: 9.31}
SHIFTS_MS = [0, 4, 8, 12, 16, 20, 24]
RUNS = [(3, "diag_p130"), (4, "diag_p130"), (5, "diag_p130"), (1, "up_p130"), (1, "up_p231")]


def one(pattern):
    f = glob.glob(pattern)
    assert len(f) == 1, (pattern, f)
    return f[0]


def phasor(m, t, y, f):
    c, _ = m._fit_sine(t, y - np.polyval(np.polyfit(t, y, 1), t), f, 0.5 * (t[0] + t[-1]))
    return c[1] - 1j * c[2]


res = {}
for j, tag0 in RUNS:
    m = importlib.import_module(f"ur5_j{j}_velocity_modulation_ident")
    fn = getattr(m, f"analyze_j{j}_velocity_modulation")
    s_k = K_AXIS[j] / m.KT_OUT
    for d, tag in ((+1, "pos"), (-1, "neg")):
        key = f"J{j}_{tag0}_{tag}"
        path = one(os.path.join(D, f"J{j}", f"j{j}velmod_{tag0}_{tag}_full_*", "*.csv"))
        with open(path, encoding="utf-8") as fh:
            rows0 = list(csv.DictReader(fh))
        t = np.array([float(r["timestamp"]) for r in rows0])
        cur = np.array([float(r[f"actual_current_{j}"]) for r in rows0])
        tcur = np.array([float(r[f"target_current_{j}"]) for r in rows0])
        tmom = np.array([float(r[f"target_moment_{j}"]) for r in rows0])

        out = {"shift": {}}
        for sh in SHIFTS_MS:
            cs = np.interp(t + sh / 1000.0, t, cur)
            rows = [dict(r, **{f"actual_current_{j}": c}) for r, c in zip(rows0, cs)]
            m.csv = types.SimpleNamespace(DictReader=lambda fh, _rows=rows: _rows)
            with contextlib.redirect_stdout(io.StringIO()):
                r = fn(path, j, m.KT_OUT, d, m.CONDITIONS, m.V_CRUISE)
            out["shift"][sh] = {str(ci): dict(f=c["freq"], A=c["amp"], J=c["J"] * s_k, se=c["se_J"] * s_k,
                                              B=c["B_local"] * s_k, Jm=r["model"].get(ci, {}).get("J_model"))
                                for ci, c in r["conditions"].items()}
            if sh == 0:
                segs = r["segments"]
        m.csv = csv

        # 相位差（每段），依條件取中位數
        lag = {}
        for sg in segs:
            ii, f = sg["idx"], sg["freq"]
            pc, pt, pm = phasor(m, t[ii], cur[ii], f), phasor(m, t[ii], tcur[ii], f), phasor(m, t[ii], tmom[ii], f)
            l_ct = -np.angle(pc / pt) / (2 * math.pi * f) * 1000
            l_cm = -np.angle(pc / pm) / (2 * math.pi * f) * 1000
            lag.setdefault(str(sg["cond"]), []).append((l_ct, l_cm, sg.get("track_lag_ms", float("nan"))))
        out["lag"] = {ci: [float(np.median([x[k] for x in v])) for k in range(3)] for ci, v in lag.items()}
        res[key] = out
        print(f"完成 {key}", flush=True)

with open(os.path.join(OUT, "wrist_delay_1119.json"), "w", encoding="utf-8") as fh:
    json.dump(res, fh, ensure_ascii=False, indent=1, default=float)

# ---------- 印出 ----------
for key, out in res.items():
    print(f"\n===== {key}")
    print("條件              電流晚 target_current  電流晚 target_moment  角度晚目標 [ms]")
    base = out["shift"][0]
    for ci in sorted(base, key=lambda c: (base[c]["f"], base[c]["A"])):
        b = base[ci]
        l = out["lag"].get(ci, [float("nan")] * 3)
        print(f"  {b['f']:.0f} Hz A={b['A']:.3f}   {l[0]:>8.1f}            {l[1]:>8.1f}            {l[2]:>6.1f}")
    print("J 對電流前移量 [ms]：" + "".join(f"{sh:>8}" for sh in SHIFTS_MS) + "     J_model")
    for ci in sorted(base, key=lambda c: (base[c]["f"], base[c]["A"])):
        b = base[ci]
        print(f"  {b['f']:.0f} Hz A={b['A']:.3f}        " + "".join(f"{out['shift'][sh][ci]['J']:>8.3f}" for sh in SHIFTS_MS) + f"   {b['Jm'] if b['Jm'] is None else round(b['Jm'], 3)}")
