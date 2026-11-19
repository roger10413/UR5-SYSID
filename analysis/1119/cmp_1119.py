import json, sys
sys.stdout.reconfigure(encoding="utf-8")
A = json.load(open(r"D:\cluade_code\ur5\analysis_1112\j1_1112.json", encoding="utf-8"))
B = json.load(open(r"D:\cluade_code\ur5\analysis_1119\analyze_1119.json", encoding="utf-8"))
s = 12.03 / 13.64
for pay in ("p130", "p231"):
    print("==", pay)
    for ci in ("0","1","2","3","4","5"):
        row = []
        for date in ("1105", "1112"):
            k = f"velmod_{date}_{pay}_pos"
            if k in A:
                P, N = A[k]["per"][ci], A[f"velmod_{date}_{pay}_neg"]["per"][ci]
                row.append(f"{date}: {(P['J']+N['J'])/2*s:.3f}")
        P, N = B[f"J1_velmod_up_{pay}_pos"]["per"][ci], B[f"J1_velmod_up_{pay}_neg"]["per"][ci]
        row.append(f"1119: {(P['J']+N['J'])/2:.3f}")
        print(f"  {P['f']:.0f} Hz A={P['A']:.3f}  " + "  ".join(row))
J2 = json.load(open(r"D:\cluade_code\ur5\analysis_1112\j2_1112.json", encoding="utf-8"))["constvel"]
s2 = 12.10 / 13.64
print("J2 UP 1112: B", J2["B"]*s2, "Tc", J2["Tc"]*s2, "H", J2["H"]*s2)
print("J2 DIAG 1119: B", B["J2_constvel_diag"]["B"], "Tc", B["J2_constvel_diag"]["Tc"], "H", B["J2_constvel_diag"]["H"])
