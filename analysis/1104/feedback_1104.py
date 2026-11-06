"""「命令減前饋」和追隨誤差的關係（10/29 速度調變資料，離線分析，不用上機）。

fb = joint_control_output − target_current（A），也就是控制器在前饋之外多加或扣掉的電流。
追隨誤差：ev = target_qd − actual_qd（速度誤差），ep = ∫ev dt（位置誤差，CSV 沒有 target_q，以每段起點為 0 積分；
起點的未知偏移由各段常數吸收）。每次、每個條件把 7 段合併做最小平方（各段常數 + 行程斜率），比較三種寫法：
  M0：fb = Ka·target_qdd                   （fb 只是扣掉一個和目標加速度成正比的項，與誤差無關）
  M1：fb = Kp·ep + Kv·ev                   （fb 是 PD 回授，全由追隨誤差決定）
  M2：fb = Kp·ep + Kv·ev + Ka·target_qdd   （兩者都有）
誤差相對 fb 的延遲（樣本數）在 −4～+4 之間取殘差最小者。Ka × 12.41（控制器的力矩/電流換算）換成 kg·m²：
若控制器在命令中直接扣掉轉子項，Ka × 12.41 應約 −1.915；若 fb 全是誤差造成的回授，Ka 約 0。
注意：調變是單一頻率的正弦，ep、ev、target_qdd 在調變頻率上只有兩個獨立分量（sin、cos），
所以 M2 的三個係數要靠調變頻率以外的成分（漣波、雜訊造成的誤差）才分得開，看標準差再判斷。
B. fb 的組成：每段依序扣掉常數與斜率、調變頻率、漣波（101·q 的諧波），看各部分占多少；
   再看調變頻率上 fb、ev 對 target_qdd 的比與相位，以及漣波各諧波的 fb 與 ev。
結果存到 OUT/feedback_1104.log 與 .json。
"""
import os, sys, io, csv, json, socket
import numpy as np

socket.socket = None                                   # 封鎖網路：只做離線分析
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "robot"))
sys.path.insert(0, os.path.join(HERE, "..", "1029"))
sys.stdout.reconfigure(encoding="utf-8")
import ur5_velocity_modulation_ident as vm
import known_mass_1029 as km

OUT = os.path.join(os.path.dirname(km.BASE), "analysis_1104"); os.makedirs(OUT, exist_ok=True)
K_CTRL = 12.41                                         # 控制器 J0 的力矩/電流換算（target_current 對 target_moment 擬合）
LAGS = range(-4, 5)
LOG = io.StringIO()


def say(*a):
    s = " ".join(str(x) for x in a)
    print(s)
    LOG.write(s + "\n")


def load(path):
    with open(path, encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) if r[k] not in ("", "None") else np.nan for r in rows]) for k in rows[0]}


def shift(x, k):
    """x 往後延 k 個樣本（k > 0：用 k 個樣本之前的誤差）；邊界補 NaN。"""
    y = np.full_like(x, np.nan)
    if k > 0:
        y[k:] = x[:-k]
    elif k < 0:
        y[:k] = x[-k:]
    else:
        y[:] = x
    return y


def fit(cs, y, cols, ng):
    """合併各段：y = Σ cols·係數 + 各段常數 + 各段行程斜率。回傳 (係數, 標準差, 殘差 RMS, R²)。"""
    Xs, ys = [], []
    for g, ii in enumerate(cs):
        one, zero = np.ones(len(ii)), np.zeros(len(ii))
        tt = np.arange(len(ii)) * vm.DT
        X = [c[ii] for c in cols] + [one if h == g else zero for h in range(ng)] + [tt if h == g else zero for h in range(ng)]
        Xs.append(np.column_stack(X)); ys.append(y[ii])
    X, Y = np.vstack(Xs), np.concatenate(ys)
    ok = np.all(np.isfinite(X), axis=1) & np.isfinite(Y)
    X, Y = X[ok], Y[ok]
    p, *_ = np.linalg.lstsq(X, Y, rcond=None)
    e = Y - X @ p
    # 每段的資料有自相關，用 HAC（沿用速度調變程式）估標準差
    XtX = X.T @ X
    cov = vm._hac_cov([X], [e], XtX, int(round(0.25 / vm.DT)))
    se = np.sqrt(np.diag(cov))
    r2 = 1 - np.sum(e ** 2) / np.sum((Y - np.mean(Y)) ** 2)
    return p[:len(cols)], se[:len(cols)], float(np.sqrt(np.mean(e ** 2))), float(r2)


def main():
    res = []
    for r in km.velmod_runs():
        D = load(r["F"])
        s = float(r["dir"])
        t = D["timestamp"]
        fb = s * (D["joint_control_output_0"] - D["target_current_0"])
        ev = s * (D["target_qd_0"] - D["actual_qd_0"])
        qdd = s * D["target_qdd_0"]
        segs, _ = vm.find_modulation_segments(t, s * D["target_qd_0"], vm.V_CRUISE, vm.CONDITIONS, vm.DT)
        say(f"\n{r['label']}（{len(segs)} 段）  fb 標準差 {np.nanstd(fb):.3f} A")
        say(f"  {'條件':<10}{'延遲':>4}{'M0 R²':>8}{'M1 R²':>8}{'M2 R²':>8}   {'M2：Ka×12.41 [kg·m²]':>22}"
            f"{'Kv [A/(rad/s)]':>18}{'Kp [A/rad]':>16}   {'M0：Ka×12.41':>14}")
        for ci, cn in enumerate(km.COND):
            cs = [sg["idx"] for sg in segs if sg["cond"] == ci]
            if len(cs) < 2:
                continue
            ng = len(cs)
            best = None
            for k in LAGS:
                evk = shift(ev, k)
                ep = np.full_like(ev, np.nan)
                for ii in cs:
                    ep[ii] = np.cumsum(np.nan_to_num(evk[ii])) * vm.DT
                m1 = fit(cs, fb, [ep, evk], ng)
                if best is None or m1[2] < best[1][2]:
                    best = (k, m1, ep, evk)
            k, m1, ep, evk = best
            m0 = fit(cs, fb, [qdd], ng)
            m2 = fit(cs, fb, [ep, evk, qdd], ng)
            (kp, kv, ka), (skp, skv, ska) = m2[0], m2[1]
            say(f"  {cn:<10}{k:>4}{m0[3]:8.3f}{m1[3]:8.3f}{m2[3]:8.3f}   {ka * K_CTRL:+10.3f} ± {ska * K_CTRL:.3f}"
                f"{kv:>11.1f} ± {skv:<5.1f}{kp:>9.0f} ± {skp:<5.0f}  {m0[0][0] * K_CTRL:+8.3f}")
            res.append(dict(label=r["label"], pose=r["pose"], w=r["w"], dir=r["dir"], cond=cn, lag=k,
                            r2=[m0[3], m1[3], m2[3]], rms=[m0[2], m1[2], m2[2]],
                            Ka_J=ka * K_CTRL, se_Ka_J=ska * K_CTRL, Kv=kv, se_Kv=skv, Kp=kp, se_Kp=skp,
                            Ka_J_m0=m0[0][0] * K_CTRL, Kv_m1=m1[0][1], Kp_m1=m1[0][0]))

    say("\n" + "=" * 90)
    J = np.array([x["Ka_J"] for x in res]); sJ = np.array([x["se_Ka_J"] for x in res])
    w = 1 / sJ ** 2
    m = float(np.sum(w * J) / np.sum(w)); sm = float(np.sqrt(1 / np.sum(w)))
    chi2 = float(np.sum(w * (J - m) ** 2))
    say(f"M2 的 Ka × 12.41：{len(J)} 個加權平均 {m:+.3f} ± {sm:.3f} kg·m²（卡方 {chi2:.1f}，自由度 {len(J) - 1}）；"
        f"範圍 {J.min():+.3f} ～ {J.max():+.3f}")
    for i in range(3):
        v = np.array([x["r2"][i] for x in res])
        say(f"M{i} 的 R²：{v.min():.3f} ～ {v.max():.3f}（中位數 {np.median(v):.3f}）")
    kv = np.array([x["Kv_m1"] for x in res]); kp = np.array([x["Kp_m1"] for x in res])
    say(f"M1 的 Kv：{kv.min():.1f} ～ {kv.max():.1f} A/(rad/s)；Kp：{kp.min():.0f} ～ {kp.max():.0f} A/rad；"
        f"延遲 {sorted(set(x['lag'] for x in res))}")

    # ---------------- B. fb 由哪些成分組成（反轉 4 次，每段分開） ----------------
    say("\n" + "=" * 90)
    say("B. fb 的組成：每段依序扣掉 常數 + 斜率 → 調變頻率 1、2 倍 → 漣波（101·q 的 24 個諧波），看剩下的標準差 [A]")
    say("   調變相量：fb、ev 在調變頻率上對 target_qdd 的比（fb 換成 kg·m²：× 12.41；ev 換成延遲 ms）")
    say("=" * 90)
    out_b = []
    for r in km.velmod_runs():
        if r["dir"] > 0:
            continue
        D = load(r["F"]); s = float(r["dir"]); t = D["timestamp"]
        fb = s * (D["joint_control_output_0"] - D["target_current_0"])
        ev = s * (D["target_qd_0"] - D["actual_qd_0"])
        qdd, q = s * D["target_qdd_0"], s * D["actual_q_0"]
        segs, _ = vm.find_modulation_segments(t, s * D["target_qd_0"], vm.V_CRUISE, vm.CONDITIONS, vm.DT)
        say(f"\n{r['label']}")
        say(f"  {'條件':<10}{'全部':>7}{'去常數斜率':>9}{'去調變':>8}{'去漣波':>8}   {'fb/qdd [kg·m²]':>16}{'ev/qdd [ms]':>14}")
        for ci, cn in enumerate(km.COND):
            sd, ph = [], []
            for sg in (sg for sg in segs if sg["cond"] == ci):
                ii = sg["idx"]; tt = t[ii] - t[ii[0]]; w = 2 * np.pi * sg["freq"]; phi = vm.GEAR_RATIO * q[ii]
                base = [np.ones(len(ii)), tt]
                mod = [np.cos(w * tt), np.sin(w * tt), np.cos(2 * w * tt), np.sin(2 * w * tt)]
                rip = [f(k * phi) for k in range(1, 25) for f in (np.cos, np.sin)]
                y = fb[ii]
                row = [np.std(y)]
                for cols in (base, base + mod, base + mod + rip):
                    X = np.column_stack(cols)
                    row.append(np.std(y - X @ np.linalg.lstsq(X, y, rcond=None)[0]))
                sd.append(row)
                X = np.column_stack(base + mod[:2])
                pz = [np.linalg.lstsq(X, v[ii], rcond=None)[0][2:4] for v in (fb, ev, qdd)]
                ph.append([complex(p[0], -p[1]) for p in pz])
            sd = np.mean(sd, axis=0); F, E, Q = np.mean(ph, axis=0)
            say(f"  {cn:<10}{sd[0]:7.3f}{sd[1]:9.3f}{sd[2]:8.3f}{sd[3]:8.3f}   "
                f"{abs(F / Q) * K_CTRL:7.2f} ∠{np.degrees(np.angle(F / Q)):+4.0f}°{abs(E / Q) * 1e3:8.2f} ∠{np.degrees(np.angle(E / Q)):+4.0f}°")
            out_b.append(dict(label=r["label"], cond=cn, std=sd.tolist(), fb_qdd_J=abs(F / Q) * K_CTRL,
                              fb_qdd_deg=float(np.degrees(np.angle(F / Q))), ev_qdd_ms=abs(E / Q) * 1e3,
                              ev_qdd_deg=float(np.degrees(np.angle(E / Q)))))

        # 漣波：整次實驗的所有段合併，fb、ev 對 101·q 的諧波
        ii = np.concatenate([sg["idx"] for sg in segs]); phi = vm.GEAR_RATIO * q[ii]
        X = np.column_stack([np.ones(len(ii))] + [f(k * phi) for k in range(1, 13) for f in (np.cos, np.sin)])
        pf = np.linalg.lstsq(X, fb[ii], rcond=None)[0]; pe = np.linalg.lstsq(X, ev[ii], rcond=None)[0]
        say("  漣波（和位置鎖定）：fb 振幅 > 0.01 A 的諧波")
        for k in range(1, 13):
            Fk = complex(pf[2 * k - 1], -pf[2 * k]); Ek = complex(pe[2 * k - 1], -pe[2 * k])
            if abs(Fk) > 0.01:
                say(f"    k = {k:2d}（巡航時約 {k * vm.GEAR_RATIO * vm.V_CRUISE / 2 / np.pi:4.1f} Hz）fb {abs(Fk):.3f} A  "
                    f"ev {abs(Ek) * 1e3:.3f} mrad/s  fb/ev {abs(Fk / Ek):7.1f} A/(rad/s) ∠{np.degrees(np.angle(Fk / Ek)):+4.0f}°")

    with open(os.path.join(OUT, "feedback_1104.log"), "w", encoding="utf-8") as fh:
        fh.write(LOG.getvalue())
    with open(os.path.join(OUT, "feedback_1104.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(rows=res, Ka_J_mean=m, se_Ka_J_mean=sm, chi2=chi2, parts=out_b), fh, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
