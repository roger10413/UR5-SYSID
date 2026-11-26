"""
UR5 六軸暖機（以 B 收斂當停止條件）—— 暖機完成後接著跑 ur5_six_axis_ident.py
================================================================

★★★ 想法 ★★★

11/26 六軸定速 4 輪：冷機第 1 輪的 B 比熱機高 8～19%，Tc 只差約 2%；而且 J1 溫度還在升時 B 已經不動，
所以「暖機完成」不能只看關節溫度，要看 B 本身有沒有停止變化。

這支程式的暖機動作就是一小段六軸定速鑑別（和 ur5_six_axis_ident.py 的 CONSTVEL 同一個姿態、同一條路徑、同一個窗口）：
  每一批 = BLOCK_SPEEDS（0.05～0.20 rad/s 先升後降）× 去程、回程，約 2 分鐘。
每批跑完、手臂停下後，用和鑑別程式相同的分析（analyze_constvel，正反合併）算出六軸的 B、Tc，
再決定要不要送下一批。暖機過程本身就是一條「B 隨時間、溫度變化」的紀錄（論文 4.2 溫度一節可用）。

結束條件（暖機完成）：
  1. 已暖機至少 MIN_MINUTES 分鐘，且
  2. 最近 N_STABLE_BLOCKS 批，CONVERGE_JOINTS（J0～J2）每一軸 B 的（最大 − 最小）／平均 < B_STABLE_FRAC
  未達成但累計超過 MAX_MINUTES 也結束（註明「未達穩定」）。結束時回到姿態 SIX，直接接著跑鑑別程式。
  手腕三軸（J3～J5）的 B 照樣每批算、印出、存檔，但不當停止條件。

門檻與判定軸的來源：11/26 四輪資料只取 0.05～0.20 的段落離線重算（analysis/1126/six_warmup_offline_1126.py）：
  - 冷機 → 熱機，暖機批的 B 降 16～21%（J0～J2）、11～31%（J3～J5）。
  - 熱機後三批之間的差異：J1 3.0%、J2 1.0%（J0 6.7%，但那是還在單調下降 27.4 → 25.9 → 25.6，不是雜訊）；
    手腕 J3 8%、J4 5%、J5 12%，半批之間 J5 甚至差到 25% —— 手腕 B 小、Tc 占大部分，2 分鐘一批估不準斜率，
    拿來當停止條件會一直判定未收斂。手腕在冷機第 1 輪後也已穩定（整輪 B 的變異 1.6～5%），和大軸同時變熱。
  - 暖機批只用 0.05～0.20，B 比鑑別程式（0.005～0.20 九檔）低 2～10%（大軸）、15～28%（手腕）（低速段斜率較大），
    所以暖機批的 B 只拿來和前後批比，不和鑑別結果比。

★★★ 使用前必讀 ★★★
1. 【新程式，尚未上機】先用 TEST_STAGE = "QUICK"（只跑 1 批、不判定），確認動作、記錄與每批分析正常。
   動作和 ur5_six_axis_ident.py CONSTVEL 完全相同（同一個 plan_constvel、build_urscript、安全檢查），只是速度檔不同。
2. 【溫度上限 TEMP_ABORT_C 未設定一律中止】任一批結束時任一軸溫度 ≥ 上限就停止、回姿態。
3. 【六軸視窗、法蘭高度等】沿用 ur5_six_axis_ident.py 的設定與 safety_check。
4. 【姿態】程式不會擺姿態：先用 ur5_goto_pose.py SIX；連上 RTDE 後以實際角度檢查六軸（±0.5°）。
5. 【教導器速度滑桿】每批送出前檢查，不是 100% 或讀不到就停止。
6. 【Ctrl+C】會送 stopj；最可靠的仍是緊急停止按鈕。
"""

import csv
import json
import math
import os
import sys
import time
from datetime import datetime

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import ur5_six_axis_ident as six   # noqa: E402  規劃、安全檢查、URScript、定速分析

vm = six.vm
rtde_receive = six.rtde_receive


# ============================================================
# 使用者設定
# ============================================================

ROBOT_IP = "192.168.50.114"   # 執行前請再次確認（第二台 UR5 只改這裡）
# QUICK：只跑 1 批，不判定收斂 —— 確認動作、記錄、每批分析
# FULL ：一批接一批，直到六軸 B 收斂或達 MAX_MINUTES
TEST_STAGE = "QUICK"

POSE = "SIX"
BLOCK_SPEEDS = [0.05, 0.08, 0.12, 0.16, 0.20, 0.20, 0.16, 0.12, 0.08, 0.05]   # rad/s，先升後降：批內的漂移在平均中抵銷

MIN_MINUTES = 6.0
MAX_MINUTES = 60.0
N_STABLE_BLOCKS = 3
CONVERGE_JOINTS = (0, 1, 2)   # 當停止條件的軸（手腕 B 每批估計雜訊 5～12%，只記錄）
B_STABLE_FRAC = 0.03          # 最近 N_STABLE_BLOCKS 批 B 的（最大 − 最小）／平均 上限（熱機雜訊 J1 3.0%、J2 1.0%）

TEMP_ABORT_C = None           # 任一軸溫度上限 [°C]，由操作者依手冊與現場決定（11/26 熱機後最高約 40 °C）

OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "outputs")


# ============================================================
# 每批分析與收斂判定（純函式，可離線測試）
# ============================================================

def block_temperatures(csv_path):
    """每軸溫度：檔案內樣本平均與最大值。"""
    sums, cnt, mx = np.zeros(6), np.zeros(6), np.full(6, -np.inf)
    with open(csv_path, encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            for j in range(6):
                x = r.get(f"joint_temp_{j}")
                if x in (None, "", "None"):
                    continue
                v = float(x)
                if math.isfinite(v):
                    sums[j] += v; cnt[j] += 1; mx[j] = max(mx[j], v)
    mean = [float(s / c) if c else float("nan") for s, c in zip(sums, cnt)]
    return mean, [float(m) if math.isfinite(m) else float("nan") for m in mx]


def block_result(csv_path, speeds=BLOCK_SPEEDS, center_deg=None):
    """一批的六軸 B、Tc（同 ur5_six_axis_ident.analyze_constvel，已乘 K_AXIS / KT_OUT）與溫度。"""
    center_deg = six.POSES[POSE]["q_deg"] if center_deg is None else center_deg
    r = six.analyze_constvel(csv_path, center_deg=center_deg, speeds=sorted(set(speeds)))
    nan = float("nan")
    T, Tmax = block_temperatures(csv_path)
    return dict(B=[r[j].get("B", nan) for j in range(6)], sd_B=[r[j].get("sd_B", nan) for j in range(6)],
                Tc=[r[j].get("Tc", nan) for j in range(6)], n_seg=[r[j]["n_pos"] + r[j]["n_neg"] for j in range(6)],
                T=T, T_max=Tmax)


def is_converged(history, elapsed_min, min_minutes=MIN_MINUTES, n_stable=N_STABLE_BLOCKS, frac=B_STABLE_FRAC,
                 joints=CONVERGE_JOINTS):
    """
    history：每批一筆 dict（至少有 B：六軸）。回傳 (是否完成, 說明, 六軸 spread)。
    spread_j = 最近 n_stable 批 B_j 的（最大 − 最小）／平均；六軸都算，只有 joints 用來判定。
    """
    spreads = [float("nan")] * 6
    if len(history) >= n_stable:
        for j in range(6):
            vals = [h["B"][j] for h in history[-n_stable:]]
            if all(math.isfinite(v) for v in vals):
                spreads[j] = (max(vals) - min(vals)) / abs(float(np.mean(vals)))
    if elapsed_min < min_minutes:
        return False, f"未滿最少暖機時間 {min_minutes:.0f} 分鐘", spreads
    if len(history) < n_stable:
        return False, f"批數不足 {n_stable}", spreads
    missing = [f"J{j}" for j in joints if not math.isfinite(spreads[j])]
    if missing:
        return False, "B 缺資料：" + "、".join(missing), spreads
    bad = [f"J{j} {100 * spreads[j]:.1f}%" for j in joints if spreads[j] >= frac]
    if bad:
        return False, f"未收斂（門檻 {100 * frac:.0f}%）：" + "、".join(bad), spreads
    return True, ("、".join(f"J{j}" for j in joints) + f" 的 B 最近 {n_stable} 批變化都 < {100 * frac:.0f}%"
                  f"（" + "、".join(f"J{j} {100 * spreads[j]:.1f}%" for j in joints) + "）"), spreads


def temperature_abort(result, limit_c):
    """回傳超過上限的軸（空清單 = 沒有）。"""
    return [j for j in range(6) if math.isfinite(result["T_max"][j]) and result["T_max"][j] >= limit_c]


def format_row(k, el, res, why):
    return (f"{k:3d}{el:7.1f} | " + " ".join(f"{b:6.2f}" for b in res["B"]) + " | "
            + " ".join(f"{t:4.1f}" for t in res["T"]) + f" | {why}")


# ============================================================
# 主流程
# ============================================================

def main():
    print("=" * 70)
    print(f" UR5 六軸暖機（B 收斂）：TEST_STAGE = {TEST_STAGE}")
    print("=" * 70)
    if TEST_STAGE not in ("QUICK", "FULL"):
        print(f"[中止] 未知的 TEST_STAGE {TEST_STAGE!r}"); return
    pose_q = list(six.POSES[POSE]["q_deg"])
    segs = six.plan_constvel(BLOCK_SPEEDS, pose_q)
    t_block = six.estimate_total_time(segs, pose_q, "CONSTVEL", pose_q)
    max_blocks = 1 if TEST_STAGE == "QUICK" else max(1, int(math.ceil(MAX_MINUTES * 60.0 / t_block)))
    print(f"每批 {len(segs)} 段（速度 {BLOCK_SPEEDS}），預估 {t_block / 60:.1f} 分鐘；最多 {max_blocks} 批")
    print(f"收斂：至少 {MIN_MINUTES:.0f} 分鐘，且最近 {N_STABLE_BLOCKS} 批 "
          + "、".join(f"J{j}" for j in CONVERGE_JOINTS) + f" 的 B 變化都 < {100 * B_STABLE_FRAC:.0f}%（手腕只記錄）")

    ok, report, (lo, hi) = six.safety_check(segs, "CONSTVEL", pose_q)
    print("\n".join(report))
    if TEMP_ABORT_C is None:
        ok = False; print("[FAIL] TEMP_ABORT_C 尚未設定，拒絕執行")
    else:
        print(f"[OK] 溫度上限 {TEMP_ABORT_C:.1f} °C")
    if not ok:
        print("\n[中止] 安全檢查未通過。"); return
    if ROBOT_IP is None or rtde_receive is None:
        print("\n[中止] ROBOT_IP 未設定或找不到 rtde_receive，僅完成離線設計檢查。"); return

    print(f"\n[注意] 六軸須先擺到 {pose_q}°（可用 ur5_goto_pose.py SIX）；六軸會同時在下列範圍內移動：")
    for j in range(6):
        print(f"   J{j}: {lo[j]:.1f}° ～ {hi[j]:.1f}°")
    if input(f"已確認周圍淨空、線材不會被拉扯，輸入 {POSE} 繼續：").strip() != POSE:
        print("已取消。"); return

    tag = f"six_warmup_{TEST_STAGE.lower()}_{datetime.now():%Y%m%d_%H%M%S}"
    out_dir = os.path.join(OUTPUT_DIR, tag)
    os.makedirs(out_dir, exist_ok=True)
    script = six.build_urscript(segs, "CONSTVEL", end_deg=pose_q)
    with open(os.path.join(out_dir, "program.script"), "w", encoding="utf-8") as fh:
        fh.write(script)

    history, result, t0, logger = [], "未完成", None, None
    print(f"\n{'批':>3}{'分鐘':>7} | " + " ".join(f"{'B' + str(j):>6}" for j in range(6)) + " | "
          + " ".join(f"{'T' + str(j):>4}" for j in range(6)) + " | 判定")
    try:
        for b in range(max_blocks):
            csv_path = os.path.join(out_dir, f"block_{b + 1:02d}.csv")
            logger = vm.RtdeLogger(ROBOT_IP, six.SAMPLE_HZ, csv_path, 0)
            logger.start()
            if not logger.started_ok.wait(timeout=10.0):
                result = f"中止（RTDE 沒有資料：{logger.error}）"; break
            q_now = [math.degrees(x) for x in logger.last_q_full]
            checks = []
            if b == 0:
                p_ok, p_lines = vm.check_pose(logger.last_q_full, POSE, six.POSES, six.POSE_TOL_DEG)
                checks.append((p_ok, "\n".join(p_lines)))
            checks.append(six.check_start_position(q_now, segs[0]["start_deg"]))
            checks.append(vm.check_speed_fraction(logger.last_speed_fraction))
            failed = [msg for good, msg in checks if not good]
            if b == 0:
                print("\n".join(msg for _, msg in checks))
            if failed:
                print(f"\n[停止] 第 {b + 1} 批送出前檢查未通過：\n" + "\n".join(failed))
                result = "停止（送出前檢查未通過）"; break
            if t0 is None:
                t0 = time.time()
            vm.send_urscript(ROBOT_IP, script)
            status, _, _ = vm.wait_for_motion_complete(logger, six.estimate_total_time(segs, q_now, "CONSTVEL", pose_q))
            time.sleep(0.3)
            logger.stop(); logger = None
            if status != "done":
                vm.send_abort(ROBOT_IP)
                result = f"中止（第 {b + 1} 批 {status}）"; break
            res = block_result(csv_path)
            el = (time.time() - t0) / 60.0
            res.update(block=b + 1, t_min=el)
            history.append(res)
            warm, why, spreads = is_converged(history, el)
            res["spread"] = spreads
            print(format_row(b + 1, el, res, why))
            hot = temperature_abort(res, TEMP_ABORT_C)
            if hot:
                result = "停止（溫度達上限：" + "、".join(f"J{j} {res['T_max'][j]:.1f} °C" for j in hot) + "）"; break
            if TEST_STAGE == "FULL" and warm:
                result = "暖機完成"; break
        else:
            result = "QUICK 完成" if TEST_STAGE == "QUICK" else f"達最長時間 {MAX_MINUTES:.0f} 分鐘，未達穩定"
    except KeyboardInterrupt:
        vm.send_abort(ROBOT_IP)
        result = "中止（Ctrl+C）"
        print("\n[中止] 使用者中斷。")
    except BaseException:
        vm.send_abort(ROBOT_IP)
        result = "中止（例外）"
        raise
    finally:
        if logger is not None:
            logger.stop()
        with open(os.path.join(out_dir, "warmup_history.json"), "w", encoding="utf-8") as fh:
            json.dump(dict(result=result, block_speeds=BLOCK_SPEEDS, b_stable_frac=B_STABLE_FRAC,
                           converge_joints=CONVERGE_JOINTS, min_minutes=MIN_MINUTES, history=history),
                      fh, ensure_ascii=False, indent=1)
        with open(os.path.join(out_dir, "summary.txt"), "w", encoding="utf-8") as fh:
            fh.write(f"模式: {TEST_STAGE}\n結果: {result}\n")
            fh.write("批,分鐘," + ",".join(f"B{j}" for j in range(6)) + "," + ",".join(f"T{j}" for j in range(6)) + "\n")
            for h in history:
                fh.write(f"{h['block']},{h['t_min']:.2f}," + ",".join(f"{x:.4f}" for x in h["B"]) + ","
                         + ",".join(f"{x:.2f}" for x in h["T"]) + "\n")

    print(f"\n[結果] {result}")
    print(f"[存檔] {out_dir}")
    if result == "暖機完成":
        print("接著請直接執行 ur5_six_axis_ident.py（MODE = CONSTVEL、TEST_STAGE = FULL），避免間隔過久又變冷。")


if __name__ == "__main__":
    main()
