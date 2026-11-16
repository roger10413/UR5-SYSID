# 本週上機計畫（暫定 2026/11/19）

這次要回答兩件事：

1. **J1 慣量隨頻率上升，是不是 4～6 Hz 彈性模態造成的？** 用低頻速度調變（`CONDITION_SET = "LOWF"`，0.5～2 Hz）加上紙杯砝碼，看 J 與 ΔJ 在低頻是不是變平、ΔJ 是不是接近 m·ρ²。
2. **摩擦往下一軸推進。** 補 J2 的 DIAG 對照，第一次跑 J3（手腕 1）定速。

另外在 DIAG 做敲擊測試，看 11/12 在 UP 量到的 4.2～4.3 Hz 會不會隨姿態改變。

J0 已知質量重做（README「J0 已知質量重做」）這次不排：J1、J2 用重力斜率已經得到換算比例約 12.1，J0 的換算比例等 J1 低頻的結果出來再決定怎麼做。

## 上機前（在實驗室電腦上）

- [ ] `git pull`，確認是包含 `ur5_j3_const_velocity_ident.py` 與 `CONDITION_SET` 的版本。
- [ ] 離線測試全部通過：`python -m unittest discover -s tests`（測試會封鎖網路，不會連到手臂）。
- [ ] 末端工具（FT300、夾爪、3D 相機）不拆；教導器負載 1.30 kg、重心 (1, 1, 74) mm。
- [ ] 紙杯砝碼（1.008 kg）與 11/12 同一個、同一種裝法；裝上後教導器負載改 2.31 kg（重心同 11/12）。
- [ ] 教導器速度滑桿 100%（每支程式都會檢查）；有人在緊急停止按鈕旁。
- [ ] 每支程式填 `ROBOT_IP`、`POSE`、`PAYLOAD_KG`、安全視窗（下表），`TEST_STAGE`、`DIRECTION` 每次跑之前改。

## 順序

依姿態排，盡量少換姿態與負載。每支新程式、新姿態、新條件組都照 QUICK → MID → FULL，每個方向跑完人工確認再改 `DIRECTION`。

| # | 姿態 | 負載 | 程式與設定 | 階段 | 預估時間 | 優先 |
|---|------|------|------------|------|----------|------|
| 1 | DIAG | 1.30 | `ur5_tap_test_record.py`：`POSE = "DIAG"`，`TAP_LABEL` = `none` → `fb` → `lr` | 各 14 s | 5 分 | 3 |
| 2 | DIAG | 1.30 | `ur5_j2_const_velocity_ident.py`：`POSE = "DIAG"`，`J2_SAFE_MIN_DEG/MAX` = 48 / 93 | QUICK±、MID±、FULL± | 8 分 | 2 |
| 3 | DIAG | 1.30 | `ur5_j3_const_velocity_ident.py`：`POSE = "DIAG"`，`J3_SAFE_MIN_DEG/MAX` = −131 / −87 | QUICK±、MID±、FULL± | 8 分 | 4 |
| 4 | UP | 1.30 | `ur5_j1_velocity_modulation_ident.py`：`CONDITION_SET = "LOWF"`，`J1_SAFE_MIN_DEG/MAX` = −114 / −70 | QUICK±、MID±、FULL± | 13 分 | **1** |
| 5 | UP | 2.31 | 同 #4（裝紙杯，`PAYLOAD_KG = 2.31`） | QUICK+、FULL± | 11 分 | **1** |
| 6 | UP, J3 = −130° | 1.30 | `ur5_j3_const_velocity_ident.py`：`POSE = "UP"`，`J3_SAFE_MIN_DEG/MAX` = −153 / −109 | QUICK±、MID±、FULL± | 8 分 | 2 |

運動時間合計約 55 分鐘；加上換姿態、換負載、每段人工確認，約 1.5～2 小時。時間不夠時依「優先」欄取捨：#4、#5 一定要做（同一天、負載相鄰才能比 ΔJ）。

### 換姿態

- **DIAG → UP（#3 → #4）**：以教導器只把 J2 從 90° 轉到 0°（J0 = 180°、J1 = −90°、J3 = −90°、J4 = −90°、J5 = 0°）。
- **#5 → #6**：拆紙杯、教導器負載改回 1.30；以教導器把 J3 從 −90° 轉到 −130°（工具朝前上方約 40°）。程式會檢查六軸 ±0.5°，忘了轉 J3 會中止。
- **結束**：J3 轉回 −90°、J2 轉回 90°，回到診斷原點。

## 各項說明

### #1 敲擊測試（DIAG）

- 做法同 11/12（README「敲擊測試」）：先 `none`（不敲）確認記錄正常，再 `fb`、`lr`；敲手腕 3 外殼或工具法蘭側面。
- 看什麼：11/12 在 UP，工具加速度的餘振 4.2～4.3 Hz（fb、lr 都一樣），加 1 kg 只降到 4.2。DIAG 時手臂的慣量分布不同；頻率跟著改變，表示是手臂本身的彈性；仍在 4.2～4.3 Hz，比較可能是底座或桌子。

### #2 J2 定速 DIAG

- 程式 11/12 在 UP 跑過（B 30.70、Tc 8.24）；DIAG 是新姿態，所以從 QUICK 開始。
- DIAG 窗口 [60.31°, 81.69°]，窗口內重力 +21.4～+23.4 N·m（前臂接近水平）；FULL 範圍 52.29°～89.71°，不低於診斷原點。
- 看什麼：Tc 與 UP 比（摩擦隨負載變多少；J1 是 9.26 vs 9.16，B 22.5 vs 19.8）；`K_mean`（target_moment 窗口平均約 −22.7，絕對值夠大）與紙杯重力斜率的 12.10、定速斜率的 11.65 比。

### #3、#6 J3 定速（新程式）

- 做法同 J1、J2 定速：窗口 6 個馬達圈，正反轉掃過同一段角度消去重力。J3 的重力小（最大約 3.6 N·m）。
- UP（主要）：J3 從 −130° 開始，窗口中心 −131.2° 為重力零點（窗口內 −0.66～+0.66 N·m），工具朝前上方 22°～60°，法蘭約 1.03 m。工具朝上，注意上方燈具與線架。
- DIAG（對照）：窗口中心 −109°（工具往前傾），窗口內重力 +2.98～+3.48 N·m；只往工具前傾一側轉，法蘭 0.43～0.51 m，不低於診斷原點。
- 法蘭離 J3 軸 0.125 m，0.20 rad/s 時 0.025 m/s；工具尖端約 0.06 m/s。
- **換算比例未知**：`KT_OUT` 暫用 J0～J2 的 13.64，J3 的 B、Tc 先當相對值。DIAG 的 `K_mean` 可作參考；UP 的重力斜率只有 J2 的 1/7，`K_slope` 誤差大（合成資料約 7%）。減速比 101 也是沿用，分析時要看漣波是否鎖定在 101 × 關節角。

### #4、#5 J1 低頻速度調變（`CONDITION_SET = "LOWF"`）

條件（巡航 0.08 rad/s，調變 4 s，每個方向 35 段、約 5.2 分鐘）：

| 條件 | f [Hz] | 速度振幅 [rad/s] | 加速度振幅 [rad/s²] | 週期數 |
|------|--------|------------------|---------------------|--------|
| 1 | 1.0 | 0.04 | 0.25 | 4 |
| 2 | 0.5 | 0.05 | 0.16 | 2 |
| 3 | 0.75 | 0.05 | 0.24 | 3 |
| 4 | 1.5 | 0.0265 | 0.25 | 6 |
| 5 | 2.0 | 0.02 | 0.25 | 8 |

- 條件 1 和 STD（11/5、11/12）的條件 1 相同，用來接上之前的資料。1～2 Hz 的加速度振幅都是 0.25，只有頻率不同。0.5、0.75 Hz 受最低速度 0.03 rad/s 限制，振幅只能到 0.05。
- 調變段長 18.33°（STD 為 13.75°）。FULL 範圍正轉 −106.73°～−79.82°、反轉 −104.18°～−77.27°，含 2° 餘裕後在 −114°～−70° 內；法蘭最低 0.96 m，線速度最高 0.14 m/s。
- 輸出資料夾 `j1velmod_up_p130_lowf_pos_full_…`（STD 的名稱不變）。
- 合成資料（J = 6，含重力）的精度：0.75～2 Hz 約 ±0.12～0.14，0.5 Hz 約 ±0.25 kg·m²。0.5 Hz 時重力造成的假慣量約 +7 kg·m²（比 J 還大），結果非常依賴重力項，只作參考；要看「重力改為各段斜率」寫法是否一致。

**判讀**：

- 單一彈性模態時，量到的 J(f) = J / (1 − (f/fn)²)。相對於低頻極限：

  | f [Hz] | 0.5 | 0.75 | 1 | 1.5 | 2 | 3 |
  |---|---|---|---|---|---|---|
  | fn = 5.6 Hz（11/12 擬合） | 1.008 | 1.018 | 1.033 | 1.077 | 1.146 | 1.40 |
  | fn = 4.3 Hz（敲擊） | 1.014 | 1.031 | 1.057 | 1.139 | 1.276 | 1.95 |

- 11/12 實測（p130）1 → 2 Hz 只多 4.4%（6.58 → 6.87），兩種 fn 預測的都是 11～21%；2 → 3 Hz 卻多 23%。1～2 Hz 這段的形狀是這次的重點：
  - 0.75～2 Hz 都平（差 < 2σ）、只有接近 3 Hz 才上升 → 低頻的 J 可以直接用，速度調變只用 ≤ 2 Hz。
  - 從 0.75 Hz 起就一路平滑上升 → 不像單一模態，要另外找原因（摩擦的動態特性、控制器）。
- 紙杯 ΔJ（同一天、同方向，p231 − p130）：低頻時應接近控制器的 ΔJ_model 0.857 ÷（K ÷ 13.64）≈ 0.94～0.97（K = 12.1～12.5，以 13.64 換算）。11/12 的 1 Hz 為 0.64 ± 0.19、0.92 ± 0.12，2 Hz 約 1.2，3 Hz 約 2.1～2.3。

## 跑完後的分析

```python
# J1 LOWF（每個方向、每種負載各一次；反轉用 -1）
import ur5_j1_velocity_modulation_ident as v1
r = v1.analyze_j1_velocity_modulation("FULL 的 CSV", 1, v1.KT_OUT, +1,
                                      v1.CONDITION_SETS["LOWF"]["conditions"], v1.V_CRUISE)

# J2 DIAG、J3 定速（姿態填實際跑的）
import ur5_j2_const_velocity_ident as j2
ta, tb = j2.window_bounds_deg(j2.POSES["DIAG"]["center_deg"], j2.N_MOTOR_REVS_WINDOW, j2.GEAR_RATIO)
rp = j2.analyze_angle_window("正轉 FULL 的 CSV", 2, j2.KT_OUT, +1, j2.SPEED_LEVELS, ta, tb)
rn = j2.analyze_angle_window("反轉 FULL 的 CSV", 2, j2.KT_OUT, -1, j2.SPEED_LEVELS, ta, tb)
j2.combine_directions(rp, rn)

import ur5_j3_const_velocity_ident as j3
ta, tb = j3.window_bounds_deg(j3.POSES["UP"]["center_deg"], j3.N_MOTOR_REVS_WINDOW, j3.GEAR_RATIO)
rp = j3.analyze_angle_window("正轉 FULL 的 CSV", 3, j3.KT_OUT, +1, j3.SPEED_LEVELS, ta, tb)
rn = j3.analyze_angle_window("反轉 FULL 的 CSV", 3, j3.KT_OUT, -1, j3.SPEED_LEVELS, ta, tb)
j3.combine_directions(rp, rn)

# 敲擊測試
import ur5_tap_test_record as tt
r = tt.analyze_tap_test("敲擊測試的 CSV")
```

每次跑完先看程式印出的：執行狀態 `done`、實測範圍在規劃範圍內、教導器速度滑桿最低 100%。任一項不符，那次資料不用，先查原因再跑下一個。
