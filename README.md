# UR5-SYSID：UR5 CB3 J0 摩擦與慣量鑑別

以 UR5（CB3）底座關節 J0 為對象，鑑別單軸摩擦與等效慣量，作為關節健康監測（PHM）的基準。

模型：

```
tau = J·qdd + B·qd + Tc·sign(qd),   tau = i × Kt × 101
```

- **定速法**：單一方向、穩態（qdd ≈ 0），多個定速值做線性回歸 → `B`、`Tc`
- **定加速度法**：同一速度下「加速段 − 減速段」消去摩擦 → `J`。10/13 實測減速段摩擦比加速段多一個與加速度無關的常數，故以多個加速度檔做含截距回歸（斜率 = `J`）
- **速度調變法**：巡航速度上疊加正弦速度（速度不歸零、不反向），力矩中與加速度同相的成分 → `J`，與速度同相的成分 → 局部摩擦斜率；以 6 種頻率、振幅檢查 `J` 是否隨頻率或振幅改變（10/22 跑前 4 種，10/29 跑 6 種）。可改在手臂外伸姿態（`POSE = "EXT45"`）執行，比較兩個姿態的 `J`；並記錄控制器的 `target_moment`，得到控制器模型的慣量 `J_model` 作對照。末端加已知質量前後各跑一次時，兩次 `J` 的差應等於 m·r²（r 為質量塊重心離 J0 軸的距離），可直接檢查力矩換算比例
- **J1 定速法**（2026/11 起）：J1 為水平軸，定速時電流含重力力矩；正反轉以定速掃過同一段角度，兩方向合併消去重力 → J1 的 `B`、`Tc`

本 repo 只含程式與文件，不含真機量測資料。

## ⚠️ 安全須知

`robot/` 內的腳本**會透過 URScript（port 30002）直接驅動真實機械臂**。

- 到 commit `735deea` 為止的上機程式都已通過安全審查（最近一次上機 2026/10/29：暖機、定速、速度調變 6 個條件、J1 靜止電流）。之後的修改**尚未經安全審查**：輸出資料夾改存 `outputs/`（commit `3e8e8ee`，六支程式各改一行 `OUTPUT_DIR`）。
- `ur5_velocity_modulation_ident.py`（2026/10/19 新增）：前 4 個條件與 `EXT45` 為 commit `72591a7`，2026/10/22 上機；條件 5、6、負載確認、教導器速度滑桿檢查為 commit `130be1d`，2026/10/29 上機（`DIAG`、`EXT45`，含夾砝碼）。新的條件或姿態仍請依 `TEST_STAGE` = QUICK → MID → FULL 執行。
  - 姿態 `POSE` 預設 `DIAG`（診斷原點）。`EXT45`（J1、J2 = −45°、45°，法蘭離 J0 軸約 0.80 m）時程式**不移動 J1~J5**：須先以教導器擺好、確認 J0 掃過的扇形範圍淨空並輸入姿態名稱；連上 RTDE 後以實際角度檢查 J1~J5（±0.5°），不符或姿態名稱未知即中止。結束後 J1~J5 仍在外伸姿態，請以教導器移回。
  - 負載 `PAYLOAD_KG`：未設定（`None`）或不是 ≥ 0 的數值即中止；送出前要輸入教導器上設定的負載質量，與 `PAYLOAD_KG` 相差超過 0.005 kg 即中止。資料夾名稱含負載標記（1.30 kg → `p130`）。加減質量後必須先在教導器更新負載質量與重心。
- `ur5_j1_static_current.py`（2026/10/27 新增，commit `735deea`，2026/10/29 上機）：J0、J2~J5 不動，J1 來回 ±0.5°（0.05 rad/s），從兩個方向靠近同一角度後各靜止 4 s。六軸姿態（含 J0 = 180°）以實際角度檢查、負載確認、速度滑桿檢查同速度調變；移動角度、速度、加速度有上限，結束時回到開始位置。
- `ur5_j1_const_velocity_ident.py`（2026/11/02 新增，**尚未審查、尚未上機**）：J1 固定角度窗口定速。J0、J2~J5 不動，J1 只在診斷原點往回抬的一側活動（FULL 為 −127.7°～−90.3°，不低於診斷原點），結束時回到 −90°。J1 安全視窗（`J1_SAFE_MIN_DEG` / `J1_SAFE_MAX_DEG`，未設定即中止）；以 UR5 官方 DH 計算規劃範圍內的法蘭高度，最低須 ≥ 0.40 m（診斷原點 0.43 m）；六軸姿態（診斷原點，含 J0 = 180°）以實際角度檢查；須輸入 `J1` 確認淨空；負載確認、速度滑桿檢查同速度調變。
- 教導器速度滑桿（2026/10/26 起，六支上機程式都有）：連上 RTDE 後讀 `target_speed_fraction`，不是 100% 或讀不到即中止（10/22 冷機定速因滑桿停在 52% 而作廢）；暖機每批送出前再檢查一次。執行中的值逐列寫進 CSV（`target_speed_fraction` 欄），最低值寫進摘要，低於 100% 時結束後警告這次的資料不可用。
- 在能連到機台的電腦上，請勿直接執行主程式；離線檢查請只用 `py_compile` 或 `tests/`（測試會封鎖網路並使用假的 RTDE）。
- 執行前務必確認 `ROBOT_IP`、J0 安全視窗（`J0_SAFE_MIN_DEG` / `J0_SAFE_MAX_DEG`；J1 定速為 `J1_SAFE_MIN_DEG` / `J1_SAFE_MAX_DEG`，未設定時程式會中止）、手臂周圍淨空，並有人在緊急停止按鈕旁。
- 正轉、反轉分開執行，每個方向結束後人工確認。

## 目錄

```
robot/                                   上機腳本（會送指令給機械臂）
  ur5_const_velocity_angle_ident.py      固定角度窗口定速實驗：所有速度、正反轉都掃過同一段 J0 角度（10 個馬達圈）
  ur5_j0_warmup.py                       暖機：J0 來回轉動，溫度與摩擦指標穩定後自動結束
  ur5_const_accel_ident_v3.py            定加速度 v3：三角形速度（定速段 0）＋起始角錯開 1/20 馬達圈，並含離線分析函式
  ur5_const_accel_ident_v3_trapezoid.py  定加速度 v3 梯形對照版：只改定速段為 4 個馬達圈與輸出前綴 constaccel_v3trap_
  ur5_velocity_modulation_ident.py       速度調變實驗：巡航 0.12 rad/s 上疊加正弦（6 種頻率、振幅 × 7 次，起點錯開 1/7 馬達圈），並含離線分析函式
  ur5_j1_static_current.py               J1 靜止電流：兩方向靠近同一角度消去靜摩擦，夾砝碼前後的電流差校正 Kt × 減速比，並含離線分析函式
  ur5_j1_const_velocity_ident.py         J1 固定角度窗口定速：正反轉掃過同一段 J1 角度（6 個馬達圈）消去重力 → J1 的 B、Tc，並含離線分析函式
  ur5_home_pose.py                       診斷原點姿態檢查（預設只檢查、不移動）
tests/                                   離線測試（假 RTDE、封鎖網路）
docs/
  EXPERIMENT_PLAN_ANGLE.md               固定角度窗口定速實驗規劃
  EXPERIMENT_PLAN_ACCEL_V3.md            定加速度 v3 規劃與合成資料驗證
analysis/
  1013/                                  10/13 實驗分析（Python）：定速 3 輪比較、定加速度 v3、報告圖
  1015/                                  10/15 實驗分析（Python）：定速溫度效應、三角形 vs 梯形 v3、卡爾曼濾波殘差分解
  1022/                                  10/22 實驗分析（Python）：速度調變兩姿態 J、控制器 J_model、力矩換電流比值
  1029/                                  10/29 實驗分析（Python）：已知質量校正 J0 換算比例、控制器 J_model 組成、J1 靜止電流
matlab/
  constvel/                              定速法分析與共用函式
  analysis_1007/                         定速重複實驗與定加速度（v2）分析
```

## J0 已知質量重做（下一次上機）

10/29 第一次做已知質量（2F-85 夾 1 kg 圓柱砝碼，`DIAG`、`EXT45`），結果見 `analysis/1029/`。當天夾砝碼時教導器負載與 `PAYLOAD_KG` 設成 1.40 kg（應為 2.30）：J0 的 `J` 由實測電流計算，不受影響；J1 靜止電流可能受影響。這次用正確的負載重做，確認換算比例隔週能否重複，並重做 J1 靜止電流。

執行前：

- 教導器速度滑桿調到 100%（程式也會檢查）。
- 末端工具（FT300、夾爪、3D 相機；教導器負載 1.30 kg、重心 (1, 1, 74) mm）不拆。
- 砝碼和 10/29 一樣橫向夾在 Robotiq 2F-85 指間（摩擦夾持額定 5 kg；夾持力設最大，砝碼表面乾淨無油），另用束帶從砝碼底下兜過、兩端固定在夾爪本體，滑落時不會掉下。末端工具加砝碼約 2.3 kg，低於 UR5 額定負載 5 kg。
- 夾好後拍照（旁邊放尺），量砝碼兩端到指尖中心線的距離：橫向夾時砝碼重心可能偏離工具軸線，會改變 r。
- 夾上、放下砝碼後，先在教導器更新負載質量（1.30 + m），再改 `PAYLOAD_KG`；程式會要求輸入教導器上的負載質量，不符即中止。重心可不改（工具垂直朝下時，重心高度不影響 J0、J1 的結果）。
- 砝碼有機會就秤重；`analysis/1029/known_mass_1029.py` 的第一個參數就是砝碼質量。

順序（速度調變只跑反轉；每一列速度調變跑完、J0 回到 180° 後，同姿態、同負載接著跑 `ur5_j1_static_current.py` FULL）：

| # | 姿態 | 法蘭負載（`PAYLOAD_KG`） | 程式 |
|---|------|------|------|
| 1 | DIAG | 末端工具（—） | `ur5_j0_warmup.py` FULL |
| 2 | DIAG | 末端工具（—） | `ur5_const_velocity_angle_ident.py` FULL 正、反轉（例行追蹤，時間不夠可省略） |
| 3 | DIAG | 末端工具（1.30） | `ur5_velocity_modulation_ident.py` FULL 反轉 → J1 靜止 |
| 4 | DIAG | 末端工具 + 砝碼（1.30 + m） | 同上 |
| 5 | EXT45 | 末端工具 + 砝碼（1.30 + m） | 同上 |
| 6 | EXT45 | 末端工具（1.30） | 同上 |

- #5 前以教導器把 J1、J2 擺到 −45°、45°，量砝碼最低點離桌面的距離，並確認 J0 掃過的扇形（反轉約 157°～207°）淨空；#6 結束後移回診斷原點（程式不動 J1~J5）。
- 加減質量的兩次在時間上相鄰（#3 與 #4、#5 與 #6），溫度變化對差值的影響較小。
- 時間夠就加速度調變正轉（`DIRECTION = +1`，J0 約 153°～203°）：先 #3、#4，再 #5、#6。
- 速度調變 FULL 每個方向約 8.6 分鐘；必做項目含暖機與換裝約 1.5～2 小時。

## J1 定速（`ur5_j1_const_velocity_ident.py`，尚未審查）

J1 是水平軸，定速時電流除了摩擦，還包含撐住手臂重力的力矩 H(q1)。沿用 J0 的固定角度窗口：所有速度、正反轉都以定速掃過同一段 J1 角度（6 個馬達圈，[−119.69°, −98.31°]），各段取窗口內的平均，重力項對每一段都相同。每個方向分開回歸 `tau = B·qd + C·d`：正轉 `C₊ = Tc + H̄`、反轉 `C₋ = Tc − H̄`，兩方向合併得 `Tc = (C₊ + C₋)/2` 與窗口平均重力保持力矩 `H̄ = (C₊ − C₋)/2`，不需要連桿質量模型（`combine_directions`）。另記錄控制器的 `target_moment`（官方連桿質量 + 教導器負載的重力模型），`H̄` 與其窗口平均之比、窗口內兩者對角度的斜率之比，可作為 J1 換算比例的參考（推測用：前提是控制器的重力模型正確）。

動作範圍（J0 = 180°、J2~J5 為診斷原點；J1 角度減少時手臂往回抬）：

- J1 在 −127.7°～−90.3° 之間活動（FULL；QUICK 為 −112.4°～−90°），不低於診斷原點，結束時回到 −90°。
- 法蘭離底座安裝面 0.43 m（−90°）→ 0.66 m（−127.7°）；手肘最多移到底座後方約 0.25 m（高度約 0.43 m）。
- 重力力矩（官方連桿質量 + 末端工具 1.30 kg）在 −119.5° 附近為零，窗口內 −0.2～+17.1 N·m。

執行前：

- 以教導器擺到診斷原點（J0 = 180°）；確認手臂上方與底座後方淨空、線材不會被拉扯。
- 填 `J1_SAFE_MIN_DEG` / `J1_SAFE_MAX_DEG`（未設定即中止）：須包含規劃範圍加 2° 餘裕，FULL 為 −129.71°～−88.00°（例如 −131、−87）。
- `PAYLOAD_KG = 1.30`（不夾砝碼）；教導器速度滑桿 100%。

順序：QUICK 正轉 → QUICK 反轉 → MID 正、反轉 → FULL 正、反轉；每個方向跑完人工確認後再改 `DIRECTION`。FULL 每個方向 12 段、約 2.4 分鐘。正、反轉 FULL 都跑完後：

```python
import ur5_j1_const_velocity_ident as j1
ta, tb = j1.window_bounds_deg(j1.WINDOW_CENTER_DEG, j1.N_MOTOR_REVS_WINDOW, j1.GEAR_RATIO)
rp = j1.analyze_angle_window("正轉 FULL 的 CSV 路徑", 1, j1.KT_OUT, +1, j1.SPEED_LEVELS, ta, tb)
rn = j1.analyze_angle_window("反轉 FULL 的 CSV 路徑", 1, j1.KT_OUT, -1, j1.SPEED_LEVELS, ta, tb)
j1.combine_directions(rp, rn)
```

- 正轉（手臂下降）時摩擦與重力方向相反，窗口內馬達力矩由正變負（減速機受力面換邊）；反轉不換邊。
- 若摩擦隨傳遞的力矩改變，兩方向的負載不同，`Tc` 與 `H̄` 會互相混入（推測）；夾砝碼前後各跑一次，`H̄` 的差只來自已知質量，可另外檢查。

## 離線測試

```bash
pip install numpy
python -m unittest discover -s tests -v
```

## 速度調變實驗的分析

```python
import ur5_velocity_modulation_ident as vm
r = vm.analyze_velocity_modulation(csv_path, 0, vm.KT_OUT, +1, vm.CONDITIONS, vm.V_CRUISE)   # 反轉用 -1
```

- 每個條件分開回歸：`J`、`B_local` ± 標準差（考慮段內殘差自相關）、實測加速度振幅、追隨增益與延遲；並以不同模型寫法（無漣波項、只取主要諧波、無 v²、行程項改為共用曲線、無行程項）重算 `J`。
- 跨條件：`J` 加權平均與卡方（各條件是否一致）、`J` 對 1/A 回歸（減速時摩擦較大時 `J` 會隨振幅變化）、`B_local` 對 ω² 回歸（電流相對編碼器的延遲）。
- 每個條件至少要 2 段才回歸，QUICK、MID 的資料只列出各段的頻率、振幅與追隨情形。
- 條件以擬合的頻率、振幅辨認，`vm.CONDITIONS`（6 個條件）也可分析 10/22 的資料（只有前 4 個條件），資料中沒有的條件會略過。
- CSV 有 `target_moment`、`target_qdd` 欄位時，另以 `target_moment` 對 `target_qdd`、`target_qd` 回歸，得到控制器模型在該姿態的 J0 慣量 `J_model`（目標值沒有量測雜訊，QUICK 也可算）。`J ÷ J_model` 不是 1 時，可能是 `Kt` 不準，也可能是控制器模型與實際不同，須比較多個姿態才能分開。

## Python 分析（`analysis/1013/`）

- `constvel_rounds.py`：固定角度窗口定速各次的 `B`、`Tc`、升降趟、力矩對角度斜率、正反轉差。
- `accel_v3_analysis.py`：定加速度 v3 標準算法、含截距回歸 `J`（bootstrap 誤差）、換階段跳過時間敏感度。
- `report_figs_1013.py`：報告圖（需 matplotlib 與中文字型 Microsoft JhengHei）。
- 資料路徑寫在 `BASE`，請改成自己的資料夾（內含各次實驗的 `constvel_angle_*`、`constaccel_v3_*` 子資料夾）。

## Python 分析（`analysis/1015/`）

需要 numpy、scipy、matplotlib。資料路徑寫在 `data_1015.py` 的 `BASE`（內含各次 `constvel_angle_*`、`constaccel_v3_*`、`constaccel_v3trap_*` 子資料夾），結果存到 `BASE/analysis_1015/`。

- `data_1015.py`：檔案清單、讀檔、逐樣本運動特徵（方向、起步後行程、加減速階段、馬達角）。
- `constvel_1015.py`：定速 3 輪的 `B`、`Tc`，以及各速度檔力矩對溫度的回歸。
- `accel_v3_1015.py`：三角形與梯形 v3 的含截距回歸（沿用 `../1013/accel_v3_analysis.py`）；摩擦基準取同日定速對溫度回歸在該次溫度的值。
- `kalman_1015.py`：卡爾曼濾波＋RTS 平滑，狀態為角度、速度、加速度與模型外力矩 `d`；雜訊參數以新息對數概似最大化調整。
- `kalman_residual_1015.py`：殘差分解——依序加入馬達角漣波、溫度、速度非線性、起步後行程、加減速階段，以運動段重抽樣估計係數標準差；並檢查 J 是否隨運動類型改變、角度與行程能否分開、行程與時間哪個較能解釋、漣波是否隨速度改變、J 對行程項節點的敏感度；再把鑑別出的各項放回濾波器重跑。
- `report_figs_1015.py`、`kalman_figs_1015.py`：報告圖。
- 執行順序：`constvel_1015.py` → `accel_v3_1015.py` → `report_figs_1015.py`；`kalman_residual_1015.py` → `kalman_figs_1015.py`（濾波結果快取在 `analysis_1015/kf_cache/`）。

## Python 分析（`analysis/1022/`）

需要 numpy、matplotlib。資料路徑寫在 `velmod_1022.py` 的 `BASE`（內含各次 `velmod_*`、`constvel_angle_*`、`warmup_*` 子資料夾），結果存到 `BASE` 的上一層 `analysis_1022/`。執行時封鎖網路。

- `velmod_1022.py`：速度調變 5 次（DIAG 正、反轉與重跑，EXT45 正、反轉）的各條件 `J`（固定用 10/22 的 4 個條件 `CONDS`）；把電流換成 `target_moment` 以同一套分析重算（檢查分析流程）；各姿態合併、`J` 對加速度振幅（`J∞ + k/A` 與直線兩種寫法）、兩姿態差 `ΔJ` 與控制器模型比較；控制器的力矩換電流比值（`target_current` 對 `target_moment` 回歸）與 J1 靜止電流；同日定速的實際速度比例與暖機各批結果。
- `report_figs_1022.py`：報告圖（`J` 對加速度振幅；實測 `J` 與 `J_model`、`ΔJ`）。
- 執行順序：`velmod_1022.py` → `report_figs_1022.py`。

## Python 分析（`analysis/1029/`）

需要 numpy、matplotlib。資料路徑寫在 `known_mass_1029.py` 的 `BASE`（內含 10/29 各次 `velmod_*`、`j1static_*`、`constvel_angle_*` 子資料夾），結果存到 `BASE` 的上一層 `analysis_1029/`。執行時封鎖網路。

- `known_mass_1029.py`：10/29 已知質量（2F-85 橫夾 1 kg 砝碼）分析——各次各條件的 `J`；官方點質量模型（UR 官方 DH 頁面的連桿質量與重心 + 末端工具 + 砝碼 m·r²）；控制器 `J_model` 扣掉官方連桿與教導器負載後剩下的常數；夾砝碼 − 不夾的 `ΔJ ÷ m·r²` 與換算比例 K；反轉 4 次的合併擬合 `J = c·X + J_d`（並與強制含 1.915 的擬合比較）；J1 靜止電流的 K；與 10/22 比較；同日定速 `B`、`Tc` 與 10/13、10/15 對照。第一個參數為砝碼質量（預設 1.000 kg，例如 `python known_mass_1029.py 1.002`）。
- 與 10/22 比較讀 `analysis_1022/velmod_1022.json`（先執行 `../1022/velmod_1022.py`，沒有就略過）；10/13、10/15 定速對照讀 `BASE` 上一層的 `constvel_angle_*_full_2026101[35]_*`。
- `report_figs_1029.py`：報告圖（已知質量的比例；實測 `J` 對官方點質量模型；與 10/22 比較）。
- 執行順序：`known_mass_1029.py` → `report_figs_1029.py`。

## MATLAB 分析（R2025a）

- `matlab/constvel/plot_constvel_analysis.m`：定速法主程式（回歸、漣波、頻譜、角度域濾波等圖），讀取 `DATA_FOLDER` 內的 `constvel_data_pos.csv`、`constvel_data_neg.csv`。
- `fit_constvel.m`：定速法統一算法——依指令速度切平台期、去掉起始 0.5 s、取尾端整數個馬達圈、6 檔平均點回歸。
- `remove_angle_ripple.m`、`ripple_vs_motor_angle.m`、`validate_ripple_model.m`：電流中鎖定在馬達角度（關節角 × 101）的漣波分析與角度域濾波。
- `sim_realistic_constvel.m`：以實測運動與漣波建立的蒙地卡羅模擬，檢驗算法偏差。
- `matlab/analysis_1007/` 的程式會自動加入 `../constvel` 路徑；資料路徑寫在各檔開頭（`base`、`file`），請改成自己的資料夾。

CSV 欄位（RTDE 125 Hz 記錄）：`timestamp`，以及各關節 `actual_q_k`、`actual_qd_k`、`actual_current_k`、`target_qd_k`、`joint_temp_k`（k = 0～5）。速度調變實驗另記錄 `target_qdd_k`、`target_moment_k`、`target_current_k`、`joint_control_output_k`（控制器或 ur_rtde 不支援時留空）。2026/10/26 起五支上機程式都另記錄 `target_speed_fraction`（教導器速度滑桿，0～1；讀不到時留空）。

上機程式的輸出資料夾（`velmod_*`、`j1static_*`、`warmup_*` 等）存在專案根目錄的 `outputs/`（2026/10/29 起，尚未審查；之前存在執行時的目前目錄），與執行位置無關，`outputs/` 不進 git。

## 備註

- `Kt = 0.135 N·m/A` 與減速比 101 借自文獻數值；`B`、`Tc`、`J` 的絕對值會隨換算比例等比例改變。10/29 以已知質量得到 J0 速度調變的等效換算比例約 10.8 N·m/A（本研究用 0.135 × 101 = 13.64），且隨調變頻率改變（各條件 9.7～11.1），原因未明，尚未用來改寫 `B`、`Tc`。
- 控制器內部的力矩換電流比值（10/22 資料，`target_current` = `target_moment` ÷ 比值 + 摩擦前饋）：J0 為 12.44、J1 為 12.05、J2 為 12.17 N·m/A。
- 控制器 `target_moment` 對應的 J0 慣量 `J_model` = 1.915 kg·m²（常數，與 Gnad 2026 表 3 UR5e 的驅動端慣量相同）+ 官方 DH 頁面連桿點質量 + 教導器負載點質量，兩個姿態、兩種負載設定都成立。實測 `J` 經已知質量校正後與後兩項一致，1.915 在 J0 電流中看不到（見 `analysis/1029/`）。
- 減速比 101 已由資料驗證：以 101 換算馬達角時，各速度檔的電流漣波對齊最好。
