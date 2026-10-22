# UR5-SYSID：UR5 CB3 J0 摩擦與慣量鑑別

以 UR5（CB3）底座關節 J0 為對象，鑑別單軸摩擦與等效慣量，作為關節健康監測（PHM）的基準。

模型：

```
tau = J·qdd + B·qd + Tc·sign(qd),   tau = i × Kt × 101
```

- **定速法**：單一方向、穩態（qdd ≈ 0），多個定速值做線性回歸 → `B`、`Tc`
- **定加速度法**：同一速度下「加速段 − 減速段」消去摩擦 → `J`。10/13 實測減速段摩擦比加速段多一個與加速度無關的常數，故以多個加速度檔做含截距回歸（斜率 = `J`）
- **速度調變法**：巡航速度上疊加正弦速度（速度不歸零、不反向），力矩中與加速度同相的成分 → `J`，與速度同相的成分 → 局部摩擦斜率；以 4 種頻率、振幅檢查 `J` 是否隨頻率或振幅改變。可改在手臂外伸姿態（`POSE = "EXT45"`）執行，比較兩個姿態的 `J`；並記錄控制器的 `target_moment`，得到控制器模型的慣量 `J_model` 作對照

本 repo 只含程式與文件，不含真機量測資料。

## ⚠️ 安全須知

`robot/` 內的腳本**會透過 URScript（port 30002）直接驅動真實機械臂**。

- `ur5_const_velocity_angle_ident.py`、`ur5_j0_warmup.py`、`ur5_const_accel_ident_v3.py` 已於 2026/10/13 上機；之後的修改（見 `docs/` 各規劃第 7 節）**尚未經安全審查**。
- `ur5_velocity_modulation_ident.py`（2026/10/19 新增）已通過安全審查（審查版本 commit `72591a7`），2026/10/22 上機（`DIAG`、`EXT45` 各 FULL 正、反轉）；新的條件或姿態仍請依 `TEST_STAGE` = QUICK → MID → FULL 執行。
  - 姿態 `POSE` 預設 `DIAG`（診斷原點）。`EXT45`（J1、J2 = −45°、45°，法蘭離 J0 軸約 0.80 m）時程式**不移動 J1~J5**：須先以示教器擺好、確認 J0 掃過的扇形範圍淨空並輸入姿態名稱；連上 RTDE 後以實際角度檢查 J1~J5（±0.5°），不符或姿態名稱未知即中止。結束後 J1~J5 仍在外伸姿態，請以示教器移回。
- 在能連到機台的電腦上，請勿直接執行主程式；離線檢查請只用 `py_compile` 或 `tests/`（測試會封鎖網路並使用假的 RTDE）。
- 執行前務必確認 `ROBOT_IP`、J0 安全視窗（`J0_SAFE_MIN_DEG` / `J0_SAFE_MAX_DEG`，未設定時程式會中止）、手臂周圍淨空，並有人在緊急停止按鈕旁。
- 正轉、反轉分開執行，每個方向結束後人工確認。

## 目錄

```
robot/                                   上機腳本（會送指令給機械臂）
  ur5_const_velocity_angle_ident.py      固定角度窗口定速實驗：所有速度、正反轉都掃過同一段 J0 角度（10 個馬達圈）
  ur5_j0_warmup.py                       暖機：J0 來回轉動，溫度與摩擦指標穩定後自動結束
  ur5_const_accel_ident_v3.py            定加速度 v3：三角形速度（定速段 0）＋起始角錯開 1/20 馬達圈，並含離線分析函式
  ur5_const_accel_ident_v3_trapezoid.py  定加速度 v3 梯形對照版：只改定速段為 4 個馬達圈與輸出前綴 constaccel_v3trap_
  ur5_velocity_modulation_ident.py       速度調變實驗：巡航 0.12 rad/s 上疊加正弦（4 種頻率、振幅 × 7 次，起點錯開 1/7 馬達圈），並含離線分析函式
  ur5_home_pose.py                       診斷原點姿態檢查（預設只檢查、不移動）
tests/                                   離線測試（假 RTDE、封鎖網路）
docs/
  EXPERIMENT_PLAN_ANGLE.md               固定角度窗口定速實驗規劃
  EXPERIMENT_PLAN_ACCEL_V3.md            定加速度 v3 規劃與合成資料驗證
analysis/
  1013/                                  10/13 實驗分析（Python）：定速 3 輪比較、定加速度 v3、報告圖
  1015/                                  10/15 實驗分析（Python）：定速溫度效應、三角形 vs 梯形 v3、卡爾曼濾波殘差分解
matlab/
  constvel/                              定速法分析與共用函式
  analysis_1007/                         定速重複實驗與定加速度（v2）分析
```

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

- `velmod_1022.py`：速度調變 5 次（DIAG 正、反轉與重跑，EXT45 正、反轉）的各條件 `J`；把電流換成 `target_moment` 以同一套分析重算（檢查分析流程）；各姿態合併、`J` 對加速度振幅（`J∞ + k/A` 與直線兩種寫法）、兩姿態差 `ΔJ` 與控制器模型比較；控制器的力矩換電流比值（`target_current` 對 `target_moment` 回歸）與 J1 靜止電流；同日定速的實際速度比例與暖機各批結果。
- `report_figs_1022.py`：報告圖（`J` 對加速度振幅；實測 `J` 與 `J_model`、`ΔJ`）。
- 執行順序：`velmod_1022.py` → `report_figs_1022.py`。

## MATLAB 分析（R2025a）

- `matlab/constvel/plot_constvel_analysis.m`：定速法主程式（回歸、漣波、頻譜、角度域濾波等圖），讀取 `DATA_FOLDER` 內的 `constvel_data_pos.csv`、`constvel_data_neg.csv`。
- `fit_constvel.m`：定速法統一算法——依指令速度切平台期、去掉起始 0.5 s、取尾端整數個馬達圈、6 檔平均點回歸。
- `remove_angle_ripple.m`、`ripple_vs_motor_angle.m`、`validate_ripple_model.m`：電流中鎖定在馬達角度（關節角 × 101）的漣波分析與角度域濾波。
- `sim_realistic_constvel.m`：以實測運動與漣波建立的蒙地卡羅模擬，檢驗算法偏差。
- `matlab/analysis_1007/` 的程式會自動加入 `../constvel` 路徑；資料路徑寫在各檔開頭（`base`、`file`），請改成自己的資料夾。

CSV 欄位（RTDE 125 Hz 記錄）：`timestamp`，以及各關節 `actual_q_k`、`actual_qd_k`、`actual_current_k`、`target_qd_k`、`joint_temp_k`（k = 0～5）。速度調變實驗另記錄 `target_qdd_k`、`target_moment_k`、`target_current_k`、`joint_control_output_k`（控制器或 ur_rtde 不支援時留空）。

## 備註

- `Kt = 0.135 N·m/A` 與減速比 101 借自文獻數值，未獨立鑑別；`B`、`Tc`、`J` 的絕對值會隨 `Kt` 等比例改變。
- 控制器內部的力矩換電流比值（10/22 資料，`target_current` = `target_moment` ÷ 比值 + 摩擦前饋）：J0 為 12.44、J1 為 12.05、J2 為 12.17 N·m/A；本研究用 0.135 × 101 = 13.64。哪個接近實際值尚未驗證。
- 減速比 101 已由資料驗證：以 101 換算馬達角時，各速度檔的電流漣波對齊最好。
