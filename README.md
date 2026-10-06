# UR5-SYSID：UR5 CB3 J0 摩擦與慣量鑑別

以 UR5（CB3）底座關節 J0 為對象，鑑別單軸摩擦與等效慣量，作為關節健康監測（PHM）的基準。

模型：

```
tau = J·qdd + B·qd + Tc·sign(qd),   tau = i × Kt × 101
```

- **定速法**：單一方向、穩態（qdd ≈ 0），多個定速值做線性回歸 → `B`、`Tc`
- **定加速度法**：梯形速度，同一速度下「加速段 − 減速段」消去摩擦 → `J`；「加速段 + 減速段」得到摩擦曲線，交叉驗證 `B`、`Tc`

本 repo 只含程式與文件，不含真機量測資料。

## ⚠️ 安全須知

`robot/` 內的腳本**會透過 URScript（port 30002）直接驅動真實機械臂**。

- `ur5_const_velocity_angle_ident.py`、`ur5_j0_warmup.py`、`ur5_const_accel_ident_v3.py` **尚未經安全審查、尚未上機**。
- 在能連到機台的電腦上，請勿直接執行主程式；離線檢查請只用 `py_compile` 或 `tests/`（測試會封鎖網路並使用假的 RTDE）。
- 執行前務必確認 `ROBOT_IP`、J0 安全視窗（`J0_SAFE_MIN_DEG` / `J0_SAFE_MAX_DEG`，未設定時程式會中止）、手臂周圍淨空，並有人在緊急停止按鈕旁。
- 正轉、反轉分開執行，每個方向結束後人工確認。

## 目錄

```
robot/                                   上機腳本（會送指令給機械臂）
  ur5_const_velocity_angle_ident.py      固定角度窗口定速實驗：所有速度、正反轉都掃過同一段 J0 角度（10 個馬達圈）
  ur5_j0_warmup.py                       暖機：J0 來回轉動，溫度與摩擦指標穩定後自動結束
  ur5_const_accel_ident_v3.py            定加速度 v3：梯形速度＋起始角錯開 1/20 馬達圈，並含離線分析函式
  ur5_home_pose.py                       診斷原點姿態檢查（預設只檢查、不移動）
tests/                                   離線測試（假 RTDE、封鎖網路）
docs/
  EXPERIMENT_PLAN_ANGLE.md               固定角度窗口定速實驗規劃
  EXPERIMENT_PLAN_ACCEL_V3.md            定加速度 v3 規劃與合成資料驗證
matlab/
  constvel/                              定速法分析與共用函式
  analysis_1007/                         定速重複實驗與定加速度（v2）分析
```

## 離線測試

```bash
pip install numpy
python -m unittest discover -s tests -v
```

## MATLAB 分析（R2025a）

- `matlab/constvel/plot_constvel_analysis.m`：定速法主程式（回歸、漣波、頻譜、角度域濾波等圖），讀取 `DATA_FOLDER` 內的 `constvel_data_pos.csv`、`constvel_data_neg.csv`。
- `fit_constvel.m`：定速法統一算法——依指令速度切平台期、去掉起始 0.5 s、取尾端整數個馬達圈、6 檔平均點回歸。
- `remove_angle_ripple.m`、`ripple_vs_motor_angle.m`、`validate_ripple_model.m`：電流中鎖定在馬達角度（關節角 × 101）的漣波分析與角度域濾波。
- `sim_realistic_constvel.m`：以實測運動與漣波建立的蒙地卡羅模擬，檢驗算法偏差。
- `matlab/analysis_1007/` 的程式會自動加入 `../constvel` 路徑；資料路徑寫在各檔開頭（`base`、`file`），請改成自己的資料夾。

CSV 欄位（RTDE 125 Hz 記錄）：`timestamp`，以及各關節 `actual_q_k`、`actual_qd_k`、`actual_current_k`、`target_qd_k`、`joint_temp_k`（k = 0～5）。

## 備註

- `Kt = 0.135 N·m/A` 與減速比 101 借自文獻數值，未獨立鑑別；`B`、`Tc`、`J` 的絕對值會隨 `Kt` 等比例改變。
- 減速比 101 已由資料驗證：以 101 換算馬達角時，各速度檔的電流漣波對齊最好。
