"""11/4 報告的計算流程圖：每個數字從哪些資料、經過哪些步驟算出來。輸出 OUT/figs/w1104_flow.png。
（Microsoft JhengHei 沒有「−」「≈」，圖中改用文字。）"""
import os, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

sys.stdout.reconfigure(encoding="utf-8")
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False})
OUT = os.path.join(r"D:\cluade_code\ur5", "analysis_1104", "figs"); os.makedirs(OUT, exist_ok=True)

INK, INK2 = "#0b0b0b", "#52514e"
C_MEAS, C_CTRL, C_REF = "#dbe8f8", "#f3e3cf", "#e4e2da"      # 量測、控制器、不靠量測的資料
C_STEP, C_RES = "#ffffff", "#2a78d6"

fig = plt.figure(figsize=(16, 11.5))
ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 100); ax.set_ylim(0, 100); ax.axis("off")
B = {}


def box(key, x, y, w, h, title, body, fc, uses=None):
    res = fc == C_RES
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h, boxstyle="round,pad=0.4,rounding_size=1.2",
                                fc=fc, ec=C_RES if res else "#8d8b83", lw=1.2, zorder=2))
    top = y + h / 2 - 1.5
    ax.text(x, top, title, ha="center", va="top", fontsize=13, fontweight="bold", color="white" if res else INK)
    if uses:
        top -= 3.3
        ax.text(x, top, "用到：" + uses, ha="center", va="top", fontsize=10, color="#d6e6fa")
    ax.text(x, top - 3.4, body, ha="center", va="top", fontsize=11, linespacing=1.5, color="white" if res else INK2)
    B[key] = (x, y, w, h)


def down(a, b):
    xa, ya, _, ha = B[a]
    xb, yb, _, hb = B[b]
    ax.annotate("", (xb, yb + hb / 2 + 0.6), (xa, ya - ha / 2 - 0.6),
                arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.4, shrinkA=0, shrinkB=0))


ax.text(1.5, 98.5, "11/4 報告的數字怎麼算出來", fontsize=18, fontweight="bold", color=INK, va="top")
ax.text(1.5, 95.0, "藍底：上機量到的　橘底：控制器自己算的　灰底：不靠量測的資料　白底：計算步驟　深藍：報告裡的結果（「用到」寫出它從哪幾步來）",
        fontsize=10.5, color=INK2, va="top")

# ---- 資料 ----
box("meas", 17, 85, 30, 10, "上機量到的（RTDE 125 Hz）", "J0 實際角度 actual_q\nJ0 電流 actual_current", C_MEAS)
box("ref", 50, 85, 30, 10, "不靠量測的資料", "UR 官網：每節連桿的質量、重心\n教導器負載 1.30 kg；砝碼 1 kg 與位置", C_REF)
box("ctrl", 83, 85, 30, 10, "控制器自己算的", "target_qdd（要的加速度）、target_moment\ntarget_current、joint_control_output", C_CTRL)

# ---- 步驟 ----
box("s1", 17, 69.5, 30, 10.5, "① 換成加速度、力矩",
    "角度擬合（多項式 + 諧波）→ 加速度 a、速度 v\n電流 × 13.64（0.135 × 101）→ 力矩 τ", C_STEP)
box("s3", 50, 69.5, 30, 10.5, "③ 官方模型 X",
    "正運動學算每節重心離 J0 軸的水平距離 r\nX = Σ m·r²（連桿 + 負載 + 砝碼）", C_STEP)
box("s4", 83, 69.5, 30, 10.5, "④ 控制器的 J_model",
    "target_moment 對 target_qdd 回歸\n斜率 = 控制器心中的 J", C_STEP)
box("s2", 17, 54, 30, 10.5, "② 每個條件回歸出 J",
    "τ = J·a + B·v + 摩擦、漣波、慢變化\n每次實驗 × 6 個條件，各得一個 J", C_STEP)
for a, b in (("meas", "s1"), ("ref", "s3"), ("ctrl", "s4"), ("s1", "s2")):
    down(a, b)
ax.text(50, 54, "下面每個結果都由上面的步驟組合而成\n（看「用到」那一行）", ha="center", va="center", fontsize=11.5, color=INK2,
        linespacing=1.5)
ax.annotate("", (50, 45.5), (50, 50.2), arrowprops=dict(arrowstyle="-|>", color=INK2, lw=1.4))

# ---- 結果 ----
W, H = 30, 15.5
box("r0", 17, 35, W, H, "第 1 節：控制器的組成", "J_model 減去 X，剩下 1.915\n兩姿態、6 次都一樣 → 固定常數\n（設定檔的轉子慣量 × 101²）",
    C_RES, uses="③、④")
box("r1", 50, 35, W, H, "圖 1（3.1）：換算比例", "夾砝碼的 J 減不夾的 J = ΔJ\nΔJ ÷ 砝碼 m·r² = 1.26\nK = 13.64 ÷ 1.26 = 10.8 N·m/A",
    C_RES, uses="② ＋ 砝碼的 m·r²")
box("r2", 83, 35, W, H, "圖 2（3.2）：校正後和官方比", "實測 J ÷ 1.26，和 X 並排\n4 個狀態都一致\n控制器（X ＋ 1.915）每次多 1.915",
    C_RES, uses="②、③、1.26")
box("r3", 17, 13, W, H, "3.2、圖 3（3.3）：截距與 c", "4 次實驗擬合 J = c·X + J_d\n截距約 0 → 電流裡沒有 1.915\nc 隨頻率 1.23 → 1.40",
    C_RES, uses="②、③")
box("r4", 50, 13, W, H, "3.5：1.915 在哪裡不見", "電流換成控制器的訊號，重做②與擬合\ntarget_moment、target_current：約 1.92\n電流命令、實測電流：約 0",
    C_RES, uses="② 的方法、③、控制器訊號")
box("r5", 83, 13, W, H, "第 1 節：完整 6 軸模型", "③ 的資料 ＋ 1.915，用牛頓–尤拉法\n算 6 軸的力矩，和 target_moment 比\nJ0 差 0.000 N·m，6 軸都重現",
    C_RES, uses="③、控制器訊號")

fig.savefig(os.path.join(OUT, "w1104_flow.png"), dpi=150)
print(os.path.join(OUT, "w1104_flow.png"))
