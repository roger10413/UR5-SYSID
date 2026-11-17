# 敲擊位置示意圖（UP 姿態；尺寸取自官方 DH，工具與紙杯為示意）
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle, Circle, FancyArrowPatch, Polygon
plt.rcParams.update({"font.family": "Microsoft JhengHei", "axes.unicode_minus": False, "font.size": 12})
INK, GRAY, ARM, RED, BLUE, TOOL = "#222222", "#9a9a9a", "#d9dde3", "#d62728", "#1f5fbf", "#f2c879"
fig, axs = plt.subplots(1, 2, figsize=(12, 8.2), gridspec_kw=dict(width_ratios=[1.25, 1]))

def platform(ax, x0, x1):
    ax.add_patch(Rectangle((x0, -0.07), x1 - x0, 0.05, fc="#c8b49a", ec=INK, lw=1))
    for xl in (x0 + 0.03, x1 - 0.06):
        ax.add_patch(Rectangle((xl, -0.28), 0.03, 0.21, fc="#c8b49a", ec=INK, lw=1))
    ax.text((x0 + x1) / 2, -0.17, "平台（會跟著手臂晃）", ha="center", va="center", fontsize=11, color="#7a5a2a")
    ax.annotate("", xy=(x1 + 0.02, -0.045), xytext=(x1 + 0.12, -0.045),
                arrowprops=dict(arrowstyle="<->", color="#7a5a2a", lw=1.6))

def arm_column(ax, y=0.0):
    ax.add_patch(Rectangle((y - 0.075, -0.02), 0.15, 0.10, fc=ARM, ec=INK, lw=1.2))                  # 底座
    ax.add_patch(Rectangle((y - 0.045, 0.09), 0.09, 0.82, fc=ARM, ec=INK, lw=1.2))                   # 上臂＋前臂
    for z, nm in ((0.089, "J1 肩"), (0.514, "J2 肘"), (0.906, "J3 腕 1")):
        ax.add_patch(Circle((y, z), 0.055, fc="white", ec=INK, lw=1.4))
        ax.text(y - 0.075, z, nm, ha="right", va="center", fontsize=10.5, color=INK)

# ---------- 左：側視（J1 轉軸指向你） ----------
ax = axs[0]
platform(ax, -0.35, 0.55)
arm_column(ax)
ax.add_patch(FancyBboxPatch((-0.05, 0.93), 0.13, 0.13, boxstyle="round,pad=0.005", fc="#b9c2cc", ec=INK, lw=1.4))   # 腕 2、腕 3 外殼
ax.text(0.015, 1.10, "手腕（腕 2、腕 3）\n金屬外殼", ha="center", va="bottom", fontsize=10.5)
ax.add_patch(Rectangle((0.082, 0.965), 0.04, 0.072, fc="#8a8f96", ec=INK, lw=1))                 # FT300
ax.add_patch(Rectangle((0.122, 0.96), 0.10, 0.082, fc=TOOL, ec=INK, lw=1))                       # 夾爪本體
ax.add_patch(Polygon([(0.222, 1.035), (0.30, 1.035), (0.30, 1.025), (0.222, 1.02)], fc=TOOL, ec=INK))   # 手指
ax.add_patch(Polygon([(0.222, 0.967), (0.30, 0.967), (0.30, 0.977), (0.222, 0.982)], fc=TOOL, ec=INK))
ax.add_patch(Rectangle((0.235, 0.977), 0.055, 0.048, fc="white", ec=INK, lw=1, ls="--"))         # 紙杯
ax.text(0.30, 0.90, "FT300 → 夾爪 → 紙杯\n（都不要敲）", ha="center", va="top", fontsize=10.5, color=GRAY)
for xx in (0.10, 0.17, 0.262):
    ax.text(xx, 1.0, "×", ha="center", va="center", fontsize=13, color=GRAY, fontweight="bold")
# fb 敲擊
ax.add_patch(FancyArrowPatch((-0.42, 1.0), (-0.06, 1.0), arrowstyle="-|>,head_width=7,head_length=12", color=RED, lw=3))
ax.text(-0.46, 1.03, "前後敲（fb）\n敲手腕外殼的背面，\n往工具伸出的方向", ha="left", va="bottom", fontsize=12, color=RED, fontweight="bold")
# J1 轉動方向
ax.add_patch(FancyArrowPatch((0.16, 0.62), (0.16, 0.82), connectionstyle="arc3,rad=0.25",
                             arrowstyle="<|-|>,head_width=4,head_length=8", color=BLUE, lw=1.6))
ax.text(0.21, 0.72, "J1 轉動時，手臂上端\n就是這樣前後擺", ha="left", va="center", fontsize=10.5, color=BLUE)
ax.set_title("側視圖（站在手臂側面，J1 轉軸指向你）", loc="left", fontsize=13)
ax.set_xlim(-0.48, 0.62); ax.set_ylim(-0.30, 1.30); ax.set_aspect("equal"); ax.axis("off")

# ---------- 右：從工具前方看 ----------
ax = axs[1]
platform(ax, -0.35, 0.45)
arm_column(ax)
ax.add_patch(Rectangle((-0.045, 0.88), 0.20, 0.06, fc=ARM, ec=INK, lw=1.2))                     # 往腕 2 的橫向段
ax.add_patch(FancyBboxPatch((0.06, 0.93), 0.10, 0.12, boxstyle="round,pad=0.005", fc="#b9c2cc", ec=INK, lw=1.4))
ax.add_patch(Circle((0.109, 1.0), 0.040, fc=TOOL, ec=INK, lw=1.2))                              # 工具朝向你
ax.text(0.109, 1.0, "×", ha="center", va="center", fontsize=12, color=GRAY, fontweight="bold")
ax.text(0.109, 1.10, "手腕外殼\n（工具朝向你）", ha="center", va="bottom", fontsize=10.5)
ax.add_patch(FancyArrowPatch((0.42, 0.99), (0.17, 0.99), arrowstyle="-|>,head_width=7,head_length=12", color=RED, lw=3))
ax.text(0.20, 0.90, "左右敲（lr）\n從側面敲手腕外殼", ha="left", va="top", fontsize=12, color=RED, fontweight="bold")
ax.set_title("從工具前方看（工具朝向你）", loc="left", fontsize=13)
ax.set_xlim(-0.40, 0.62); ax.set_ylim(-0.30, 1.30); ax.set_aspect("equal"); ax.axis("off")

fig.suptitle("敲擊位置（手臂朝上 UP）：橡膠槌輕敲，手臂看得出輕微晃動就夠；敲完手馬上離開", x=0.01, ha="left", fontsize=14.5, fontweight="bold")
fig.text(0.01, 0.01, "輕敲那一次：同一個位置、同一個方向，改用手掌隔一條折起來的毛巾輕拍。"
         "尺寸取自官方 DH，工具、紙杯為示意。", ha="left", va="bottom", fontsize=10.5, color="#555555")
fig.tight_layout(rect=(0, 0.03, 1, 0.95))
fig.savefig(r"D:\cluade_code\ur5\analysis_1112\w1112_tap_where.png", dpi=150)
print("ok")
