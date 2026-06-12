import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 定義要讀取的檔案與標題
cases = [
    {"file": "test/convection_Ra500.csv", "title": "No Convection (Ra ~ 500)"},
    {"file": "test/convection_Ra2000.csv", "title": "Moderate Convection (Ra ~ 2,000)"},
    {"file": "test/convection_Ra10000.csv", "title": "Strong Convection (Ra ~ 10,000)"}
]

# 建立 1x3 的子圖
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, case in zip(axes, cases):
    # 讀取 CSV
    df = pd.read_csv(case["file"])
    
    # 轉回二維網格
    X = df.pivot(index='y', columns='x', values='x').values
    Y = df.pivot(index='y', columns='x', values='y').values
    U = df.pivot(index='y', columns='x', values='ux').values
    V = df.pivot(index='y', columns='x', values='uy').values
    T = df.pivot(index='y', columns='x', values='T').values

    # 畫等溫線與流線
    contour = ax.contourf(X, Y, T, levels=50, cmap='inferno')
    ax.streamplot(X, Y, U, V, color='white', density=1.2, linewidth=0.8, arrowsize=1)
    
    # 設定標題與標籤
    ax.set_title(case["title"])
    ax.set_xlabel('x')
    if ax == axes[0]:
        ax.set_ylabel('y')

# 共用一個 Colorbar
#fig.colorbar(contour, ax=axes.ravel().tolist(), label='Temperature')
plt.tight_layout()

output_path = 'test/convection_compare.png'
plt.savefig(output_path, dpi=200)
print(f"Plot saved to {output_path}")