import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 定義要讀取的檔案與標題 (固定 Ra = 10000, 改變 Pr)
cases = [
    {"file": "test/convection_Ra10000_Pr0_16.csv", "title": "Low Prandtl (Pr = 0.16)"},
    {"file": "test/convection_Ra10000_Pr1.csv",    "title": "Base Case (Pr = 1.0)"},
    {"file": "test/convection_Ra10000_Pr6_25.csv", "title": "High Prandtl (Pr = 6.25)"}
]

# 建立 1x3 的子圖
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for ax, case in zip(axes, cases):
    try:
        # 讀取 CSV
        df = pd.read_csv(case["file"])
        
        # --- 新增：計算每個格點的速度大小，並抓出最大值 ---
        v_mag = np.sqrt(df['ux']**2 + df['uy']**2)
        v_max = v_mag.max()
        
        # 轉回二維網格
        X = df.pivot(index='y', columns='x', values='x').values
        Y = df.pivot(index='y', columns='x', values='y').values
        U = df.pivot(index='y', columns='x', values='ux').values
        V = df.pivot(index='y', columns='x', values='uy').values
        T = df.pivot(index='y', columns='x', values='T').values

        # 畫等溫線與流線
        contour = ax.contourf(X, Y, T, levels=50, cmap='inferno', extend='both')
        ax.streamplot(X, Y, U, V, color='white', density=1.2, linewidth=0.8, arrowsize=1)
        
        # --- 修改標題：把 v_max 顯示出來 ---
        ax.set_title(f"{case['title']}\n$V_{{max}}$ = {v_max:.6f}", fontsize=14)
        ax.set_xlabel('x')
        if ax == axes[0]:
            ax.set_ylabel('y')
            
    except FileNotFoundError:
        ax.set_title(f"File not found:\n{case['file']}")
        ax.axis('off')
        print(f"警告：找不到檔案 {case['file']}")

plt.tight_layout()

# 輸出圖檔
output_path = 'test/convection_Pr_compare_with_Vmax.png'
plt.savefig(output_path, dpi=200)
print(f"Plot saved to {output_path}")