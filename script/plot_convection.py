import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 讀取剛剛生成的最後一步 CSV
df = pd.read_csv("test/convection_field.csv")

# 將一維資料轉回二維網格形狀
X = df.pivot(index='y', columns='x', values='x').values
Y = df.pivot(index='y', columns='x', values='y').values
U = df.pivot(index='y', columns='x', values='ux').values
V = df.pivot(index='y', columns='x', values='uy').values
T = df.pivot(index='y', columns='x', values='T').values

plt.figure(figsize=(10, 8))
# 畫出溫度等高線 (暖色系)
contour = plt.contourf(X, Y, T, levels=50, cmap='inferno')
plt.colorbar(contour, label='Temperature')

# 畫出速度流線圖 (白色)
plt.streamplot(X, Y, U, V, color='white', density=1.5, linewidth=0.8, arrowsize=1)

plt.title('Rayleigh-Bénard Convection')
plt.xlabel('x')
plt.ylabel('y')
plt.tight_layout()

output_path = 'test/convection_result.png'
plt.savefig(output_path, dpi=200)
print(f"Plot saved to {output_path}")