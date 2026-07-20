import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

W = np.load("W_matrix.npy")
H = np.load("H_matrix.npy")

fig, axs = plt.subplots(1, 2, figsize=(14, 5))
sns.heatmap(W, annot=False, fmt=".2f", cmap = "coolwarm", cbar = True, ax=axs[0])
axs[0].set_title("W")
sns.heatmap(H, annot=False, fmt=".2f", cmap = "coolwarm", cbar = True, ax=axs[1])
axs[1].set_title("H")
plt.tight_layout()
plt.show()




