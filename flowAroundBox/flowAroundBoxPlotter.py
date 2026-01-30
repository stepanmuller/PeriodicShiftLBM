import os
import numpy as np
import matplotlib.pyplot as plt

# 1. Read binary data
file_path = "/dev/shm/sim_data.bin"
with open(file_path, "rb") as f:
    plotNumber = np.fromfile(f, dtype=np.int32, count=1)[0]
    ny, nz, _ = np.fromfile(f, dtype=np.int32, count=3)
    data = np.fromfile(f, dtype=np.float32).reshape((ny, nz, -1))

# 2. Extract and calculate uMag
uMag = np.sqrt(data[:, :, 2]**2 + data[:, :, 3]**2)

# 3. Minimalist Plot
fig, ax = plt.subplots()
ax.imshow(uMag, origin="lower", cmap="viridis")
ax.axis("off") # Removes axes, ticks, and labels

# 4. Save results (matching your original path logic)
os.makedirs("results", exist_ok=True)
save_path = os.path.join("results", f"{plotNumber}.png")
plt.savefig(save_path, dpi=500, bbox_inches="tight", pad_inches=0)
plt.close(fig)
