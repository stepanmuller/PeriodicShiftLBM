import os
import numpy as np
import matplotlib.pyplot as plt

# 1. Read binary data
file_path = "/dev/shm/sim_data.bin"
with open(file_path, "rb") as f:
    plotNumber = np.fromfile(f, dtype=np.int32, count=1)[0]
    ny, nz, _ = np.fromfile(f, dtype=np.int32, count=3)
    data = np.fromfile(f, dtype=np.float32).reshape((ny, nz, -1))

rho = data[:, :, 0]
u1 = data[:, :, 1]
u2 = data[:, :, 2]
u3 = data[:, :, 3]

# 2. Calculate uMag
uMag = np.sqrt(u1**2 + u2**2)

# 3. Minimalist Plot
fig, ax = plt.subplots()
ax.imshow(uMag, origin="lower", cmap="viridis")
ax.axis("off") # Removes axes, ticks, and labels

# 4. Save results (matching your original path logic)
os.makedirs("results", exist_ok=True)
save_path = os.path.join("results", f"{plotNumber}.png")
plt.savefig(save_path, dpi=500, bbox_inches="tight", pad_inches=0)
plt.close(fig)
