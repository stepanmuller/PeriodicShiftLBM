import os
import numpy as np
import matplotlib.pyplot as plt

# 1. Read dimensions and data directly from RAM-disk
with open("/dev/shm/sim_data.bin", "rb") as f:
    dims = np.fromfile(f, dtype=np.int32, count=2)
    # Read the rest as float32
    data = np.fromfile(f, dtype=np.float32).reshape(dims)

# Use LaTeX-style fonts and smaller text
plt.rcParams.update({
    "text.usetex": True,                  # enable LaTeX
    "font.family": "serif",               # serif font
    "font.size": 8,                       # global font size
    "axes.labelsize": 8,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 7,
    "figure.titlesize": 9
})

# Create results directory if it doesn't exist
results_dir = "results"
os.makedirs(results_dir, exist_ok=True)

plt.figure(figsize=(6, 4))
img = plt.imshow(data, origin="lower", cmap="viridis", aspect="equal")

# Save image into results directory
plt.savefig(os.path.join(results_dir, "result.png"),
            dpi=1000, bbox_inches="tight")

