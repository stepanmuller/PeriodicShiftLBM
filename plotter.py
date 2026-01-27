import os
import numpy as np
import matplotlib.pyplot as plt

# 1. Read dimensions and data directly from RAM-disk
with open("/dev/shm/sim_data.bin", "rb") as f:
	plotNumber = np.fromfile(f, dtype=np.int32, count=1)[0]
	dims = np.fromfile(f, dtype=np.int32, count=3)
	# Read the rest as float32
	ny, nz, variables = dims
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

rho = data[:, :, 0]
ux = data[:, :, 1]
uy = data[:, :, 2]
uz = data[:, :, 3]
mask = data[:, :, 4]

uMag = (uy**2 + uz**2)**0.5

plt.imshow(uMag, origin="lower", cmap="viridis", interpolation="nearest", zorder=3)
ax = plt.gca()
for spine in ax.spines.values():
	spine.set_zorder(1)

# Save image into results directory
filename = str(plotNumber) + ".png"
plt.savefig(os.path.join(results_dir, filename),
			dpi=1000, bbox_inches="tight")
