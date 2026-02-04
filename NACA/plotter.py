import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy.ma as ma

# --- Configuration & Styling ---
plt.rcParams.update({
    "text.usetex": True, # Set to True if you have working LaTeX on your system
    "font.family": "serif",
    "font.size": 10
})

file_path = "/dev/shm/sim_data.bin"

# 1. Read binary data
with open(file_path, "rb") as f:
	plotNumber = np.fromfile(f, dtype=np.int32, count=1)[0]
	dims = np.fromfile(f, dtype=np.int32, count=3)
	ny, nx, n_vars = dims
	data = np.fromfile(f, dtype=np.float32).reshape((ny, nx, n_vars))

# 2. Extract variables
# Data is saved as {p, ux, uy, uz, mask}
p  = data[:, :, 0]
ux   = data[:, :, 1]
uy   = data[:, :, 2]
uz   = data[:, :, 3]
mask = data[:, :, 4]  # 1.0 = solid, 0.0 = fluid

uMag = np.sqrt(ux**2 + uy**2)

# 3. Setup Figure (2 plots side-by-side)
fig = plt.figure(figsize=(10, 4))
gs = fig.add_gridspec(1, 2, width_ratios=[1, 1], wspace=0.1, left=0.05, right=0.95, top=0.95, bottom=0.15)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])

# Create masks for plotting (Mask where mask == 1)
# We use a boolean array where True means "hide this"
is_solid = mask > 0.5 
uMag_masked = ma.array(uMag, mask=is_solid)
p_masked  = ma.array(p, mask=is_solid)

# --- Velocity plot ---
imgv = ax1.imshow(uMag_masked, origin="lower", cmap="viridis", aspect="equal")
imgv.cmap.set_bad(color="black")
divider1 = make_axes_locatable(ax1)
caxv = divider1.append_axes("bottom", size="5%", pad=0.4)  # moved below
cbarv = fig.colorbar(imgv, cax=caxv, orientation="horizontal")
cbarv.set_label("In-surface velocity magnitude [m/s]", labelpad=5)
imgv.set_clim(0.0, 120)

# --- streamlines ---
s_vals = np.arange(ny) # Vertical axis (y)
x_vals = np.arange(nx) # Horizontal axis (x)
ux_masked_stream = np.where(mask > 0.5, np.nan, ux)
uy_masked_stream = np.where(mask > 0.5, np.nan, uy)

streamlines = ax1.streamplot(
    x_vals, s_vals, ux_masked_stream, uy_masked_stream,
    density=1.5, color="white", linewidth=0.3, arrowsize=0.4
)
ax1.set_xlim(x_vals[0], x_vals[-1])
ax1.set_ylim(s_vals[0], s_vals[-1])
ax1.set_aspect("equal")


# --- Pressure plot ---
imgp = ax2.imshow(p_masked, origin="lower", cmap="viridis", aspect="equal")
imgp.cmap.set_bad(color="black")
divider2 = make_axes_locatable(ax2)
caxp = divider2.append_axes("bottom", size="5%", pad=0.4)  # moved below
cbarp = fig.colorbar(imgp, cax=caxp, orientation="horizontal")
cbarp.set_label("Static pressure [Pa]", labelpad=5)
imgp.set_clim(np.min(p), np.max(p))

# 4. Save results
os.makedirs("results", exist_ok=True)
save_path = os.path.join("results", f"{plotNumber}.png")
plt.savefig(save_path, dpi=1000, bbox_inches="tight")
print(f"Saved plot to {save_path}")
plt.close(fig)
