import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
})

# 1. Read binary data
file_path = "/dev/shm/sim_data.bin"
with open(file_path, "rb") as f:
    plotNumber = np.fromfile(f, dtype=np.int32, count=1)[0]
    ny, nz, _ = np.fromfile(f, dtype=np.int32, count=3)
    data = np.fromfile(f, dtype=np.float32).reshape((ny, nz, -1))

rho = data[:, :, 0]
ux = data[:, :, 1]
uy = data[:, :, 2]
uz = data[:, :, 3]
mask = data[:, :, 4]

# 2. Calculate uMag
uMag = np.sqrt(ux**2 + uy**2 + uz**2)

# Mask points where mask == 1
uMag_masked = np.ma.masked_where(mask == 1, uMag)

# 3. Plotting
fig, ax = plt.subplots(figsize=(8, 6))

#cmap = plt.cm.jet.copy()
#cmap.set_bad(color='grey')
cmap = plt.cm.viridis.copy()
cmap.set_bad(color='black')

im = ax.imshow(uMag_masked, origin="lower", cmap=cmap)
ax.axis("off")


# 4. Add Inset Colorbar
# loc=8 is "lower center". width/height are percentages of the parent axes.
# borderpad provides a tiny bit of breathing room from the bottom edge.
axins = inset_axes(ax,
                   width="30%",  
                   height="4%",  
                   loc='upper left',
                   borderpad=1)

# Create the colorbar and style it
cbar = fig.colorbar(im, cax=axins, orientation="horizontal")
cbar.ax.tick_params(labelsize=8, colors='black')
cbar.set_label('LBM velocity magnitude', color='black', fontsize=8)

# 5. Save results
os.makedirs("results", exist_ok=True)
save_path = os.path.join("results", f"{plotNumber}.png")

# Use bbox_inches="tight" and pad_inches=0 to ensure the image fills the entire file
plt.savefig(save_path, dpi=500, bbox_inches="tight", pad_inches=0)
plt.close(fig)
