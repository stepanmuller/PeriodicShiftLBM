"""
import numpy as np
import matplotlib.pyplot as plt

# --- Load CSV ---
data = np.loadtxt("result.csv", delimiter=",")

# --- Create figure ---
plt.figure(figsize=(6, 4))
img = plt.imshow(data, origin="lower", cmap="viridis", aspect="equal")

# --- Add colorbar ---
#cbar = plt.colorbar(img)
#cbar.set_label("Velocity magnitude")

# --- Save image ---
plt.savefig("velocity_slice.png", dpi=500, bbox_inches="tight")

# --- Optionally show ---
#plt.show()
"""
import numpy as np
import matplotlib.pyplot as plt

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

# --- Load CSV ---
data = np.loadtxt("result3.csv", delimiter=",")

# --- Create figure ---
plt.figure(figsize=(6, 4))
img = plt.imshow(data, origin="lower", cmap="viridis", aspect="equal")

# --- Add colorbar ---
# cbar = plt.colorbar(img)
# cbar.set_label(r"Velocity magnitude")

# --- Save image ---
plt.savefig("velocity_slice.png", dpi=1000, bbox_inches="tight")

# --- Optionally show ---
# plt.show()

