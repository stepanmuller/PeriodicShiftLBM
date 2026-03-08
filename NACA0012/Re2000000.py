import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

ladson = np.array([
    [-4.25, -0.4300],
    [-2.10, -0.2150],
    [ 0.00,  0.0000],
    [ 1.85,  0.1940],
    [ 4.25,  0.4450],
    [ 6.05,  0.6250],
    [ 8.15,  0.8550],
    [10.15,  1.0450],
    [11.15,  1.1350],
    [12.10,  1.2180],
    [13.08,  1.2900],
    [14.25,  1.3620],
    [15.25,  1.4080],
    [16.25,  0.7530]
])
    
# ==========================================
# 1. DATA DEFINITION
# ==========================================
filename = "resultsAll"
supTitle = r"\textbf{NACA 0012 Lift coefficient vs Angle of attack at $Re = 6 \times 10^6$}"
xLabel = r"$\alpha$ [deg]"
yLabel = r"$C_L$ [1]"

xMin = -5
xMax = 20
yMin = -0.5
yMax = 2

# Grouping for easy iteration
xList = []
yList = []
yLegendList = []
yColorList = []
yStyleList = []

# ladson
xList.append(ladson[:, 0])
yList.append(ladson[:, 1])
yLegendList.append("Ladson")
yColorList.append("gray")
yStyleList.append(".-")

# Cumulant AllOne, D/res = 2000, LES Sm 0.1
xList.append(np.array([10]))
yList.append(np.array([0.848]))
yLegendList.append("Cumulant AllOne, c/res = 2000, LES Sm = 0.1")
yColorList.append("purple")
yStyleList.append("X") 

# ==========================================
# 2. PLOT STYLING (LaTeX Enabled)
# ==========================================
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "font.size": 16,
    "axes.titlesize": 16,
    "axes.labelsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "legend.fontsize": 16,
})


# ==========================================
# 3. PLOT GENERATION
# ==========================================
# Single subplot, maintaining the 16:9 aspect ratio
fig, ax = plt.subplots(1, 1, figsize=(16, 9))

# Loop through all datasets and plot them
for i in range(len(yList)):
    ax.plot(xList[i], yList[i], yStyleList[i], color=yColorList[i], label=yLegendList[i], linewidth=2.0, markersize=8, alpha=0.8)

# Apply limits
ax.set_xlim(xMin, xMax)
ax.set_ylim(yMin, yMax)

# Set labels and title
ax.set_xlabel(xLabel)
ax.set_ylabel(yLabel)
ax.set_title(supTitle)

# Visual grid and locators
# 1. Y-axis Locators: Increase major bins slightly and enable minor ticks
ax.yaxis.set_major_locator(MaxNLocator(nbins=10)) 
ax.yaxis.set_minor_locator(AutoMinorLocator(5)) # Divides each major y-interval into 5 minor intervals
# 2. Denser, highly visible grid
# Major grid lines (solid, darker)
ax.grid(which='major', color='black', linestyle='-', linewidth=0.8, alpha=0.6)
# Minor grid lines (dashed, lighter, fills in the log scale and new y-axis minor ticks)
ax.grid(which='minor', color='dimgray', linestyle='--', linewidth=0.5, alpha=0.4)

# Add a legend, mimicking the style of the old square bbox
ax.legend(loc='upper left', frameon=True, fancybox=False, edgecolor='black', framealpha=0.8, borderpad=0.5)

# Clean layout
plt.tight_layout()

# Save the plot
plt.savefig(filename, dpi=300)
plt.close()

print(f"Successfully plotted and saved to: {filename}")
