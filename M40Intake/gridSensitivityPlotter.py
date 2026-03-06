import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

# ==========================================
# 1. DATA DEFINITION
# ==========================================
filename = "resultsGridSensitivity"
supTitle = r"\textbf{Grid Sensitivity: Mass Flow vs Resolution}"
xLabel = r"res [mm]"
yLabel = r"$m_p$ [kg/s]"

xMin = 0.1
xMax = 1.2
yMin = 0
yMax = 5

# Grouping for easy iteration
xList = []
yList = []
yLegendList = []
yColorList = []
yStyleList = []

# achenbach 1
xList.append(np.array([1, 0.7, 0.6, 0.5, 0.4, 0.25]))
yList.append(np.array([1.427, 1.946, 3.038, 2.870, 3.426, 4.126]))
yLegendList.append("NACA")
yColorList.append("black")
yStyleList.append(".")

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

#ax.set_xscale('log')

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
ax.legend(loc='upper right', frameon=True, fancybox=False, edgecolor='black', framealpha=0.8, borderpad=0.5)

# Clean layout
plt.tight_layout()

# Save the plot
plt.savefig(filename, dpi=300)
plt.close()

print(f"Successfully plotted and saved to: {filename}")
