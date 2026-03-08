import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, AutoMinorLocator

ladson = np.array([
	[-4.05, -0.4280],
    [-2.00, -0.2150],
    [ 0.05,  0.0040],
    [ 1.98,  0.2080],
    [ 4.18,  0.4520],
    [ 6.20,  0.6630],
    [ 8.22,  0.8800],
    [10.18,  1.0880],
    [11.08,  1.1800],
    [12.25,  1.2920],
    [13.10,  1.3680],
    [14.28,  1.4580],
    [15.20,  1.5280],
    [16.18,  1.5900],
    [16.90,  1.6180],
    [17.35,  1.6600],
    [17.65,  1.6450],
    [18.65,  1.0050]])
    
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

# Cumulant AllOne, D/res = 500, LES Sm 0.1
xList.append(np.array([10]))
yList.append(np.array([0.863]))
yLegendList.append("Cumulant AllOne, c/res = 500, LES Sm = 0.1")
yColorList.append("black")
yStyleList.append("X") 

# Cumulant AllOne, D/res = 1000, LES Sm 0.1
xList.append(np.array([10]))
yList.append(np.array([0.765]))
yLegendList.append("Cumulant AllOne, c/res = 1000, LES Sm = 0.1")
yColorList.append("teal")
yStyleList.append("X") 

# Cumulant AllOne, D/res = 2000, LES Sm 0.1
xList.append(np.array([10]))
yList.append(np.array([0.847]))
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
