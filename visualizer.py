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

data = np.loadtxt("result.csv", delimiter=",")

plt.figure(figsize=(6, 4))
img = plt.imshow(data, origin="lower", cmap="viridis", aspect="equal")

plt.savefig("result.png", dpi=1000, bbox_inches="tight")
