import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# Enable LaTeX rendering
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

def set_smart_ylim(ax, data):
    """Sets y-limits based on the second half of data with a 5% margin."""
    if len(data) < 2:
        return
    mid = len(data) // 2
    recent_data = data[mid:]
    y_min, y_max = np.min(recent_data), np.max(recent_data)
    padding = (y_max - y_min) * 0.05
    if padding == 0: padding = 1e-7
    ax.set_ylim(y_min - padding, y_max + padding)

def add_average_diagnostics(ax, iterations, data, unit_str="", fmt=".3f"):
    """Calculates last 20% average, draws the line, and returns the formatted string."""
    count = len(data)
    window = max(1, count // 5)
    avg_val = np.mean(data[-window:])
    
    ax.hlines(y=avg_val, xmin=iterations[-window], xmax=iterations[-1], 
              colors='black', linestyles='--', linewidth=2.0, alpha=0.9, zorder=3)
    
    return f"{avg_val:{fmt}}{unit_str}"

def plot_history(file_number):
    try:
        with open("/dev/shm/historyData.bin", "rb") as f:
            count = np.fromfile(f, dtype=np.int32, count=1)[0]
            # Read mass flow first, then ETA
            mass_flow = np.fromfile(f, dtype=np.float32, count=count)
            eta = np.fromfile(f, dtype=np.float32, count=count)
    except Exception as e:
        print(f"Error reading binary data: {e}")
        return

    # Two subplots stacked vertically, sharing the x-axis
    fig, axs = plt.subplots(2, 1, figsize=(16, 9), sharex=True)
    iterations = np.arange(count)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1, alpha=0.8)

    # Set the main figure title
    fig.suptitle(r"\textbf{Mass Flow, Intake Efficiency History}", fontsize=20)

    # --- Mass Flow Plot (Top) ---
    ax_mf = axs[0]
    ax_mf.plot(iterations, mass_flow, color='black', linewidth=1.2, alpha=0.7)
    ax_mf.set_ylabel(r"$m_p$ [kg/s]")
    ax_mf.grid(True, linestyle='--', alpha=0.4)
    ax_mf.yaxis.set_major_locator(MaxNLocator(nbins=8))
    
    set_smart_ylim(ax_mf, mass_flow)
    label_mf = add_average_diagnostics(ax_mf, iterations, mass_flow)
    ax_mf.text(0.98, 0.95, f'\\textbf{{Avg $m_p$: {label_mf}}}', 
               transform=ax_mf.transAxes, ha='right', va='top', bbox=bbox_props)

    # --- Intake Efficiency (ETA) Plot (Bottom) ---
    ax_eta = axs[1]
    ax_eta.plot(iterations, eta, color='black', linewidth=1.2, alpha=0.7)
    ax_eta.set_ylabel(r"$\eta$ [1]")
    ax_eta.set_xlabel(r"Iteration")
    ax_eta.grid(True, linestyle='--', alpha=0.4)
    ax_eta.yaxis.set_major_locator(MaxNLocator(nbins=8))
    
    set_smart_ylim(ax_eta, eta)
    label_eta = add_average_diagnostics(ax_eta, iterations, eta)
    ax_eta.text(0.98, 0.95, rf'\textbf{{Avg $\eta$: {label_eta}}}', 
                transform=ax_eta.transAxes, ha='right', va='top', bbox=bbox_props)

    # Adjust layout to accommodate the suptitle nicely
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    os.makedirs("results/History", exist_ok=True)
    filename = str(file_number) + "History.png"
    plt.savefig(os.path.join("results/History", filename), dpi=300)
    plt.close()

if __name__ == "__main__":
    file_num = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    plot_history(file_num)
