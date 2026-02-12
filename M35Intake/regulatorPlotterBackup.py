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
    """Calculates last 10% average, draws the line, and returns the formatted string."""
    count = len(data)
    window = max(1, count // 10)
    avg_val = np.mean(data[-window:])
    
    # Draw horizontal line for the last 10%
    # We use iterations[-window:] to ensure the line starts exactly at the 10% mark
    ax.hlines(y=avg_val, xmin=iterations[-window], xmax=iterations[-1], 
              colors='black', linestyles='--', linewidth=2.0, alpha=0.9, zorder=3)
    
    return f"{avg_val:{fmt}}{unit_str}"

def plot_regulator():
    try:
        with open("/dev/shm/regulator_history.bin", "rb") as f:
            count = np.fromfile(f, dtype=np.int32, count=1)[0]
            i_reg = np.fromfile(f, dtype=np.float32, count=count)
            m_flow = np.fromfile(f, dtype=np.float32, count=count)
            eta = np.fromfile(f, dtype=np.float32, count=count)
    except Exception:
        return

    fig, axs = plt.subplots(3, 1, figsize=(16, 9), sharex=True)
    iterations = np.arange(count)
    bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1, alpha=0.8)

    # --- Subplot 1: Regulator Value ---
    axs[0].plot(iterations, i_reg, color='black', linewidth=1.2, alpha=0.7)
    axs[0].set_ylabel(r"Regulator $\Delta \rho$ [1]")
    axs[0].set_title(r"\textbf{History}")
    axs[0].grid(True, linestyle='--', alpha=0.4)
    axs[0].yaxis.set_major_locator(MaxNLocator(nbins=8))
    set_smart_ylim(axs[0], i_reg)
    
    label_i = add_average_diagnostics(axs[0], iterations, i_reg, fmt=".6f")
    axs[0].text(0.98, 0.90, f'\\textbf{{{label_i}}}', 
                transform=axs[0].transAxes, ha='right', va='top', bbox=bbox_props)

    # --- Subplot 2: Mass Flow ---
    axs[1].plot(iterations, m_flow, color='black', linewidth=1.2, alpha=0.7)
    axs[1].set_ylabel(r"$\dot{m}$ [kg/s]")
    axs[1].grid(True, linestyle='--', alpha=0.4)
    axs[1].yaxis.set_major_locator(MaxNLocator(nbins=8))
    set_smart_ylim(axs[1], m_flow)
    
    label_m = add_average_diagnostics(axs[1], iterations, m_flow, unit_str=" kg/s")
    axs[1].text(0.98, 0.90, f'\\textbf{{{label_m}}}', 
                transform=axs[1].transAxes, ha='right', va='top', bbox=bbox_props)

    # --- Subplot 3: Intake Efficiency ---
    axs[2].plot(iterations, eta, color='black', linewidth=1.2, alpha=0.7)
    axs[2].set_ylabel(r"$\eta_{Intake}$ [1]")
    axs[2].set_xlabel(r"Iteration")
    axs[2].grid(True, linestyle='--', alpha=0.4)
    axs[2].yaxis.set_major_locator(MaxNLocator(nbins=8))
    set_smart_ylim(axs[2], eta)
    
    label_e = add_average_diagnostics(axs[2], iterations, eta)
    axs[2].text(0.98, 0.90, f'\\textbf{{{label_e}}}', 
                transform=axs[2].transAxes, ha='right', va='top', bbox=bbox_props)

    plt.tight_layout()
    os.makedirs("results", exist_ok=True)
    plt.savefig("results/000History.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_regulator()
