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

def plot_regulator():
    try:
        # Read the binary data from shared memory
        with open("/dev/shm/regulator_history.bin", "rb") as f:
            count = np.fromfile(f, dtype=np.int32, count=1)[0]
            i_reg = np.fromfile(f, dtype=np.float32, count=count)
            m_flow = np.fromfile(f, dtype=np.float32, count=count)
            eta = np.fromfile(f, dtype=np.float32, count=count)
    except Exception:
        # Silently exit if file is being written to or doesn't exist yet
        return

    # Create figure with 16:9 aspect ratio
    fig, axs = plt.subplots(3, 1, figsize=(16, 9), sharex=True)
    iterations = np.arange(count)
    
    # Text box properties for the top-right status values
    bbox_props = dict(boxstyle="square,pad=0.3", fc="white", ec="black", lw=1, alpha=0.8)

    # --- Subplot 1: Regulator Value (Density Offset) ---
    axs[0].plot(iterations, i_reg, color='black', linewidth=1.2)
    axs[0].set_ylabel(r"i Regulator $\Delta \rho$ [1]")
    axs[0].set_title(r"\textbf{History}")
    
    # Denser grid: Force ~10 divisions on Y
    axs[0].yaxis.set_major_locator(MaxNLocator(nbins=8))
    axs[0].grid(True, linestyle='--', alpha=0.6)
    
    # Current value annotation
    axs[0].text(0.98, 0.90, f'\\textbf{{{i_reg[-1]:.6f}}}', 
                transform=axs[0].transAxes, ha='right', va='top', bbox=bbox_props)

    # --- Subplot 2: Mass Flow ---
    axs[1].plot(iterations, m_flow, color='black', linewidth=1.2)
    axs[1].set_ylabel(r"$\dot{m}$ [kg/s]")
    
    # Denser grid: Force ~10 divisions on Y
    axs[1].yaxis.set_major_locator(MaxNLocator(nbins=8))
    axs[1].grid(True, linestyle='--', alpha=0.6)
    
    # Current value annotation
    axs[1].text(0.98, 0.90, f'\\textbf{{{m_flow[-1]:.3f} kg/s}}', 
                transform=axs[1].transAxes, ha='right', va='top', bbox=bbox_props)

    # --- Subplot 3: Intake Efficiency ---
    axs[2].plot(iterations, eta, color='black', linewidth=1.2)
    axs[2].set_ylabel(r"$\eta_{Intake}$ [1]")
    axs[2].set_xlabel(r"Iteration")
    axs[2].set_ylim(0, 1)
    
    # Denser grid: Force ~10 divisions on Y
    axs[2].yaxis.set_major_locator(MaxNLocator(nbins=8))
    axs[2].grid(True, linestyle='--', alpha=0.6)
    
    # Current value annotation
    axs[2].text(0.98, 0.90, f'\\textbf{{{eta[-1]:.3f}}}', 
                transform=axs[2].transAxes, ha='right', va='top', bbox=bbox_props)

    # Clean up and export
    plt.tight_layout()
    os.makedirs("results", exist_ok=True)
    plt.savefig("results/000History.png", dpi=300)
    plt.close()

if __name__ == "__main__":
    plot_regulator()
