import os
import numpy as np
import pyvista as pv

# --- Configuration ---
FILE_PATH = "/dev/shm/sim_data.bin"

# Visualization Settings
N_STREAMLINES = 1000     # Density of lines
LINE_WIDTH = 1.5         # Thickness
MAX_LEN = 1000           # Length of streamlines (increased to ensure they span domain)
OPACITY = 1.0            

# View Controls (Applies AFTER auto-centering)
VIEW_AZIMUTH = -150.0     # Rotate camera left/right (Degrees)
VIEW_ELEVATION = 10.0    # Rotate camera up/down (Degrees)
ZOOM_FACTOR = 1.5        # 1.0 = Fit to screen. 1.2 = Zoom in 20%. 0.8 = Zoom out.

# 1. Read Binary Data
with open(FILE_PATH, "rb") as f:
    plotNumber = np.fromfile(f, dtype=np.int32, count=1)[0]
    dims = np.fromfile(f, dtype=np.int32, count=4)
    nx, ny, nz, n_vars = dims
    data = np.fromfile(f, dtype=np.float32).reshape((nx, ny, nz, n_vars))

# 2. Extract Variables (Flattened for PyVista/VTK)
ux = data[:, :, :, 1].flatten(order='F')
uy = data[:, :, :, 2].flatten(order='F')
uz = data[:, :, :, 3].flatten(order='F')
mask = data[:, :, :, 4].flatten(order='F')

# 3. Create Grid
grid = pv.ImageData(dimensions=(nx, ny, nz))
vectors = np.column_stack((ux, uy, uz))
grid.point_data["velocity"] = vectors
grid.point_data["mask"] = mask
grid["mag"] = np.linalg.norm(vectors, axis=1)

# 4. Filter Solid Parts
fluid = grid.threshold(value=[0.0, 0.5], scalars="mask")

# 5. Generate Streamlines
# source_radius checks points within the bounding box of the fluid
streams = fluid.streamlines(
    vectors="velocity", 
    n_points=N_STREAMLINES, 
    max_steps=MAX_LEN,
    integration_direction="both"
)

# 6. Plotting
# off_screen=True prevents the window from popping up
pl = pv.Plotter(off_screen=True, window_size=[1920, 1080])

pl.add_mesh(
    streams, 
    scalars="mag", 
    cmap="viridis", 
    line_width=LINE_WIDTH, 
    opacity=OPACITY,
    lighting=True,
    show_scalar_bar=False
)

# Custom Colorbar
pl.add_scalar_bar(
    title="Velocity Magnitude [m/s]",
    title_font_size=20,
    label_font_size=16,
    shadow=True,
    n_labels=5,
    fmt="%.1f",
    font_family="arial",
    position_x=0.85,  # Position on right side
    position_y=0.05,
    width=0.1,
    height=0.7
)

pl.remove_bounds_axes()
pl.set_background("white") 

# --- CAMERA LOGIC (The Fix) ---

# 1. Force "Y" to be the Up direction
pl.camera.up = (0, 1, 0)

# 2. Auto-center the camera on the data
pl.reset_camera()

# 3. Apply User Rotations (Azimuth/Elevation)
# Note: We orbit around the focal point (the data center)
pl.camera.azimuth += VIEW_AZIMUTH
pl.camera.elevation += VIEW_ELEVATION

# 4. Apply Zoom
pl.camera.zoom(ZOOM_FACTOR)

# 7. Save Result
os.makedirs("results/3D", exist_ok=True)
save_path = os.path.join("results/3D", f"{plotNumber}.png")
pl.screenshot(save_path)
print(f"3D plot saved to {save_path}")

pl.close()
