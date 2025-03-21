import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # needed for 3D plotting

#############################################
# Functions for loading VOF data from HDF5
#############################################

def load_F(file_path):
    """Load the volume fraction field F from the /VOF/F dataset."""
    with h5py.File(file_path, "r") as f:
        F = f["/VOF/F"][:]  # Expect F to be a 3D array (NZ, NY, NX)
    return F

def load_time(file_path):
    """Load the simulation time from the /time dataset."""
    with h5py.File(file_path, "r") as f:
        t = f["/time"][0]
    return t

#############################################
# Main Processing Code
#############################################

# Create a list of files from 0 to 20.
files = [f"Data_{i}.h5" for i in range(21)]

# Containers for storing data
times = []
F_all_cells = []
center_F_vals = []

# Loop over files, load data and extract information.
for file_path in files:
    F = load_F(file_path)
    t = load_time(file_path)
    times.append(t)
    F_all_cells.append(F.flatten())  # Store flattened F values
    
    # Extract center cell volume fraction using grid center:
    NZ, NY, NX = F.shape
    center_F = F[NZ//2, NY//2, NX//2]
    center_F_vals.append(center_F)

# Convert to numpy arrays
times = np.array(times)
F_all_cells = np.array(F_all_cells)  # Shape: (num_files, total_cells)
center_F_vals = np.array(center_F_vals)

# ---------------------------
# Plot F vs. time for all cells
# ---------------------------
plt.figure(figsize=(8, 6))
for i in range(F_all_cells.shape[1]):  # Loop over all cells
    plt.plot(times, F_all_cells[:, i], alpha=0.3, linewidth=0.8)

plt.xlabel("Time")
plt.ylabel("Volume Fraction (F)")
plt.title("Volume Fraction vs. Time for All Cells")
plt.grid(True)
plt.savefig("all_cells_F_vs_time.png", dpi=300)
plt.close()

# ---------------------------
# Plot center cell F vs. time
# ---------------------------
plt.figure(figsize=(6, 4))
plt.plot(times, center_F_vals, 'o-', markersize=8, linewidth=2)
plt.xlabel("Time")
plt.ylabel("Center Cell Volume Fraction")
plt.title("Center Cell Volume Fraction vs. Time")
plt.grid(True)
plt.savefig("center_cell_F_vs_time.png", dpi=300)
plt.close()

# ---------------------------
# Plot 3D colormap of F from the last file
# ---------------------------
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Using the F from the last file processed in the loop
NZ, NY, NX = F.shape
# Create a grid of indices that match the array dimensions.
# Note: F is indexed as F[z, y, x]
z, y, x = np.indices((NZ, NY, NX))
# Flatten the arrays and plot a scatter plot with color mapping based on F values.
sc = ax.scatter(x.flatten(), y.flatten(), z.flatten(), c=F.flatten(), cmap='viridis')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.title("3D Colormap of Volume Fraction F (last file)")
plt.colorbar(sc, shrink=0.5)
plt.savefig("F_3D_colormap.png", dpi=300)
plt.close(fig)
print("Saved 3D colormap plot as F_3D_colormap.png")

# ---------------------------
# Collect detailed debug info from all files and save to a text file.
# ---------------------------
def get_vof_debug(file_path):
    """Collects debugging information from VOF-related datasets."""
    out = []
    with h5py.File(file_path, "r") as f:
        out.append(f"File: {file_path}\n")
        out.append(f"Time = {f['/time'][0]:.3f}\n")
        grid_group = f["/grid"]
        NX, NY, NZ = int(grid_group["NX"][0]), int(grid_group["NY"][0]), int(grid_group["NZ"][0])
        out.append(f"Grid dimensions: NX={NX}, NY={NY}, NZ={NZ}\n")
        F = f["/VOF/F"][:]
        out.append(f"F min={np.min(F):.4f}, max={np.max(F):.4f}, sum={np.sum(F):.4f}\n")
    return "".join(out)

debug_output = "".join([get_vof_debug(file_path) for file_path in files])
with open("vof_debug.txt", "w") as fout:
    fout.write(debug_output)
print("Debug output saved to vof_debug.txt")

# ---------------------------
# Write all F values for all cells through time to a new text file.
# ---------------------------
with open("all_cells_F_values.txt", "w") as fout:
    # Create a header listing time and each cell index
    header = "Time," + ",".join([f"F_cell_{i}" for i in range(F_all_cells.shape[1])]) + "\n"
    fout.write(header)
    for i in range(len(times)):
        line = f"{times[i]:.3f}," + ",".join([f"{val:.4f}" for val in F_all_cells[i]]) + "\n"
        fout.write(line)
print("All cell F values saved to all_cells_F_values.txt")

# ---------------------------
# Extract and save normal and flux for cell (3,3,3) and its direct neighbors,
# and add the value of F for each cell, in normal_fluxes.txt.
# ---------------------------
# Define the cell indices for (3,3,3) and its six direct neighbors.
# The indexing order is (z, y, x).
neighbors = {
    "center": (3, 3, 3),
    "x_minus": (3, 3, 2),
    "x_plus":  (3, 3, 4),
    "y_minus": (3, 2, 3),
    "y_plus":  (3, 4, 3),
    "z_minus": (2, 3, 3),
    "z_plus":  (4, 3, 3)
}

with open("normal_fluxes.txt", "w") as fout:
    fout.write("Time, Cell, normal_x, normal_y, normal_z, flux_x, flux_y, flux_z, F\n")
    for file_path in files:
        t = load_time(file_path)
        with h5py.File(file_path, "r") as f:
            normal_x = f["/VOF/normal_x"][:]
            normal_y = f["/VOF/normal_y"][:]
            normal_z = f["/VOF/normal_z"][:]
            flux_x = f["/VOF/flux_x"][:]
            flux_y = f["/VOF/flux_y"][:]
            flux_z = f["/VOF/flux_z"][:]
            F_data = f["/VOF/F"][:]
            # Loop over each cell of interest (center and its neighbors)
            for cell, (k, j, i) in neighbors.items():
                n_x = normal_x[k, j, i]
                n_y = normal_y[k, j, i]
                n_z = normal_z[k, j, i]
                f_x = flux_x[k, j, i]
                f_y = flux_y[k, j, i]
                f_z = flux_z[k, j, i]
                F_val = F_data[k, j, i]
                fout.write(f"{t:.3f}, {cell}, {n_x:.4f}, {n_y:.4f}, {n_z:.4f}, {f_x:.4f}, {f_y:.4f}, {f_z:.4f}, {F_val:.4f}\n")
print("Normal, flux, and F values for (3,3,3) and its neighbors saved to normal_fluxes.txt")
