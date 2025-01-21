import h5py
import numpy as np
import matplotlib.pyplot as plt
import glob
import re
import os

def find_latest_data_file(pattern="Data_*.h5"):
    """
    Finds the latest Data_*.h5 file by sorting filenames numerically.
    Assumes that the filenames have numbers, e.g., Data_1.h5, Data_2.h5, etc.
    """
    data_files = glob.glob(pattern)
    
    if not data_files:
        raise FileNotFoundError(f"No files matching pattern '{pattern}' found.")
    
    # Extract the numerical part from the filenames and sort by that
    data_files.sort(key=lambda x: int(re.findall(r'\d+', x)[0]))  # Sort numerically by the number in the filename
    
    latest_file = data_files[-1]
    print(f"Latest data file found: {latest_file}")
    return latest_file

def extract_velocity_profile(h5_file):
    """
    Extracts the u-velocity profile from the specified .h5 file.
    
    Parameters:
        h5_file (str): Path to the .h5 file.
        
    Returns:
        x, y (numpy.ndarray): Arrays of x and y coordinates.
        u (numpy.ndarray): 2D velocity field (velocity in the x-direction).
    """
    with h5py.File(h5_file, 'r') as f:
        # Access grid information
        grid = f['grid']
        xc = grid['xc'][:]  # Cell center x-coordinates
        yc = grid['yc'][:]  # Cell center y-coordinates
        
        # Access u-velocity data
        u = f['u'][:]  # Shape: [z, y, x]
    
    # Remove the last y and x indices to account for the staggered grid
    u = u[:, :-1, :-1]  # Shape: [z, y, x]
    
    # Compute the slice at the center of the domain in the z-direction
    z_center_idx = u.shape[0] // 2  # Get the index for the center in the z-direction
    u_slice = u[z_center_idx, :, :]  # Take the slice in the x-y plane at z_center
    
    return xc[:-1], yc[:-1], u_slice

def plot_velocity_field_with_particle(x, y, u, particle_center, particle_radius, output_file=None):
    """
    Plots the 2D fluid velocity field with a particle at the center.
    
    Parameters:
        x (numpy.ndarray): X-coordinates of the grid.
        y (numpy.ndarray): Y-coordinates of the grid.
        u (numpy.ndarray): 2D velocity field.
        particle_center (tuple): (x, y) coordinates of the particle center.
        particle_radius (float): Radius of the particle.
        output_file (str): File path to save the plot (optional).
    """
    # Set up the figure and axis
    fig, ax = plt.subplots(figsize=(8, 6))

    # Create a pseudocolor plot for the velocity field
    cmap = plt.get_cmap('jet')  # Color map similar to the example
    norm = plt.Normalize(vmin=0, vmax=1.5)  # Normalize to velocity range

    # Plot the velocity field
    c = ax.pcolormesh(x, y, u, cmap=cmap, norm=norm, shading='auto')

    # Add a color bar to show the velocity magnitude
    cbar = fig.colorbar(c, ax=ax)
    cbar.set_label(r'$u/u_{ref}$', fontsize=12)

    # Plot the particle as a filled circle
    particle = plt.Circle(particle_center, particle_radius, color='gray', zorder=10)
    ax.add_patch(particle)

    # Add labels and set aspect ratio
    ax.set_xlabel(r'$x$', fontsize=12)
    ax.set_ylabel(r'$y$', fontsize=12)
    ax.set_aspect('equal', 'box')

    # Optional: Save the plot if an output path is provided
    if output_file:
        plt.savefig(output_file, dpi=300)
        print(f"Plot saved as: {output_file}")

    # Show the plot
    plt.show()

def main():
    # Particle parameters
    particle_center = (1, 1)  # Adjust to your particle's center in the domain
    particle_radius = 0.5         # Adjust to your particle's radius

    # Find the latest data file
    latest_file = find_latest_data_file()

    # Extract the velocity profile from the file
    x, y, u_slice = extract_velocity_profile(latest_file)

    # Plot the velocity field with the particle
    plot_velocity_field_with_particle(x, y, u_slice, particle_center, particle_radius, output_file="velocity_field.png")

if __name__ == "__main__":
    main()
