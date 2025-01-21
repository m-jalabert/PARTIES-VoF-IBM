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
        y (numpy.ndarray): Array of adimensional y-coordinates .
        u_mean (numpy.ndarray): Averaged adimensional u-velocity profile .
    """
    with h5py.File(h5_file, 'r') as f:
        # Access grid information
        grid = f['grid']
        xc = grid['xc'][:]  # Cell center x-coordinates
        yc = grid['yc'][:]  # Cell center y-coordinates
        zc = grid['zc'][:]  # Cell center z-coordinates
        xu = grid['xu'][:]  # Cell wall x-coordinates
        yv = grid['yv'][:]  # Cell wall y-coordinates
        zw = grid['zw'][:]  # Cell wall z-coordinates
        
        # Access u-velocity data
        u = f['u'][:]  # Shape: [z, y, x]
    
    # Remove the last y and x indices to account for the staggered grid
    u = u[:, :-1, :-1]  # Now shape: [z, y, x]
    
    # Extract y-coordinates corresponding to the u-velocity
    y = yc[:-1]  # Exclude the last y-coordinate
    
    # Compute the mean velocity by averaging over z and x directions
    u_mean = np.mean(u, axis=(0, 2))  # Shape: [y]
    
    return y, u_mean

def compute_analytical_profile(Re, dpdx, mu):
    """
    Computes the analytical velocity profile for Plane Poiseuille Flow.
    
    Parameters:
        V_bulk (float): Mean (bulk) velocity (m/s).
        h (float): Gap height between plates (m).
        mu (float): Dynamic viscosity (Pa·s).
        
    Returns:
        u_analytical (numpy.ndarray): Analytical velocity profile (m/s).
    """

    u_analytical = Re/2 * (-dpdx) * (y-y**2)
    return u_analytical

def plot_velocity_profiles(y, u_simulated, u_analytical, h, save_plot=True, output_dir="plots"):
    """
    Plots the simulated and analytical velocity profiles.
    
    Parameters:
        y (numpy.ndarray): Array of y-coordinates (meters).
        u_simulated (numpy.ndarray): Simulated mean u-velocity profile (m/s).
        u_analytical (numpy.ndarray): Analytical velocity profile (m/s).
        h (float): Gap height between plates (m).
        save_plot (bool): Whether to save the plot as an image file.
        output_dir (str): Directory where the plot image will be saved.
    """
    plt.figure(figsize=(8, 6))
    plt.plot(u_simulated, y, label='Simulated Bulk Velocity', marker='o', linestyle='-', markersize=4)
    plt.plot(u_analytical, y, label='Analytical Velocity Profile', linestyle='--')
    plt.xlabel('Adimensional Velocity')
    plt.ylabel('Adimensional y')
    plt.title('Plane Poiseuille Flow: Simulated vs. Analytical Velocity Profile')
    plt.legend()
    plt.grid(True)
    plt.gca().invert_yaxis()  # Optional: Invert y-axis to match typical flow diagrams
    plt.tight_layout()
    
    if save_plot:
        # Create output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        # Define the filename based on parameters or timestamp
        filename = os.path.join(output_dir, "velocity_profile.png")
        plt.savefig(filename, dpi=300)
        print(f"Plot saved as: {filename}")
    
    plt.show()

def check_mean_velocity(u_simulated, target_velocity=1):
    """
    Checks if the mean of the simulated velocity matches the expected bulk velocity.

    Parameters:
        u_simulated (numpy.ndarray): Simulated velocity profile.
        target_velocity (float): Expected bulk velocity (default is 10 m/s).
    """
    mean_velocity = np.mean(u_simulated)
    if np.isclose(mean_velocity, target_velocity, atol=0.1):
        print(f"Mean adimensional velocity is correct: {mean_velocity:.2f} ")
    else:
        print(f"Warning: Mean adimensional velocity is {mean_velocity:.2f} , expected {target_velocity:.2f} ")    

def main():
    # Simulation Parameters
    V_bulk = 1       # Adimensional Mean velocity
    rho_f = 1000      # Fluid density (kg/m^3)
    mu = 0.05         # Dynamic viscosity (Pa·s)
    h = 0.02          # Gap height between plates (m)
    Re = 4            # Reynolds number
    dpdx = 3          # Constant Pressure gradient adimensional
    
    # Find the latest Data_*.h5 file
    try:
        latest_data_file = find_latest_data_file()
    except FileNotFoundError as e:
        print(e)
        return
    
    # Extract the velocity profile from the latest data file
    try:
        y, u_simulated = extract_velocity_profile(latest_data_file)
    except KeyError as e:
        print(f"Error accessing datasets in the file: {e}")
        return
    
    # Compute the analytical velocity profile
    u_analytical = Re/2 * (-dpdx) * (y-y**2)
    
    # Plot the results
    plot_velocity_profiles(y, u_simulated, u_analytical, h)

    # Check if the mean velocity is close to the target value
    check_mean_velocity(-u_simulated, target_velocity=V_bulk)


    # Optionally, save data to a CSV file for further analysis
    save_data = True  # Set to True if you want to save the data
    if save_data:
        import pandas as pd
        data = pd.DataFrame({
            'y (m)': y,
            'u_simulated (m/s)': u_simulated,
            'u_analytical (m/s)': u_analytical
        })
        csv_filename = os.path.join("plots", "velocity_profiles.csv")
        data.to_csv(csv_filename, index=False)
        print(f"Data saved as: {csv_filename}")

if __name__ == "__main__":
    main()
