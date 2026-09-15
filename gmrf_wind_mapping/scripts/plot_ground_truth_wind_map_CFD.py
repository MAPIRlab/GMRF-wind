# This script visualizes the Experimental scenarios together with the ground truth wind data (CFD)
# For each scenario, the occupancy and the wind_at_cell_centers is necessary
# ------------------------------------------------------------
import os
import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import cv2
from scipy.interpolate import griddata

# --- CONFIGURATION ---

DATA_FOLDER = "./install/test_env/share/test_env/scenarios/10x6_maze"
WIND_FILE = os.path.join(DATA_FOLDER, "wind_simulations/1ms/wind_at_cell_centers_0.csv")
MAP_YAML = os.path.join(DATA_FOLDER, "occupancy.yaml")

def load_ros_map(yaml_path):
    """
    Loads a ROS2 occupancy grid map using the YAML metadata and OpenCV.
    Returns:
        occ_array: The image as a numpy array.
        extent: [xmin, xmax, ymin, ymax] for matplotlib plotting.
    """
    # 1. Load YAML metadata
    if not os.path.exists(yaml_path):
        raise FileNotFoundError(f"Could not find YAML file at: {yaml_path}")
    
    with open(yaml_path, 'r') as f:
        config = yaml.safe_load(f)
    
    # ROS2 YAMLs usually provide a relative path to the image
    image_filename = config['image']
    dir_path = os.path.dirname(os.path.abspath(yaml_path))
    image_path = os.path.join(dir_path, image_filename)
    
    # Load the image in grayscale
    # cv2.IMREAD_GRAYSCALE is much more robust for PGM files
    occ_array = cv2.imread(image_path, cv2.IMREAD_GRAYSCALE)
    
    if occ_array is None:
        raise FileNotFoundError(f"OpenCV could not open PGM file at: {image_path}")

    occ_array = cv2.flip(occ_array, 0) # Flip vertically to match ROS map orientation

    # Extract ROS metadata
    resolution = config['resolution']
    origin = config['origin']  # [x, y, yaw]
    
    # Calculate map dimensions in meters
    height, width = occ_array.shape
    xmin = origin[0]
    ymin = origin[1]
    xmax = xmin + (width * resolution)
    ymax = ymin + (height * resolution)
    
    # ROS map convention:
    # 255/white is usually "free" (0% occupancy)
    # 0/black is "occupied" (100% occupancy)
    # In some ROS versions, 254 is free and 0 is occupied.
    # We return the array and the extent for plt.imshow(origin='lower')
    extent = [xmin, xmax, ymin, ymax]
    
    return occ_array, extent

def plot_wind(mode='quiver'):
    # 1. Load Occupancy Map
    occ_array, extent = load_ros_map(MAP_YAML)
    
    # 2. Load Wind Data (points from CFD)
    df = pd.read_csv(WIND_FILE)
    x = df['Points:0'].values
    y = df['Points:1'].values
    u = df['U:0'].values
    v = df['U:1'].values
    
    # 3. Plotting
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plot obstacles (inverted grayscale: high values/free = white, low/obstacle = black)
    ax.imshow(occ_array, cmap='gray', extent=extent, origin='lower', zorder=1)
    
    # 2. Create Interpolation Grid
    # Match the grid density to something reasonable (e.g., 100-200 points)
    xi = np.linspace(extent[0], extent[1], 150)
    yi = np.linspace(extent[2], extent[3], 100)
    X, Y = np.meshgrid(xi, yi)
    
    # 3. Interpolate
    Ui = griddata((x, y), u, (X, Y), method='linear')
    Vi = griddata((x, y), v, (X, Y), method='linear')
    
    # --- THE FIX: MASKING ---
    # Convert real-world (X, Y) coordinates back to pixel indices to check obstacles
    res = (extent[1] - extent[0]) / occ_array.shape[1]
    
    # Calculate which pixel index corresponds to each grid point
    px = ((X - extent[0]) / res).astype(int)
    py = ((Y - extent[2]) / res).astype(int)
    
    # Clip indices to stay within array bounds
    px = np.clip(px, 0, occ_array.shape[1] - 1)
    py = np.clip(py, 0, occ_array.shape[0] - 1)
    
    # Create a mask where the map is "occupied" (usually 0 in ROS PGM)
    # If your obstacles are black (0), mask where occ_array <= threshold
    mask = occ_array[py, px] < 128 
    
    # Set wind to NaN in those areas so streamplot won't draw there
    Ui[mask] = np.nan
    Vi[mask] = np.nan
    # -----------------------

    if mode == 'quiver':
        # Subsample if there are too many points for quiver
        #skip = (slice(None, None, 10)) 
        #ax.quiver(x[skip], y[skip], u[skip], v[skip], color='green', scale=20, width=0.002)
        ax.quiver(x, y, u, v, color='red', zorder=2)
        ax.set_title("Wind Map: Quiver Plot")
        
    elif mode == 'streamplot':
        speed = np.sqrt(Ui**2 + Vi**2)
        speed_clipped = np.clip(speed, 0, 1.0)
        # Use zorder=2 to ensure wind is drawn over the map
        strm =ax.streamplot(X, Y, Ui, Vi, color=speed_clipped, cmap='jet', linewidth=1, density=2.0, zorder=2)
        cbar = fig.colorbar(strm.lines, label='Wind Speed (m/s)')
        cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])
        cbar.ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '≥1.0'])

        ## Create a grid for interpolation
        #xi = np.linspace(extent[0], extent[1], 100)
        #yi = np.linspace(extent[2], extent[3], 100)
        #X, Y = np.meshgrid(xi, yi)
        
        # Interpolate unstructured data onto the grid
        #Ui = griddata((x, y), u, (X, Y), method='linear')
        #Vi = griddata((x, y), v, (X, Y), method='linear')
        
        #speed = np.sqrt(Ui**2 + Vi**2)
        #strm = ax.streamplot(X, Y, Ui, Vi, color=speed, cmap='jet', linewidth=1)
        #fig.colorbar(strm.lines, label='Wind Speed (m/s)')
        #ax.set_title("Wind Map: Streamplot")

    ax.set_xlim(extent[0], extent[1])
    ax.set_ylim(extent[2], extent[3])
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    plt.tight_layout()
    plt.show()

# ====================================================
if __name__ == "__main__":
    try:
        #print("Generating Quiver Plot...")
        #plot_wind(mode='quiver')
        print("Generating Streamplot...")
        plot_wind(mode='streamplot')
    except Exception as e:
        print(f"Error: {e}")