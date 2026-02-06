# This script visualizes the Experimental scenarios together with the ground truth wind data (CFD)
# ------------------------------------------------------------
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import glob
import os


# ----------------------------------
# 1. CONFIGURATION
# ----------------------------------
folders = ['NLPD_optimization_results_central_obstacle_01ms_phase1', 
           'NLPD_optimization_results_central_obstacle_05ms_phase1', 
           'NLPD_optimization_results_central_obstacle_1ms_phase1',
           'NLPD_optimization_results_snake_01ms_phase1',
           'NLPD_optimization_results_snake_05ms_phase1',
           'NLPD_optimization_results_snake_1ms_phase1',
           ]

csv_file = "gmrf_opt_estimation_obs_10.csv" # Any file has the Ground-Truth


def main():

    # Create subplots for the 6 scenarios
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    cmap_error = 'YlOrRd' # Yellow to Red for errors
    cmap_uncert = 'Purples'  # Distinct color for uncertainty

    # Helper to apply black to obstacles
    def apply_obstacles(ax, mask):
        # Draw black obstacles
        ax.imshow(np.where(mask, 0, np.nan), cmap='gray', vmin=0, vmax=1, origin='lower')
        # ax.imshow(np.where(mask, 1, np.nan), cmap='gray', vmin=0, vmax=1, origin='lower')

    for k,folder in enumerate(folders):
        # Load a sample CSVs in the folder
        test_file = glob.glob(os.path.join(folder, csv_file))
        test_file = test_file[0]
        if not os.path.isfile(test_file):
            raise FileNotFoundError(f"CSV file not found: {test_file}")

        # 1. Load Metadata manually from the first few lines
        # ------------------------------------------------------------
        meta = {}
        with open(test_file, 'r') as f:
            for _ in range(3):
                line = f.readline().strip().split(',')
                meta[line[0]] = float(line[1])
        nx, ny = int(meta['Dimensions_x']), int(meta['Dimensions_y'])
        
        # Load the actual data (skipping the metadata lines)
        df = pd.read_csv(test_file, skiprows=3)

        # 2. Reshape flat columns into 2D grids (Y, X) for plotting
        # We use order='C' or 'F' depending on how your cell_index was generated
        def to_grid(series):
            return series.values.reshape((ny, nx))
        
        gt_x = to_grid(df['gt_wind_x'])
        gt_y = to_grid(df['gt_wind_y'])

        # Create Mask for Obstacles (where data is NaN)
        # We'll use this to color cells black
        mask = np.isnan(gt_x)

        # 3. Visualization
        # ------------------------------------------------------------  
        # Select subplot
        ax = axs[k//3, k%3]
        # Ground Truth (Quiver)
        ax.quiver(gt_x, gt_y, color='green')
        apply_obstacles(ax, mask)
        ax.set_title("Ground Truth Wind")

    # Set aspect ratio and remove ticks for all subplots
    for ax in axs.flat:
        ax.set_aspect('equal')
        ax.set_xticks([]) # Optional: hide ticks for cleaner look
        ax.set_yticks([])

    plt.tight_layout()
    plt.savefig('scenarios_ground_truth.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()