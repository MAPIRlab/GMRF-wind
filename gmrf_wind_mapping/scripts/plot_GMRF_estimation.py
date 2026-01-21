import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# Configuration
csv_file = "NLPD_optimization_results_central_obstacle_1ms_negative/gmrf_opt_estimation_obs_51.csv"
base_name = os.path.splitext(csv_file)[0]
out_file = base_name + '.png'


def main():
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")

    # 1. Load Metadata manually from the first few lines
    meta = {}
    with open(csv_file, 'r') as f:
        for _ in range(3):
            line = f.readline().strip().split(',')
            meta[line[0]] = float(line[1])

    nx, ny = int(meta['Dimensions_x']), int(meta['Dimensions_y'])
    
    # Load the actual data (skipping the metadata lines)
    df = pd.read_csv(csv_file, skiprows=3)

    # 2. Reshape flat columns into 2D grids (Y, X) for plotting
    # We use order='C' or 'F' depending on how your cell_index was generated
    def to_grid(series):
        return series.values.reshape((ny, nx))

    est_x = to_grid(df['gmrf_wind_x'])
    est_y = to_grid(df['gmrf_wind_y'])
    std_x = to_grid(df['gmrf_std_x'])
    std_y = to_grid(df['gmrf_std_y'])
    gt_x = to_grid(df['gt_wind_x'])
    gt_y = to_grid(df['gt_wind_y'])

    # Create Mask for Obstacles (where data is NaN)
    # We'll use this to color cells black
    mask = np.isnan(est_x)

    # 3. Compute Metrics Cell-by-Cell
    # RMSE = sqrt( (u_est - u_gt)^2 + (v_est - v_gt)^2 )
    rmse_grid = np.sqrt((est_x - gt_x)**2 + (est_y - gt_y)**2)

    # Combined Uncertainty (Total Standard Deviation)
    uncertainty_grid = np.sqrt(std_x**2 + std_y**2)

    # NLPD for 2D Independent Gaussian: 
    # NLPD = 0.5 * [ ln(2π * σx²) + (x-μx)²/σx² + ln(2π * σy²) + (y-μy)²/σy² ]
    def calc_nlpd(mu_x, mu_y, sx, sy, tx, ty):
        term_x = np.log(2 * np.pi * sx**2) + ((tx - mu_x)**2 / (sx**2))
        term_y = np.log(2 * np.pi * sy**2) + ((ty - mu_y)**2 / (sy**2))
        return 0.5 * (term_x + term_y)

    nlpd_grid = calc_nlpd(est_x, est_y, std_x, std_y, gt_x, gt_y)

    # 4. Visualization
    fig, axs = plt.subplots(2, 3, figsize=(20, 10))
    cmap_error = 'YlOrRd' # Yellow to Red for errors
    cmap_uncert = 'Purples'  # Distinct color for uncertainty

    # Helper to apply black obstacles
    def apply_obstacles(ax):
        ax.imshow(np.where(mask, 1, np.nan), cmap='gray', vmin=0, vmax=1, origin='lower')

    # --- Row 1: The Wind Fields ---
    # [0,0] GMRF Estimation (Quiver) ---
    apply_obstacles(axs[0,0])
    axs[0,0].quiver(est_x, est_y, color='blue')
    axs[0,0].set_title("GMRF Wind Estimation")

    # [0,1] Ground Truth (Quiver) ---
    apply_obstacles(axs[0,1])
    axs[0,1].quiver(gt_x, gt_y, color='green')
    axs[0,1].set_title("Ground Truth Wind")

    # [0,2] Combined Uncertainty (The new plot)
    apply_obstacles(axs[0,2])
    im_u = axs[0,2].imshow(uncertainty_grid, cmap=cmap_uncert, origin='lower')
    plt.colorbar(im_u, ax=axs[0,2], label='Combined Sigma')
    axs[0,2].set_title("Total Uncertainty ($\sqrt{\sigma_x^2 + \sigma_y^2}$)")

    # --- Row 2: Performance Metrics ---
    # [1,0] RMSE Heatmap ---
    apply_obstacles(axs[1,0])
    im3 = axs[1,0].imshow(rmse_grid, cmap=cmap_error, origin='lower')
    plt.colorbar(im3, ax=axs[1,0], label='RMSE')
    axs[1,0].set_title("RMSE (Velocity Error)")

    # [1,1] NLPD Heatmap ---
    apply_obstacles(axs[1,1])
    im4 = axs[1,1].imshow(nlpd_grid, cmap=cmap_error, origin='lower')
    plt.colorbar(im4, ax=axs[1,1], label='NLPD')
    axs[1,1].set_title("NLPD (Uncertainty-Aware Error)")

    # [1,2] Empty or Comparison Plot (Let's show Error vs Uncertainty ratio)
    # This helps see where the model is "confidently wrong"
    apply_obstacles(axs[1,2])
    conf_wrong = rmse_grid / (uncertainty_grid + 1e-6)
    im_c = axs[1,2].imshow(conf_wrong, cmap='coolwarm', origin='lower', vmin=0, vmax=5)
    plt.colorbar(im_c, ax=axs[1,2], label='RMSE / Sigma')
    axs[1,2].set_title("Calibration (RMSE/$\sigma$) ratio")

    for ax in axs.flat:
        ax.set_aspect('equal')
        ax.set_xticks([]) # Optional: hide ticks for cleaner look
        ax.set_yticks([])

    plt.tight_layout()
    plt.savefig('wind_analysis.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()