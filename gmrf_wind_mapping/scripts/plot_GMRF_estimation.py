# Plots different maps to visualize GMRF Wind Estimation Results 
#   GT: Quiver plot
#   GMRF estimation: Quiver + Uncertainty heatmap
#   Performance Metrics: AAE (error in angle), Module Error (error in magnitude), NSP (normalized scalar product)
#   NLPD heatmaps: can also include  the Negative Log Predictive Density to see uncertainty-aware errors
# ===============================================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os


# Configuration
csv_file = "NSP_optimization_results_central_obstacle_1ms/gmrf_opt_estimation_obs_60.csv"
#csv_file ="gmrf_opt_estimation_obs_50.csv"

base_name = os.path.splitext(csv_file)[0]
out_file = base_name + '.png'


def cartesian_to_polar_with_uncertainty(df):
    # Prepare new columns
    df['gmrf_wind_r'] = 0.0
    df['gmrf_wind_theta'] = 0.0
    df['gmrf_std_r'] = 0.0
    df['gmrf_std_theta'] = 0.0

    for index, row in df.iterrows():
        # Cartesian components and variances
        x = row['gmrf_wind_x']
        y = row['gmrf_wind_y']
        vx = row['gmrf_std_x'] ** 2  # Variance in x
        vy = row['gmrf_std_y'] ** 2  # Variance in y

        r, theta, var_r, var_theta = propagate_uncertainty_cartesian_to_polar(x, y, vx, vy)

        # Store results
        df.at[index, 'gmrf_wind_r'] = r
        df.at[index, 'gmrf_wind_theta'] = theta
        df.at[index, 'gmrf_std_r'] = np.sqrt(var_r)
        df.at[index, 'gmrf_std_theta'] = np.sqrt(var_theta)
    return df

def propagate_uncertainty_cartesian_to_polar(x, y, vx, vy):
    # Convert (x, y) to (r, theta) and propagate uncertainties
    # using first-order error propagation (Jacobian method)

    # 1. Compute r and theta
    r = np.sqrt(pow(x, 2) + pow(y, 2))
    theta = np.atan2(y, x)

    # Jacobian elements
    dr_dx = x / r if r > 0 else 0
    dr_dy = y / r if r > 0 else 0
    dth_dx = -y / (r*r) if r > 0 else 0
    dth_dy = x / (r*r) if r > 0 else 0

    # Uncertainty propagation: Sigma_polar = J * Sigma_cartesian * J^T
    # Assuming independence between x and y (no covariance) --> Sigma_cartesian = [vx 0; 0 vy]
    # Resulting Sigma_r^2 = (dr/dx)^2 * vx + (dr/dy)^2 * vy
    var_r = (dr_dx * dr_dx * vx) + (dr_dy * dr_dy * vy)
    # Resulting Sigma_theta^2 = (dth/dx)^2 * vx + (dth/dy)^2 * vy
    var_theta = (dth_dx * dth_dx * vx) + (dth_dy * dth_dy * vy)

    return r, theta, var_r, var_theta



def main():
    if not os.path.isfile(csv_file):
        raise FileNotFoundError(f"CSV file not found: {csv_file}")

    # 1. Load Metadata and Data
    # ====================================
    meta = {}
    with open(csv_file, 'r') as f:
        for _ in range(3):
            line = f.readline().strip().split(',')
            meta[line[0]] = float(line[1])

    nx, ny = int(meta['Dimensions_x']), int(meta['Dimensions_y'])
    
    # Load the actual data (skipping the metadata lines)
    df = pd.read_csv(csv_file, skiprows=3)

    # 2. Convert to polar coordinates (module, direction) and propagate uncertainties
    # ====================================
    df = cartesian_to_polar_with_uncertainty(df)


    # 3. Reshape columns into 2D grids (Y, X) for plotting
    # ====================================
    def to_grid(series):
        return series.values.reshape((ny, nx))

    est_x = to_grid(df['gmrf_wind_x'])
    est_y = to_grid(df['gmrf_wind_y'])
    std_x = to_grid(df['gmrf_std_x'])
    std_y = to_grid(df['gmrf_std_y'])

    est_r = to_grid(df['gmrf_wind_r'])
    est_theta = to_grid(df['gmrf_wind_theta'])
    std_r = to_grid(df['gmrf_std_r'])
    std_theta = to_grid(df['gmrf_std_theta'])

    gt_x = to_grid(df['gt_wind_x'])
    gt_y = to_grid(df['gt_wind_y'])

    # Create Mask for Obstacles (where data is NaN)
    # We'll use this to color cells black
    mask = np.isnan(est_x)


    # ====================================
    # 3. Compute Metrics Cell-by-Cell
    # ====================================

    # 3.1 AAE = Average Angular Error (in degrees)
    def angular_error(u1, v1, u2, v2):
        dot = u1 * u2 + v1 * v2
        mag1 = np.sqrt(u1**2 + v1**2)
        mag2 = np.sqrt(u2**2 + v2**2)
        cos_angle = dot / (mag1 * mag2 + 1e-6)  # Avoid division by zero
        cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Numerical stability
        angle = np.arccos(cos_angle)
        return np.degrees(angle)
    aae_grid = angular_error(est_x, est_y, gt_x, gt_y)

    # 3.2 Module Error
    module_error_grid = np.abs(np.sqrt(est_x**2 + est_y**2) - np.sqrt(gt_x**2 + gt_y**2))


    # 3.3 RMSE = sqrt( (u_est - u_gt)^2 + (v_est - v_gt)^2 )
    rmse_grid = np.sqrt((est_x - gt_x)**2 + (est_y - gt_y)**2)

    # 3.4 ANSP = Average Normalized Scalar Product
    def normalized_scalar_product(u1, v1, u2, v2):
        dot = u1 * u2 + v1 * v2
        mag1_sq = (u1**2 + v1**2)
        mag2_sq = (u2**2 + v2**2)
        c = 0.01 #1e-6  # Small constant to avoid division by zero
        nsp = (2*dot + c) / (mag1_sq + mag2_sq + c)  # Avoid division by zero
        return nsp
    ansp_grid = normalized_scalar_product(est_x, est_y, gt_x, gt_y)

    # 3.5 NLPD = Negative Log Predictive Density (polars)
    def calc_nlpd_polar(mu_r, mu_theta, sr, stheta, gtr, gttheta):

        # Magnitude term
        term_r = np.log(2 * np.pi * sr**2) + ((gtr - mu_r)**2 / (sr**2))
        
        # Direction term (circular)
        diff_theta = gttheta - mu_theta
        diff_theta = (diff_theta + np.pi) % (2 * np.pi) - np.pi  # Wrap-around to [-π, π]
        term_theta = np.log(2 * np.pi * stheta**2) + ((diff_theta)**2 / (stheta**2))
        
        # Weights for magnitude and direction
        w_mag = 0.3
        w_dir = 0.7
        return 0.5 * (w_mag * term_r + w_dir * term_theta)
    nlpd_polar_grid = calc_nlpd_polar(est_r, est_theta, std_r, std_theta, gt_x, gt_y)

    # 3.6 Cartesian NLPD  = 0.5 * [ ln(2π * σx²) + (x-μx)²/σx² + ln(2π * σy²) + (y-μy)²/σy² ]
    def calc_nlpd_cartesian(mu_x, mu_y, sx, sy, tx, ty):
        term_x = np.log(2 * np.pi * sx**2) + ((tx - mu_x)**2 / (sx**2))
        term_y = np.log(2 * np.pi * sy**2) + ((ty - mu_y)**2 / (sy**2))
        return 0.5 * (term_x + term_y)

    nlpd_cart_grid = calc_nlpd_cartesian(est_x, est_y, std_x, std_y, gt_x, gt_y)

    # 3.7 Combined Uncertainty (Total Standard Deviation)
    uncertainty_grid = np.sqrt(std_x**2 + std_y**2)



    # ====================================
    # 4. Visualization
    # ====================================
    fig, axs = plt.subplots(2, 3, figsize=(18, 8))
    cmap_error = 'YlOrRd' # Yellow to Red for errors
    cmap_uncert = 'Purples'  # Distinct color for uncertainty
    w_scale = np.nanmax(np.sqrt(gt_x**2 + gt_y**2))
    norm_colors = plt.Normalize(vmin=0, vmax=w_scale)

    # Helper to apply black obstacles
    def apply_obstacles(ax):
        ax.imshow(np.where(mask, 0, np.nan), cmap='gray', vmin=0, vmax=1, origin='lower')

    # --- Row 1: ---
    # ------------------------------------------------------------------------------------------
    # [0,0] Ground Truth (Quiver) ---
    apply_obstacles(axs[0,0])
    gt_r = np.sqrt(gt_x**2 + gt_y**2)
    q = axs[0,0].quiver(gt_x, gt_y, gt_r, norm=norm_colors, cmap='viridis', pivot='mid', scale=w_scale, scale_units='xy')
    plt.colorbar(q, label='Wind Speed (m/s)')
    axs[0,0].set_title("Ground Truth Wind")

    # [0,1] GMRF Estimation (Quiver) ---
    apply_obstacles(axs[0,1])
    q = axs[0,1].quiver(est_x, est_y, est_r, norm=norm_colors, cmap='viridis', pivot='mid', scale=w_scale, scale_units='xy')
    plt.colorbar(q, label='Wind Speed (m/s)')
    axs[0,1].set_title("GMRF Wind Estimation")

    # [0,2] Combined Uncertainty Heatmap ---
    apply_obstacles(axs[0,2])
    im_u = axs[0,2].imshow(uncertainty_grid, cmap=cmap_uncert, origin='lower', vmin=0, vmax=3*w_scale)
    plt.colorbar(im_u, ax=axs[0,2], label='Combined Sigma')
    axs[0,2].set_title("Total Uncertainty ($\sqrt{\sigma_x^2 + \sigma_y^2}$)")
    # Add average uncertainty text
    avg_uncertainty = np.nanmean(uncertainty_grid)
    axs[0,2].text(0.5, 0.1, f'Avg Uncertainty: {avg_uncertainty:.2f} m/s', transform=axs[0,2].transAxes, ha='center', va='top')

    # --- Row 2: Performance Metrics ---
    # ------------------------------------------------------------------------------------------
    # [1,0] AAE Heatmap ---
    apply_obstacles(axs[1,0])
    im1 = axs[1,0].imshow(aae_grid, cmap=cmap_error, origin='lower', vmin=0, vmax=180)
    plt.colorbar(im1, ax=axs[1,0], label='AAE (degrees)')
    axs[1,0].set_title("Average Angular Error (AAE)")
    # Add average AAE text
    avg_aae = np.nanmean(aae_grid)
    axs[1,0].text(0.5, 0.1, f'Avg AAE: {avg_aae:.2f}°', transform=axs[1,0].transAxes, ha='center', va='top')

    # [1,1] Module Error Heatmap ---
    apply_obstacles(axs[1,1])
    im2 = axs[1,1].imshow(module_error_grid, cmap=cmap_error, origin='lower')
    plt.colorbar(im2, ax=axs[1,1], label='Module Error (m/s)')
    axs[1,1].set_title("Module Error")
    # Add average Module Error text
    avg_mod_err = np.nanmean(module_error_grid)
    axs[1,1].text(0.5, 0.1, f'Avg Module Error: {avg_mod_err:.2f} m/s', transform=axs[1,1].transAxes, ha='center', va='top')

    # [1,2] Normalized Scalar Product (module and angle)---
    apply_obstacles(axs[1,2])
    im3 = axs[1,2].imshow(ansp_grid, cmap='viridis', origin='lower', vmin=-1, vmax=1)
    plt.colorbar(im3, ax=axs[1,2], label='Normalized Scalar Product')
    axs[1,2].set_title("Normalized Scalar Product (Performance Metric)")
    # Add average ANSP text
    avg_ansp = np.nanmean(ansp_grid)
    axs[1,2].text(0.5, 0.1, f'Avg ANSP: {avg_ansp:.2f}', transform=axs[1,2].transAxes, ha='center', va='top')

    # [1,3] NRMSE Heatmap ---
    #apply_obstacles(axs[1,3])
    #im4 = axs[1,3].imshow(rmse_grid, cmap=cmap_error, origin='lower')
    #plt.colorbar(im4, ax=axs[1,3], label='RMSE (m/s)')
    #axs[1,3].set_title("Root Mean Square Error (RMSE)")
    # Add average RMSE text
    #avg_rmse = np.nanmean(rmse_grid)
    #axs[1,3].text(0.5, 0.1, f'Avg RMSE: {avg_rmse:.2f} m/s', transform=axs[1,3].transAxes, ha='center', va='top')

    # ------------------------------------------------------------------------------------------
    # [1,1] NLPD polar Heatmap ---
    #apply_obstacles(axs[1,1])
    #im4 = axs[1,1].imshow(nlpd_polar_grid, cmap=cmap_error, origin='lower', vmin=0, vmax=5)
    #plt.colorbar(im4, ax=axs[1,1], label='NLPD')
    #axs[1,1].set_title("Weighted Polar NLPD (Uncertainty-Aware Error)")

    # [1,2] Error vs Uncertainty ratio
    # This helps see where the model is "confidently wrong"
    #apply_obstacles(axs[1,2])
    #conf_wrong = rmse_grid / (uncertainty_grid + 1e-6)
    #im_c = axs[1,2].imshow(aae_grid, cmap='coolwarm', origin='lower', vmin=0, vmax=180)
    #plt.colorbar(im_c, ax=axs[1,2], label='AAE / Sigma')
    #axs[1,2].set_title("Calibration (AAE/$\sigma$) ratio")


    # Adjust layout
    for ax in axs.flat:
        ax.set_aspect('equal')
        ax.set_xticks([]) # Optional: hide ticks for cleaner look
        ax.set_yticks([])

    plt.tight_layout()
    plt.savefig('wind_analysis.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()