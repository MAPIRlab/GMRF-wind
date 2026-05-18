# Plots different maps to visualize GMRF Wind Estimation Results 
#   GT: Quiver plot
#   GMRF estimation: Quiver + Uncertainty heatmap
#   Performance Metrics: AAE (error in angle), Module Error (error in magnitude), NSP (normalized scalar product)
#   NLPD heatmaps: can also include the Negative Log Predictive Density to see uncertainty-aware errors
# ===============================================================================================
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

from pyrsistent import ny

# Configuration
#csv_files = ("gmrf_estimation_1_iters.csv", "gmrf_estimation_2_iters.csv", "gmrf_estimation_10_iters.csv", "gmrf_estimation_100_iters.csv", "gmrf_estimation_1000_iters.csv")
csv_files = ("experiment2/repetition_0/gmrf_estimation_expC_45_obs.csv",)
csv_files = ("gmrf_estimation.csv",)
threshold = 0.001 # Threshold to consider a cell as "no wind" and avoid plotting streamlines there

def cartesian_to_polar_with_uncertainty(df):
    # Prepare new columns
    df['gmrf_wind_r'] = 0.0
    df['gmrf_wind_theta'] = 0.0
    df['gmrf_var_r'] = 0.0
    df['gmrf_var_theta'] = 0.0
    df['gmrf_cov_rtheta'] = 0.0

    for index, row in df.iterrows():
        # Cartesian components and variances
        x = row['gmrf_wind_x']
        y = row['gmrf_wind_y']
        vx = row['gmrf_var_x']   # Variance in x
        vy = row['gmrf_var_y']   # Variance in y
        cov_xy = row['gmrf_cov_xy']  # Covariance between x and y

        r, theta, var_r, var_theta, cov_rtheta = propagate_uncertainty_cartesian_to_polar(x, y, vx, vy, cov_xy)

        # Store results
        df.at[index, 'gmrf_wind_r'] = r
        df.at[index, 'gmrf_wind_theta'] = theta
        df.at[index, 'gmrf_var_r'] = var_r
        df.at[index, 'gmrf_var_theta'] = var_theta
        df.at[index, 'gmrf_cov_rtheta'] = cov_rtheta
    return df

def propagate_uncertainty_cartesian_to_polar(x, y, vx, vy, cov_xy):
    # Convert (x, y) to (r, theta) and propagate uncertainties
    # using first-order error propagation (Jacobian method)
    # and teh full 2x2 covariance matrix

    # 1. Compute r and theta
    r = np.sqrt(pow(x, 2) + pow(y, 2))
    theta = np.arctan2(y, x)

    # Jacobian elements
    dr_dx = x / r if r > 0 else 0
    dr_dy = y / r if r > 0 else 0
    dth_dx = -y / (r*r) if r > 0 else 0
    dth_dy = x / (r*r) if r > 0 else 0

    # Uncertainty propagation: Sigma_polar = J * Sigma_cartesian * J^T
    var_r = dr_dx**2 * vx + dr_dy**2 * vy + 2 * dr_dx * dr_dy * cov_xy
    var_theta = dth_dx**2 * vx + dth_dy**2 * vy + 2 * dth_dx * dth_dy * cov_xy
    cov_rtheta = dr_dx * dth_dx * vx + dr_dy * dth_dy * vy + (dr_dx * dth_dy + dr_dy * dth_dx) * cov_xy

    return r, theta, var_r, var_theta, cov_rtheta



def main():
    for csv_file in csv_files:
        # Load CSV file
        if not os.path.isfile(csv_file):
            raise FileNotFoundError(f"CSV file not found: {csv_file}")

        # 1. Load Metadata (Header)
        # ==========================
        meta = {}
        with open(csv_file, 'r') as f:
            for _ in range(5):
                line = f.readline().strip().split(',')
                # line 5 contains the indices of the observed cells
                if line[0] == "Obs_idx":
                    meta['Obs_idx'] = [int(idx) for idx in line[1:]]
                else:
                    meta[line[0]] = float(line[1])

        nx, ny = int(meta['Dimensions_x']), int(meta['Dimensions_y'])
        
        # Load the actual data (skipping the metadata lines)
        df = pd.read_csv(csv_file, skiprows=5)

        # 2. Convert to polar coordinates (module, direction) and propagate uncertainties
        # ====================================
        df = cartesian_to_polar_with_uncertainty(df)


        # 3. Reshape columns into 2D grids (Y, X) for plotting
        # ====================================
        def to_grid(series):
            return series.values.reshape((ny, nx))
        # Cartesian
        est_x = to_grid(df['gmrf_wind_x'])
        est_y = to_grid(df['gmrf_wind_y'])
        var_x = to_grid(df['gmrf_var_x'])
        var_y = to_grid(df['gmrf_var_y'])
        cov_xy = to_grid(df['gmrf_cov_xy'])
        # Polar
        est_r = to_grid(df['gmrf_wind_r'])
        est_theta = to_grid(df['gmrf_wind_theta'])
        var_r = to_grid(df['gmrf_var_r'])
        var_theta = to_grid(df['gmrf_var_theta'])
        cov_rtheta = to_grid(df['gmrf_cov_rtheta'])
        # Ground Truth
        gt_x = to_grid(df['gt_wind_x'])
        gt_y = to_grid(df['gt_wind_y'])

        # Create Mask for Obstacles (where data is NaN)
        # We'll use this to color cells black
        mask = np.isnan(est_x)


        # ====================================
        # 3. Compute Metrics Cell-by-Cell
        # ====================================

        # 3.1 AE = Angular Error (in degrees)
        def angular_error(u1, v1, u2, v2):
            dot = u1 * u2 + v1 * v2
            mag1 = np.sqrt(u1**2 + v1**2)
            mag2 = np.sqrt(u2**2 + v2**2)
            cos_angle = dot / (mag1 * mag2 + 1e-6)  # Avoid division by zero
            cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Numerical stability
            angle = np.arccos(cos_angle)
            return np.degrees(angle)
        ae_grid = angular_error(est_x, est_y, gt_x, gt_y)

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

        # 3.5 Cartesian NLPD = Negative Log Predictive Density
        def calc_nlpd_cartesian(mu_x, mu_y, vx, vy, cov_xy, tx, ty):
            # Determinant of the covariance matrix
            det = vx * vy - cov_xy**2
            det = np.maximum(det, 1e-12)  # Numerical stability
                
            # Explicit 2x2 inverse formula
            inv_vx = vy / det
            inv_vy = vx / det
            inv_vxy = -cov_xy / det

            inv_cov = np.array([
                [inv_vx, inv_vxy],
                [inv_vxy, inv_vy]
            ])
            
            # Differences (Error vector)
            diffx = tx - mu_x
            diffy = ty - mu_y
            
            # NLPD Calculation
            # 0.5 * (ln|Sigma| + diff.T * inv_Sigma * diff + 2 * ln(2*pi))
            term_log = np.log(det)
            term_mahalanobis = (diffx**2 * inv_vx) + (2 * diffx * diffy * inv_vxy) + (diffy**2 * inv_vy)
            term_const = 2 * np.log(2 * np.pi)
            
            nlpd = 0.5 * (term_log + term_mahalanobis + term_const)
            return nlpd

            # Construct the covariance matrix for the prediction
            #cov_matrix = np.array([[vx, cov_xy], [cov_xy, vy]])
            #inv_cov = np.linalg.inv(cov_matrix + 1e-6 * np.eye(2))  # Add small value for numerical stability
            #diff = np.array([tx - mu_x, ty - mu_y])
            #nlpd = 0.5 * (np.log(np.linalg.det(cov_matrix) + 1e-6) + diff.T @ inv_cov @ diff + 2 * np.log(2 * np.pi))
            #return nlpd
        nlpd_cart_grid = calc_nlpd_cartesian(est_x, est_y, var_x, var_y, cov_xy, gt_x, gt_y)

        # 3.6 Polar NLPD = Negative Log Predictive Density 
        def calc_nlpd_polar(mu_r, mu_theta, vr, vtheta, cov_rtheta, gtr, gttheta):
            # Construct the covariance matrix for the prediction
            cov_matrix = np.array([[vr, cov_rtheta], [cov_rtheta, vtheta]])
            inv_cov = np.linalg.inv(cov_matrix + 1e-6 * np.eye(2))  # Add small value for numerical stability
            diff = np.array([gtr - mu_r, gttheta - mu_theta])
            nlpd = 0.5 * (np.log(np.linalg.det(cov_matrix) + 1e-6) + diff.T @ inv_cov @ diff + 2 * np.log(2 * np.pi))
            return nlpd
        #nlpd_polar_grid = calc_nlpd_polar(est_r, est_theta, var_r, var_theta, cov_rtheta, gt_x, gt_y)

        # 3.7 Uncertainty using the covariance matrix in Cartesian coordinates:    
        max_eig_grid = np.sqrt( 0.5 * (var_x + var_y + np.sqrt((var_x - var_y)**2 + 4 * cov_xy**2)) ) # Largest eigenvalue (m/s)
        det_uncertainty_grid = var_x * var_y - cov_xy**2  # Geometric mean of the eigenvalues
        trace_uncertainty_grid = var_x + var_y  # Sum of variances (upper bound on total uncertainty)


        # ====================================
        # 4. Visualization
        # ====================================
        fig, axs = plt.subplots(3, 3, figsize=(18, 8))
        cmap_error = 'jet' #'YlOrRd' # Yellow to Red for errors
        cmap_uncert = 'Purples'  # Distinct color for uncertainty
        w_scale = np.nanmax(np.sqrt(gt_x**2 + gt_y**2))
        norm_colors = plt.Normalize(vmin=0, vmax=w_scale)
        # Meshgrid for streamlines
        x_coords = np.arange(nx)
        y_coords = np.arange(ny)
        X, Y = np.meshgrid(x_coords, y_coords)

        # Helper to apply black obstacles
        def apply_obstacles(ax):
            ax.imshow(np.where(mask, 0, np.nan), cmap='gray', vmin=0, vmax=1, origin='lower')

        # --- Row 1:  GT 
        # ------------------------------------------------------------------------------------------
        # [0,0] Ground Truth (Quiver) ---
        apply_obstacles(axs[0,0])
        gt_r = np.sqrt(gt_x**2 + gt_y**2)
        q = axs[0,0].quiver(gt_x, gt_y, gt_r, norm=norm_colors, cmap='jet', pivot='mid', scale=w_scale, scale_units='xy')
        plt.colorbar(q, label='Wind Speed (m/s)')
        axs[0,0].set_title("Ground Truth Wind")

        # [0,1] GT Streamlines 
        apply_obstacles(axs[0,1])
        wind_module = np.sqrt(gt_x**2 + gt_y**2)
        # Sustituyes los componentes U y V por NaN donde el viento sea menor al umbral        
        gt_x_clean = np.where(wind_module < threshold, np.nan, gt_x)
        gt_y_clean = np.where(wind_module < threshold, np.nan, gt_y)
        strm = axs[0,1].streamplot(X, Y, gt_x_clean, gt_y_clean, color=wind_module, linewidth=1, density=1.5, cmap='jet', norm=norm_colors)
        plt.colorbar(strm.lines, ax=axs[0,1], label='Wind Speed (m/s)')
        axs[0,1].set_title("Streamlines of Ground Truth Wind")

        # [0,2] Uncertainty Heatmap
        apply_obstacles(axs[0,2])
        im_u = axs[0,2].imshow(max_eig_grid, cmap=cmap_uncert, origin='lower', vmin=0, vmax=w_scale)
        plt.colorbar(im_u, ax=axs[0,2], label='Max Eigenvalue of Covariance (m/s)')
        axs[0,2].set_title("GMRF Estimation Uncertainty (Max Eigenvalue)")
        # Add average uncertainty text
        avg_uncertainty = np.nanmean(max_eig_grid)
        axs[0,2].text(0.5, 0.1, f'Avg Max Eigenvalue: {avg_uncertainty:.2f} m/s', transform=axs[0,2].transAxes, ha='center', va='top')

        # Add observed cells in the plot
        obs_idx = meta['Obs_idx']
        obs_x = [idx % nx for idx in obs_idx]
        obs_y = [idx // nx for idx in obs_idx]
        axs[0,1].scatter(obs_x, obs_y, color='red', marker='x', label='Observations')
        #axs[0,1].legend(loc='upper right')

        # --- Row 2: GMRF Estimation ---
        # ------------------------------------------------------------------------------------------
        # [1,0] GMRF Estimation (Quiver)
        apply_obstacles(axs[1,0])
        q = axs[1,0].quiver(est_x, est_y, est_r, norm=norm_colors, cmap='jet', pivot='mid', scale=w_scale, scale_units='xy')
        plt.colorbar(q, label='Wind Speed (m/s)')
        axs[1,0].set_title("GMRF Wind Estimation")

        # [1,1] GMRF Streamlines 
        apply_obstacles(axs[1,1])
        wind_module = np.sqrt(est_x**2 + est_y**2)
        # Sustituyes los componentes U y V por NaN donde el viento sea menor al umbral
        est_x_clean = np.where(wind_module < threshold, np.nan, est_x)
        est_y_clean = np.where(wind_module < threshold, np.nan, est_y)
        strm = axs[1,1].streamplot(X, Y, est_x_clean, est_y_clean, color=wind_module, linewidth=1, density=1.5, cmap='jet', norm=norm_colors)
        plt.colorbar(strm.lines, ax=axs[1,1], label='Wind Speed (m/s)')
        axs[1,1].set_title("Streamlines of GMRF Wind Estimation")

        # [1,2] NLPD cartesian Heatmap ---
        apply_obstacles(axs[1,2])
        im4 = axs[1,2].imshow(nlpd_cart_grid, cmap=cmap_error, origin='lower', vmin=-5, vmax=5)
        plt.colorbar(im4, ax=axs[1,2], label='NLPD')
        axs[1,2].set_title("Cartesian NLPD (Uncertainty-Aware Error)")
        # Add average NLPD text
        avg_nlpd_cart = np.nanmean(nlpd_cart_grid)
        axs[1,2].text(0.5, 0.1, f'Avg NLPD: {avg_nlpd_cart:.2f}', transform=axs[1,2].transAxes, ha='center', va='top')

        # --- Row 3: Performance Metrics ---
        # ------------------------------------------------------------------------------------------
        # [2,0] Angular Error (AE) Heatmap ---
        apply_obstacles(axs[2,0])
        im1 = axs[2,0].imshow(ae_grid, cmap=cmap_error, origin='lower', vmin=0, vmax=180)
        plt.colorbar(im1, ax=axs[2,0], label='AE (degrees)')
        axs[2,0].set_title("Angular Error (AE)")
        # Add average (AAE) text
        avg_ae = np.nanmean(ae_grid)
        axs[2,0].text(0.5, -0.1, f'AAE: {avg_ae:.2f}°', transform=axs[2,0].transAxes, ha='center', va='top')

        # [2,1] Module Error (ME) Heatmap ---
        apply_obstacles(axs[2,1])
        im2 = axs[2,1].imshow(module_error_grid, cmap=cmap_error, origin='lower')
        plt.colorbar(im2, ax=axs[2,1], label='ME (m/s)')
        axs[2,1].set_title("Module Error (ME)")
        # Add average Module Error text
        avg_mod_err = np.nanmean(module_error_grid)
        axs[2,1].text(0.5, -0.1, f'AME: {avg_mod_err:.2f} m/s', transform=axs[2,1].transAxes, ha='center', va='top')

        # [2,2] Normalized Scalar Product (module and angle)---
        apply_obstacles(axs[2,2])
        im3 = axs[2,2].imshow(ansp_grid, cmap='jet', origin='lower', vmin=-1, vmax=1)
        plt.colorbar(im3, ax=axs[2,2], label='Normalized Scalar Product')
        axs[2,2].set_title("Normalized Scalar Product (Performance Metric)")
        # Add average ANSP text
        avg_ansp = np.nanmean(ansp_grid)
        axs[2,2].text(0.5, 0.1, f'Avg ANSP: {avg_ansp:.2f}', transform=axs[2,2].transAxes, ha='center', va='top')

        

        # Adjust layout
        for ax in axs.flat:
            ax.set_aspect('equal')
            ax.set_xticks([]) # Optional: hide ticks for cleaner look
            ax.set_yticks([])

        plt.tight_layout()
        base_name = os.path.splitext(csv_file)[0]
        out_file = base_name + '.png'
        plt.savefig(out_file, dpi=300)
        plt.show()


        # ====================================
        # 5. PAPER Visualizations
        # ====================================
        fig, axs = plt.subplots(1, 2, figsize=(12, 6))
        cmap_error = 'YlOrRd' # Yellow to Red for errors
        cmap_uncert = 'Purples'  # Distinct color for uncertainty
        w_scale = np.nanmax(np.sqrt(gt_x**2 + gt_y**2))
        norm_colors = plt.Normalize(vmin=0, vmax=w_scale)
        # Meshgrid for streamlines
        x_coords = np.arange(nx)
        y_coords = np.arange(ny)
        X, Y = np.meshgrid(x_coords, y_coords)

        # Helper to apply black obstacles
        def apply_obstacles(ax):
            ax.imshow(np.where(mask, 0, np.nan), cmap='gray', vmin=0, vmax=1, origin='lower')

        plot_gt = True
        if (plot_gt):
            # [0,0] Ground Truth (Quiver) ---
            apply_obstacles(axs[0])
            gt_r = np.sqrt(gt_x**2 + gt_y**2)
            q = axs[0].quiver(gt_x, gt_y, gt_r, norm=norm_colors, cmap='jet', pivot='mid', scale=w_scale, scale_units='xy')
            plt.colorbar(q, label='Wind Speed (m/s)')
            axs[0].set_title("Ground Truth Wind")

            # [0,1] GT Streamlines 
            apply_obstacles(axs[1])
            wind_module = np.sqrt(gt_x**2 + gt_y**2)
            # Sustituyes los componentes U y V por NaN donde el viento sea menor al umbral        
            gt_x_clean = np.where(wind_module < threshold, np.nan, gt_x)
            gt_y_clean = np.where(wind_module < threshold, np.nan, gt_y)
            strm = axs[1].streamplot(X, Y, gt_x_clean, gt_y_clean, color=wind_module, linewidth=1, density=1.5, cmap='jet', norm=norm_colors)
            plt.colorbar(strm.lines, ax=axs[1], label='Wind Speed (m/s)')
            axs[1].set_title("Streamlines of Ground Truth Wind")
        else:
            # [1,0] GMRF Estimation (Quiver)
            apply_obstacles(axs[0])
            q = axs[0].quiver(est_x, est_y, est_r, norm=norm_colors, cmap='jet', pivot='mid', scale=w_scale, scale_units='xy')
            plt.colorbar(q, label='Wind Speed (m/s)')
            axs[0].set_title("GMRF Wind Estimation")

            # [1,1] GMRF Streamlines 
            apply_obstacles(axs[1])
            wind_module = np.sqrt(est_x**2 + est_y**2)
            # Sustituyes los componentes U y V por NaN donde el viento sea menor al umbral
            est_x_clean = np.where(wind_module < threshold, np.nan, est_x)
            est_y_clean = np.where(wind_module < threshold, np.nan, est_y)
            strm = axs[1].streamplot(X, Y, est_x_clean, est_y_clean, color=wind_module, linewidth=1, density=1.5, cmap='jet', norm=norm_colors)
            plt.colorbar(strm.lines, ax=axs[1], label='Wind Speed (m/s)')
            axs[1].set_title("Streamlines of GMRF Wind Estimation")

        # Adjust layout
        for ax in axs.flat:
            ax.set_aspect('equal')
            ax.set_xticks([]) # Optional: hide ticks for cleaner look
            ax.set_yticks([])

        plt.tight_layout()        
        plt.show()


if __name__ == "__main__":
    main()