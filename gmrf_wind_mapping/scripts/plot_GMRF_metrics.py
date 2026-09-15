#!/usr/bin/env python
# -------------------------------------------------------------------------------
# Plots the evolution of AAE and AME with increasing number of observations (N_Obs) using the results from the GMRF optimization experiments.
# -------------------------------------------------------------------------------
import os
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

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
    theta = np.atan2(y, x)

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
    # load results from CSV files
    results_dir = "experiment2"
    num_repetitions = 100
    num_observations_list = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 100, 125, 150, 175, 200]
    metrics_data = {"N_Obs": [], "AAE_mean": [], "AAE_std": [], "AME_mean": [], "AME_std": []}

    for num_observations in num_observations_list:
        aae_values = []
        ame_values = []

        for repeat in range(num_repetitions):
            
            # Load CSV file
            filename_csv = os.path.join(results_dir, f"repetition_{repeat}", f"gmrf_estimation_expC_{num_observations}_obs.csv")
            if not os.path.isfile(filename_csv):
                raise FileNotFoundError(f"CSV file not found: {filename_csv}")
            
            # ==========================
            # 1. Load Metadata (Header)
            # ==========================
            meta = {}
            with open(filename_csv, 'r') as f:
                for _ in range(5):
                    line = f.readline().strip().split(',')
                    # line 5 contains the indices of the observed cells
                    if line[0] == "Obs_idx":
                        meta['Obs_idx'] = [int(idx) for idx in line[1:]]
                    else:
                        meta[line[0]] = float(line[1])

            nx, ny = int(meta['Dimensions_x']), int(meta['Dimensions_y'])
            
            # Load the actual data (skipping the metadata lines)
            df = pd.read_csv(filename_csv, skiprows=5)

            # 2. Convert to polar coordinates (module, direction) and propagate uncertainties
            df = cartesian_to_polar_with_uncertainty(df)

            # 3. Reshape columns into 2D grids (Y, X) for plotting
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
            # 2. Compute Metrics Cell-by-Cell
            # ====================================

            # 2.1 AE = Angular Error (in degrees)
            def angular_error(u1, v1, u2, v2):
                dot = u1 * u2 + v1 * v2
                mag1 = np.sqrt(u1**2 + v1**2)
                mag2 = np.sqrt(u2**2 + v2**2)
                cos_angle = dot / (mag1 * mag2 + 1e-6)  # Avoid division by zero
                cos_angle = np.clip(cos_angle, -1.0, 1.0)  # Numerical stability
                angle = np.arccos(cos_angle)
                return np.degrees(angle)
            ae_grid = angular_error(est_x, est_y, gt_x, gt_y)

            # Average AE across all non-obstacle cells
            aae = np.mean(ae_grid[~mask])
            aae_values.append(aae)


            # 2.2 ME = Module Error (in m/s)
            me_grid = np.abs(est_r - np.sqrt(gt_x**2 + gt_y**2))
            ame = np.mean(me_grid[~mask])
            ame_values.append(ame)

        # Keep data for metrics evolution plot
        metrics_data["N_Obs"].append(num_observations)
        metrics_data["AAE_mean"].append(np.mean(aae_values))
        metrics_data["AAE_std"].append(np.std(aae_values))
        metrics_data["AME_mean"].append(np.mean(ame_values))
        metrics_data["AME_std"].append(np.std(ame_values))

    # ====================================
    # 3. Plot Metrics Evolution with N_Obs (AAE and AME with error bars on different y-axis)
    # ====================================
    fig, ax1 = plt.subplots(figsize=(12, 6))
    # --- Graficar AAE en el primer eje (Izquierdo) ---
    color_aae = 'tab:blue'
    ax1.set_xlabel("Number of Observations (N_Obs)")
    ax1.set_ylabel("AAE (degrees)", color=color_aae)
    line1 = ax1.errorbar(metrics_data["N_Obs"], metrics_data["AAE_mean"], 
                        yerr=metrics_data["AAE_std"], label="AAE (degrees)", 
                        marker='o', color=color_aae)
    ax1.tick_params(axis='y', labelcolor=color_aae)
    ax1.grid(True, linestyle='--', alpha=0.7)

    # --- Crear el segundo eje (Derecho) ---
    ax2 = ax1.twinx() 

    # --- Graficar AME en el segundo eje (Derecho) ---
    color_ame = 'tab:red'
    ax2.set_ylabel("AME (m/s)", color=color_ame)
    line2 = ax2.errorbar(metrics_data["N_Obs"], metrics_data["AME_mean"], 
                        yerr=metrics_data["AME_std"], label="AME (m/s)", 
                        marker='o', color=color_ame)
    ax2.tick_params(axis='y', labelcolor=color_ame)

    # --- Ajustes finales ---
    plt.title("Evolution of AAE and AME with Increasing Number of Observations")

    # Unificar las leyendas de ambos ejes en una sola
    lines = [line1, line2]
    labels = [l.get_label() for l in lines]
    ax1.legend(lines, labels, loc='upper right')
    plt.savefig("metrics_evolution.png")
    plt.show()

if __name__ == "__main__":
    main()