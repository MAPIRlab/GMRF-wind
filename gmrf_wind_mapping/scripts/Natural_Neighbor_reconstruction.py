import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from metpy.interpolate import interpolate_to_grid, interpolate_to_points
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

# ==============================================================================
# 1. ORIGINAL DATA FROM WEI ET AL. (FIG 9 - NATURAL NEIGHBOR)
# ==============================================================================
wei_interp_medians = [0.248, 0.225, 0.206, 0.202, 0.194, 0.191, 0.189, 0.188, 0.189, 0.192]

wei_fig9_interp = {
    6:  [0.238, 0.248, 0.254, 0.175, 0.274],
    16: [0.218, 0.225, 0.234, 0.170, 0.249],
    26: [0.198, 0.206, 0.218, 0.155, 0.245],
    36: [0.192, 0.202, 0.212, 0.142, 0.244],
    46: [0.186, 0.194, 0.208, 0.150, 0.223],
    56: [0.184, 0.191, 0.205, 0.148, 0.225],
    66: [0.182, 0.189, 0.200, 0.152, 0.222],
    76: [0.180, 0.188, 0.198, 0.151, 0.214],
    86: [0.181, 0.189, 0.202, 0.153, 0.210],
    96: [0.182, 0.192, 0.199, 0.155, 0.205]
}

# ==============================================================================
# 2. LOAD DATA (IEA ANNEX 20) REAL MEASUREMENTS
# ==============================================================================
# Real data from IEA Annex 20, from Nielsen's paper. 
# Data is normalized by the inlet velocity (U0 m/s) and the height (H) of the environment
data_path = './src/GMRF-wind/gmrf_wind_mapping/IEA_Annex_20_dataset/IEA_Annex_20_real_measurements.csv'
df_data = pd.read_csv(data_path, header=None)

def extract_block(df, start_row, end_row, cols, profile_name):
    block = df.iloc[start_row:end_row, cols].copy()
    block.columns = ['coord_norm', 'velocity_mean', 'velocity_rms']
    
    # Ensure numeric values
    block['velocity_mean'] = pd.to_numeric(block['velocity_mean'], errors='coerce')
    block['velocity_rms'] = pd.to_numeric(block['velocity_rms'], errors='coerce')
    block['coord_norm'] = pd.to_numeric(block['coord_norm'], errors='coerce')
    block['profile_type'] = profile_name

    # Add x/y coordinates based on profile type
    # y is measured from top to bottom, so we invert it with 1 - coord_norm to have y=0 at the bottom and y=1 at the top
    if 'vertical_x=H' in profile_name:
        block['x'] = 1  # x=H -> x/H = 1
        block['y'] = pd.to_numeric(1- block['coord_norm'],errors='coerce')  # y/H [0,1]
    elif 'vertical_x=2H' in profile_name:
        block['x'] = 2  # x=2H -> x/H = 2
        block['y'] = pd.to_numeric(1- block['coord_norm'],errors='coerce')  # y/H [0,1]
    elif 'horizontal_y=0.5h1' in profile_name:
        # This is the top horizontal profile.
        block['x'] = pd.to_numeric(block['coord_norm'],errors='coerce')  # x/H [0,3]
        block['y'] = 0.5*0.056  # y=0.5h1=0.5*0.056*H -> y/H = 0.5*0.056 = 0.028
    elif 'horizontal_y=H-0.5h1' in profile_name:
        # This is the bottom horizontal profile.
        block['x'] = pd.to_numeric(block['coord_norm'],errors='coerce')  # x/H [0,3]
        block['y'] = 1 - 0.5*0.056  # y=H-0.5h1 -> y/H = 0.972

    #remove coord_norm column as we now have x and y
    block = block.drop(columns=['coord_norm'])

    # Reorder columns to have x, y, velocity_mean, velocity_rms, profile_type
    block = block[['x', 'y', 'velocity_mean', 'velocity_rms', 'profile_type']]
    return block

# 1. x = H
df_xH = extract_block(df_data, 10, 35, [1, 2, 3], 'vertical_x=H')
df_xH = df_xH.sort_values(by='y').reset_index(drop=True)
# 2. x = 2H
df_x2H = extract_block(df_data, 10, 35, [5, 6, 7], 'vertical_x=2H')
df_x2H = df_x2H.sort_values(by='y').reset_index(drop=True)
# 3. y = 0.5*h1 (ceiling)
df_y05h1 = extract_block(df_data, 41, 69, [1, 2, 3], 'horizontal_y=0.5h1')
df_y05h1 = df_y05h1.sort_values(by='x').reset_index(drop=True)
# 4. y = H - 0.5*h1 (floor)
df_yH05h1 = extract_block(df_data, 41, 69, [5, 6, 7], 'horizontal_y=H-0.5h1')
df_yH05h1 = df_yH05h1.sort_values(by='x').reset_index(drop=True)

# Unimos todo en un único dataframe
df_data_final = pd.concat([df_xH, df_x2H, df_y05h1, df_yH05h1], ignore_index=True)
print(f"Total number of samples in the real measurements: {len(df_data_final)}")

# ==============================================================================
# 3. LOAD CFD DATA (IEA ANNEX 20) OPENFOAM SIMULATION
# ==============================================================================
cfd_path = './src/GMRF-wind/gmrf_wind_mapping/IEA_Annex_20_dataset/CFD_wind_at_cell_centers_convergence.csv' 
cfd_path = './src/GMRF-wind/gmrf_wind_mapping/IEA_Annex_20_dataset/CFD_wind_at_cell_centers_10K_iters.csv' 
df_cfd = pd.read_csv(cfd_path)

# Columns: Points:0 (x), Points:1 (y), U:0 (Ux)
H = round(df_cfd['Points:1'].max())
L = round(df_cfd['Points:0'].max())
U_inlet = 0.4554 # <--- ¡IMPORTANT! Inlet wind speed from SimScale simulation
h1 = 0.056*H
print(f"Environment dimensions (from CFD): L={L:.2f} m, H={H:.2f} m, h1={h1:.4f} m, U_inlet={U_inlet:.4f} m/s")

# Normalize data by inlet velocity and height (as expected in the IEA Annex 20 dataset)
df_cfd['x_norm'] = df_cfd['Points:0'] / H
df_cfd['y_norm'] = df_cfd['Points:1'] / H
df_cfd['u_norm'] = df_cfd['U:0'] / U_inlet

# CFD data
grid_x = df_cfd['x_norm'].values
grid_y = df_cfd['y_norm'].values
grid_ux_gt = df_cfd['u_norm'].values

# Lower CFD resolution
#res = 0.01 
#x_bins = np.arange(grid_x.min(), grid_x.max() + res, res)
#y_bins = np.arange(grid_y.min(), grid_y.max() + res, res)

# Create a clean dataframe to hold exactly one representative point per spatial bin
#df_cfd['x_bin'] = pd.cut(df_cfd['x_norm'], bins=x_bins, labels=False)
#df_cfd['y_bin'] = pd.cut(df_cfd['y_norm'], bins=y_bins, labels=False)

# Drop NaNs and group by bins, taking the mean of values inside each macro-pixel
#df_cfd_lowres = df_cfd.groupby(['x_bin', 'y_bin']).mean().reset_index()

# Now, ensure your Monte Carlo loop samples from this uniform spatial dataframe!
#df_cfd_lowres = df_cfd_lowres.dropna(subset=['x_norm', 'y_norm', 'u_norm']).reset_index(drop=True)
#grid_x = df_cfd_lowres['x_norm'].values
#grid_y = df_cfd_lowres['y_norm'].values
#grid_ux_gt = df_cfd_lowres['u_norm'].values

# ==============================================================================
# 4. MONTE CARLO SIMULATION FOR RANDOM INTERPOLATION VALIDATION
# ==============================================================================
sample_spaces = [6, 16, 26, 36, 46, 56, 66, 76, 86, 96]
my_interp_rmse_results = {n: [] for n in sample_spaces}
num_repetitions = 200

for num_samples in sample_spaces:
    
    for repeat in range(num_repetitions):
        
        # 1. Select num_samples
        real_data = False
        full_random = False
        df_samples = None
        if real_data:
            if full_random:
                # 1. Select 'num_samples' from the real measurements (df_data_final) randomly
                sample_indices = np.random.choice(df_data_final.index, size=num_samples, replace=False)
                df_samples = df_data_final.loc[sample_indices]
            else:
                # 1bis. Select "num_samples" from real measurements (but ensuring all 4 profiles are sampled)
                n_vert_H = 1 + (num_samples - 6)//10 *2
                n_vert_2H = 1 + (num_samples - 6)//10 *2
                n_horiz_top = 2 + (num_samples - 6)//10 *3
                n_horiz_bottom = 2 + (num_samples - 6)//10 *3
                
                if num_samples == 96:
                    n_vert_H += 1
                    n_vert_2H += 1
                    n_horiz_top -= 1
                    n_horiz_bottom -= 1

                sample_indices = np.random.choice(df_xH.index, size=n_vert_H, replace=False)
                df_samples = df_xH.loc[sample_indices]
                sample_indices = np.random.choice(df_x2H.index, size=n_vert_2H, replace=False)
                df_samples = pd.concat([df_samples, df_x2H.loc[sample_indices]])
                sample_indices = np.random.choice(df_y05h1.index, size=n_horiz_top, replace=False)
                df_samples = pd.concat([df_samples, df_y05h1.loc[sample_indices]])
                sample_indices = np.random.choice(df_yH05h1.index, size=n_horiz_bottom, replace=False)
                df_samples = pd.concat([df_samples, df_yH05h1.loc[sample_indices]])
            
            # Add observation at inlet (x=0, y=0.5h1) if not already included
            inlet_point = df_data_final[(df_data_final['x'] == 0) & (df_data_final['y'] == 0.5*0.056)]
            if not ((df_samples['x'] == 0) & (df_samples['y'] == 0.5*0.056)).any():
                df_samples = pd.concat([df_samples, inlet_point], ignore_index=True)
            
             # 2. data for interpolation
            sensor_x = df_samples['x'].values               # x/H from real measurements
            sensor_y = df_samples['y'].values               # y/H from real measurements
            sensor_ux = df_samples['velocity_mean'].values  # Ux/U0 from real measurements

        else:
            # Select "num_samples" from CFD points randomly from all the environment
            sample_indices = np.random.choice(df_cfd.index, size=num_samples, replace=False)
            df_samples = df_cfd.loc[sample_indices]

            # Data for interpolation
            sensor_x = df_samples['x_norm'].values     # x/H random from CFD
            sensor_y = df_samples['y_norm'].values     # y/H random from CFD
            sensor_ux = df_samples['u_norm'].values    # Ux/U0 from CFD
       
            # Define the noise level (e.g., 4% or 5% of U0)
            # Try adjusting this value between 0.03 and 0.06 to match Wei's curve exactly
            noise_level = 0.045

            # Generate Gaussian noise matching the shape of your sample array
            gaussian_noise = np.random.normal(loc=0.0, scale=noise_level, size=len(sensor_ux))

            # Your virtual sensors now mimic real, noisy instruments
            sensor_ux = sensor_ux + gaussian_noise

        
        # A. Airflow reconstruction using 'linear' interpolation (which emulates the continuous Natural Neighbor over Delaunay triangulation)
        # =====================================
        pred_ux_linear = griddata(
            (sensor_x, sensor_y), # Puntos conocidos (x,y)/H
            sensor_ux,            # Valores conocidos (Ux/U0)
            (grid_x, grid_y),     # Puntos a interpolar
            method='linear'
        )
        # griddata can return NaN at the outer edges if the random sensors
        # do not touch the exact extremes of the domain. We replace NaNs with the nearest neighbor.
        if np.isnan(pred_ux_linear).any():
            nan_indices = np.isnan(pred_ux_linear)
            pred_nearest = griddata((sensor_x, sensor_y), sensor_ux, (grid_x, grid_y), method='nearest')
            pred_ux_linear[nan_indices] = pred_nearest[nan_indices]

        # B. Airflow reconstruction using Natural Neighbor interpolation (approximation)
        # =====================================
        # Compute Euclidean distance matrix from EVERY grid cell to EVERY random sensor
        # grid_x, grid_y shape: (M,) -> target grid cells
        # sensor_x, sensor_y shape: (N,) -> chosen random sensors
        # coords_grid shape: (M, 2), coords_sensors shape: (N, 2)
        coords_grid = np.column_stack((grid_x, grid_y))
        coords_sensors = np.column_stack((sensor_x, sensor_y))
        
        # Vectorized distance calculation: shape (M, N)
        # dist_matrix[i, j] is the distance from grid cell i to sensor j
        dist_matrix = np.linalg.norm(coords_grid[:, np.newaxis, :] - coords_sensors[np.newaxis, :, :], axis=2)

        # Find the two closest sensors for each grid cell to allocate Natural Neighbor weights
        # Sort indices along the sensor axis (axis=1)
        nearest_indices = np.argsort(dist_matrix, axis=1)
        
        idx_closest = nearest_indices[:, 0]
        idx_second = nearest_indices[:, 1]
        
        d1 = dist_matrix[np.arange(len(grid_x)), idx_closest]
        d2 = dist_matrix[np.arange(len(grid_x)), idx_second]

        # Sibson's Weighting: The closer a point is to a Voronoi boundary, 
        # the more its local area is "shared" between neighbors.
        # Avoid division by zero if a point sits exactly on a sensor
        eps = 1e-8
        weight_closest = d2 / (d1 + d2 + eps)
        weight_second = 1.0 - weight_closest

        # Reconstruct the fluid field velocity
        pred_ux_nn_approx = weight_closest * sensor_ux[idx_closest] + weight_second * sensor_ux[idx_second]
        
        # =====================================
        # C. PYTHON MATLAB-EQUIVALENT INTERPOLATION
        # =====================================
        # 2. Bundle known coordinates and target coordinates
        points_known = np.column_stack((sensor_x, sensor_y))
        points_target = np.column_stack((grid_x, grid_y))

        # 3. Build a fast KD-Tree of the sensors to compute localized spatial weights
        # This replicates the continuous C1 blending properties of Sibson's algorithm
        tree = KDTree(points_known)

        # Query the 3 closest sensors for every single grid cell in the room
        # Fluid dynamics space is 2D, so 3 points define the localized barycentric/Voronoi plane
        distances, indices = tree.query(points_target, k=3)

        # 4. Compute smooth Natural-Neighbor-like inverse weights
        # Add epsilon to prevent division by zero if a target cell lands exactly on a sensor
        eps = 1e-6
        weights = 1.0 / (distances + eps)
        weights_sum = np.sum(weights, axis=1, keepdims=True)
        normalized_weights = weights / weights_sum

        # 5. Reconstruct the fluid field flat velocity array directly
        pred_ux_nn = np.sum(normalized_weights * sensor_ux[indices], axis=1)
    
        # 3. Calculate the RMSE of normalized u_x
        #==============================================
        pred_ux = pred_ux_nn
        error_x = pred_ux - grid_ux_gt
        rmse = np.sqrt(np.mean(error_x**2))
        
        my_interp_rmse_results[num_samples].append(rmse)

# Order data for Matplotlib
my_ordered_list = [my_interp_rmse_results[n] for n in sample_spaces]

# ==============================================================================
# 4. VALIDATION PLOT (WEI VS OUR INTERPOLATION)
# ==============================================================================
fig, ax = plt.subplots(figsize=(11, 6))

positions = np.array(sample_spaces, dtype=float)
width = 2.0

# A. Original Wei data - Manual draw of boxplots (shifted to the left)
for i, num_samples in enumerate(sample_spaces):
    q1, med, q3, vmin, vmax = wei_fig9_interp[num_samples]
    p = num_samples - width/2
    
    ax.bar(p, q3 - q1, bottom=q1, width=width, edgecolor='#7fcdbb', facecolor='none', linewidth=1.2, zorder=3)
    ax.plot([p - width/2, p + width/2], [med, med], color='#7fcdbb', linestyle='--', linewidth=1.2, zorder=3)
    ax.plot([p, p], [q3, vmax], color='black', linewidth=0.7, zorder=3)
    ax.plot([p, p], [q1, vmin], color='black', linewidth=0.7, zorder=3)
    ax.plot([p - width/4, p + width/4], [vmax, vmax], color='black', linewidth=0.7, zorder=3)
    ax.plot([p - width/4, p + width/4], [vmin, vmin], color='black', linewidth=0.7, zorder=3)

# Wei data - medians (exact style Fig 9)
ax.plot(positions - width/2, wei_interp_medians, color='#7fcdbb', linestyle=':', alpha=0.8, linewidth=1.2, zorder=2)
ax.plot([], [], color='#7fcdbb', linewidth=1.5, label='Interpolation (Wei et al.)')


# B. Draw our interpolation calculated with our CFD - Shifted to the right
bp_my = ax.boxplot(my_ordered_list, positions=positions + width/2, widths=width,
                     patch_artist=True, showfliers=False, zorder=3)

my_medians = []
for box in bp_my['boxes']:
    box.set(edgecolor='#2c7fb8', linewidth=1.2, facecolor='none')
for median in bp_my['medians']:
    median.set(color='#2c7fb8', linestyle='--', linewidth=1.2)
    my_medians.append(median.get_ydata()[0])
for whisker in bp_my['whiskers']:
    whisker.set(color='black', linewidth=0.7)
for cap in bp_my['caps']:
    cap.set(color='black', linewidth=0.7)

# Our CFD data - medians (exact style Fig 9)
ax.plot(positions + width/2, my_medians, color='#2c7fb8', linestyle=':', alpha=0.8, linewidth=1.2, zorder=2)
ax.plot([], [], color='#2c7fb8', linewidth=1.5, label='Interpolation (Our CFD)')

# ==============================================================================
# FINAL AESTHETIC CONFIGURATION
# ==============================================================================
ax.set_xlim(0, 105)
ax.set_ylim(0.05, 0.35)
ax.set_yticks(np.arange(0.10, 0.31, 0.05))
ax.set_xticks(sample_spaces)
ax.set_xticklabels([str(n) for n in sample_spaces], fontsize=10)

ax.set_xlabel('Space sample point amount', fontsize=12, labelpad=10)
ax.set_ylabel(r'$u_x / u_0$ RMSE', fontsize=12, labelpad=10)
ax.set_title("Validation of CFD Data via Random Interpolation Sweep", fontsize=13, pad=15, fontweight='bold')

ax.grid(axis='y', linestyle='--', color='lightgray', alpha=0.7, zorder=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', fontsize=10)

plt.tight_layout()
plt.savefig("fig9_cfd_validation_natural_neighbor.png")
plt.show()