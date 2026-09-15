import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np


def estimate_complexity_order(N_vals_df, T_vals_df):
    # Transform to np arrays
    N_vals = np.array(N_vals_df)
    T_vals = np.array(T_vals_df)

    # Fit a power law T = a * N^b to the data
    log_N = np.log(N_vals)
    log_T = np.log(T_vals)

    # Perform linear regression on the log-log data
    coeffs = np.polyfit(log_N, log_T, 1)
    b = coeffs[0]  # The exponent in the power law

    # Anchor the constant a to the last/largest data point
    N_last = N_vals[-1]
    T_last = T_vals[-1]
    a = T_last / (N_last ** b)

    # Return a, b
    return a, b


# 1. Load the data
column_names = ['Phase', 'NumCells', 'Iteration', 'Time_ms']
df = pd.read_csv('gmrfw_timing_log.csv', names=column_names, header=0)

# Remove NumCells = 10000
df = df[df['NumCells'] != 10000]

# 2. Separate data types
df_total_map = df[df['Phase'] == 'MAP_Total'].sort_values(by='NumCells')
df_total_unc = df[df['Phase'] == 'Uncertainty_Total'].sort_values(by='NumCells')
df_iterations = df[df['Phase'] == 'MAP_Iteration'].sort_values(by='NumCells')

# For each Numcells there are 100 repetitions of the MAP estimation, which includes multiple Picard iterations. 
# We want to analyze both the total time for the MAP estimation (mean +/- std) and the time per Picard iteration (mean +/- std)
# as a function  of the number of cells in the grid.
avg_total_map = df_total_map.groupby('NumCells')['Time_ms'].agg(['mean', 'std']).reset_index()
avg_total_unc = df_total_unc.groupby('NumCells')['Time_ms'].agg(['mean', 'std']).reset_index()
avg_iterations = df_iterations.groupby('NumCells')['Time_ms'].agg(['mean', 'std']).reset_index()

# Estimate the complexity order
a_map, b_map = estimate_complexity_order(avg_total_map['NumCells'], avg_total_map['mean'])
map_order_curve = a_map * avg_total_map['NumCells']**b_map

a_unc, b_unc = estimate_complexity_order(avg_total_unc['NumCells'], avg_total_unc['mean'])
unc_order_curve = a_unc * avg_total_unc['NumCells']**b_unc

a_iter, b_iter = estimate_complexity_order(avg_iterations['NumCells'], avg_iterations['mean'])
iter_order_curve = a_iter * avg_iterations['NumCells']**b_iter


# 3. Create the plots
# ---------------------------
fig, axes = plt.subplots(1, 1, figsize=(15, 6))

# --- MAP vs Number of Cells ---
sns.lineplot(data=avg_total_map, x='NumCells', y='mean', marker='o', ax=axes, label='MAP Time')
#axes.fill_between(avg_total_map['NumCells'], avg_total_map['mean'] - avg_total_map['std'], avg_total_map['mean'] + avg_total_map['std'], alpha=0.2)
# use same color for the order curve
axes.plot(avg_total_map['NumCells'], map_order_curve, label=f'MAP: {a_map:.4f}·O(N^{b_map:.2f})', linestyle='--', color=sns.color_palette()[0])

sns.lineplot(data=avg_total_unc, x='NumCells', y='mean', marker='s', ax=axes, label='Uncertainty Time')
#axes.fill_between(avg_total_unc['NumCells'], avg_total_unc['mean'] - avg_total_unc['std'], avg_total_unc['mean'] + avg_total_unc['std'], alpha=0.2)
axes.plot(avg_total_unc['NumCells'], unc_order_curve, label=f'Uncertainty: {a_unc:.4f}·O(N^{b_unc:.2f})', linestyle='--', color=sns.color_palette()[1])

# --- Average Time per Picard Iteration vs Number of Cells ---
sns.lineplot(data=avg_iterations, x='NumCells', y='mean', marker='^', ax=axes, label='Single Picard Iteration')
#axes.fill_between(avg_iterations['NumCells'], avg_iterations['mean'] - avg_iterations['std'], avg_iterations['mean'] + avg_iterations['std'], alpha=0.2)
axes.plot(avg_iterations['NumCells'], iter_order_curve, label=f'Iteration: {a_iter:.4f}·O(N^{b_iter:.2f})', linestyle='--', color=sns.color_palette()[2])

# Add reference N^(1.5) curve for comparison
#N_values = avg_iterations['NumCells']
#reference_curve = 0.0225 * N_values**1.5  # Adjust the scaling
#axes.plot(N_values, reference_curve, label='Reference N^(1.5)', linestyle='--')

axes.set_title('Computational Time vs Number of Grid Cells')
axes.set_xlabel('Number of Cells (N)')
axes.set_ylabel('Time (ms)')
axes.grid(True)
axes.legend()

plt.tight_layout()
plt.show()



