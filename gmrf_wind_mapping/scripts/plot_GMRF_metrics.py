#!/usr/bin/env python
import os
import pandas as pd
import numpy as np
import matplotlib
# Switch to the 'Agg' backend. This is CRITICAL for remote servers without X11.
# It allows Matplotlib to render to files (PNG, PDF, etc.) without needing a display.
# This call must be made BEFORE importing matplotlib.pyplot.
matplotlib.use('Agg')
import matplotlib.pyplot as plt


# Configuration
metrics_file = "NLPD_optimization_results_central_obstacle_05ms/gmrf_opt_metrics_lambda_Nobs.csv"
base_name = os.path.splitext(metrics_file)[0]
out_file = base_name + '.png'
output_dir = "metrics_plots"

def main():
    if not os.path.isfile(metrics_file):
        raise FileNotFoundError(f"Metrics file not found: {metrics_file}")

    # Attempt to load CSV
    try:
        df = pd.read_csv(metrics_file)
    except Exception as e:
        raise RuntimeError(f"Failed to read CSV: {e}")

    print("CSV loaded. Columns:", df.columns.tolist())
    if df.empty:
        raise RuntimeError("CSV is empty")

    # Output
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

   
    # Sort by N_Obs to ensure lines are drawn correctly
    df = df.sort_values('N_Obs')

    # --- Calculation: Global Mean of Lambdas ---
    # We select the columns that end in '_mean'
    mean_cols = [c for c in df.columns if c.endswith('_mean')]
    global_averages = df[mean_cols].mean()

    print("-" * 30)
    print("GLOBAL LAMBDA AVERAGES (Across all N_Obs):")
    for col, val in global_averages.items():
        print(f"{col.replace('_mean', ''):<20}: {val:.6f}")
    print("-" * 30)


    # Create a figure with two subplots (side by side)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 7))

    # --- Plot 1: Metrics (AAE, RMSE, NLPD) ---
    ax1.plot(df['N_Obs'], df['AAE'], marker='o', label='AAE', linewidth=2)
    ax1.plot(df['N_Obs'], df['RMSE'], marker='s', label='RMSE', linewidth=2)
    ax1.plot(df['N_Obs'], df['NLPD'], marker='^', label='NLPD', linewidth=2)
    
    ax1.set_title('Evolution of Error Metrics', fontsize=14)
    ax1.set_xlabel('Number of Observations (N_Obs)', fontsize=12)
    ax1.set_ylabel('Metric Value', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()

    # --- Plot 2: Lambda Parameters with Uncertainty (Std) ---
    lambda_pairs = [        
        ('Lambda_Flux_mean', 'Lambda_Flux_std'),
        ('Lambda_Obstacles_mean', 'Lambda_Obstacles_std'),
        ('Lambda_Reg_mean', 'Lambda_Reg_std')
    ]
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    markers = ['o', 's', '^', 'd']
    
    # Track min/max of MEANS ONLY to set the Y-limits
    all_means = df[mean_cols]
    y_min = all_means.min().min()
    y_max = all_means.max().max()

    # Add 10% padding to the range for better visibility
    padding = (y_max - y_min) * 0.1
    if padding == 0: padding = 1.0 # Avoid error if all values are identical
    ax2.set_ylim(y_min - padding, y_max + padding)

    for (m_col, s_col), color, m in zip(lambda_pairs, colors, markers):
        # Plot the mean line
        ax2.plot(df['N_Obs'], df[m_col], marker=m, label=m_col.replace('_mean', ''), 
                 linewidth=2, color=color)
        
        # Plot the uncertainty area (Mean +/- Std)
        # alpha controls transparency of the shaded region
        ax2.fill_between(df['N_Obs'], 
                         df[m_col] - df[s_col], 
                         df[m_col] + df[s_col], 
                         color=color, alpha=0.2)

    ax2.set_title('Lambda Parameters with Uncertainty ($\sigma$)', fontsize=14, fontweight='bold')
    ax2.set_xlabel('Number of Observations (N_Obs)', fontsize=12)
    ax2.set_ylabel('Lambda Value', fontsize=12)
    
    # Given the high variance in your example (std ~ 791), log scale is usually better
    #ax2.set_yscale('log') 
    
    ax2.grid(True, which="both", linestyle='--', alpha=0.4)
    ax2.legend(loc='best', fontsize='small')

    plt.tight_layout()
    plt.savefig(out_file, dpi=300)
    print(f"Plot saved to: {out_file}")
    plt.show()


if __name__ == "__main__":
    main()