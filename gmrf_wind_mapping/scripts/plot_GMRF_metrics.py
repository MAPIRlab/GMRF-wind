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
metrics_file = "NLPD_optimization_results_central_obstacle_1ms_negative/gmrf_opt_metrics_lambda_Nobs.csv"
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

    # Create a figure with two subplots (side by side)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # --- Plot 1: Metrics (AAE, RMSE, NLPD) ---
    ax1.plot(df['N_Obs'], df['AAE'], marker='o', label='AAE', linewidth=2)
    ax1.plot(df['N_Obs'], df['RMSE'], marker='s', label='RMSE', linewidth=2)
    ax1.plot(df['N_Obs'], df['NLPD'], marker='^', label='NLPD', linewidth=2)
    
    ax1.set_title('Evolution of Error Metrics', fontsize=14)
    ax1.set_xlabel('Number of Observations (N_Obs)', fontsize=12)
    ax1.set_ylabel('Metric Value', fontsize=12)
    ax1.grid(True, linestyle='--', alpha=0.7)
    ax1.legend()

    # --- Plot 2: Lambda Parameters ---
    lambdas = ['Lambda_Reg', 'Lambda_Flux', 'Lambda_Obstacles', 'Lambda_Observations']
    markers = ['o', 's', '^', 'd']
    
    for lambda_col, m in zip(lambdas, markers):
        ax2.plot(df['N_Obs'], df[lambda_col], marker=m, label=lambda_col, linewidth=2)

    ax2.set_title('Evolution of Lambda Parameters', fontsize=14)
    ax2.set_xlabel('Number of Observations (N_Obs)', fontsize=12)
    ax2.set_ylabel('Lambda Value', fontsize=12)
    
    # If your Lambdas vary by orders of magnitude, uncomment the line below:
    # ax2.set_yscale('log') 
    
    ax2.grid(True, linestyle='--', alpha=0.7)
    ax2.legend()

    plt.tight_layout()
    plt.savefig('optimization_metrics_Nobs.png', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()