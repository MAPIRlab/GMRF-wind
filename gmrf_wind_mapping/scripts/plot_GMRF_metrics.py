#!/usr/bin/env python
"""
    Plots the GMRF factor-graph with the relations between the nodes
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

lambda_suf = "flux"
metrics_file = "gmrf_validation_metrics_lambda_" + lambda_suf + ".csv"
# output filename: 'data' + '.png'
base_name = os.path.splitext(metrics_file)[0]
out_file = base_name + '.png'

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

    # Define columns to plot
    x_col = 'GMRF_lambdaPrior_' + lambda_suf
    y_cols = ['NRMSE', 'NCosSim', 'NLL']

    # Ensure required columns exist
    missing_cols = [c for c in ([x_col] + y_cols) if c not in df.columns]
    if missing_cols:
        raise KeyError(f"Missing required columns: {missing_cols}")

    # Convert x and y columns to numeric (fails early if non-numeric)
    df[x_col] = pd.to_numeric(df[x_col], errors='raise')
    for yc in y_cols:
        df[yc] = pd.to_numeric(df[yc], errors='coerce')

    # Drop rows with NaN in y columns, but keep track
    before = len(df)
    df = df.dropna(subset=y_cols)
    after = len(df)
    if after < before:
        print(f"Dropped {before - after} rows with NaN in y-columns")

    # Sort by x for sensible line plots
    df = df.sort_values(by=x_col)

    # Create a figure with subplots (robust for single or many y_cols)
    n = len(y_cols)
    fig, axes = plt.subplots(n, 1, figsize=(8, 3*n), sharex=True)
    if n == 1:
        axes = [axes]

    for i, y_col in enumerate(y_cols):
        ax = axes[i]
        ax.plot(df[x_col].values, df[y_col].values, marker='o', linestyle='-', color='b')
        ax.set_title(f'Plot of {y_col} vs {x_col}', fontsize=12)
        ax.set_ylabel(y_col, fontsize=10)
        ax.grid(True, linestyle='--', alpha=0.6)

    # Limit xticks to avoid overcrowding
    x_unique = np.unique(df[x_col].values)
    if x_unique.size <= 10:
        xticks = x_unique
    else:
        xticks = np.linspace(x_unique.min(), x_unique.max(), 10)
    axes[-1].set_xticks(xticks)
    axes[-1].set_xlabel(x_col, fontsize=12)

    # Auto layout and save
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.suptitle('Metrics as a function of GMRF_lambdaPrior', fontsize=14, y=0.99)
    fig.savefig(out_file, dpi=150)
    print("Plot saved to", out_file)

if __name__ == "__main__":
    main()