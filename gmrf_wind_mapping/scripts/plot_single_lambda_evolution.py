# This scripts plots the evolution of optimization metrics (RMSE, NLPD) and lambda parameters
# over multiple repetitions of a single experiment (Nobs fixed), comparing optimal values from Ceres solver against
# baseline fixed values.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os


def plot_optimization_results(df):
    # 1. Identify Baseline vs Optimal rows
    # The first occurrence of a value gets 0, the second gets 1, etc.
    df['occurrence'] = df.groupby('Repetition').cumcount()

    # Separate into two different DataFrames
    df_optimal = df[df['occurrence'] == 0].drop(columns=['occurrence'])
    df_baseline = df[df['occurrence'] == 1].drop(columns=['occurrence'])

    # 2. Create subplots for the metrics
    metrics = ['AAE', 'RMSE', 'NLPD']
    fig, axes = plt.subplots(1, 4, figsize=(22, 5))

    # 3. Plot metrics
    for i, metric in enumerate(metrics):
        # Plot Optimal results
        axes[i].plot(df_optimal['Repetition'], df_optimal[metric], 
                     marker='o', linestyle='-', label='Optimal (Ceres)', color='#1f77b4')
        
        # Plot Baseline results
        axes[i].plot(df_baseline['Repetition'], df_baseline[metric], 
                     marker='x', linestyle='--', label='Baseline (Fixed)', color='#ff7f0e')
        
        # Formatting
        axes[i].set_title(f'{metric} Comparison', fontsize=14, fontweight='bold')
        axes[i].set_xlabel('Repetition ID')
        axes[i].set_ylabel(metric)
        axes[i].legend()
        axes[i].grid(True, linestyle=':', alpha=0.6)

    # 4. Plot Lambda parameters evolution for Optimal only
    lambdas = ['Lambda_Reg', 'Lambda_Flux', 'Lambda_Obstacles', 'Lambda_Observations']
    for l_col in lambdas:
        line = axes[3].plot(df_optimal['Repetition'], df_optimal[l_col], 
                     marker='s', markersize=4, label=f'Opt: {l_col}')
        current_color = line[0].get_color()

        # Plot Baseline results
        axes[3].plot(df_baseline['Repetition'], df_baseline[l_col], 
                     marker='x', linestyle='--', color=current_color, alpha=0.7, label=f'GlobalOpt: {l_col}')
        

    #axes[3].set_yscale('log') # Log scale is vital because values range from 0.0002 to 3600
    axes[3].set_title('$\lambda$ Parameter Values (Log Scale)', fontweight='bold')
    axes[3].set_xlabel('Repetition')
    axes[3].set_ylabel('Value')
    axes[3].grid(True, which="both", ls="-", alpha=0.2)
    axes[3].legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize='small')

    plt.tight_layout()
    plt.savefig('optimization_comparison.png', dpi=300)
    plt.show()



if __name__ == "__main__":
    
    # 1. Load your data
    file_name = "NLPD_optimization_results_central_obstacle_05ms_phase1/Lambda_values_obs_50.csv"
    df = pd.read_csv(file_name)

    # Call the plotting function
    plot_optimization_results(df)