# This scripts plots the evolution of optimization metrics (RMSE, NLPD) and lambda parameters
# over multiple repetitions of a single experiment (Nobs fixed), comparing optimal values from Ceres solver against
# baseline fixed values.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os


def plot_optimization_results(df):
    
    # Data
    metrics = ['AAE', 'AME', 'NRMSE', 'ANSP']
    lambdas = ['Lambda_Advection', 'Lambda_Mass', 'Lambda_Diffusion', 'Lambda_Obstacles']
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 1. Boxplot for the optimization metrics
    boxplot_data = df.melt(id_vars=['Repetition'], value_vars=metrics, var_name='Metric', value_name='Value')
    sns.boxplot(x='Metric', y='Value', data=boxplot_data, ax=axes[0],showfliers=True)
    axes[0].set_title('Optimization Metrics Across Repetitions', fontweight='bold')
    axes[0].set_xlabel('Metric')
    axes[0].set_ylabel('Value')
    axes[0].grid(True, which="both", ls="-", alpha=0.2)

    
    # boxplot for all lambda values
    lambda_data = df.melt(id_vars=['Repetition'], value_vars=lambdas, var_name='Lambda', value_name='Value')
    sns.boxplot(x='Lambda', y='Value', data=lambda_data, ax=axes[1],showfliers=True)
    axes[1].set_title('Lambda Values Across Repetitions', fontweight='bold')
    axes[1].set_xlabel('Lambda Parameter')
    axes[1].set_ylabel('Value')
    axes[1].grid(True, which="both", ls="-", alpha=0.2)

    plt.tight_layout()
    plt.savefig('optimization_comparison.png', dpi=300)
    plt.show()



if __name__ == "__main__":
    
    # 1. Load your data
    file_name = "Lambda_values_6_obs_ME&AE_centralObs_01ms.csv"
    df = pd.read_csv(file_name)

    # Call the plotting function
    plot_optimization_results(df)