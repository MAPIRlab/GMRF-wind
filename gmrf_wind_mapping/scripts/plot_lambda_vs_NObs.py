# -------------------------------------------------------------------------------
# Script to visualize the evolution and distribution of lambda parameters
# across multiple scenarios and varying the number of observations
# -------------------------------------------------------------------------------
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import re

# ---------------------------------------
# 1. CONFIGURATION (fodlers to consider)
# ---------------------------------------
folders = [#'NLPD_optimization_results_central_obstacle_01ms', 
           #'NLPD_optimization_results_central_obstacle_05ms', 
           'NLPD_optimization_results_central_obstacle_1ms',
           #'NLPD_optimization_results_snake_01ms',
           #'NLPD_optimization_results_snake_05ms',
           #'NLPD_optimization_results_snake_1ms',
           ]
lambda_cols = ['Lambda_Mass', 'Lambda_Vorticity', 'Lambda_Obstacles', 'Lambda_Reg']
cleaned_data_list = []


def extract_Number(filename):
    """Extrae el primer número encontrado en el nombre del archivo."""
    numbers = re.findall(r'\d+', filename)
    return int(numbers[0]) if numbers else filename


# 2. Process each folder/scenario
# ----------------------------------
for folder in folders:
    # Load all CSVs in the folder
    files = glob.glob(os.path.join(folder, "Lambda_values_obs_*.csv"))
    folder_dfs = []    

    for file in files:
        df = pd.read_csv(file)
        
        # Use only the first line of each Repetition (the second is a baseline)
        experiment_results = df.groupby('Repetition').nth(0).reset_index()
        
        # Add metadata so we can distinguish between folders in the plot
        experiment_results['Source_Folder'] = folder
        experiment_results['File_Name'] = os.path.basename(file)
        folder_dfs.append(experiment_results)             

    if not folder_dfs:
        continue
        
    folder_combined = pd.concat(folder_dfs, ignore_index=True)          # Combine all experiments (Num Obs) in this folder
    
    # 4. Per-Folder Outlier Removal
    # --------------------------------------------------
    # Remove outliers based on IQR for each lambda column
    mask = pd.Series(True, index=folder_combined.index)
    for col in lambda_cols:
        Q1 = folder_combined[col].quantile(0.25)
        Q3 = folder_combined[col].quantile(0.75)
        IQR = Q3 - Q1
        lower = Q1 - 1.5 * IQR
        upper = Q3 + 1.5 * IQR
        mask &= (folder_combined[col] >= lower) & (folder_combined[col] <= upper)
    
    df_folder_clean = folder_combined[mask]
    cleaned_data_list.append(df_folder_clean)
    print(f"Folder '{folder}': Kept {len(df_folder_clean)}/{len(folder_combined)} rows.")

# ALL FOLDERS PROCESSED
# ----------------------------------

# 5. Combine everything into one big DataFrame
final_df_clean = pd.concat(cleaned_data_list, ignore_index=True)

# 6. Calculate the "General Hyperparameters" (Centroid)
general_mean = final_df_clean[lambda_cols].mean()
general_median = final_df_clean[lambda_cols].median()

print("\n--- SUGGESTED GLOBAL VALUES (Cleaned Data) ---")
summary_table = pd.DataFrame({'Mean': general_mean, 'Median': general_median})
print(summary_table)


# 7. Visualization: Pair Plot
plt.figure(figsize=(12, 10))
sns.set_theme(style="whitegrid")

# PLOT1: Boxplot without outliers
df_melted = final_df_clean.melt(id_vars=['Source_Folder'], value_vars=lambda_cols, 
                                var_name='Parameter', value_name='Value')

# 'showfliers=False' is an extra safety, though we already removed them manually
sns.boxplot(data=df_melted, x='Parameter', y='Value', hue='Source_Folder', showfliers=False)
plt.title("Distribution of Lambda Parameters (Outliers Removed Per Folder)")
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig('lambda_boxplot_all_folders.png')
plt.show()


# Plot 2: Pairplot with Global Median marked
g = sns.pairplot(final_df_clean, vars=lambda_cols, hue='Source_Folder', 
                 diag_kind='kde', plot_kws={'alpha': 0.4})

for i, y_var in enumerate(lambda_cols):
    for j, x_var in enumerate(lambda_cols):
        if i != j:
            # Mark the global median with a red star
            g.axes[i, j].scatter(general_median[x_var], general_median[y_var], 
                                 marker='*', color='red', s=250, label='Global Median', zorder=5)

g.fig.suptitle("Cleaned Global Clusters and Suggested Medians (*)", y=1.02, fontsize=16)
plt.savefig('lambda_pairplot_all_folders.png')
plt.show()



# 8. Calculate "General" Values (Summary Statistics)
summary = final_df_clean[lambda_cols].describe()
print("Summary Statistics for Lambda Parameters:")
print(summary)