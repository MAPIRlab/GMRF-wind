# Script to visualize the evolution and distribution of lambda parameters
# across multiple experimental runs, varying the scenario, wind speed and number of observations.

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import glob
import os
import numpy as np
import re

# ----------------------------------
# 1. CONFIGURATION
# ----------------------------------
folders = [#'',
           #'optimization_results_central_obstacle_01ms',
           #'optimization_results_central_obstacle_1ms',
           #'optimization_results_snake_01ms',
           #'optimization_results_snake_1ms',
           #'optimization_results_expC_2-3,5',
           'optimization_results_expC_1-6',
           ]
lambda_cols = ['Lambda_Advection', 'Lambda_Mass', 'Lambda_Diffusion', 'Lambda_Obstacles']
normalized_cols = [f'Normalized_{col}' for col in lambda_cols]
metric_cols = ['AAE', 'AME']
cleaned_data_list = []
folder_general_medians = []
remove_outliers = False

def extract_Number(filename):
    """Extrae el primer número encontrado en el nombre del archivo."""
    numbers = re.findall(r'\d+', filename)
    return int(numbers[0]) if numbers else filename


# 2. Process each folder/scenario
# ----------------------------------
for folder in folders:
    # Load all CSVs in the folder
    files = glob.glob(os.path.join(folder, "Lambda_values_*.csv"))
    folder_dfs = []

    for file in files:
        # Load CSV
        df = pd.read_csv(file)
        
        # Add metadata so we can distinguish between folders in the plot
        experiment_results = df
        experiment_results['Source_Folder'] = folder
        experiment_results['File_Name'] = os.path.basename(file)

        # Set normalized values in new columns (relative to sum of lambdas)
        lambda_sum = experiment_results[lambda_cols].sum(axis=1)
        for col in lambda_cols:
            experiment_results[f'Normalized_{col}'] = experiment_results[col] / lambda_sum

        # Outlier Removal (per file)
        if remove_outliers:
            # For each lambda column, remove rows where the value is outside 1.5*IQR
            mask = pd.Series(True, index=experiment_results.index)
            for col in lambda_cols:
                Q1 = experiment_results[col].quantile(0.25)
                Q3 = experiment_results[col].quantile(0.75)
                IQR = Q3 - Q1
                lower = Q1 - 1.5 * IQR
                upper = Q3 + 1.5 * IQR
                mask &= (experiment_results[col] >= lower) & (experiment_results[col] <= upper)
            experiment_results = experiment_results[mask]
            print(f"File '{file}': Kept {len(experiment_results)}/{len(df)} rows after outlier removal.")

        # Append cleaned experiment results to the folder list
        folder_dfs.append(experiment_results)

    if not folder_dfs:
        continue
        
    # Combine all experiments (NumObs) in this folder into one DataFrame
    folder_combined = pd.concat(folder_dfs, ignore_index=True)          
   
    # 4. Per-Folder Outlier Removal
    # --------------------------------------------------
    # Remove outliers based on IQR for each lambda column
    if remove_outliers:
        mask = pd.Series(True, index=folder_combined.index)
        for col in lambda_cols:
            Q1 = folder_combined[col].quantile(0.25)
            Q3 = folder_combined[col].quantile(0.75)
            IQR = Q3 - Q1
            lower = Q1 - 1.5 * IQR
            upper = Q3 + 1.5 * IQR
            mask &= (folder_combined[col] >= lower) & (folder_combined[col] <= upper)
        
        folder_combined = folder_combined[mask]
        print(f"Folder '{folder}': Kept {len(folder_combined)}/{len(folder_combined)} rows.")

    # Store the cleaned data for this folder
    cleaned_data_list.append(folder_combined)
    
    # Calculate the "per-folder"=scenario general hyperparameters (median) and store
    folder_general_median = folder_combined[lambda_cols].median()
    folder_general_median['Source_Folder'] = folder
    folder_general_medians.append(folder_general_median)    

# ALL FOLDERS PROCESSED
# ----------------------------------

# 5. Combine everything into one big DataFrame
final_df_clean = pd.concat(cleaned_data_list, ignore_index=True)

# 6. Calculate the "General Hyperparameters" (Centroid)
general_mean = final_df_clean[normalized_cols].mean()
general_median = final_df_clean[normalized_cols].median()

print("\n--- SUGGESTED GLOBAL VALUES (Cleaned Data among all experiments) ---")
summary_table = pd.DataFrame({'Mean': general_mean, 'Median': general_median})
print(summary_table)


# 7. Visualization:
plt.figure(figsize=(12, 10))
sns.set_theme(style="whitegrid")

# PLOT1: Lambda Boxplot
df_melted = final_df_clean.melt(id_vars=['Source_Folder'], value_vars=lambda_cols, 
                                var_name='Parameter', value_name='Value')
sns.boxplot(data=df_melted, x='Parameter', y='Value', hue='Source_Folder', 
            showfliers=True, 
            medianprops={'color': 'red', 'linewidth': 2}, 
            showmeans=True, 
            meanprops={
                "marker": "D",              # "D" para diamante, "o" para círculo, "^" para triángulo
                "markerfacecolor": "white", # Interior blanco para que resalte
                "markeredgecolor": "black", # Borde negro
                "markersize": "6"           # Tamaño del marcador
            })
plt.title("Distribution of Lambda Parameters Across All Folders", fontsize=16)
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig('lambda_boxplot_all_folders.png')
plt.show()

# PLOT 2: Boxplot (normalized values)
df_normalized_melted = final_df_clean.melt(id_vars=['Source_Folder'], value_vars=normalized_cols, 
                                var_name='Parameter', value_name='Value')
sns.boxplot(data=df_normalized_melted, x='Parameter', y='Value', hue='Source_Folder', 
            showfliers=True, 
            medianprops={'color': 'red', 'linewidth': 2}, 
            showmeans=True, 
            meanprops={
                "marker": "D",              # "D" para diamante, "o" para círculo, "^" para triángulo
                "markerfacecolor": "white", # Interior blanco para que resalte
                "markeredgecolor": "black", # Borde negro
                "markersize": "6"           # Tamaño del marcador
            })
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Optimal Normalized Lambda Parameters", fontsize=16)
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig('normalized_lambda_boxplot_all_folders.png')
plt.show()

# PLOT: Metrics Boxplot
df_melted = final_df_clean.melt(id_vars=['Source_Folder'], value_vars=metric_cols, 
                                var_name='Parameter', value_name='Value')
sns.boxplot(data=df_melted, x='Parameter', y='Value', hue='Source_Folder', 
            showfliers=True, 
            medianprops={'color': 'red', 'linewidth': 2}, 
            showmeans=True, 
            meanprops={
                "marker": "D",              # "D" para diamante, "o" para círculo, "^" para triángulo
                "markerfacecolor": "white", # Interior blanco para que resalte
                "markeredgecolor": "black", # Borde negro
                "markersize": "6"           # Tamaño del marcador
            })
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.title("Optimized Metrics", fontsize=16)
plt.xticks(rotation=15)
plt.tight_layout()
plt.savefig('metrics_boxplot_all_folders.png')
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



# 9. Per-Folder General Medians Summary and Plot
folder_general_medians_df = pd.DataFrame(folder_general_medians)
print("\n--- PER-FOLDER GENERAL MEDIANS ---")
#folder_general_medians_df[lambda_cols] = folder_general_medians_df[lambda_cols].div(folder_general_medians_df[lambda_cols].sum(axis=1), axis=0) * 100
print(folder_general_medians_df)

# Plot per-folder general medians
plt.figure(figsize=(10, 6))
folder_general_medians_df.set_index('Source_Folder')[lambda_cols].T.plot(kind='bar', width=0.8)
plt.title("Per-Folder General Lambda Medians")
plt.xlabel("Source Folder")
plt.ylabel("Lambda Value")
plt.xticks(rotation=45)
plt.legend(title="Lambda Parameter")
plt.tight_layout()
plt.savefig('per_folder_lambda_medians.png')
plt.show()