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
folders = [#'NLPD_optimization_results_central_obstacle_01ms_phase1', 
           #'NLPD_optimization_results_central_obstacle_05ms_phase1', 
           #'NLPD_optimization_results_central_obstacle_1ms_phase1',
           #'NLPD_optimization_results_snake_01ms_phase1',
           #'NLPD_optimization_results_snake_05ms_phase1',
           #'NLPD_optimization_results_snake_1ms_phase1',
           'NSP_optimization_results_central_obstacle_01ms', 
           'NSP_optimization_results_central_obstacle_05ms', 
           'NSP_optimization_results_central_obstacle_1ms',
           #'NSP_optimization_results_snake_01ms',
           #'NSP_optimization_results_snake_05ms',
           #'NSP_optimization_results_snake_1ms',
           ]
lambda_cols = ['Lambda_Reg', 'Lambda_Flux', 'Lambda_Obstacles', 'Lambda_Observations']
cleaned_data_list = []
folder_general_medians = []

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
    folder_medians = []

    for file in files:

        # Exclude files with less than 50 observations
        n_obs = extract_Number(os.path.basename(file))
        if n_obs < 50:
            continue

        # Load CSV
        df = pd.read_csv(file)
        
        # Use only the first line of each Repetition (the second is a baseline)
        experiment_results = df.groupby('Repetition').nth(0).reset_index()
        
        # Add metadata so we can distinguish between folders in the plot
        experiment_results['Source_Folder'] = folder
        experiment_results['File_Name'] = os.path.basename(file)
        folder_dfs.append(experiment_results)

        # Calculamos la mediana de este archivo específico (de las 20 repeticiones)
        file_median = experiment_results[lambda_cols].median().to_frame().T
        file_median['Samples'] = extract_Number(os.path.basename(file))
        file_median['Source_Folder'] = folder
        folder_medians.append(file_median)        

    if not folder_dfs:
        continue
        
    folder_combined = pd.concat(folder_dfs, ignore_index=True)          # Combine all experiments (Num Obs) in this folder
    folder_medians_df = pd.concat(folder_medians, ignore_index=True)    # Combine all medians
    # reorder rows by Samples
    folder_medians_df = folder_medians_df.sort_values(by='Samples')
    
    # 3. Plot Median Lambda values vs Nobs (samples) in this folder
    # ------------------------------------------------
    filename_medians = folder.replace("NLPD_optimization_results_", "").replace("_phase1", "")
    plt.figure(figsize=(10, 6))
    for col in lambda_cols:
        plt.plot(folder_medians_df['Samples'], folder_medians_df[col], marker='o', label=col)#, linestyle='None')
    plt.title(f'Evolución de Lambdas según muestras - Carpeta: {filename_medians}')
    plt.xlabel('Número de Muestras')
    plt.ylabel('Valor de Lambda (Mediana)')
    #plt.xscale('log')  # Escala logarítmica suele ser útil si las muestras crecen mucho
    plt.grid(True, which="both", ls="-", alpha=0.5)
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'lambda_evolution_{filename_medians}.png', dpi=300)
    plt.show()


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

    # 5. Calculate the "per-folder" general hyperparameters (median) and store
    folder_general_median = df_folder_clean[lambda_cols].median()
    folder_general_median['Source_Folder'] = folder
    folder_general_medians.append(folder_general_median)
    

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



# 9. Per-Folder General Medians Summary and Plot (normalizing sum = 100)
folder_general_medians_df = pd.DataFrame(folder_general_medians)
print("\n--- PER-FOLDER GENERAL MEDIANS (NORMALIZED)---")
folder_general_medians_df[lambda_cols] = folder_general_medians_df[lambda_cols].div(folder_general_medians_df[lambda_cols].sum(axis=1), axis=0) * 100
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