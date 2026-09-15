import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ==============================================================================
# 1. VALORES APROXIMADOS EXTRAÍDOS DE LA FIGURA 7 DE WEI ET AL. (TEST A)
# ==============================================================================
# Los diccionarios contienen aproximaciones de: [Q1, Mediana, Q3, Mínimo, Máximo]
wei_interp_data = {
    16: [0.030, 0.075, 0.145, 0.000, 0.320],
    32: [0.030, 0.072, 0.148, 0.000, 0.330],
    98: [0.030, 0.063, 0.140, 0.000, 0.300]
}

wei_ann_data = {
    16: [0.032, 0.078, 0.152, 0.000, 0.335],
    32: [0.020, 0.053, 0.115, 0.000, 0.255],
    98: [0.025, 0.057, 0.106, 0.000, 0.230]
}

wei_pinn_data = {
    16: [0.018, 0.054, 0.098, 0.000, 0.216],
    32: [0.012, 0.033, 0.070, 0.000, 0.160],
    98: [0.020, 0.051, 0.094, 0.000, 0.205]
}

# ==============================================================================
# 2. PROCESAMIENTO DE TUS EXPERIMENTOS GMRF-W (20 ITERACIONES)
# ==============================================================================
sample_spaces = [98, 32, 16]
gmrf_boxplot_data = {16: [], 32: [], 98: []}

for num_samples in [16, 32, 98]:
    for repeat in range(100):
        # Construir el nombre del archivo exactamente como se guardó en C++
        filename = f"IEA_Annex_20_results_MAE_NoInlet/gmrf_estimation_IAE_annex20_{num_samples}obs_iter{repeat}.csv"
        
        try:
            # Load CSV file
            if not os.path.isfile(filename):
                raise FileNotFoundError(f"CSV file not found: {filename}")

            # 1. Load Metadata (Header)
            # ==========================
            meta = {}
            with open(filename, 'r') as f:
                for _ in range(5):
                    line = f.readline().strip().split(',')
                    # line 5 contains the indices of the observed cells
                    if line[0] == "Obs_idx":
                        meta['Obs_idx'] = [int(idx) for idx in line[1:]]
                    else:
                        meta[line[0]] = float(line[1])

            nx, ny = int(meta['Dimensions_x']), int(meta['Dimensions_y'])
            
            # Load the actual data (skipping the metadata lines)
            df = pd.read_csv(filename, skiprows=5)
            
            # Avoid NaN (obstacles)
            mask = np.isnan(df['gmrf_wind_x'])
            free = ~mask

            # Errors
            error_x = df.loc[free, 'gmrf_wind_x'] - df.loc[free, 'gt_wind_x']
            error_y = df.loc[free, 'gmrf_wind_y'] - df.loc[free, 'gt_wind_y']
            
            # MAE (Mean Absolute Error) normalized
            error_magnitude = np.sqrt(error_x**2 + error_y**2)
            u0 = 0.4554 # wind speed at inlet
            
            normalized_mae_x = np.mean(np.abs(error_x))/u0
            normalized_mae_y = np.mean(np.abs(error_y))/u0
            normalized_mae = np.mean(error_magnitude)/u0

            # save data for boxplot
            gmrf_boxplot_data[num_samples].append(normalized_mae_x)            
            
        except FileNotFoundError:
            print(f"[Error] File not found: {filename}.")
            continue

    # End repeat loop
    gmrf_ordered_list = [gmrf_boxplot_data[98], gmrf_boxplot_data[32], gmrf_boxplot_data[16]]
    print(f"Num samples: {num_samples} - GMRF-W Normalized MAE: {gmrf_boxplot_data[num_samples]}")

# ==============================================================================
# 3. BOXPLOT COMPARISON (FIGURE 7 from Wei et al.) - 
# ==============================================================================
fig, ax = plt.subplots(figsize=(10, 6))

# Config
positions = np.array([1, 2, 3])
width = 0.16

# Helper to plot manual boxplots ("manual" boxplot for Wei Data on Figu 7)
def draw_manual_boxplot(ax, pos, stats_dict, color, label):
    for i, num_samples in enumerate(sample_spaces):
        q1, med, q3, vmin, vmax = stats_dict[num_samples]
        p = pos[i]
        # Cuerpo de la caja (Q1 a Q3)
        ax.bar(p, q3 - q1, bottom=q1, width=width, edgecolor=color, facecolor='none', linewidth=1.5, zorder=3)
        # Línea de la Mediana
        ax.plot([p - width/2, p + width/2], [med, med], color=color, linewidth=2, zorder=3)
        # Bigotes (Whiskers)
        ax.plot([p, p], [q3, vmax], color=color, linestyle='--', linewidth=1.2, zorder=3)
        ax.plot([p, p], [q1, vmin], color=color, linestyle='--', linewidth=1.2, zorder=3)
        # Topes de los bigotes
        ax.plot([p - width/4, p + width/4], [vmax, vmax], color=color, linewidth=1.2, zorder=3)
        ax.plot([p - width/4, p + width/4], [vmin, vmin], color=color, linewidth=1.2, zorder=3)
    # Dummy plot para la leyenda
    ax.plot([], [], color=color, linewidth=1.5, label=label)

# 1. Interpolation 
draw_manual_boxplot(ax, positions - 1.5*width, wei_interp_data, color='#7fcdbb', label='Interpolation')

# 2. ANN 
draw_manual_boxplot(ax, positions - 0.5*width, wei_ann_data, color='#fec44f', label='ANN')

# 3. PINN 
draw_manual_boxplot(ax, positions + 0.5*width, wei_pinn_data, color='#f46d43', label='PINN')

# 4. GMRF-W reales
if any(gmrf_ordered_list):
    bp_gmrf = ax.boxplot(gmrf_ordered_list, positions=positions + 1.5*width, widths=width,
                         patch_artist=True, showfliers=True, zorder=3)
    
    for box in bp_gmrf['boxes']:
        box.set(edgecolor='#1f78b4', linewidth=1.2, facecolor='none')
    for median in bp_gmrf['medians']:
        median.set(color='#1f78b4', linestyle='--', linewidth=1.2)
    for whisker in bp_gmrf['whiskers']:
        whisker.set(color='black', linewidth=0.8)
    for cap in bp_gmrf['caps']:
        cap.set(color='black', linewidth=0.8)
    for flier in bp_gmrf['fliers']:
        flier.set(marker='o', markerfacecolor='#1f78b4', markeredgecolor='none', alpha=0.5, markersize=4)
        
    ax.plot([], [], color='#1f78b4', linewidth=1.5, label='GMRF-W (Ours)')

# ==============================================================================
# 4. AXES AESTHETIC CONFIGURATION
# ==============================================================================
ax.set_ylim(-0.01, 0.35) # Rango lineal que cubre holgadamente el máximo de 0.30 del paper
ax.set_yticks(np.arange(0.00, 0.31, 0.05)) # Ticks idénticos a la gráfica de Wei

ax.set_xticks(positions)
ax.set_xticklabels(['98 points', '32 points', '16 points'], fontsize=11)
ax.set_xlabel('Space sample point amount', fontsize=12, labelpad=10)
ax.set_ylabel(r'$u_x$ absolute error/$u_0$', fontsize=12, labelpad=10)

# Rejilla horizontal gris discontinua únicamente en el eje Y (Estilo exacto de la imagen)
ax.grid(axis='y', linestyle='--', color='lightgray', alpha=0.7, zorder=0)
# Desactivar las líneas de contorno superior y derecha para limpiar el gráfico
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Leyenda idéntica a la distribución lateral
ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', fontsize=10)

plt.tight_layout()
plt.savefig("fig7_perfect_replication.png", dpi=300)
plt.show()