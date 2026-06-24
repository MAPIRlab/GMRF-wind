import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ==============================================================================
# DATOS COMPLETOS EXTRAÍDOS DE LA FIGURA 9 DE WEI ET AL. (BARRIDO 6:10:96)
# Estructura: [Q1, Mediana/Media, Q3, Mínimo, Maksimo]
# ==============================================================================
sample_spaces = [6, 16, 26, 36, 46, 56, 66, 76, 86, 96]

wei_fig9_interp = {
    6:  [0.238, 0.248, 0.254, 0.175, 0.274],
    16: [0.218, 0.225, 0.234, 0.170, 0.249],
    26: [0.198, 0.206, 0.218, 0.155, 0.245],
    36: [0.192, 0.202, 0.212, 0.142, 0.244],
    46: [0.186, 0.194, 0.208, 0.150, 0.223],
    56: [0.184, 0.191, 0.205, 0.148, 0.225],
    66: [0.182, 0.189, 0.200, 0.152, 0.222],
    76: [0.180, 0.188, 0.198, 0.151, 0.214],
    86: [0.181, 0.189, 0.202, 0.153, 0.210],
    96: [0.182, 0.192, 0.199, 0.155, 0.205]
}

wei_fig9_ann = {
    6:  [0.195, 0.225, 0.252, 0.115, 0.298],
    16: [0.160, 0.201, 0.233, 0.110, 0.290],
    26: [0.182, 0.211, 0.231, 0.138, 0.270],
    36: [0.128, 0.156, 0.181, 0.090, 0.220],
    46: [0.120, 0.144, 0.167, 0.088, 0.205],
    56: [0.115, 0.131, 0.146, 0.085, 0.181],
    66: [0.111, 0.126, 0.135, 0.085, 0.165],
    76: [0.106, 0.122, 0.132, 0.085, 0.155],
    86: [0.103, 0.120, 0.126, 0.085, 0.145],
    96: [0.101, 0.112, 0.122, 0.085, 0.135]
}

wei_fig9_pinn = {
    6:  [0.207, 0.226, 0.247, 0.157, 0.303],
    16: [0.147, 0.166, 0.179, 0.120, 0.226],
    26: [0.135, 0.149, 0.160, 0.110, 0.195],
    36: [0.125, 0.138, 0.146, 0.105, 0.170],
    46: [0.120, 0.129, 0.135, 0.100, 0.165],
    56: [0.115, 0.126, 0.131, 0.098, 0.150],
    66: [0.112, 0.123, 0.128, 0.095, 0.153],
    76: [0.110, 0.120, 0.125, 0.095, 0.145],
    86: [0.110, 0.119, 0.124, 0.095, 0.144],
    96: [0.110, 0.115, 0.121, 0.095, 0.136]
}

# ==============================================================================
# 2. PROCESAMIENTO DE TUS EXPERIMENTOS GMRF-W
# ==============================================================================
gmrf_rmse_data = {n: [] for n in sample_spaces}

for num_samples in sample_spaces:
    for repeat in range(200):
        # Construir el nombre del archivo exactamente como se guardó en C++
        filename = f"IEA_Annex_20_results_RMSE_CFD_highNoise/gmrf_estimation_IAE_annex20_{num_samples}obs_iter{repeat}.csv"
        
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
            u0 = 0.4554 # wind speed at inlet
            
            normalized_rmse_x = np.sqrt(np.mean(error_x**2))/u0
            normalized_rmse_y = np.sqrt(np.mean(error_y**2))/u0
            normalized_rmse = np.sqrt(np.mean(error_x**2 + error_y**2))/u0

            # save data for boxplot
            gmrf_rmse_data[num_samples].append(normalized_rmse_x)            
            
        except FileNotFoundError:
            print(f"[Error] File not found: {filename}.")
            continue

    # End repeat loop
    gmrf_ordered_list = [gmrf_rmse_data[n] for n in sample_spaces]
    
# ==============================================================================
# 3. RMSE COMPARISON (FIGURE 9 from Wei et al.) - 
# ==============================================================================
fig, ax = plt.subplots(figsize=(12, 6))

# Config
positions = np.array(sample_spaces, dtype=float)
width = 1.6

# Helper to plot manual boxplots ("manual" boxplot for Wei Data on Figu 7)
def draw_manual_boxplot(ax, pos_offsets, stats_dict, color, label):
    medians_y = []
    for num_samples in sample_spaces:
        q1, med, q3, vmin, vmax = stats_dict[num_samples]
        p = num_samples + pos_offsets
        medians_y.append(med)
        
        # Caja, bigotes y topes
        ax.bar(p, q3 - q1, bottom=q1, width=width, edgecolor=color, facecolor='none', linewidth=1.1, zorder=3)
        ax.plot([p - width/2, p + width/2], [med, med], color=color, linestyle='--', linewidth=1.1, zorder=3)
        ax.plot([p, p], [q3, vmax], color='black', linewidth=0.7, zorder=3)
        ax.plot([p, p], [q1, vmin], color='black', linewidth=0.7, zorder=3)
        ax.plot([p - width/4, p + width/4], [vmax, vmax], color='black', linewidth=0.7, zorder=3)
        ax.plot([p - width/4, p + width/4], [vmin, vmin], color='black', linewidth=0.7, zorder=3)
        
    # Dibujar la línea de tendencia discontinua que une las medias (Estilo exacto Fig 9)
    ax.plot(positions + pos_offsets, medians_y, color=color, linestyle=':', alpha=0.7, linewidth=1.2, zorder=2)
    ax.plot([], [], color=color, linewidth=1.5, label=label)

# Comparison (in the same order as the original figure: Interpolation, ANN, PINN, GMRF-W)
draw_manual_boxplot(ax, -1.5*width, wei_fig9_interp, color='#7fcdbb', label='Interpolation')
draw_manual_boxplot(ax, -0.5*width, wei_fig9_ann, color='#fec44f', label='ANN')
draw_manual_boxplot(ax, 0.5*width, wei_fig9_pinn, color='#f46d43', label='PINN')

# 4. GMRF-W reales
if any(gmrf_ordered_list):
    # Boxplot nativo para el GMRF-W desviado a la derecha extrema
    bp_gmrf = ax.boxplot(gmrf_ordered_list, positions=positions + 1.5*width, widths=width,
                         patch_artist=True, showfliers=True, zorder=3)
    
    gmrf_medians = []
    for box, num_samples in zip(bp_gmrf['boxes'], sample_spaces):
        box.set(edgecolor='#1f78b4', linewidth=1.1, facecolor='none')
    for median in bp_gmrf['medians']:
        median.set(color='#1f78b4', linestyle='--', linewidth=1.1)
        gmrf_medians.append(median.get_ydata()[0])
    for whisker in bp_gmrf['whiskers']:
        whisker.set(color='black', linewidth=0.7)
    for cap in bp_gmrf['caps']:
        cap.set(color='black', linewidth=0.7)
    for flier in bp_gmrf['fliers']:
        flier.set(marker='o', markerfacecolor='#1f78b4', markeredgecolor='none', alpha=0.4, markersize=3.5)
        
    # Línea de tendencia para el GMRF-W
    ax.plot(positions + 1.5*width, [np.mean(gmrf_rmse_data[n]) for n in sample_spaces], 
            color='#1f78b4', linestyle=':', alpha=0.7, linewidth=1.2, zorder=2)
    
    ax.plot([], [], color='#1f78b4', linewidth=1.5, label='GMRF-W (Ours)')

# ==============================================================================
# 4. AXES AESTHETIC CONFIGURATION
# ==============================================================================
ax.set_xlim(0, 105)
ax.set_ylim(0.05, 0.35)
ax.set_yticks(np.arange(0.10, 0.31, 0.05))

# Forzamos los ticks del eje X a mostrar exactamente las etiquetas del paper
ax.set_xticks(sample_spaces)
ax.set_xticklabels([str(n) for n in sample_spaces], fontsize=10)

ax.set_xlabel('Space sample point amount', fontsize=12, labelpad=10)
ax.set_ylabel(r'$u_x / u_0$ RMSE', fontsize=12, labelpad=10)

# Rejilla horizontal idéntica
ax.grid(axis='y', linestyle='--', color='lightgray', alpha=0.7, zorder=0)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.legend(loc='upper right', frameon=True, facecolor='white', edgecolor='none', fontsize=10)

plt.tight_layout()
plt.savefig("fig9_replicated_full_sweep.png", dpi=300)
plt.show()