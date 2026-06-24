import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

# ==================================================
# 1. CARGAR Y NORMALIZAR DATOS 
# ==================================================
# (A) Real data from IEA Annex 20 (Ux), from Nielsen's paper. 
# We will extract the 4 profiles (2 vertical and 2 horizontal) based on the structure of the CSV file.
data_path = './src/GMRF-wind/gmrf_wind_mapping/IEA_Annex_20_dataset/IEA_Annex_20_real_measurements.csv'
df_data = pd.read_csv(data_path, header=None)

# (B) CFD data (Simscale simulation of the same case)
cfd_path = './src/GMRF-wind/gmrf_wind_mapping/IEA_Annex_20_dataset/CFD_wind_at_cell_centers_convergence.csv'
df_cfd = pd.read_csv(cfd_path)

# Asumiendo columnas: Points:0 (x), Points:1 (y), U:0 (Ux)
# Detectar dimensiones para normalizar
H = round(df_cfd['Points:1'].max())
L = round(df_cfd['Points:0'].max())
U_inlet = 0.4554 # <--- ¡IMPORTANTE! de SimScale
h1 = 0.056*H # Altura del jet de entrada, según la descripción del caso
print(f"Dimensiones detectadas en CFD: L={L:.2f} m, H={H:.2f} m, h1={h1:.4f} m, U_inlet={U_inlet:.4f} m/s")

df_cfd['x_norm'] = df_cfd['Points:0'] / H
df_cfd['y_norm'] = df_cfd['Points:1'] / H
df_cfd['u_norm'] = df_cfd['U:0'] / U_inlet

# ==================================================
# 2. PROCESAR DATOS EXPERIMENTALES
# =================================================
# Función para extraer bloques específicos de un CSV
def extract_block(df, start_row, end_row, cols, profile_name):
    block = df.iloc[start_row:end_row, cols].copy()
    block.columns = ['coord_norm', 'velocity_mean', 'velocity_rms']
    
    # Ensure numeric values
    block['velocity_mean'] = pd.to_numeric(block['velocity_mean'], errors='coerce')
    block['velocity_rms'] = pd.to_numeric(block['velocity_rms'], errors='coerce')
    block['coord_norm'] = pd.to_numeric(block['coord_norm'], errors='coerce')
    # Limpiamos filas que no sean números (metadatos o celdas vacías)
    #block['coord_norm'] = pd.to_numeric(block['coord_norm'], errors='coerce')
    #block = block.dropna(subset=['coord_norm']).reset_index(drop=True)
    block['profile_type'] = profile_name

    # Add x/y coordinates based on profile type
    # y is measured from top to bottom, so we invert it with 1 - coord_norm to have y=0 at the bottom and y=1 at the top
    if 'vertical_x=H' in profile_name:
        block['x'] = 1  # x=H -> x/H = 1
        block['y'] = pd.to_numeric(1- block['coord_norm'],errors='coerce')  # y/H [0,1]
    elif 'vertical_x=2H' in profile_name:
        block['x'] = 2  # x=2H -> x/H = 2
        block['y'] = pd.to_numeric(1- block['coord_norm'],errors='coerce')  # y/H [0,1]
    elif 'horizontal_y=0.5h1' in profile_name:
        # This is the top horizontal profile.
        block['x'] = pd.to_numeric(block['coord_norm'],errors='coerce')  # x/H [0,3]
        block['y'] = 0.5*0.056  # y=0.5h1=0.5*0.056*H -> y/H = 0.5*0.056 = 0.028
    elif 'horizontal_y=H-0.5h1' in profile_name:
        # This is the bottom horizontal profile.
        block['x'] = pd.to_numeric(block['coord_norm'],errors='coerce')  # x/H [0,3]
        block['y'] = 1 - 0.5*0.056  # y=H-0.5h1 -> y/H = 0.972

    #remove coord_norm column as we now have x and y
    block = block.drop(columns=['coord_norm'])

    # Reorder columns to have x, y, velocity_mean, velocity_rms, profile_type
    block = block[['x', 'y', 'velocity_mean', 'velocity_rms', 'profile_type']]
    return block


# Extraemos las 4 tablas de datos, basadas en la estructura de Nielsen:
# 1. Perfil Vertical x = H
df_xH = extract_block(df_data, 10, 35, [1, 2, 3], 'vertical_x=H')
df_xH = df_xH.sort_values(by='y').reset_index(drop=True)
# 2. Perfil Vertical x = 2H
df_x2H = extract_block(df_data, 10, 35, [5, 6, 7], 'vertical_x=2H')
df_x2H = df_x2H.sort_values(by='y').reset_index(drop=True)
# 3. Perfil Horizontal y = 0.5*h1 (cerca del suelo)
df_y05h1 = extract_block(df_data, 41, 69, [1, 2, 3], 'horizontal_y=0.5h1')
df_y05h1 = df_y05h1.sort_values(by='x').reset_index(drop=True)
# 4. Perfil Horizontal y = H - 0.5*h1 (línea del jet de entrada, en el techo)
df_yH05h1 = extract_block(df_data, 41, 69, [5, 6, 7], 'horizontal_y=H-0.5h1')
df_yH05h1 = df_yH05h1.sort_values(by='x').reset_index(drop=True)

# Unimos todo en un único dataframe
df_final = pd.concat([df_xH, df_x2H, df_y05h1, df_yH05h1], ignore_index=True)


# ==================================================
# 3. PLOTTING COMPARATIVO
# ==================================================
fig = plt.figure(figsize=(12, 8))
tol = 0.01 # Tolerancia para considerar puntos cercanos a las líneas de interés
gs = gridspec.GridSpec(2, 3, width_ratios=[2, 2, 1], wspace=0.3, hspace=0.5)

# Row 0, Col 0: Top-left plot -> Perfil Horizontal y = 0.5*h1 (top, inlet jet)
ax1 = fig.add_subplot(gs[0, 0:2])
#plt.subplot(2, 2, 1)
subset = df_y05h1
ax1.plot(subset['x'], subset['velocity_mean'], label='Mean Velocity', marker='o')
ax1.fill_between(subset['x'], subset['velocity_mean'] - subset['velocity_rms'], subset['velocity_mean'] + subset['velocity_rms'], color='gray', alpha=0.3, label='RMS')
# Filtramos celdas cerca del techo (y/H ~ 0.97)
cfd_slice = df_cfd[df_cfd['y_norm'] > 1-tol].sort_values('x_norm')
ax1.plot(cfd_slice['x_norm'], cfd_slice['u_norm'], 'r-', linewidth=2, label='CFD (OpenFOAM)')
ax1.set_title('Horizontal Profile at y=0.5*h1 (top, inlet jet)')
ax1.set_xlabel('$x/H$')
ax1.set_ylabel('$U_x/U_0$')
ax1.grid(True, linestyle='--', alpha=0.6)
ax1.legend()
# set y-ticks (wind speed) to be between -0.5 and 1.5 with a step of 0.5
ax1.set_ylim(-0.5, 1.5)
ax1.set_yticks(np.arange(-0.5, 1.6, 0.5))


# Row 0, Col 2: Top-right plot -> Perfil Vertical x = 2H
ax2 = fig.add_subplot(gs[0, 2])
#plt.subplot(2, 2, 2)
subset = df_x2H
ax2.plot(subset['velocity_mean'], subset['y'], label='Mean Velocity', marker='o')
ax2.fill_betweenx(subset['y'], subset['velocity_mean'] - subset['velocity_rms'], subset['velocity_mean'] + subset['velocity_rms'], color='gray', alpha=0.3, label='RMS')
cfd_slice = df_cfd[(df_cfd['x_norm'] > 2-tol) & (df_cfd['x_norm'] < 2+tol)].sort_values('y_norm')
ax2.plot(cfd_slice['u_norm'], cfd_slice['y_norm'], 'r-', linewidth=2, label='CFD (OpenFOAM)')
ax2.set_title('Vertical Profile at x=2H')
ax2.set_xlabel('$U_x/U_0$')
ax2.set_ylabel('$y/H$')
ax2.grid(True, linestyle='--', alpha=0.6)
#ax2.legend()
# set x-ticks (wind speed) to be between -0.5 and 1 with a step of 0.5
ax2.set_xlim(-0.5, 1)
ax2.set_xticks(np.arange(-0.5, 1.1, 0.5))


# Row 1, Col 0: Bottom-left plot -> Perfil Horizontal y = H - 0.5*h1 (bottom, near the ground)
ax3 = fig.add_subplot(gs[1, 0:2]) 
#plt.subplot(2, 2, 3)
subset = df_yH05h1
ax3.plot(subset['x'], subset['velocity_mean'], label='Mean Velocity', marker='o')
ax3.fill_between(subset['x'], subset['velocity_mean'] - subset['velocity_rms'], subset['velocity_mean'] + subset['velocity_rms'], color='gray', alpha=0.3, label='RMS')
# Filtramos celdas cerca del suelo (y/H ~ 0.03)
cfd_slice = df_cfd[df_cfd['y_norm'] < tol].sort_values('x_norm')
ax3.plot(cfd_slice['x_norm'], cfd_slice['u_norm'], 'r-', linewidth=2, label='CFD (OpenFOAM)')
ax3.set_title('Horizontal Profile at y=H-0.5*h1 (bottom, near the ground)')
ax3.set_xlabel('$x/H$')
ax3.set_ylabel('$U_x/U_0$')
ax3.grid(True, linestyle='--', alpha=0.6)
ax3.legend()
# set y-ticks (wind speed)
ax3.set_ylim(-0.5, 0.5)
ax3.set_yticks(np.arange(-0.5, 0.6, 0.2))


# Row 1, Col 2: Bottom-right plot -> Perfil Vertical x = H
ax4 = fig.add_subplot(gs[1, 2])
#plt.subplot(2, 2, 4)
subset = df_xH
# Experimental data: mean velocity with RMS as shaded area
ax4.plot(subset['velocity_mean'], subset['y'], label='Mean Velocity', marker='o')
ax4.fill_betweenx(subset['y'], subset['velocity_mean'] - subset['velocity_rms'], subset['velocity_mean'] + subset['velocity_rms'], color='gray', alpha=0.3, label='RMS')
# CFD (cells close to x/H = 1)
cfd_slice = df_cfd[(df_cfd['x_norm'] > 1-tol) & (df_cfd['x_norm'] < 1+tol)].sort_values('y_norm')
ax4.plot(cfd_slice['u_norm'], cfd_slice['y_norm'], 'r-', linewidth=2, label='CFD (OpenFOAM)')
ax4.set_title('Vertical Profile at x=H')
ax4.set_xlabel('$U_x/U_0$')
ax4.set_ylabel('$y/H$')
ax4.grid(True, linestyle='--', alpha=0.6)
#ax4.legend()
# x margins linear between -0.5 and 1. Ticks with a step of 0.5
ax4.set_xlim(-0.5, 1)
ax4.set_xticks(np.arange(-0.5, 1.1, 0.5))


plt.tight_layout()
plt.show()