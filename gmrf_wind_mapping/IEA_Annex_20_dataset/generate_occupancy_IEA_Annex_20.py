# Generates the occupancy grid map for the IEA Annex 20 dataset, following ROS 2 conventions.
# The map is a simple rectangle with two openings (inlet and outlet) based on the dimensions of the experimental setup. 
# The output is a PGM image and a YAML file for ROS 2.
import numpy as np
import cv2
import yaml

# ==========================================
# 1. PARAMS (METERS)
# ==========================================
H = 3.0       # Room height
L = 9.0       # Room length (Ratio L/H = 3)
h1 = 0.168    # Inlet height (0.056 * H)
h2 = 0.48     # Outlet height (0.16 * H)

resolution = 0.01  # 0.01 meters per pixel (1 cm/pixel) -> 100 pixels/meter
wall_thickness_m = 0.10 #m

# Calculate dimensions in pixels
free_width_px = int(L / resolution)
free_height_px = int(H / resolution)
thick_px = int(wall_thickness_m / resolution)
width_px = free_width_px + 2 * thick_px
height_px = free_height_px + 2 * thick_px

h1_px = int(h1 / resolution)
h2_px = int(h2 / resolution)

# ==========================================
# 2. CREATE MAP MATRIX (ROS Style: 0 = occupied, 255 = free)
# ==========================================
# In ROS: 0 is free (white in standard image, but ROS reads it inversely)
# For PNG image: 255 is free (white), 0 is occupied/wall (black)
map_img = np.zeros((height_px, width_px), dtype=np.uint8)

# Fill the interior with "free" (255 = white)
# Leave the borders as 0 (black = wall)
map_img[thick_px : thick_px + free_height_px, thick_px : thick_px + free_width_px] = 255

# --- VENTILATION HOLES (Open the walls) ---
# Inlet: Top-left corner (attached to the ceiling)
# The Y-axis of the image starts at the top (0) and goes down to height_px
map_img[thick_px : thick_px + h1_px, 0 : thick_px] = 255

# Outlet: Bottom-right corner (attached to the floor)
map_img[thick_px + free_height_px - h2_px : thick_px + free_height_px, thick_px + free_width_px : width_px] = 255

# Save the image in PGM format (preferred by ROS 2)
map_name = "IEA_Annex_20_occupancy"
cv2.imwrite(f"{map_name}.pgm", map_img)

# ==========================================
# 3. GENERATE YAML FILE FOR ROS 2
# ==========================================
origin_x = -thick_px * resolution
origin_y = -thick_px * resolution

yaml_data = {
    'image': f"{map_name}.pgm",
    'resolution': resolution,
   'origin': [float(origin_x), float(origin_y), 0.0], # Origin at the bottom-left corner (X, Y, Theta)
    'negate': 0,
    'occupied_thresh': 0.65,
    'free_thresh': 0.196
}

with open(f"{map_name}.yaml", 'w') as f:
    yaml.dump(yaml_data, f, default_flow_style=False)

print(f"Map successfully generated!")
print(f"Image dimensions: {width_px}x{height_px} pixels.")
print(f"Generated files: {map_name}.pgm and {map_name}.yaml")