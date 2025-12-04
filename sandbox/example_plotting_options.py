#!/usr/bin/env python
"""
Plotting Options Example

This script demonstrates various plotting options available in pycartogram.
"""

import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from pycartogram import WardCartogram

# Create a simple 3x3 grid of wards with varying densities
wards = []
densities = []

for i in range(3):
    for j in range(3):
        ward = Polygon([
            (i, j), (i+1, j), (i+1, j+1), (i, j+1)
        ])
        wards.append(ward)
        # Density increases toward bottom-right
        densities.append((i + 1) * (j + 1) * 10)

# Create and compute cartogram
carto = WardCartogram(
    wards=wards,
    ward_density=densities,
    norm_density=True,
    x_raster_size=128,
    y_raster_size=128,
    margin_ratio=0.3,
)
carto.compute(verbose=True)

# Create figure with different plotting options
fig, axes = plt.subplots(2, 2, figsize=(10, 10))

# 1. Original wards only
carto.plot(
    show_new_wards=False,
    show_density_matrix=False,
    ward_colors='lightblue',
    edge_colors='darkblue',
    ax=axes[0, 0],
)
axes[0, 0].set_title('Original Wards', fontsize=12)

# 2. Density matrix visualization
carto.plot(
    show_new_wards=False,
    show_density_matrix=True,
    ax=axes[0, 1],
)
axes[0, 1].set_title('Density Matrix', fontsize=12)

# 3. Cartogram with custom colors
colors = plt.cm.viridis([d / max(densities) for d in densities])
carto.plot(
    show_new_wards=True,
    ward_colors=colors,
    edge_colors='white',
    ax=axes[1, 0],
)
axes[1, 0].set_title('Cartogram (Colored by Density)', fontsize=12)

# 4. Cartogram with outline
carto.plot(
    show_new_wards=True,
    ward_colors='#ffcccc',
    edge_colors='darkred',
    outline_whole_shape=True,
    ax=axes[1, 1],
)
axes[1, 1].set_title('Cartogram with Outline', fontsize=12)

fig.tight_layout()

# Save figure
output_path = 'img/example_plotting_options.png'
fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

plt.show()
