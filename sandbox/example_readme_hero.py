#!/usr/bin/env python
"""
README Hero Image - Appealing cartogram comparison

Creates a side-by-side comparison of original wards vs cartogram
with creative ward shapes and density-based coloring.
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.cm as cm
import numpy as np
from shapely.geometry import Polygon
from pycartogram import WardCartogram

mpl.rcParams['figure.dpi'] = 72

# Define five adjacent wards with creative, organic shapes
# Inspired by real city districts with winding boundaries

# Ward A: Northwestern - rural, low density (farmland)
A = Polygon([
    (0.0, 1.8), (0.2, 2.5), (0.6, 2.8), (1.2, 2.7),
    (1.5, 2.3), (1.3, 1.9), (1.4, 1.5), (1.0, 1.2),
    (0.5, 1.3), (0.1, 1.5), (0.0, 1.8)
])

# Ward B: Northeastern - suburban, medium-low density
B = Polygon([
    (1.5, 2.3), (1.2, 2.7), (1.8, 3.0), (2.5, 2.9),
    (2.9, 2.5), (2.7, 2.0), (2.3, 1.7), (1.8, 1.6),
    (1.4, 1.5), (1.3, 1.9), (1.5, 2.3)
])

# Ward C: Central - urban core, highest density
C = Polygon([
    (1.0, 1.2), (1.4, 1.5), (1.8, 1.6), (2.3, 1.7),
    (2.5, 1.3), (2.3, 0.9), (1.9, 0.7), (1.4, 0.6),
    (1.0, 0.8), (0.8, 1.0), (1.0, 1.2)
])

# Ward D: Southwestern - industrial, medium density
D = Polygon([
    (0.1, 1.5), (0.5, 1.3), (1.0, 1.2), (0.8, 1.0),
    (1.0, 0.8), (0.8, 0.4), (0.4, 0.2), (0.0, 0.5),
    (0.0, 1.0), (0.1, 1.5)
])

# Ward E: Southeastern - mixed use, medium-high density
E = Polygon([
    (1.0, 0.8), (1.4, 0.6), (1.9, 0.7), (2.3, 0.9),
    (2.5, 1.3), (2.8, 1.0), (2.9, 0.5), (2.5, 0.1),
    (1.8, 0.0), (1.2, 0.2), (0.8, 0.4), (1.0, 0.8)
])

wards = [A, B, C, D, E]
# Densities reflecting urban structure: rural < suburban < industrial < mixed < urban core
densities = [3., 12., 95., 25., 45.]

# Normalize densities for colormap
norm_densities = np.array(densities)
norm_densities = (norm_densities - norm_densities.min()) / (norm_densities.max() - norm_densities.min())

# Use a perceptually uniform colormap (viridis or plasma)
cmap = cm.YlOrRd  # Yellow-Orange-Red shows density intuitively
colors = [cmap(d) for d in norm_densities]

# Create and compute cartogram
carto = WardCartogram(
    wards=wards,
    ward_density=densities,
    norm_density=True,
    margin_ratio=0.25,
    x_raster_size=256,
    y_raster_size=256,
)
carto.compute(verbose=True)

# Create comparison figure
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
fig.patch.set_alpha(0.0)  # Transparent figure background

for ax in axes:
    ax.set_aspect('equal')
    ax.axis('off')
    ax.patch.set_alpha(0.0)  # Transparent axes background

# Plot original wards
ax1 = axes[0]
for ward, color in zip(wards, colors):
    x, y = ward.exterior.xy
    ax1.fill(x, y, facecolor=color, edgecolor='#444444', linewidth=1.2)
ax1.set_title('Original Geography', fontsize=13, fontweight='medium', color='#333333', pad=10)

# Plot cartogram
ax2 = axes[1]
for ward, color in zip(carto.new_wards, colors):
    x, y = ward.exterior.xy
    ax2.fill(x, y, facecolor=color, edgecolor='#444444', linewidth=1.2)
ax2.set_title('Density-Equalizing Cartogram', fontsize=13, fontweight='medium', color='#333333', pad=10)

# Match axis limits across both plots
all_wards = wards + carto.new_wards
all_x = []
all_y = []
for w in all_wards:
    x, y = w.exterior.xy
    all_x.extend(x)
    all_y.extend(y)

padding = 0.2
x_min, x_max = min(all_x) - padding, max(all_x) + padding
y_min, y_max = min(all_y) - padding, max(all_y) + padding

for ax in axes:
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)

plt.tight_layout()

# Save with transparent background
output_path = 'img/example_readme_hero.png'
fig.savefig(output_path, dpi=150, bbox_inches='tight', transparent=True)
print(f"Saved: {output_path}")

plt.show()
