#!/usr/bin/env python
"""
Point Cartogram Example

This script generates the output figure for the PointCartogram usage section.
"""

import numpy as np
import matplotlib.pyplot as plt
from pycartogram import PointCartogram

# Generate random points (clustered in one area)
np.random.seed(42)
points = np.random.randn(1000, 2)
points[:500] *= 0.3  # Cluster first half tightly
points[:500] += [0.5, 0.5]  # Offset the cluster

# Create cartogram from points
carto = PointCartogram(
    points=points,
    x_raster_size=128,
    y_raster_size=128,
)

# Compute density and cartogram
carto.cast_density_to_matrix(verbose=True)
carto.compute(verbose=True)

# Plot original and transformed points
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Original points with density
carto.plot_points(
    show_density_matrix=True,
    show_new_points=False,
    ax=axes[0],
    ms=2,
    mfc='blue',
    alpha=0.5,
)
axes[0].set_title('Original Points with Density', fontsize=12)

# Transformed points
carto.plot_points(
    show_density_matrix=False,
    show_new_points=True,
    ax=axes[1],
    ms=2,
    mfc='red',
    alpha=0.5,
)
axes[1].set_title('Transformed Points (Density Equalized)', fontsize=12)

fig.tight_layout()

# Save figure
output_path = 'img/example_point_cartogram.png'
fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

plt.show()
