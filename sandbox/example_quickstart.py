#!/usr/bin/env python
"""
Quick Start Example - Basic WardCartogram

This script generates the output figure for the README Quick Start section.
"""

import matplotlib.pyplot as plt
from shapely.geometry import Polygon
from pycartogram import WardCartogram

# Define three adjacent rectangular wards
A = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
B = Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
C = Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])

wards = [A, B, C]
ward_population = [2., 20., 200.]

# Create and compute cartogram
carto = WardCartogram(
    wards=wards,
    ward_density=ward_population,
    norm_density=True,
    margin_ratio=0.5,
    x_raster_size=128,
    y_raster_size=64,
)
carto.compute(verbose=True)

# Plot result
fig, ax = carto.plot(show_new_wards=True, edge_colors='k')
ax.set_title('Cartogram: Areas Proportional to Population', fontsize=12)
fig.tight_layout()

# Save figure
output_path = 'img/example_quickstart.png'
fig.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
print(f"Saved: {output_path}")

plt.show()
