import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
from shapely.geometry import LineString, MultiPolygon
from shapely.ops import polygonize, unary_union

x = np.array([ 0.38517325,  0.40859912,  0.43296919,  0.4583215 ,  0.4583215 ,
               0.43296919,  0.40859912,  0.38517325,  0.36265506,  0.34100929])
y = np.array([ 62.5       ,  56.17977528,  39.39698492,   0.        ,
               0.        ,  17.34605377,  39.13341671,  60.4180932 ,
               76.02574417,  85.47008547])

# Create the original (potentially self-intersecting) polygon
original_polygon = Polygon(zip(x, y))

# Fix self-intersections using LineString + polygonize approach
ls = LineString(np.c_[x, y])
lr = LineString(ls.coords[:] + ls.coords[:1])
mls = unary_union(lr)
fixed_polygons = list(polygonize(mls))
mp = MultiPolygon(fixed_polygons)

for polygon in fixed_polygons:
    print(polygon)

# Plot before and after
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

# Before: original self-intersecting polygon
ax1 = axes[0]
ax1.set_title('Before: Self-intersecting polygon')
ax1.plot(x, y, 'b-', linewidth=2, label='Boundary')
ax1.plot(x, y, 'ro', markersize=6)
# Close the polygon for visualization
ax1.plot([x[-1], x[0]], [y[-1], y[0]], 'b-', linewidth=2)
# Mark the vertices with numbers
for i, (xi, yi) in enumerate(zip(x, y)):
    ax1.annotate(str(i), (xi, yi), textcoords="offset points",
                 xytext=(5, 5), fontsize=8)
ax1.set_xlabel('x')
ax1.set_ylabel('y')
ax1.grid(True, alpha=0.3)

# After: fixed polygons
ax2 = axes[1]
ax2.set_title(f'After: {len(fixed_polygons)} valid polygon(s)')
colors = ['blue', 'green', 'red', 'orange', 'purple']
for i, poly in enumerate(fixed_polygons):
    px, py = poly.exterior.xy
    ax2.fill(px, py, alpha=0.4, color=colors[i % len(colors)],
             edgecolor=colors[i % len(colors)], linewidth=2,
             label=f'Polygon {i+1} (area={poly.area:.2f})')
ax2.set_xlabel('x')
ax2.set_ylabel('y')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig('self_intersections.png', dpi=150)
plt.show()

