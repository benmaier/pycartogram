import cartopy.crs as ccrs
import numpy as np
from cartopy.io.shapereader import Reader
import shapely.geometry as sgeom
from shapely.geometry import Polygon
import pycartogram

import matplotlib.pyplot as pl
import matplotlib.colors as mcolors
from colorsys import hsv_to_rgb

def random_harmonious_colors(n, saturation=0.5, value=0.85, seed=None):
    """Generate n random colors with consistent saturation and brightness."""
    if seed is not None:
        np.random.seed(seed)
    hues = np.random.random(n)
    return [hsv_to_rgb(h, saturation, value) for h in hues]

def get_random_exponential(carto):

    points = np.random.random((10000,2))

    ymin = carto.orig_bbox.bounds[1]
    xmin = carto.orig_bbox.bounds[0]

    mean = (carto.orig_height*0.2)
    new_y = np.log(1.-points[:,1])/(-1./mean)
    indices = np.where(new_y>carto.orig_height)[0]
    new_y[indices] = carto.orig_height * np.random.random((len(indices),))
    points[:,1] = new_y + ymin
    points[:,0] = carto.orig_width * points[:,0] + xmin

    return points

def load_wards(shapefile):


    shape_data = [
            geom[0] if hasattr(geom, '__getitem__') else geom \
            for geom in Reader(shapefile).geometries()\
          ]

    target_proj = ccrs.UTM(zone=33, southern_hemisphere=False)  # used to be '33N'
    source_proj = ccrs.PlateCarree()

    new_shape_data = []
    for geom in shape_data:
        poly_coords = []
        for x,y in geom.exterior.coords:
            poly_coords.append(
                      target_proj.transform_point(
                               x,
                               y,
                               source_proj
                            )
            )
        new_shape_data.append(Polygon(poly_coords))

    return new_shape_data

berlin_wards = load_wards('../data/berlin_postleitzahlen.shp')

carto = pycartogram.WardCartogram(
            wards = berlin_wards,
            ward_density = None,
            x_raster_size = 1024/8,
            y_raster_size = 786/8,
            margin_ratio = 0.2, # added margin will be this ratio of original width
            map_orientation = 'landscape',
        )

berlin_whole = carto.whole_shape

points = get_random_exponential(carto)
carto.compute_ward_density_from_locations(points,verbose=True)
carto.compute(verbose=True)

fig, ax = pl.subplots(2,3,figsize=(13,6))

carto.plot(show_new_wards=False,
           show_density_matrix=False,
           ward_colors = 'w',
           ax=ax[0,0])
ax[0,0].set_title('location data')

ax[0,0].plot(points[:,0],points[:,1],'.',c='g',ms=0.5)

carto.plot(show_new_wards=False,
           show_density_matrix=True,
           ax=ax[0,1])
ax[0,1].set_title('density matrix binned by wards')

carto.plot(show_new_wards=False,
           ax=ax[1,0])
ax[1,0].set_title('ward bins as polygons')

# Random but harmonious colors (pastel-ish, same saturation/brightness)
ward_colors = random_harmonious_colors(len(berlin_wards), saturation=0.45, value=0.9, seed=42)

carto.plot(
        ward_colors = ward_colors,
        edge_colors = 'w',
        ax=ax[1,1])
ax[1,1].set_title('cartogram after ward binning')

carto_point = pycartogram.PointCartogram(
            points = points,
            wards = [berlin_whole],
            x_raster_size = 1024/8,
            y_raster_size = 786/8,
            margin_ratio = 0.2, # added margin will be this ratio of original width
            )

carto_point.cast_density_to_matrix(verbose=True)

carto_point.plot_points(
           show_density_matrix=True,
           show_new_points=False,
           ax=ax[0,2],
           ms = 0.2,
           mfc = 'g'
           )
ax[0,2].set_title('density matrix binned in pixels')

carto_point.compute(verbose=True)
carto_point.transform_wards(verbose=True)

carto_point.plot(
        ward_colors = random_harmonious_colors(1, saturation=0.45, value=0.9, seed=7),
        edge_colors = 'w',
        ax=ax[1,2])
ax[1,2].set_title('cartogram after pixel binning')

fig.tight_layout()
fig.savefig('./berlin.png')

pl.show()


