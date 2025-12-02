Usage Guide
===========

This guide covers the main use cases for pycartogram.

Basic Cartogram from Polygons
-----------------------------

The most common use case is creating a cartogram from a list of polygons
with associated density values:

.. code-block:: python

   from shapely.geometry import Polygon
   from pycartogram import WardCartogram

   # Define ward polygons
   wards = [
       Polygon([(0, 0), (1, 0), (1, 1), (0, 1)]),
       Polygon([(1, 0), (2, 0), (2, 1), (1, 1)]),
   ]

   # Define density for each ward
   densities = [10.0, 100.0]

   # Create and compute cartogram
   carto = WardCartogram(
       wards=wards,
       ward_density=densities,
       norm_density=True,  # Normalize by area
       x_raster_size=256,
       y_raster_size=256,
   )
   carto.compute(verbose=True)

   # Access transformed wards
   new_wards = carto.new_wards

Cartogram from Point Data
-------------------------

If you have point location data instead of predefined densities,
use :class:`~pycartogram.PointCartogram`:

.. code-block:: python

   import numpy as np
   from pycartogram import PointCartogram

   # Generate random points (clustered in one area)
   np.random.seed(42)
   points = np.random.randn(1000, 2)
   points[:500] *= 0.3  # Cluster first half

   # Create cartogram from points
   carto = PointCartogram(
       points=points,
       x_raster_size=128,
       y_raster_size=128,
   )
   carto.compute(verbose=True)

   # Plot original and transformed points
   fig, ax = carto.plot_points()

Using with GeoPandas
--------------------

For working with shapefiles and GeoDataFrames, use
:class:`~pycartogram.geopandas_cartogram.GeoDataFrameWardCartogram`:

.. code-block:: python

   import geopandas as gpd
   from pycartogram.geopandas_cartogram import GeoDataFrameWardCartogram

   # Load geographic data
   gdf = gpd.read_file("regions.geojson")
   gdf['density'] = gdf['population'] / gdf.geometry.area

   # Create cartogram
   carto = GeoDataFrameWardCartogram(
       gdf,
       ward_density_column='density',
       x_raster_size=1024,
       y_raster_size=768,
   )
   carto.compute(verbose=True)

   # Get result as GeoDataFrame
   result_gdf = carto.get_cartogram_geo_df()

   # Create animated interpolation
   for t in [0.0, 0.25, 0.5, 0.75, 1.0]:
       frame = carto.get_interpolated_geo_df(t)
       frame.plot()

Plotting Options
----------------

The :meth:`~pycartogram.WardCartogram.plot` method supports various options:

.. code-block:: python

   # Plot with density matrix background
   fig, ax = carto.plot(
       show_density_matrix=True,
       show_new_wards=True,
       ward_colors='viridis',  # Colormap or list of colors
       edge_colors='white',
       bg_color='black',
       outline_whole_shape=True,
   )

   # Plot original wards instead of transformed
   fig, ax = carto.plot(show_new_wards=False)

Parameters Explained
--------------------

Key parameters for cartogram generation:

``x_raster_size``, ``y_raster_size``
   Resolution of the density grid. Higher values give more accurate results
   but take longer to compute. Typical values: 256-2048.

``margin_ratio``
   Fraction of the map size to add as margin around the wards.
   Default: 0.2 (20% margin).

``norm_density``
   If True, divide density values by ward area. Use this when you have
   absolute values (like population counts) rather than densities.

``map_orientation``
   Either 'landscape' or 'portrait'. Determines aspect ratio handling.
