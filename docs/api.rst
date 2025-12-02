API Reference
=============

Main Classes
------------

WardCartogram
~~~~~~~~~~~~~

.. autoclass:: pycartogram.WardCartogram
   :members:
   :undoc-members:
   :show-inheritance:

PointCartogram
~~~~~~~~~~~~~~

.. autoclass:: pycartogram.PointCartogram
   :members:
   :undoc-members:
   :show-inheritance:

VoronoiCartogram
~~~~~~~~~~~~~~~~

.. autoclass:: pycartogram.VoronoiCartogram
   :members:
   :undoc-members:
   :show-inheritance:

Data Loaders
------------

GoogleShapeProject
~~~~~~~~~~~~~~~~~~

.. autoclass:: pycartogram.GoogleShapeProject
   :members:
   :undoc-members:
   :show-inheritance:

GeoPandas Integration
---------------------

GeoDataFrameWardCartogram
~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: pycartogram.geopandas_cartogram.GeoDataFrameWardCartogram
   :members:
   :undoc-members:
   :show-inheritance:

Utility Functions
-----------------

.. autofunction:: pycartogram.polygon_patch

.. autofunction:: pycartogram.enrich_polygon_with_points

.. autofunction:: pycartogram.enrich_polygon_to_n_points

.. autofunction:: pycartogram.match_vertex_count

.. autofunction:: pycartogram.coarse_grain_wards

.. autofunction:: pycartogram.get_json

.. autofunction:: pycartogram.get_polygon_network

.. autofunction:: pycartogram.voronoi_finite_polygons_2d
