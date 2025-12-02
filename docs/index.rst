pycartogram
============

A Python package for generating cartograms using the Gastner-Newman diffusion algorithm.

Cartograms are maps where region sizes are distorted to be proportional to a variable
of interest (e.g., population), while maintaining the topology and shape of the regions
as much as possible.

.. note::

   This package requires `cCartogram <https://github.com/benmaier/cCartogram>`_,
   the Python port of Mark Newman's original C implementation.

Installation
------------

.. code-block:: bash

   # Install core package
   pip install pycartogram

   # Install with geographic extras (cartopy, geopandas)
   pip install pycartogram[geo]

   # Install with all extras
   pip install pycartogram[all]

Quick Start
-----------

.. code-block:: python

   import matplotlib.pyplot as plt
   from shapely.geometry import Polygon
   from pycartogram import WardCartogram

   # Define three adjacent rectangular wards
   A = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
   B = Polygon([(1, 0), (2, 0), (2, 1), (1, 1)])
   C = Polygon([(2, 0), (3, 0), (3, 1), (2, 1)])

   wards = [A, B, C]
   ward_population = [2., 20., 200.]

   # Create cartogram
   carto = WardCartogram(
       wards=wards,
       ward_density=ward_population,
       norm_density=True,
       margin_ratio=0.5,
       x_raster_size=128,
       y_raster_size=64,
   )

   # Compute the transformation
   carto.compute(verbose=True)

   # Plot result
   fig, ax = carto.plot(show_new_wards=True, edge_colors='k')
   plt.show()

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   usage
   api

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
