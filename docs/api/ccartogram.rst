C Cartogram Backend
===================

.. module:: cCartogram

Low-level C extension for computing density-equalizing cartograms using
the diffusion method of Gastner and Newman (PNAS, 2004).

This module provides a Python wrapper around Mark Newman's C implementation
of the diffusion-based cartogram algorithm. Cartograms are map projections
where geographic regions are rescaled according to some variable (e.g.,
population), distorting the geometry while preserving topology.

References
----------
- Gastner & Newman (2004): https://doi.org/10.1073/pnas.0400280101
- Original C code: http://www-personal.umich.edu/~mejn/cart/

Functions
---------

.. function:: compute_cartogram(density, offset=0.005, blur=0.0, show_progress=False)

   Compute a density-equalizing cartogram transformation grid.

   Uses the diffusion method of Gastner and Newman to compute a cartogram
   from a 2D density matrix. The result is a displacement grid that can be
   used with :func:`remap_coordinates` to transform arbitrary coordinates.

   :param density: A 2D matrix of density values (e.g., population density).
       All values must be positive. The matrix shape determines the cartogram
       grid size.
   :type density: list of lists of float
   :param offset: Small positive value added to density to avoid division by
       zero and improve numerical stability.
   :type offset: float, optional
   :param blur: Gaussian blur radius applied to the density before computation.
       Useful for smoothing noisy input data.
   :type blur: float, optional
   :param show_progress: If True, print progress information during computation.
   :type show_progress: bool, optional
   :returns: A flattened list of (x, y) displacement pairs representing the
       transformed grid coordinates. The list has length ``xsize * ysize``,
       stored in row-major order.
   :rtype: list of tuple of (float, float)

   **Example**::

       import numpy as np
       import cCartogram as cart

       # Create a density matrix with higher values in one region
       density = np.ones((128, 128))
       density[32:96, 32:96] = 4.0  # Higher density in center
       cartogram = cart.compute_cartogram(density.tolist(), blur=2.0)


.. function:: remap_coordinates(coordinates, cartogram, xsize, ysize)

   Transform coordinates according to a precomputed cartogram.

   Takes a list of (x, y) coordinates and transforms them using the
   displacement grid computed by :func:`compute_cartogram`. Coordinates are
   interpolated bilinearly between grid points.

   :param coordinates: List of (x, y) coordinate pairs to transform.
       Coordinates should be in the range [0, xsize) and [0, ysize) respectively.
   :type coordinates: list of tuple of (float, float)
   :param cartogram: The cartogram displacement grid returned by
       :func:`compute_cartogram`.
   :type cartogram: list of tuple of (float, float)
   :param xsize: Width of the original density matrix used to compute the cartogram.
   :type xsize: int
   :param ysize: Height of the original density matrix used to compute the cartogram.
   :type ysize: int
   :returns: The transformed coordinates in the same order as the input.
   :rtype: list of tuple of (float, float)

   **Example**::

       import cCartogram as cart

       # After computing a cartogram...
       old_coords = [(10.5, 20.0), (50.0, 75.5)]
       new_coords = cart.remap_coordinates(old_coords, cartogram, 128, 128)
       for old, new in zip(old_coords, new_coords):
           print(f"{old} -> {new}")
