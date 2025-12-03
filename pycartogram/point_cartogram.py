"""
Point-based cartogram generation.

This module provides the PointCartogram class for creating cartograms
directly from point location data, without requiring predefined ward densities.
"""

from __future__ import annotations

from typing import Any
import numpy as np
from numpy.typing import ArrayLike, NDArray
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.geometry import Polygon
from pycartogram.tools import polygon_patch
import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib.figure import Figure
from matplotlib.axes import Axes
import cCartogram as cart
from pycartogram import WardCartogram
import visvalingamwyatt as vw
import scipy.sparse as sprs
from scipy.spatial import ConvexHull


class PointCartogram(WardCartogram):
    """
    Create cartograms from point location data.

    Extends WardCartogram to work directly with point data. Density is
    computed from point counts in grid cells rather than predefined values.

    Parameters
    ----------
    points : array-like of shape (n, 2)
        Point coordinates as (x, y) pairs.
    wards : list of shapely.geometry.Polygon, optional
        Ward boundaries. If None, uses convex hull of points.
    norm_density : bool, optional
        Normalize density by area (default: False).
    margin_ratio : float, optional
        Margin as fraction of map size (default: 0.2).
    map_orientation : str, optional
        'landscape' or 'portrait' (default: 'landscape').
    x_raster_size : int, optional
        Grid width in pixels (default: 1024).
    y_raster_size : int, optional
        Grid height in pixels (default: 768).

    Attributes
    ----------
    points : numpy.ndarray
        Original point coordinates.
    new_points : numpy.ndarray
        Transformed point coordinates (after compute()).

    Examples
    --------
    >>> points = np.random.rand(1000, 2)
    >>> carto = PointCartogram(points)
    >>> carto.compute(verbose=True)
    >>> carto.plot_points()
    """

    def __init__(
        self,
        points: ArrayLike,
        wards: list[Polygon] | None = None,
        norm_density: bool = False,
        margin_ratio: float = 0.2,
        map_orientation: str = 'landscape',
        x_raster_size: int = 1024,
        y_raster_size: int = 768,
    ) -> None:
        """Initialize point cartogram with point data and optional wards."""
        self.no_wards_given = wards is None
        ward_density = None

        points = np.array(points,dtype=float)
        self.points = points

        if self.no_wards_given:
            hull = ConvexHull(self.points)
            whole_shape = Polygon(list(zip(self.points[hull.vertices,0],
                                      self.points[hull.vertices,1])))
            wards = [whole_shape]
            self.no_wards_given = False

        WardCartogram.__init__(self,
                               wards = wards,
                               ward_density = None,
                               margin_ratio = margin_ratio,
                               map_orientation = map_orientation,
                               x_raster_size = x_raster_size,
                               y_raster_size = y_raster_size,
                              )

    def _get_whole_shape_matrix(self, verbose: bool = False) -> NDArray[np.floating]:
        """Create binary matrix marking the whole shape area."""
        A = np.zeros((self.xsize, self.ysize))
        self._mark_matrix_with_shape(A, self.whole_shape)
        return A

    def cast_density_to_matrix(self, verbose: bool = False) -> NDArray[np.floating]:
        """
        Rasterize point locations to density matrix.

        Parameters
        ----------
        verbose : bool, optional
            Show progress (default: False).

        Returns
        -------
        numpy.ndarray
            Density matrix from point counts.
        """
        self.density_matrix = self.cast_points_to_matrix(self.points, verbose)
        return self.density_matrix

    def fast_density_to_matrix(self, verbose: bool = False):
        return cast_density_to_matrix(verbose=verbose)

    def compute(self, verbose: bool = False) -> None:
        """
        Compute the complete point cartogram transformation.

        Runs all steps: cast density, compute cartogram, transform wards,
        and transform points.

        Parameters
        ----------
        verbose : bool, optional
            Show progress (default: False).
        """
        self.cast_density_to_matrix(verbose)
        self.compute_cartogram(verbose=verbose)
        self.transform_wards(verbose)
        # Transform the point coordinates
        x, y = self.transform_coords(self.points[:, 0], self.points[:, 1])
        n_ = len(x)
        self.new_points = np.concatenate((x.reshape(n_, 1), y.reshape(n_, 1)), axis=1)

    def plot_points(
        self,
        ax: Axes | None = None,
        show_density_matrix: bool = False,
        show_new_points: bool = True,
        plot_wards: bool = True,
        ward_colors: Any = 'None',
        bg_color: Any = 'w',
        edge_colors: Any = None,
        outline_whole_shape: bool = True,
        use_new_density: bool = False,
        **kwargs: Any,
    ) -> tuple[Figure, Axes] | Axes:
        """
        Plot points on the cartogram.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates new figure.
        show_density_matrix : bool, optional
            Show density matrix as background (default: False).
        show_new_points : bool, optional
            If True, show transformed points; else original (default: True).
        plot_wards : bool, optional
            Show ward boundaries (default: True).
        ward_colors : color-like, optional
            Ward face colors (default: 'None' for transparent).
        bg_color : color-like, optional
            Background color (default: 'w').
        edge_colors : color-like, optional
            Ward edge colors.
        outline_whole_shape : bool, optional
            Draw outline around all wards (default: True).
        use_new_density : bool, optional
            Use transformed density (default: False).
        **kwargs : dict
            Additional arguments passed to ax.plot() for points.

        Returns
        -------
        fig, ax : matplotlib Figure and Axes
            Only if ax was None; otherwise returns just ax.
        """

        generate_figure = ax is None

        if (not plot_wards) or\
           (self.no_wards_given):
            #color = list(mpl.colors.to_rgba(bg_color))
            #color[-1] = 0.
            #print(color)
            temp = self.plot(
                 ax = ax,
                 show_density_matrix = show_density_matrix,
                 show_new_wards = show_new_points,
                 ward_colors = 'None',
                 bg_color = 'None',
                 edge_colors = 'None',
                 outline_whole_shape = False,
                 use_new_density = use_new_density,
             )
        else:
            temp = self.plot(
                 ax = ax,
                 show_density_matrix = show_density_matrix,
                 show_new_wards = show_new_points,
                 ward_colors = ward_colors,
                 bg_color = bg_color,
                 edge_colors = edge_colors,
                 outline_whole_shape = outline_whole_shape,
                 use_new_density = use_new_density,
             )

        if generate_figure:
            fig, ax = temp
        else:
            ax = temp

        if show_new_points:
            x, y = self.new_points[:,0], self.new_points[:,1]
        else:
            x, y = self.points[:,0], self.points[:,1]

        if 'marker' not in kwargs:
            kwargs['marker'] = '.'
        if ('markerfacecolor' not in kwargs) and\
           ('mfc' not in kwargs):
            kwargs['mfc'] = 'k'
        if ('markeredgecolor' not in kwargs) and\
           ('mec' not in kwargs):
            kwargs['mec'] = 'k'

        kwargs['ls'] = 'None'

        ax.plot(x,y,**kwargs)

        if generate_figure:
            return fig, ax
        else:
            return ax

if __name__ == "__main__":
        A = Polygon([
                      (0.,0.),
                      (1.,0.),
                      (1.,1.),
                      (0.,1.),
                      (0.,0.),
                    ])
        B = Polygon([
                      (1.,0.),
                      (2.,0.),
                      (2.,1.),
                      (1.,1.),
                      (1.,0.),
                    ])
        C = Polygon([
                      (2.,0.),
                      (3.,0.),
                      (3.,1.),
                      (2.,1.),
                      (2.,0.),
                    ])

        wards = [A,B,C]
        points = np.random.random((1000,2))

        transform = 'y'

        if transform == 'x':
            for i, x in enumerate(points[:,0]):
                u = np.log(1-x)/(-1./1.5)
                if u > 3.:
                    u = 3*np.random.rand()
                points[i,0] = u
        else:
            for i, y in enumerate(points[:,1]):
                u = np.log(1-y)/(-1./.1)
                if u > 1.:
                    u = 1*np.random.rand()
                points[i,1] = u
                points[i,0] *= 3

        carto = PointCartogram(points,
                               #wards=wards,
                               margin_ratio=0.5,
                               x_raster_size=128,
                               y_raster_size=64,
                              )

        density_matrix = carto.cast_density_to_matrix(True)
        #fig, ax = carto.plot_points(show_density_matrix=True,show_new_points=False,show_wards=True)
        fig, ax = carto.plot_points(show_new_points=False)
        #x,y = carto.get_ward_bounds(carto.big_bbox)
        #pl.imshow(density_matrix.T,extent=x+y,origin='lower')
        #carto.compute_cartogram((True))
        #pts = carto.points
        #x, y = pts[:,0], pts[:,1]
        #pl.plot(x,y,'.')
        carto.compute(verbose=True)
        fig, ax = carto.plot_points(show_density_matrix=True)
        #print density_matrix
        #print carto.big_bbox
        #fig = pl.figure()
        #pts = carto.new_points
        #x, y = pts[:,0], pts[:,1]
        #pl.plot(x,y,'.')
        #print(carto.points)
        #print(carto.new_points)
        #fig, ax = carto.plot_points()
        #ax.plot()

        """
        carto.plot_points(
                            show_new_points=False,
                            show_density_matrix = True,
                            plot_wards = False,
                         )
        fig, ax = carto.plot_points(
                    show_new_points=True,
                  )
        """
        pl.show()



