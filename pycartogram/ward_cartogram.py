"""
Ward-based cartogram generation using the Gastner-Newman algorithm.

This module provides the WardCartogram class for creating density-equalizing
cartograms from polygon regions (wards) with associated density values.

A cartogram is a map where regions are distorted so their area becomes
proportional to a variable of interest (e.g., population, GDP).
"""

import numpy as np
import shapely.geometry as sgeom
from shapely.ops import unary_union, polygonize
from shapely.geometry import Polygon, LineString, MultiPolygon, Point
from pycartogram.tools import (
    polygon_patch,
    add_intersection_points_to_wards,
    coarse_grain_wards,
    enrich_polygon_with_points,
    scale,
    logify_and_scale,
    is_iter,
)
import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib.path import Path
import progressbar
import cCartogram as cart
from tqdm import tqdm

class WardCartogram:
    """
    Create density-equalizing cartograms from polygon regions.

    This class implements the Gastner-Newman diffusion-based cartogram
    algorithm. Regions (wards) are distorted so their areas become
    proportional to specified density values.

    Parameters
    ----------
    wards : list of shapely.geometry.Polygon
        Polygon regions to transform.
    ward_density : list of float
        Density value for each ward (e.g., population, GDP).
    norm_density : bool, optional
        If True, normalize densities by ward area (default: False).
    margin_ratio : float, optional
        Margin around wards as fraction of map size (default: 0.2).
    map_orientation : str, optional
        'landscape' or 'portrait' (default: 'landscape').
    x_raster_size : int, optional
        Width of computation grid in pixels (default: 1024).
    y_raster_size : int, optional
        Height of computation grid in pixels (default: 768).
    threshold_for_polygon_coarse_graining : float, optional
        If set, simplify polygons using Visvalingam-Whyatt algorithm.

    Attributes
    ----------
    wards : list of Polygon
        Input ward polygons.
    ward_density : list of float
        Density values per ward.
    new_wards : list of Polygon
        Transformed ward polygons (after calling compute()).
    density_matrix : numpy.ndarray
        Rasterized density grid.
    cartogram : object
        Computed cartogram transformation.

    Examples
    --------
    >>> from shapely.geometry import Polygon
    >>> wards = [Polygon([(0,0), (1,0), (1,1), (0,1)]),
    ...          Polygon([(1,0), (2,0), (2,1), (1,1)])]
    >>> density = [100, 200]  # second ward has twice the density
    >>> carto = WardCartogram(wards, density, norm_density=True)
    >>> carto.compute(verbose=True)
    >>> carto.plot()
    """

    def __init__(self,
                 wards,
                 ward_density,
                 norm_density=False,
                 margin_ratio=0.2,
                 map_orientation='landscape',
                 x_raster_size=1024,
                 y_raster_size=768,
                 threshold_for_polygon_coarse_graining=None):
        """Initialize the cartogram with wards and density values."""

        self.ward_density = ward_density

        add_intersection_points_to_wards(wards)

        if threshold_for_polygon_coarse_graining is not None:
            wards = coarse_grain_wards(wards, threshold_for_polygon_coarse_graining)

        self.wards = wards

        if norm_density:
            self.ward_density = list(ward_density)
            for i, w in enumerate(self.wards):
                self.ward_density[i] /= w.area

        self._process_wards()

        self.compute_bounding_box(
            margin_ratio = margin_ratio,
            map_orientation = map_orientation,
            x_raster_size = int(x_raster_size),
            y_raster_size = int(y_raster_size),
            )

        self.density_matrix = None
        self.cartogram = None
        self.ward_matrix_coordinates = None

    def enrich_wards_with_points(self, delta_for_enrichment=None):
        """
        Add interpolated points to ward boundaries.

        Increases vertex density for smoother cartogram transformations.

        Parameters
        ----------
        delta_for_enrichment : float, optional
            Maximum distance between vertices. Defaults to tile_size.
        """
        if delta_for_enrichment is None:
            delta_for_enrichment = self.tile_size

        self.wards = [enrich_polygon_with_points(ward, delta_for_enrichment) for ward in self.wards]

    def _mark_matrix_with_shape(self, A_, shape, new_val=1., old_val=None):
        """
        Rasterize a polygon onto a matrix.

        Uses matplotlib.path.contains_points for fast vectorized computation.

        Parameters
        ----------
        A_ : numpy.ndarray
            Matrix to mark (modified in place).
        shape : shapely.geometry.Polygon
            Polygon to rasterize.
        new_val : float, optional
            Value to set inside polygon (default: 1.0).
        old_val : float, optional
            Unused, kept for API compatibility.
        """
        # Create grid of all matrix coordinates
        i_coords = np.arange(A_.shape[0])
        j_coords = np.arange(A_.shape[1])
        ii, jj = np.meshgrid(i_coords, j_coords, indexing='ij')

        # Convert matrix indices to world coordinates
        x_coords = self.x_from_i(ii.ravel())
        y_coords = self.y_from_j(jj.ravel())
        points = np.column_stack([x_coords, y_coords])

        # Create matplotlib Path from polygon exterior
        poly_coords = np.array(shape.exterior.coords)
        path = Path(poly_coords)

        # Check which points are inside (vectorized, fast)
        inside = path.contains_points(points).reshape(A_.shape)

        # Set values
        A_[inside] = new_val

    def compute_ward_density_from_locations(self, locations, verbose=False):
        """
        Compute ward densities from point locations.

        Counts points within each ward and normalizes by area.

        Parameters
        ----------
        locations : array-like or list of shapely.geometry.Point
            Point coordinates as (x, y) pairs or Point objects.
        verbose : bool, optional
            Show progress bar (default: False).

        Returns
        -------
        numpy.ndarray
            Computed density for each ward.
        """
        coarse_matrix = self.cast_points_to_matrix(locations, replace_value_zero=False, verbose=verbose)
        i, j = np.nonzero(coarse_matrix)

        hist = np.zeros((len(self.wards),),dtype=float)

        # those are the logarithmically distributed time points
        # at which the wards are sorted for relevance
        sort_events = np.unique(\
                               np.array(
                               np.logspace(2,
                                           np.log(len(i)-1)/np.log(10),
                                           100,
                                          ),
                               dtype=int
                               )
                              )

        # start with the first sort event (at point_id == 100)
        next_sort_event = 0

        # in the beginning each ward is equally important
        ward_order = np.arange(len(self.wards),dtype=int)

        # iterate through locations
        if verbose:
            bar = progressbar.ProgressBar(
                        max_value = len(i)-1,
                        widgets = [
                            progressbar.SimpleProgress()," ",
                            progressbar.ETA()," computing ward densities ..."
                        ]
                )

        for point_id, (i_,j_) in enumerate(zip(i,j)):

            P = Point(self.x_from_i(i_),self.y_from_j(j_))
            visit_count = coarse_matrix[i_,j_]
            # this would've been a check whether or not the location
            # is actually in Berlin but ultimately this was too expensive
            #if not berlin_polygon.contains(P):
            #    continue

            # go trough the ward_order containing the most relevant wards
            for trial, ward_id in enumerate(ward_order):
                # get the relevant polygon
                geom = self.wards[ward_id]

                # check whether the point is in the polygon
                if geom.contains(P):
                    hist[ward_id] += visit_count
                    break # if that's the case stop the search for the point

            #if trial == len(ward_order)-1:
            #    print "found outlier"

            # check if we should sort the ward order for relevance
            if next_sort_event < len(sort_events) and point_id == sort_events[next_sort_event]:
                next_sort_event += 1
                # the ward with the highest probability to contain a location will
                # now be at the beginning of the list
                ward_order = np.argsort(-hist,kind='mergesort')
                # print "sort_event at point_id =", point_id

            if verbose:
                bar.update(point_id)
        # NORM
        self.ward_density = np.array(hist)
        for i, w in enumerate(self.wards):
            self.ward_density[i] /= w.area

        return self.ward_density


    def compute_bounding_box(self,
                             margin_ratio=0.1,
                             map_orientation='landscape',
                             x_raster_size=1024,
                             y_raster_size=768):
        """
        Compute the bounding box and coordinate transforms for the raster grid.

        Sets up the mapping between world coordinates and matrix indices.

        Parameters
        ----------
        margin_ratio : float, optional
            Margin as fraction of map extent (default: 0.1).
        map_orientation : str, optional
            'landscape' or 'portrait' (default: 'landscape').
        x_raster_size : int, optional
            Grid width in pixels (default: 1024).
        y_raster_size : int, optional
            Grid height in pixels (default: 768).
        """
        self.xsize = x_raster_size
        self.ysize = y_raster_size

        use_width = map_orientation == 'landscape'

        x_, y_ = self.get_ward_bounds()

        if use_width:
            margin_x = self.orig_width * margin_ratio
            self.new_width = self.orig_width + 2*margin_x
            self.tile_size = self.new_width / self.xsize

            self.new_height = self.tile_size * self.ysize
            margin_y = 0.5 * (self.new_height - self.orig_height)

            if margin_y < 0:
                raise ValueError("New height is smaller than old height, please increase y_raster_size")
        else:
            margin_y = self.orig_height * margin_ratio
            self.new_height = self.orig_height + 2*margin_y

            self.tile_size = self.new_height / self.ysize

            self.new_width = self.tile_size * self.xsize
            margin_x = 0.5 * (self.new_width - self.orig_width)

            if margin_x < 0:
                raise ValueError("New width is smaller than old width, please increase x_raster_size")

        x_ = x_[0] - margin_x, x_[1] + margin_x
        y_ = y_[0] - margin_y, y_[1] + margin_y

        self.big_bbox = Polygon([(x_[0],y_[0]),
                                 (x_[1],y_[0]),
                                 (x_[1],y_[1]),
                                 (x_[0],y_[1]),
                                ])

        size_x = self.tile_size
        size_y = self.tile_size

        self.x_from_i = lambda i: i * size_x + size_x / 2. + x_[0]
        self.y_from_j = lambda j: j * size_y + size_y / 2. + y_[0]
        self.i_from_x = lambda x: (x - size_x / 2. - x_[0] ) / size_x
        self.j_from_y = lambda y: (y - size_y / 2. - y_[0] ) / size_y

    def _get_matrix_coordinate_bounds(self,ward):
        x0 = self.big_bbox.bounds[0]
        y0 = self.big_bbox.bounds[1]
        x_ = (ward.bounds[0]-x0,ward.bounds[2]-x0)
        y_ = (ward.bounds[1]-y0,ward.bounds[3]-y0)
        imin = int(np.floor(x_[0]/self.tile_size))
        imax = int(np.ceil(x_[1]/self.tile_size))
        jmin = int(np.floor(y_[0]/self.tile_size))
        jmax = int(np.ceil(y_[1]/self.tile_size))
        return imin, imax, jmin, jmax

    def get_ward_bounds(self,shape=None):

        if shape is None:
            shape = self.whole_shape

        x_ = (
              shape.bounds[0],
              shape.bounds[2]
             )
        y_ = (shape.bounds[1],
              shape.bounds[3]
             )

        return x_, y_

    def _process_wards(self):
        self.whole_shape = unary_union(self.wards)
        x_, y_ = self.get_ward_bounds()

        self.orig_bbox = Polygon([(x_[0],y_[0]),
                                  (x_[1],y_[0]),
                                  (x_[1],y_[1]),
                                  (x_[0],y_[1]),
                                 ])
        self.orig_width = x_[1] - x_[0]
        self.orig_height = y_[1] - y_[0]

    def fast_density_to_matrix(self, verbose=False, set_boundary_to='mean', **kwargs):
        """
        Fast rasterization of ward densities to matrix using polygon filling.

        Parameters
        ----------
        verbose : bool, optional
            Show progress (default: False).
        set_boundary_to : str or float, optional
            Value for areas outside wards: 'mean', 'min', or numeric.
        **kwargs : dict
            Additional arguments (unused).

        Returns
        -------
        numpy.ndarray
            Density matrix of shape (xsize, ysize).
        """
        ward_dens = np.array(self.ward_density)

        min_density = np.min(ward_dens[ward_dens>0.])
        mean_density = np.mean(ward_dens[ward_dens>0.])
        offset_density = min_density / 10.
        #print("offset density", offset_density)
        ward_dens[ward_dens==0.] = offset_density
        #density = mean_density * np.ones((self.xsize,self.ysize),dtype=float)
        density = np.zeros((self.xsize,self.ysize),dtype=float)

        for iward, ward in enumerate(self.wards):
            self._mark_matrix_with_shape(density,ward,new_val=ward_dens[iward],old_val=0.)

        density[np.where(density==0.)] = mean_density
        self.density_matrix = density
        return self.density_matrix

    def cast_density_to_matrix(self,verbose=False,set_boundary_to='mean',ward_matrix_coordinates=[],**kwargs):

        ward_dens =  np.array(self.ward_density)
        nnz = ward_dens[ward_dens>0.]
        min_density = np.min(nnz)
        mean_density = np.mean(nnz)
        ratio_of_nonzero_values = len(nnz) / float(len(ward_dens))
        offset_density = ratio_of_nonzero_values * min_density
        #print("number of nnz density values", len(nnz))
        #print("ratio of nnz density values", ratio_of_nonzero_values)
        #print( offset_density
        #if len(nnz) < len(ward_dens) / 10:
        #    offset_density = 0.005*min_density
        #else:
        #    offset_density = min_density
        #print("offset density", offset_density)
        ward_dens[ward_dens==0.] = offset_density
        density = np.zeros((self.xsize,self.ysize),dtype=float)

        if verbose:
            #print("casting ward density to discrete matrix values")
            bar = progressbar.ProgressBar(
                max_value = len(self.wards) - 1,
                widgets = [
                    progressbar.SimpleProgress()," ",
                    progressbar.ETA(),
                    " cast density to matrix ...",
                ]
            )
        #offset_i = []
        #offset_j = []
        if len(ward_matrix_coordinates) == len(self.wards):
            self.ward_matrix_coordinates = ward_matrix_coordinates

        no_matrix_coordinates_given = self.ward_matrix_coordinates is None

        if no_matrix_coordinates_given:

            self.ward_matrix_coordinates = [[] for ward in self.wards ]

            for iward, ward in enumerate(self.wards):

                imin,imax,jmin,jmax = self._get_matrix_coordinate_bounds(ward)

                #if verbose:
                #    print "ward", iward, "has matrix bounds", imin, imax, jmin, jmax

                for i in range(imin,imax):
                    for j in range(jmin,jmax):

                        this_point = sgeom.Point(self.x_from_i(i),
                                                 self.y_from_j(j))

                        if ward.contains(this_point):
                            density[i,j] = ward_dens[iward]
                            self.ward_matrix_coordinates[iward].append((i,j))

                            #if self.ward_density[iward] == 0.:
                            #    offset_i.append(i)
                            #    offset_j.append(j)
                            #    #print("set offset density", offset_density)
                if verbose:
                    bar.update(iward)
        else:
            for iward, ward in enumerate(self.wards):
                for i, j in self.ward_matrix_coordinates[iward]:
                    density[i,j] = ward_dens[iward]

        #offset_i = np.array(offset_i)
        #offset_j = np.array(offset_j)

        #mean_density = np.mean(density[density>0.])
        #print("mean_density",mean_density)
        #density[(offset_i,offset_j)] = offset_density
        #print("density at first offset", density[offset_i[0],offset_j[0]])

        #print("min_density", np.amin(density))
        if isinstance(set_boundary_to,str):
            if set_boundary_to == 'mean':
                density[np.where(density==0.)] = mean_density
            elif set_boundary_to == 'min':
                density[np.where(density==0.)] = min_density
            elif set_boundary_to.endswith('percent_mean'):
                factor = float(set_boundary_to.replace('percent_mean','')) / 100
                density[np.where(density==0.)] = factor * mean_density
            elif set_boundary_to.endswith('percent_min'):
                factor = float(set_boundary_to.replace('percent_min','')) / 100
                density[np.where(density==0.)] = factor * min_density
        else:
            density[np.where(density==0.)] = set_boundary_to
        #print("min_density", np.amin(density))

        self.density_matrix = density
        return self.density_matrix

    def cast_points_to_matrix(self,
                              locations,
                              verbose = False,
                              replace_value_zero = True,
                              ):

        if type(locations[0]) == Point:
            # convert to array
            locations = np.array([ (P.x, P.y) for P in locations])
        else:
            # get copy
            locations = np.array(locations)

        size_x = size_y = self.tile_size
        xmin = self.big_bbox.bounds[0]
        ymin = self.big_bbox.bounds[1]

        # sparse matrix containing a 1 if a grid was visited
        density = np.zeros((self.xsize,self.ysize))

        # fill matrix with points visited
        if verbose:
            bar = progressbar.ProgressBar(
                max_value = locations.shape[0] - 1,
                widgets = [ progressbar.SimpleProgress()," ",
                            progressbar.ETA(),
                            " casting points to matrix..."
                          ]
            )

        for point_id in range(locations.shape[0]):
            # show progress
            x, y = locations[point_id]
            x_id = int(np.floor((x - xmin) / size_x))
            y_id = int(np.floor((y - ymin) / size_y))
            if x_id>=0 and x_id < self.xsize and\
               y_id>=0 and y_id < self.ysize:
                density[x_id, y_id] += 1
                if verbose:
                    bar.update(point_id)

        if replace_value_zero:
            mean_density = np.mean(density[density>0.])

            min_density = np.min(density[density>0.])
            points_within_shape = self._get_whole_shape_matrix()
            offset_density = min_density / 10.
            density[np.where(np.logical_and(density==0.,points_within_shape==1.))] = offset_density

            density[density==0.] = mean_density

        return density

    def transform_coords(self,x,y):
        old_x = x
        old_y = y
        i_ = self.i_from_x(np.array(old_x))
        j_ = self.j_from_y(np.array(old_y))
        new_ij = cart.remap_coordinates(list(zip(i_,j_)),
                                        self.cartogram,
                                        self.xsize,
                                        self.ysize,
                                        )
        new_x = np.array([self.x_from_i(ij[0]) for ij in new_ij ])
        new_y = np.array([self.y_from_j(ij[1]) for ij in new_ij ])

        return new_x, new_y

    def compute_cartogram(self,offset=0.005,blur=0.,verbose=False,**kwargs):

        if verbose:
            print()
            print("computing cartogram...")
        self.cartogram = cart.compute_cartogram(
                                            self.density_matrix.tolist(),
                                            offset = offset,
                                            blur = blur,
                                            show_progress = verbose,
                                            )
        return self.cartogram

    def transform_wards(self,
                        verbose=False,
                        ignore_self_intersection=True,
                        enrich_wards_with_points=True,
                        delta_for_enrichment=None):

        if delta_for_enrichment is None:
            delta_for_enrichment = self.tile_size

        new_wards = []
        new_ward_density = []

        self.new_ward_coords = []
        self.old_ward_coords = []
        self.new_old_wards = []

        if verbose:
            bar = progressbar.ProgressBar(
                max_value = len(self.wards),
                widgets = [progressbar.SimpleProgress()," transforming wards to new coordinates..."]
            )

        for iward,ward in enumerate(self.wards):
            if enrich_wards_with_points:
                temp_ward = enrich_polygon_with_points(ward,delta_for_enrichment)
                old_x, old_y = temp_ward.exterior.coords.xy
                self.new_old_wards.append(temp_ward)
            else:
                old_x, old_y = ward.exterior.coords.xy
                self.new_old_wards.append(ward)

            self.old_ward_coords.append(list(zip(old_x,old_y)))
            #number_of_points_old = len(old_x)
            #new_number_of_points = 3*number_of_points_old
            #old_x = []
            #old_y = []
            #old_boundary = ward.boundary
            #for val in np.linspace(0,1,new_number_of_points):
            #    p = ward.boundary.interpolate(val, normalized=True)
            #    old_x.append(p.x)
            #    old_y.append(p.y)
            #new_boundary = old_boundary.union(LineString(zip(old_x,old_y)))
            #new_ward = polygonize(new_boundary)
            #print([g for g in new_ward])
            #old_x, old_y = new_ward.coords.xy

            i_ = self.i_from_x(np.array(old_x))
            j_ = self.j_from_y(np.array(old_y))

            new_ij = cart.remap_coordinates(list(zip(i_,j_)),self.cartogram,self.xsize,self.ysize)
            new_coords = [ (self.x_from_i(i), self.y_from_j(j)) for i,j in new_ij]
            self.new_ward_coords.append(new_coords)
            """
            if iward == 2:
                x = np.array([ 0.38517325,  0.40859912,  0.43296919,  0.4583215 ,  0.4583215 ,
                               0.43296919,  0.40859912,  0.38517325,  0.36265506,  0.34100929])
                x /= x.max()*3
                y = np.array([ 62.5       ,  56.17977528,  39.39698492,   0.        ,
                               0.        ,  17.34605377,  39.13341671,  60.4180932 ,
                               76.02574417,  85.47008547])
                y /= y.max()
                new_coords = zip(x,y)
            """
            ls = LineString(new_coords)
            lr = LineString(ls.coords[:] + ls.coords[:1])
            new_ward = Polygon(new_coords)
            #if lr.is_simple:
            #    new_ward = Polygon(new_coords)
            #    #new_ward = new_ward.buffer(0)
            #else:
            #    mls = unary_union(lr)
            #    mp = MultiPolygon(list(polygonize(mls)))
            #    new_ward = max(multipolygon, key=lambda a: a.area)
            #    #new_ward = unary_union(mp)
            #    ##new_ward = new_ward.buffer(0)
            #    ##for new_ward in polygonize(mls):

            #if type(new_ward) == MultiPolygon:
            #    new_ward = new_ward.buffer(self.tile_size/250.)

            new_wards.append(new_ward)


            if verbose:
                bar.update(iward)

        self.new_ward_density = np.array(new_ward_density)
        self.new_wards = new_wards
        #self.new_whole_shape = cascaded_union(new_wards)

        # this is such a dirty hack
        # (this is to prevent a really strange bug in unary_union (or cascaded union)
        # where putting all wards at once raises a TopologyException
        try:
            size = 300
            chunks = [self.new_wards[size*i:size*(i+1)] for i in range(len(self.new_wards)//size + 1)]
            poly1 = unary_union(chunks[0])
            for i,poly2 in enumerate(chunks[1:]):
                poly1 = unary_union([poly1]+poly2)
        except ValueError as e:
            if verbose:
                print('method produced error')
                print(e)
                print('will fall back to slower method')
                bar = progressbar.ProgressBar(
                        max_value = len(self.new_wards),
                        widgets = [
                                progressbar.SimpleProgress()," ",
                                progressbar.ETA()," joining wards to whole shape ...",
                            ]
                        )

            poly1 = self.new_wards[0]
            for i, poly2 in enumerate(self.new_wards[1:]):
                poly1 = poly1.union(poly2)
                if verbose:
                    bar.update(i)
        # dirty hack ends here

        # delete holes (construct new Polygon from exterior)
        if isinstance(poly1, Polygon):
            poly1 = Polygon(poly1.exterior.coords)

        self.new_whole_shape = poly1
        x_, y_ = self.get_ward_bounds(self.new_whole_shape)
        self.new_bbox = Polygon([(x_[0],y_[0]),
                                 (x_[1],y_[0]),
                                 (x_[1],y_[1]),
                                 (x_[0],y_[1]),
                                 ])

        return self.new_wards

    def compute(self, verbose=False, **kwargs):
        """
        Compute the complete cartogram transformation.

        Convenience method that runs all steps: cast density to matrix,
        compute cartogram, and transform wards.

        Parameters
        ----------
        verbose : bool, optional
            Show progress (default: False).
        **kwargs : dict
            Additional arguments passed to sub-methods.

        Returns
        -------
        list of shapely.geometry.Polygon
            Transformed ward polygons.
        """
        if self.density_matrix is None:
            self.cast_density_to_matrix(verbose)
        if self.cartogram is None:
            self.compute_cartogram(verbose=verbose)
        return self.transform_wards(verbose)

    def _convert_wards_to_list(self,wards):
        """return list of wards if a dictionary was supplied"""
        pass

    def plot(self,
             ax=None,
             show_density_matrix=False,
             show_new_wards=True,
             ward_colors=None,
             bg_color='w',
             edge_colors=None,
             outline_whole_shape=True,
             whole_shape_color=[0.05, 0., 0.],
             whole_shape_linewidth=1,
             use_new_density=False,
             intensity_range=[0.1, 0.9]):
        """
        Plot the cartogram.

        Parameters
        ----------
        ax : matplotlib.axes.Axes, optional
            Axes to plot on. If None, creates new figure.
        show_density_matrix : bool, optional
            Show density matrix as background image (default: False).
        show_new_wards : bool, optional
            If True, show transformed wards; if False, show original (default: True).
        ward_colors : str, list, or color-like, optional
            Ward face colors. Can be:
            - 'density': color by density values
            - 'log_density': color by log of density
            - list of colors: one per ward
            - single color: all wards same color
        bg_color : color-like, optional
            Background color (default: 'w').
        edge_colors : color-like or list, optional
            Ward edge colors (default: same as face).
        outline_whole_shape : bool, optional
            Draw outline around all wards (default: True).
        whole_shape_color : color-like, optional
            Color for whole shape outline.
        whole_shape_linewidth : float, optional
            Line width for outline (default: 1).
        use_new_density : bool, optional
            Use transformed density values (default: False).
        intensity_range : list of float, optional
            [min, max] for color intensity mapping (default: [0.1, 0.9]).

        Returns
        -------
        fig, ax : matplotlib Figure and Axes
            Only if ax was None; otherwise returns just ax.
        """

        generate_figure = ax is None

        if generate_figure:
            fig, ax = pl.subplots(1,1)

        if show_new_wards:
            wards = self.new_wards
            bbox = self.new_bbox
            whole = self.new_whole_shape
            density = self.new_ward_density
        else:
            wards = self.wards
            bbox = self.orig_bbox
            whole = self.whole_shape
            density = self.ward_density


        if ward_colors is None and not show_density_matrix:
            ward_colors = 'log_density'
        elif show_density_matrix:
            ward_colors = list(mpl.colors.to_rgba(bg_color))
            ward_colors[-1] = 0.
            bg_color = list(mpl.colors.to_rgba(bg_color))
            bg_color[-1] = 0.

        if mpl.colors.is_color_like(ward_colors):
            color = lambda iward: ward_colors
        elif isinstance(ward_colors, str) and ward_colors in ['density', 'log_density']:
            if use_new_density:
                density = self.new_ward_density
            else:
                density = self.ward_density
            if ward_colors == 'density':
                values, intensity = scale(self.ward_density,new_min=intensity_range[0],new_max=intensity_range[1])
            else:
                values, intensity = logify_and_scale(self.ward_density,new_min=intensity_range[0],new_max=intensity_range[1])

            if bg_color in ('w','white'):
                color = lambda iward: [1, 1-intensity(values[iward]), 1-intensity(values[iward])]
            else:
                color = lambda iward: [intensity(values[iward]), 0, 0 ]
        elif is_iter(ward_colors):
            color = lambda iward: ward_colors[iward]

        if edge_colors is None:
            edge_color = color
        elif mpl.colors.is_color_like(edge_colors):
            edge_color = lambda iward: edge_colors
        elif is_iter(edge_colors):
            edge_color = lambda iward: edge_colors[iward]

        # set background patch
        patch = polygon_patch(self.big_bbox,
                             facecolor = bg_color,
                             edgecolor = 'None',
                             #alpha = 1,
                             lw = 0,
                            )
        ax.add_patch(patch)

        # get bounds of map
        x_, y_ = self.get_ward_bounds(bbox)

        if show_density_matrix:
            x_b, y_b = self.get_ward_bounds(self.big_bbox)
            ax.imshow(self.density_matrix.T,
                      extent=x_b+y_b,
                      origin='lower')
        else:
            ax.set_aspect('equal')

        # plot every ward
        for ward_id, ward in enumerate(wards):
            fc = color(ward_id)
            patch = polygon_patch(ward,
                                 facecolor = fc,
                                 edgecolor = edge_color(ward_id),
                                 #alpha = 1,
                                 lw = 0.5,
                                )
            ax.add_patch(patch)

        # set whole shape background
        if outline_whole_shape:
            patch = polygon_patch(whole,
                                 facecolor = 'None',
                                 edgecolor = whole_shape_color,
                                 alpha = 1,
                                 lw = whole_shape_linewidth,
                                )
            ax.add_patch(patch)

        ax.set_xlim(x_)
        ax.set_ylim(y_)

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
        ward_population = [0.,2.,200.]

        points = np.random.random((1000,2))
        points[:,0] *= 3

        carto = WardCartogram(wards=wards,
                              ward_density=ward_population,
                              norm_density=True, # needs to be normed by area
                              margin_ratio=0.5,
                              x_raster_size=128,
                              y_raster_size=64)

        carto.compute_ward_density_from_locations(points,True)

        density_matrix = carto.cast_density_to_matrix(verbose=True)
        carto.compute_cartogram((True))
        carto.transform_wards((True))
        #carto.compute(verbose=True)
        #print density_matrix
        #print carto.big_bbox


        carto.plot(show_new_wards=False,
                          show_density_matrix=True
                   #ward_colors = [0.3,0.5,0.8],
                   #edge_colors = 'k',
                   #bg_color = 'k',
                          )
        fig, ax = carto.plot(show_new_wards=True,
                   edge_colors = 'k'
                  #show_density_matrix=True
                  )

        pl.show()



