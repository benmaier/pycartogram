from __future__ import print_function
import numpy as np
import shapely.geometry as sgeom
from shapely.ops import cascaded_union
from shapely.geometry import Polygon
from descartes.patch import PolygonPatch
import matplotlib as mpl
import matplotlib.pyplot as pl
import progressbar
import cCartogram as cart
from pycartogram.tools import *
import visvalingamwyatt as vw

class WardCartogram():

    def __init__(self,
                 wards,
                 ward_density,
                 norm_density=False,
                 margin_ratio = 0.2,
                 dominant_dimension = 'x',
                 x_raster_size = 1024,
                 y_raster_size = 768,
                 threshold_for_polygon_coarse_graining = None,
                 ):
        """
        wards:        list of wards as shapely.Polygon
        ward_density: list of density values for the wards
        norm_density: if this value is 'True' the ward_density will be normed by area (default: False)
        """

        if threshold_for_polygon_coarse_graining is not None:
            self.wards = self._coarse_grain_wards(wards,threshold_for_polygon_coarse_graining)
        else:
            self.wards = wards
        self.ward_density = ward_density

        if norm_density:
            self.ward_density = list(ward_density)
            for i, w in enumerate(self.wards):
                self.ward_density[i] /= w.area

        self._process_wards()

        self.compute_bounding_box(
            margin_ratio = margin_ratio,
            dominant_dimension = dominant_dimension,
            x_raster_size = x_raster_size,
            y_raster_size = y_raster_size,
            )

        self.density_matrix = None

    def _coarse_grain_wards(self,wards,th):
        new_wards = []
        for ward in wards:
            coo = ward.exterior.coords.xy
            new_coo = vw.Simplifier(zip(*coo)).simplify(threshold=th)
            new_ward = Polygon(new_coo)
            new_wards.append(new_ward)
        return new_wards

    def compute_bounding_box(
            self,
            margin_ratio = 0.1,
            dominant_dimension = 'x',
            x_raster_size = 1024,
            y_raster_size = 768,
            ):

        self.xsize = x_raster_size
        self.ysize = y_raster_size

        use_width = dominant_dimension == 'x'

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

            self.tile_size = self.new_width / self.ysize

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
        self.whole_shape = cascaded_union(self.wards)
        x_, y_ = self.get_ward_bounds()

        self.orig_bbox = Polygon([(x_[0],y_[0]),
                                  (x_[1],y_[0]),
                                  (x_[1],y_[1]),
                                  (x_[0],y_[1]),
                                 ])
        self.orig_width = x_[1] - x_[0]
        self.orig_height = y_[1] - y_[0]
        
    def cast_density_to_matrix(self,verbose=False):

        min_density = np.min(self.ward_density[self.ward_density>0.])
        offset_density = min_density / 10.
        density = np.zeros((self.xsize,self.ysize),dtype=float)

        if verbose:
            print("casting ward density to discrete matrix values")
            bar = progressbar.ProgressBar(
                max_value = len(self.wards) - 1,
                widgets = [
                    progressbar.SimpleProgress()," ",
                    progressbar.ETA(),
                ]
            )

        for iward, ward in enumerate(self.wards):

            imin,imax,jmin,jmax = self._get_matrix_coordinate_bounds(ward)

            #if verbose:
            #    print "ward", iward, "has matrix bounds", imin, imax, jmin, jmax

            for i in range(imin,imax):
                for j in range(jmin,jmax):

                    this_point = sgeom.Point(self.x_from_i(i),
                                             self.y_from_j(j))

                    if ward.contains(this_point):
                        density[i,j] = self.ward_density[iward]

                        if self.ward_density[iward] == 0.:
                            density[i,j] = offset_density                      
                            #print("set offset density", offset_density)
            if verbose:
                bar.update(iward)

        mean_density = np.mean(density[density>0.])
        density[density==0.] = mean_density

        self.density_matrix = density
        return self.density_matrix

    def transform_coords(self,x,y):
        old_x = x
        old_y = y
        i_ = self.i_from_x(np.array(old_x))
        j_ = self.j_from_y(np.array(old_y)) 
        new_ij = cart.remap_coordinates(zip(i_,j_),
                                        self.cartogram,
                                        self.xsize,
                                        self.ysize,
                                        )
        new_x = np.array([self.x_from_i(ij[0]) for ij in new_ij ])
        new_y = np.array([self.y_from_j(ij[1]) for ij in new_ij ])
            
        return new_x, new_y

    def compute_cartogram(self,offset=0.005,blur=0.,verbose=False):

        if verbose:
            print("computing cartogram...")
        self.cartogram = cart.compute_cartogram(
                                            self.density_matrix.tolist(),
                                            offset = offset,
                                            blur = blur,
                                            show_progress = verbose,
                                            )
        return self.cartogram

    def transform_wards(self,verbose=False,ignore_self_intersection=True):
        new_wards = []
        new_ward_density = [] 

        if verbose:
            print("transforming wards to new coordinates...")
            bar = progressbar.ProgressBar(
                max_value = len(self.wards),
                widgets = [progressbar.SimpleProgress()]
            )

        for iward,ward in enumerate(self.wards):
            old_x, old_y = ward.exterior.coords.xy
            i_ = self.i_from_x(np.array(old_x))
            j_ = self.j_from_y(np.array(old_y)) 
            new_ij = cart.remap_coordinates(zip(i_,j_),self.cartogram,self.xsize,self.ysize)
            new_coords = [ (self.x_from_i(i), self.y_from_j(j)) for i,j in new_ij]
            new_ward = Polygon(new_coords)
            #try:
            #except TopologyException:
            if ignore_self_intersection and not new_ward.is_valid:
                new_ward = new_ward.buffer(0)

            new_wards.append(new_ward)
            if self.ward_density is not None:
                new_ward_density.append(self.ward_density[iward] * ward.area / new_ward.area)

            if verbose:
                bar.update(iward)
            
        self.new_ward_density = np.array(new_ward_density)
        self.new_wards = new_wards
        self.new_whole_shape = cascaded_union(new_wards)
        x_, y_ = self.get_ward_bounds(self.new_whole_shape)
        self.new_bbox = Polygon([(x_[0],y_[0]),
                                 (x_[1],y_[0]),
                                 (x_[1],y_[1]),
                                 (x_[0],y_[1]),
                                 ])

        return self.new_wards

    def compute(self,verbose=False):
        """do everything after init"""
        self.cast_density_to_matrix(verbose)
        self.compute_cartogram(verbose=verbose)
        return self.transform_wards(verbose)

    def _convert_wards_to_list(self,wards):
        """return list of wards if a dictionary was supplied"""
        pass

    def plot(self,
             ax = None,                    
             show_density_matrix = False,
             show_new_wards = True,
             ward_colors = None,
             bg_color = 'w',
             edge_colors = None,
             outline_whole_shape = True,
             use_new_density = False,
            ):
        """
            ward_colors can be
                - an iterable container of colors - plot ward face colors according to a provided list
                - "density" - plot ward face colors proportional to density
                - "log_density" - plot ward face colors proportional to density
        """

        generate_figure = ax is None

        if generate_figure:
            fig, ax = pl.subplots(1,1)

        if show_new_wards:
            wards = self.new_wards
            bbox = self.new_bbox
            whole = self.new_whole_shape
        else:
            wards = self.wards
            bbox = self.orig_bbox
            whole = self.whole_shape

        if ward_colors is None and not show_density_matrix:
            ward_colors = 'log_density'
        elif show_density_matrix:
            ward_colors = list(mpl.colors.to_rgba(bg_color))
            ward_colors[-1] = 0.
            bg_color = list(mpl.colors.to_rgba(bg_color))
            bg_color[-1] = 0.

        if mpl.colors.is_color_like(ward_colors):
            color = lambda iward: ward_colors
        elif ward_colors in ['density', 'log_density']:
            if use_new_density:
                density = self.new_ward_density
            else:
                density = self.ward_density
            if ward_colors == 'density':
                values, intensity = scale(self.ward_density)
            else:
                values, intensity = logify_and_scale(self.ward_density)

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
        patch = PolygonPatch(bbox,
                             facecolor = bg_color,
                             edgecolor = 'None',
                             #alpha = 1,
                             lw = 0,
                            )
        ax.add_patch(patch)
            
        # set whole shape background
        if outline_whole_shape:
            patch = PolygonPatch(whole,
                                 facecolor = 'None',
                                 edgecolor = [0.25,0,0],
                                 alpha = 1,
                                 lw = 2,
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
            patch = PolygonPatch(ward,
                                 facecolor = fc,
                                 edgecolor = edge_color(ward_id),
                                 #alpha = 1,
                                 lw = 0.5,
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
        ward_population = [0.,2.,20.]

        carto = WardCartogram(wards=wards,
                              ward_density=ward_population,
                              norm_density=True, # needs to be normed by area
                              margin_ratio=0.5,
                              x_raster_size=128,
                              y_raster_size=64,)

        #density_matrix = carto.cast_density_to_matrix(True)
        #carto.compute_cartogram((True))
        carto.compute(verbose=True)
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


    
