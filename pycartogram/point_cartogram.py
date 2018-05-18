
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
from pycartogram import WardCartogram
import visvalingamwyatt as vw
import scipy.sparse as sprs
from scipy.spatial import ConvexHull

class PointCartogram(WardCartogram):

    def __init__(self,
                 points, # this is a list of (x,y)-coordinates
                 wards = None,
                 norm_density=False,
                 margin_ratio = 0.2,
                 map_orientation = 'landscape',
                 x_raster_size = 1024,
                 y_raster_size = 768,
                 ):
        """
        wards:        list of wards as shapely.Polygon
        ward_density: list of density values for the wards
        norm_density: if this value is 'True' the ward_density will be normed by area (default: False)
        """
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

    def _get_whole_shape_matrix(self,verbose=False):
        A = np.zeros((self.xsize,self.ysize))
        self._mark_matrix_with_shape(A,self.whole_shape)
        return A

    def cast_density_to_matrix(self,verbose=False):
        self.density_matrix = self.cast_points_to_matrix(self.points,verbose)
        return self.density_matrix

    def compute(self,verbose=False):
        """do everything after init"""
        self.cast_density_to_matrix(verbose)
        self.compute_cartogram(verbose=verbose)
        self.transform_wards(verbose)
        x, y = self.transform_coords(self.points[:,0],self.points[:,1])
        n_ = len(x)
        self.new_points = np.concatenate((x.reshape(n_,1),y.reshape(n_,1)),axis=1)

    def plot_points(self,
             ax = None,                    
             show_density_matrix = False,
             show_new_points = True,
             plot_wards = True,
             ward_colors = 'None',
             bg_color = 'w',
             edge_colors = None,
             outline_whole_shape = True,
             use_new_density = False,
             **kwargs # matplotlib arguments for points
            ):
        """
            ward_colors can be
                - an iterable container of colors - plot ward face colors according to a provided list
                - "density" - plot ward face colors proportional to density
                - "log_density" - plot ward face colors proportional to density
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
        kwargs['linestyle'] = 'None'

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


    
