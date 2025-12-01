from __future__ import print_function
import numpy as np
import shapely.geometry as sgeom
from shapely.ops import unary_union
from shapely.geometry import Polygon
from pycartogram.tools import polygon_patch, voronoi_finite_polygons_2d
import matplotlib as mpl
import matplotlib.pyplot as pl
from matplotlib import collections  as mplcoll
import progressbar
import cCartogram as cart
from pycartogram import WardCartogram
from pycartogram import PointCartogram
import scipy.sparse as sprs
from itertools import combinations
from scipy.spatial import ConvexHull

from scipy.spatial import Voronoi
from scipy.spatial import Delaunay

class VoronoiCartogram(PointCartogram):

    def __init__(self,
                 points, # this is a list of (x,y)-coordinates
                 wards = None,
                 margin_ratio = 0,
                 map_orientation = 'landscape',
                 x_raster_size = 1024,
                 y_raster_size = 768,
                 ):
        """
        wards:        list of wards as shapely.Polygon
        ward_density: list of density values for the wards
        """

        PointCartogram.__init__(self,
                               points,
                               wards = wards,
                               margin_ratio = margin_ratio,
                               map_orientation = map_orientation,
                               x_raster_size = x_raster_size,
                               y_raster_size = y_raster_size,
                              )

        self.compute()


    def cast_density_to_matrix(self,verbose=False):
        self.density_matrix = self.cast_points_to_matrix(self.points,verbose,replace_value_zero=False)
        return self.density_matrix

    def compute(self,verbose=False):
        """do everything after init"""
        self.cast_density_to_matrix(verbose)
        self.compute_voronoi(verbose=verbose)

    def compute_voronoi(self,verbose=False):
        x_indices, y_indices = np.nonzero(self.density_matrix)
        x_ = np.array(x_indices,dtype=float)
        y_ = np.array(y_indices,dtype=float)
        x_coords = self.x_from_i(x_)
        y_coords = self.y_from_j(y_)
        points_for_voronoi = np.concatenate((
                                    self.x_from_i(x_).reshape(len(x_),1),
                                    self.y_from_j(y_).reshape(len(y_),1)
                                   ),axis=1)
        vor = Voronoi(points_for_voronoi)
        width = self.orig_width
        point_regions, vertices = voronoi_finite_polygons_2d(vor,radius=width)

        self.new_wards = []
        self.new_ward_density = []

        for ipoint, reg in enumerate(point_regions):
            # get closed polygon
            geom = Polygon(vertices[reg+[reg[0]]].tolist())
            # only take the part which is actually in Berlin
            geom = geom.intersection(self.orig_bbox)
            A = geom.area
            if A > 0:
                self.new_ward_density.append(1./A)
                self.new_wards.append(geom)
            else:
                print("got rid of a polygon")

        self.new_bbox = self.orig_bbox
        self.new_whole_shape = self.whole_shape
        self.ward_density = self.new_ward_density

        edges = set()
        tri = Delaunay(points_for_voronoi)

        for p in tri.vertices:

            for i,j in combinations(p,2):
                if self.new_wards[i].touches(self.new_wards[j]):
                    if i<j:
                        edges.add((i,j))
                    else:
                        edges.add((j,i))

        self.network_lines = [ (points_for_voronoi[i,:].tolist(),
                                points_for_voronoi[j,:].tolist(),) for i,j in edges ]



    def plot_voronoi(self,
                     draw_point_network=True,
                     network_color = [1,1,1],
                     network_alpha = 0.5,
                     network_linewidth = 1,
                     **kwargs):

        fig, ax = self.plot(**kwargs)

        """
        
        if 'bg_color' in kwargs:
            bg_color = kwargs['bg_color']            
        else:
            bg_color = 'w'

        if bg_color == 'w':
            col = [0,0,0,0.25]
        elif bg_color == 'k':
            col = [1,1,1,0.25]
        else:
            col = [0,0,0,0.25]
        """

        col = list(mpl.colors.to_rgba(network_color))
        col[-1] = network_alpha


        if draw_point_network:
            lc = mplcoll.LineCollection(self.network_lines,
                                        colors=col,
                                        linewidths=network_linewidth)
            ax.add_collection(lc)

        else:
            ax.plot(self.points_for_voronoi[:,0],
                    self.points_for_voronoi[:,1],
                    'o',
                    markersize=0.4,
                    mew = 0,
                    mfc = col,
                   )  

        return fig, ax



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

        carto = VoronoiCartogram(points,
                               x_raster_size=128,
                               y_raster_size=64,
                              )

        #density_matrix = carto.cast_density_to_matrix(True)
        #fig, ax = carto.plot_points(show_density_matrix=True,show_new_points=False,show_wards=True)
        #fig, ax = carto.plot_points(show_new_points=False)
        #x,y = carto.get_ward_bounds(carto.big_bbox)
        #pl.imshow(density_matrix.T,extent=x+y,origin='lower')
        #carto.compute_cartogram((True))
        #pts = carto.points
        #x, y = pts[:,0], pts[:,1]
        #pl.plot(x,y,'.')
        #carto.compute(verbose=True)
        fig, ax = carto.plot_voronoi(bg_color='k',intensity_range=[0,0.6])
        print(sorted(carto.ward_density))
        print(len(carto.ward_density))
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


    
