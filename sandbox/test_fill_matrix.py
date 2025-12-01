from __future__ import print_function
import numpy as np
import shapely.geometry as sgeom
from shapely.ops import unary_union, polygonize
from shapely.geometry import Polygon, LineString, MultiPolygon, Point
from pycartogram.tools import polygon_patch
import matplotlib as mpl
import matplotlib.pyplot as pl
import progressbar
import cCartogram as cart
from pycartogram import WardCartogram
import visvalingamwyatt as vw

x,y = ([389612.8319566059, 389853.54164727655, 390584.8584439841, 391486.60920954373, 391518.6112988488, 391311.3541827696, 391468.9528185115, 390058.6043322335, 389632.0783692974, 389161.581598248, 389612.8319566059],
 [5822146.7455827175, 5822149.384505943, 5821580.726679202, 5822252.992637813, 5821901.652236092, 5820928.236719232, 5820445.920328324, 5820517.711397066, 5821025.231596797, 5821930.263471162, 5822146.7455827175])

#A = Polygon([
#              (0.,0.),
#              (1.,1.5),
#              (1.,0.),
#            ])
A = Polygon(zip(x,y))

carto = WardCartogram([A],None,x_raster_size=5,y_raster_size=10)

mat = 0.5*np.ones((carto.xsize,carto.ysize))

carto._mark_matrix_with_shape(mat,A,new_val=0.1)

pl.imshow(mat)
pl.tight_layout()
pl.show()
