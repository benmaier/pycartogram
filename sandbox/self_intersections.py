import numpy as np
from shapely.geometry import Polygon, MultiPolygon
import matplotlib.pyplot as plt
from shapely.geometry import LineString, MultiPolygon
from shapely.ops import polygonize, unary_union, cascaded_union

x = np.array([ 0.38517325,  0.40859912,  0.43296919,  0.4583215 ,  0.4583215 ,
               0.43296919,  0.40859912,  0.38517325,  0.36265506,  0.34100929])
y = np.array([ 62.5       ,  56.17977528,  39.39698492,   0.        ,
               0.        ,  17.34605377,  39.13341671,  60.4180932 ,
               76.02574417,  85.47008547])
ls = LineString(np.c_[x, y])
lr = LineString(ls.coords[:] + ls.coords[:1])
mls = unary_union(lr)
mp = MultiPolygon(list(polygonize(mls)))
print 

for polygon in polygonize(mls):
    print(polygon)

