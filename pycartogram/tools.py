import numpy as np
import matplotlib.pyplot as pl
from shapely.geometry import Polygon, LineString, MultiPolygon, Point
from shapely.ops import polygonize
import itertools
import visvalingamwyatt as vw

def pair_iterate(l):
    for i in range(1,len(l)):
        yield l[i-1], l[i]

def enrich_polygon_with_points(geom,delta):
    coords = []
    for seg_start, seg_end in pair_iterate(geom.exterior.coords):
        #line_start = Point(seg_start)
        #line_end = Point(seg_end)
        segment = LineString([seg_start,seg_end])
        n_vals = segment.length / delta
        if n_vals > 1:
            interpol_vals = np.linspace(0,1,n_vals+1)
            for val in interpol_vals[:-1]:
                P = segment.interpolate(val,normalized=True)
                coords.append(P.coords[0])
        else:
            coords.append(seg_start)
    new_geom = Polygon(coords)

    return new_geom

def fill_matrix(A,i,j,new_val=1.,old_val=None):
    index_list = [(i,j)]
    if old_val is None:
        old_val = A[i,j]
    x_, y_ = A.shape

    while len(index_list) > 0:
        i,j = index_list.pop()
        A[i,j] = new_val
        if i+1 < x_ and A[i+1,j] == old_val:
            index_list.append((i+1,j))
        if j+1 < y_ and A[i,j+1] == old_val:
            index_list.append((i,j+1))
        if j-1 >= 0 and A[i,j-1] == old_val:
            index_list.append((i,j-1))
        if i-1 >= 0 and A[i-1,j] == old_val:
            index_list.append((i-1,j))            
        
def savefig_marginless(fn,fig,ax,**kwargs):
    ax.set_axis_off()
    fig.subplots_adjust(top = 1, bottom = 0, right = 1, left = 0, 
            hspace = 0, wspace = 0)
    ax.margins(0,0)
    ax.xaxis.set_major_locator(pl.NullLocator())
    ax.yaxis.set_major_locator(pl.NullLocator())
    fig.savefig(fn, bbox_inches = 'tight',
    pad_inches = 0,**kwargs)

def scale(hist,
                     new_min=0.2,
                     new_max=0.9,
                     take_inverse=False,
                     replace_nan = True,
                     get_nan_min_max = False
                    ):
    if take_inverse:
        factor = -1
    else:
        factor = 1
    log_hist = factor * np.array(hist)
    nan_min = np.nanmin(log_hist)
    nan_max = np.nanmax(log_hist)
    if replace_nan:
        log_hist[np.isnan(log_hist)] = nan_min - 1
    min_ = nan_min - 1
    max_ = nan_max
    intensity = lambda x: (x-min_)/ (max_-min_) * (new_max-new_min) + new_min
    if not get_nan_min_max:
        return log_hist, intensity
    else:
        return log_hist, intensity, nan_min, nan_max

def logify_and_scale(hist,
                     new_min=.2,
                     new_max=.9,
                     take_inverse=False,
                     replace_nan = True,
                     get_nan_min_max = False
                    ):
    if take_inverse:
        factor = -1
    else:
        factor = 1
    hist = np.array(hist)
    log_hist = np.array(hist)
    log_hist[hist>0.] = factor * np.log(hist[hist>0])
    log_hist[log_hist==0] = np.nan
    nan_min = np.nanmin(log_hist)
    nan_max = np.nanmax(log_hist)
    if replace_nan:
        log_hist[np.isnan(log_hist)] = nan_min - 1
    min_ = nan_min - 1
    max_ = nan_max
    intensity = lambda x: (x-min_)/ (max_-min_) * (new_max-new_min) + new_min
    if not get_nan_min_max:
        return log_hist, intensity
    else:
        return log_hist, intensity, nan_min, nan_max

def coarse_grain_wards(self,wards,th):
    new_wards = []
    for ward in wards:
        coo = ward.exterior.coords.xy
        new_coo = vw.Simplifier(zip(*coo)).simplify(threshold=th)
        new_ward = Polygon(new_coo)
        new_wards.append(new_ward)
    return new_wards

def is_iter(obj):
    try:
        _ = iter(obj)
        return True
    except TypeError as e:
        return False


def add_intersection_points_to_wards(wards):

    for g0, g1 in itertools.combinations(wards,2):
        
        if g0.touches(g1):

            g0_ring = LineString(list(g0.exterior.coords))
            g1_ring = LineString(list(g1.exterior.coords))

            # union of linestring
            union = g0_ring.union(g1_ring)

            # now if you polygonize, the resulting polygons
            # carry the coordinates of the intersections
            result = [ g for g in polygonize(union) ]
            g0 = result[0]
            g1 = result[1]

if __name__ == "__main__":

    A = Polygon([ 
                  (0.,0.),
                  (1.,0.),
                  (1.,1.),
                  (0.,1.), 
                  (0.,0.),
                ])

    B = enrich_polygon_with_points(A,0.1)

    fig, ax = pl.subplots(1,2)

    x,y = A.exterior.coords.xy
    ax[0].plot(x,y,'o')
    x,y = B.exterior.coords.xy
    ax[1].plot(x,y,'s')



    pl.show()

