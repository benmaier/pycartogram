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

def get_polygon_network(wards):

    edges = []

    for ig0, ig1 in itertools.combinations(range(len(wards)),2):

        g0 = wards[ig0]
        g1 = wards[ig1]
        
        if g0.touches(g1):
            edges.append((ig0,ig1))

    return edges


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

def voronoi_finite_polygons_2d(vor, radius=None):
    """
    copied from https://gist.github.com/pv/8036995

    Reconstruct infinite voronoi regions in a 2D diagram to finite
    regions.

    Parameters
    ----------
    vor : Voronoi
        Input diagram
    radius : float, optional
        Distance to 'points at infinity'.

    Returns
    -------
    regions : list of tuples
        Indices of vertices in each revised Voronoi regions.
    vertices : list of tuples
        Coordinates for revised Voronoi vertices. Same as coordinates
        of input vertices, with 'points at infinity' appended to the
        end.

    """

    if vor.points.shape[1] != 2:
        raise ValueError("Requires 2D input")

    new_regions = [ None for i in range(len(vor.point_region)) ]
    new_vertices = vor.vertices.tolist()

    center = vor.points.mean(axis=0)
    if radius is None:
        radius = vor.points.ptp().max()

    # Construct a map containing all ridges for a given point
    all_ridges = {}
    for (p1, p2), (v1, v2) in zip(vor.ridge_points, vor.ridge_vertices):
        all_ridges.setdefault(p1, []).append((p2, v1, v2))
        all_ridges.setdefault(p2, []).append((p1, v1, v2))

    # Reconstruct infinite regions
    for p1, region in enumerate(vor.point_region):

        vertices = vor.regions[region]

        if all(v >= 0 for v in vertices):
            # finite region
            new_regions[p1] = vertices
            continue

        # reconstruct a non-finite region
        ridges = all_ridges[p1]
        new_region = [v for v in vertices if v >= 0]

        for p2, v1, v2 in ridges:
            if v2 < 0:
                v1, v2 = v2, v1
            if v1 >= 0:
                # finite ridge: already in the region
                continue

            # Compute the missing endpoint of an infinite ridge

            t = vor.points[p2] - vor.points[p1] # tangent
            t /= np.linalg.norm(t)
            n = np.array([-t[1], t[0]])  # normal

            midpoint = vor.points[[p1, p2]].mean(axis=0)
            direction = np.sign(np.dot(midpoint - center, n)) * n
            far_point = vor.vertices[v2] + direction * radius

            new_region.append(len(new_vertices))
            new_vertices.append(far_point.tolist())

        # sort region counterclockwise
        vs = np.asarray([new_vertices[v] for v in new_region])
        c = vs.mean(axis=0)
        angles = np.arctan2(vs[:,1] - c[1], vs[:,0] - c[0])
        new_region = np.array(new_region)[np.argsort(angles)]

        # finish
        new_regions[p1] = new_region.tolist()
    return new_regions, np.asarray(new_vertices)

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

