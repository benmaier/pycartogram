"""
Utility functions for cartogram generation and polygon manipulation.

This module provides helper functions for:
- Converting shapely geometries to matplotlib patches
- Enriching polygons with additional vertices for smoother transformations
- Scaling and normalizing density values
- Simplifying polygon boundaries
- Computing polygon adjacency networks
"""

from __future__ import annotations

from typing import Any, Callable, Generator, Sequence, Union
import numpy as np
from numpy.typing import ArrayLike, NDArray
import matplotlib.pyplot as pl
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from matplotlib.path import Path
from matplotlib.patches import PathPatch
from shapely.geometry import Polygon, LineString, MultiPolygon, Point
from shapely.validation import make_valid
from shapely.ops import polygonize
from scipy.spatial import Voronoi
import itertools
import visvalingamwyatt as vw


def polygon_patch(
    polygon: Union[Polygon, MultiPolygon],
    **kwargs: Any
) -> PathPatch:
    """
    Create a matplotlib PathPatch from a shapely Polygon or MultiPolygon.

    Drop-in replacement for the deprecated descartes.PolygonPatch.

    Parameters
    ----------
    polygon : shapely.geometry.Polygon or shapely.geometry.MultiPolygon
        The polygon geometry to convert to a matplotlib patch.
    **kwargs : dict
        Additional keyword arguments passed to matplotlib.patches.PathPatch
        (e.g., facecolor, edgecolor, alpha, linewidth).

    Returns
    -------
    matplotlib.patches.PathPatch
        A patch that can be added to a matplotlib axes via ax.add_patch().

    Raises
    ------
    TypeError
        If polygon is not a Polygon or MultiPolygon.

    Examples
    --------
    >>> from shapely.geometry import Polygon
    >>> poly = Polygon([(0, 0), (1, 0), (1, 1), (0, 1)])
    >>> patch = polygon_patch(poly, facecolor='blue', edgecolor='black')
    >>> ax.add_patch(patch)
    """
    def ring_to_codes(n):
        """Generate path codes for a ring with n points."""
        codes = [Path.LINETO] * n
        codes[0] = Path.MOVETO
        codes[-1] = Path.CLOSEPOLY
        return codes

    def polygon_to_path(poly):
        """Convert a single Polygon to vertices and codes."""
        vertices = []
        codes = []

        # Exterior ring
        ext_coords = list(poly.exterior.coords)
        vertices.extend(ext_coords)
        codes.extend(ring_to_codes(len(ext_coords)))

        # Interior rings (holes)
        for interior in poly.interiors:
            int_coords = list(interior.coords)
            vertices.extend(int_coords)
            codes.extend(ring_to_codes(len(int_coords)))

        return vertices, codes

    vertices = []
    codes = []

    if isinstance(polygon, MultiPolygon):
        for poly in polygon.geoms:
            v, c = polygon_to_path(poly)
            vertices.extend(v)
            codes.extend(c)
    elif isinstance(polygon, Polygon):
        vertices, codes = polygon_to_path(polygon)
    else:
        raise TypeError(f"Expected Polygon or MultiPolygon, got {type(polygon)}")

    path = Path(vertices, codes)
    return PathPatch(path, **kwargs)

def pair_iterate(l: Sequence[Any]) -> Generator[tuple[Any, Any], None, None]:
    """
    Iterate over consecutive pairs of elements in a sequence.

    Parameters
    ----------
    l : sequence
        Input sequence to iterate over.

    Yields
    ------
    tuple
        Consecutive pairs (l[i-1], l[i]) for i in range(1, len(l)).

    Examples
    --------
    >>> list(pair_iterate([1, 2, 3, 4]))
    [(1, 2), (2, 3), (3, 4)]
    """
    for i in range(1, len(l)):
        yield l[i-1], l[i]


def enrich_polygon_with_points(geom: Polygon, delta: float) -> Polygon:
    """
    Add interpolated points along polygon edges at regular intervals.

    This increases the vertex count of a polygon by adding points along
    each edge. Useful for ensuring smooth cartogram transformations.

    Parameters
    ----------
    geom : shapely.geometry.Polygon
        Input polygon to enrich with additional vertices.
    delta : float
        Maximum distance between consecutive points along edges.

    Returns
    -------
    shapely.geometry.Polygon
        New polygon with additional interpolated vertices.
    """
    coords = []
    for seg_start, seg_end in pair_iterate(geom.exterior.coords):
        segment = LineString([seg_start, seg_end])
        n_vals = int(segment.length / delta)
        if n_vals > 1:
            # Interpolate points along the segment
            interpol_vals = np.linspace(0, 1, n_vals + 1)
            for val in interpol_vals[:-1]:
                P = segment.interpolate(val, normalized=True)
                coords.append(P.coords[0])
        else:
            coords.append(seg_start)
    return Polygon(coords)


def enrich_polygon_to_n_points(geom: Polygon, n_total: int) -> Polygon:
    """
    Enrich a polygon to have exactly n_total vertices.

    Distributes additional points evenly across all segments to reach
    the target vertex count. Used for matching vertex counts between
    original and transformed polygons for animation.

    Parameters
    ----------
    geom : shapely.geometry.Polygon
        Input polygon to enrich.
    n_total : int
        Target total number of vertices.

    Returns
    -------
    shapely.geometry.Polygon
        New polygon with exactly n_total vertices, or original if
        it already has >= n_total vertices.
    """
    n_points = len(geom.exterior.coords.xy[0])
    n_segments = n_points - 1
    new_points = n_total - n_points

    if new_points <= 0:
        return geom

    # Distribute new points evenly across segments (round-robin)
    points_per_segment = [0 for i in range(n_segments)]
    current_segment = 0
    while new_points > 0:
        points_per_segment[current_segment % n_segments] += 1
        new_points -= 1
        current_segment += 1

    # Interpolate points along each segment
    coords = []
    current_segment = 0
    for seg_start, seg_end in pair_iterate(geom.exterior.coords):
        segment = LineString([seg_start, seg_end])
        n_vals = points_per_segment[current_segment] + 1
        if n_vals > 1:
            interpol_vals = np.linspace(0, 1, n_vals + 1)
            for val in interpol_vals[:-1]:
                P = segment.interpolate(val, normalized=True)
                coords.append(P.coords[0])
        else:
            coords.append(seg_start)
        current_segment += 1

    return Polygon(coords)


def get_cumulative_relative_exterior_length(geom: Polygon) -> NDArray[np.floating]:
    """
    Compute cumulative relative length along polygon exterior.

    Returns an array where each element represents the fraction of
    total perimeter length up to that vertex.

    Parameters
    ----------
    geom : shapely.geometry.Polygon
        Input polygon.

    Returns
    -------
    numpy.ndarray
        Array of cumulative relative lengths, normalized to [0, 1].
    """
    x, y = geom.exterior.xy
    length = np.zeros_like(np.array(x))
    for i in range(1, len(x)):
        r1 = np.array([x[i], y[i]])
        r0 = np.array([x[i-1], y[i-1]])
        length[i] = np.linalg.norm(r1 - r0)
    length = np.cumsum(length)
    return length / length[-1]


def match_vertex_count(geom0: Polygon, geom1: Polygon) -> tuple[Polygon, Polygon]:
    """
    Match vertex counts between two polygons by enriching the shorter one.

    Adds vertices to the polygon with fewer points so both have the same
    count. Vertices are added at positions that correspond to the same
    relative perimeter positions as in the longer polygon.

    Parameters
    ----------
    geom0 : shapely.geometry.Polygon
        First polygon.
    geom1 : shapely.geometry.Polygon
        Second polygon.

    Returns
    -------
    tuple of shapely.geometry.Polygon
        (geom0, geom1) with matched vertex counts. The polygon that
        originally had fewer vertices is enriched.
    """
    n0 = len(geom0.exterior.coords.xy[0])
    n1 = len(geom1.exterior.coords.xy[0])

    if n0 == n1:
        return geom0, geom1
    elif n0 > n1:
        short_geom = geom1
        long_geom = geom0
    else:
        short_geom = geom0
        long_geom = geom1

    # Match vertices by relative perimeter position
    short_length = get_cumulative_relative_exterior_length(short_geom)
    long_length = get_cumulative_relative_exterior_length(long_geom)
    matchings = [np.argmin(np.abs(_short - long_length)) for _short in short_length]

    # Interpolate to add vertices at matching positions
    coords = []
    for i, (seg_start, seg_end) in enumerate(pair_iterate(short_geom.exterior.coords)):
        dpoints = matchings[i+1] - matchings[i]
        segment = LineString([seg_start, seg_end])
        n_vals = dpoints + 1
        interpol_vals = np.linspace(0, 1, n_vals)
        if i < len(matchings) - 2:
            interpol_vals = interpol_vals[:-1]

        for val in interpol_vals:
            P = segment.interpolate(val, normalized=True)
            coords.append(P.coords[0])

    new_geom = Polygon(coords)

    if n0 > n1:
        return long_geom, new_geom
    else:
        return new_geom, long_geom


def fix_invalid_geometry(geom: Polygon) -> Polygon:
    """
    Fix invalid geometry while preserving Polygon type.

    Uses shapely's make_valid to fix self-intersections and other
    topology issues. If make_valid produces a MultiPolygon (from
    splitting a self-intersecting polygon), returns only the largest
    polygon by area.

    Parameters
    ----------
    geom : shapely.geometry.Polygon
        Input polygon that may be invalid.

    Returns
    -------
    shapely.geometry.Polygon
        Valid polygon geometry.
    """
    if geom.is_valid:
        return geom
    fixed = make_valid(geom)
    if isinstance(fixed, MultiPolygon):
        return max(fixed.geoms, key=lambda g: g.area)
    return fixed


def fix_geodataframe_geometries(gdf: "gpd.GeoDataFrame") -> "gpd.GeoDataFrame":
    """
    Return a copy of a GeoDataFrame with all geometries made valid.

    Useful before operations like dissolve() that require valid geometries.
    Note: This fixes geometries independently, which may break topological
    consistency between neighboring polygons. Use only when necessary
    (e.g., before dissolve), not for rendering where seamless boundaries matter.

    Parameters
    ----------
    gdf : geopandas.GeoDataFrame
        Input GeoDataFrame with potentially invalid geometries.

    Returns
    -------
    geopandas.GeoDataFrame
        Copy of input with all geometries validated.

    Examples
    --------
    >>> # Before dissolve, fix geometries to avoid TopologyException
    >>> gdf_valid = fix_geodataframe_geometries(gdf)
    >>> states = gdf_valid.dissolve(by="STATE_ID")
    """
    result = gdf.copy()
    result["geometry"] = result["geometry"].apply(fix_invalid_geometry)
    return result


def savefig_marginless(fn: str, fig: Figure, ax: Axes, **kwargs: Any) -> None:
    """
    Save a figure with no margins or whitespace.

    Removes all axes, labels, and padding for clean map exports.

    Parameters
    ----------
    fn : str
        Output filename.
    fig : matplotlib.figure.Figure
        Figure to save.
    ax : matplotlib.axes.Axes
        Axes to configure.
    **kwargs : dict
        Additional arguments passed to fig.savefig().
    """
    ax.set_axis_off()
    fig.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
    ax.margins(0, 0)
    ax.xaxis.set_major_locator(pl.NullLocator())
    ax.yaxis.set_major_locator(pl.NullLocator())
    fig.savefig(fn, bbox_inches='tight', pad_inches=0, **kwargs)


def scale(
    hist: ArrayLike,
    new_min: float = 0.2,
    new_max: float = 0.9,
    take_inverse: bool = False,
    replace_nan: bool = True,
    get_nan_min_max: bool = False,
) -> Union[
    tuple[NDArray[np.floating], Callable[[float], float]],
    tuple[NDArray[np.floating], Callable[[float], float], float, float]
]:
    """
    Scale values to a new range for visualization.

    Parameters
    ----------
    hist : array-like
        Input values to scale.
    new_min : float, optional
        Minimum of output range (default: 0.2).
    new_max : float, optional
        Maximum of output range (default: 0.9).
    take_inverse : bool, optional
        If True, invert values before scaling (default: False).
    replace_nan : bool, optional
        If True, replace NaN values with min - 1 (default: True).
    get_nan_min_max : bool, optional
        If True, also return original min/max before NaN replacement.

    Returns
    -------
    scaled_values : numpy.ndarray
        Scaled values.
    intensity : callable
        Function to map original values to scaled range.
    nan_min, nan_max : float, optional
        Original min/max (only if get_nan_min_max=True).
    """
    factor = -1 if take_inverse else 1
    scaled = factor * np.array(hist)
    nan_min = np.nanmin(scaled)
    nan_max = np.nanmax(scaled)
    if replace_nan:
        scaled[np.isnan(scaled)] = nan_min - 1
    min_ = nan_min - 1
    max_ = nan_max
    intensity = lambda x: (x - min_) / (max_ - min_) * (new_max - new_min) + new_min
    if not get_nan_min_max:
        return scaled, intensity
    else:
        return scaled, intensity, nan_min, nan_max


def logify_and_scale(
    hist: ArrayLike,
    new_min: float = 0.2,
    new_max: float = 0.9,
    take_inverse: bool = False,
    replace_nan: bool = True,
    get_nan_min_max: bool = False,
) -> Union[
    tuple[NDArray[np.floating], Callable[[float], float]],
    tuple[NDArray[np.floating], Callable[[float], float], float, float]
]:
    """
    Apply log transform and scale values to a new range.

    Similar to scale() but applies logarithm first. Useful for
    density data that spans many orders of magnitude.

    Parameters
    ----------
    hist : array-like
        Input values (should be positive for log transform).
    new_min : float, optional
        Minimum of output range (default: 0.2).
    new_max : float, optional
        Maximum of output range (default: 0.9).
    take_inverse : bool, optional
        If True, invert values before scaling (default: False).
    replace_nan : bool, optional
        If True, replace NaN values with min - 1 (default: True).
    get_nan_min_max : bool, optional
        If True, also return original min/max before NaN replacement.

    Returns
    -------
    log_scaled : numpy.ndarray
        Log-transformed and scaled values.
    intensity : callable
        Function to map log values to scaled range.
    nan_min, nan_max : float, optional
        Original min/max (only if get_nan_min_max=True).
    """
    factor = -1 if take_inverse else 1
    hist = np.array(hist)
    log_hist = np.array(hist)
    log_hist[hist > 0.] = factor * np.log(hist[hist > 0])
    log_hist[log_hist == 0] = np.nan
    nan_min = np.nanmin(log_hist)
    nan_max = np.nanmax(log_hist)
    if replace_nan:
        log_hist[np.isnan(log_hist)] = nan_min - 1
    min_ = nan_min - 1
    max_ = nan_max
    intensity = lambda x: (x - min_) / (max_ - min_) * (new_max - new_min) + new_min
    if not get_nan_min_max:
        return log_hist, intensity
    else:
        return log_hist, intensity, nan_min, nan_max


def coarse_grain_wards(wards: list[Polygon], th: float) -> list[Polygon]:
    """
    Simplify ward polygons using Visvalingam-Whyatt algorithm.

    Reduces vertex count while preserving shape characteristics.
    Useful for faster cartogram computation with complex boundaries.

    Parameters
    ----------
    wards : list of shapely.geometry.Polygon
        Input polygons to simplify.
    th : float
        Simplification threshold. Higher values = more simplification.

    Returns
    -------
    list of shapely.geometry.Polygon
        Simplified polygons.
    """
    new_wards = []
    for ward in wards:
        coo = ward.exterior.coords.xy
        new_coo = vw.Simplifier(list(zip(*coo))).simplify(threshold=th)
        new_ward = Polygon(new_coo)
        new_wards.append(new_ward)
    return new_wards


def is_iter(obj: Any) -> bool:
    """
    Check if an object is iterable (but not a string).

    Parameters
    ----------
    obj : any
        Object to check.

    Returns
    -------
    bool
        True if obj is iterable, False otherwise.
    """
    try:
        _ = iter(obj)
        return True
    except TypeError:
        return False


def get_json(wards: list[Polygon]) -> list[list[list[float]]]:
    """
    Convert ward polygons to JSON-serializable format.

    Parameters
    ----------
    wards : list of shapely.geometry.Polygon
        Input polygons.

    Returns
    -------
    list of list
        List of polygons, each as a list of [x, y] coordinate pairs.
    """
    polygons = []
    for ward in wards:
        this_poly = [[x, y] for x, y in ward.exterior.coords]
        polygons.append(this_poly)
    return polygons


def get_polygon_network(wards: list[Polygon]) -> list[tuple[int, int]]:
    """
    Build adjacency network of touching polygons.

    Parameters
    ----------
    wards : list of shapely.geometry.Polygon
        Input polygons.

    Returns
    -------
    list of tuple
        List of (i, j) tuples indicating which polygon pairs touch.
    """
    edges = []
    for ig0, ig1 in itertools.combinations(range(len(wards)), 2):
        if wards[ig0].touches(wards[ig1]):
            edges.append((ig0, ig1))
    return edges


def add_intersection_points_to_wards(wards: list[Polygon]) -> None:
    """
    Add shared boundary points between adjacent wards.

    Ensures that touching wards share exactly the same vertices
    along their common boundary. This prevents gaps in the
    transformed cartogram.

    Parameters
    ----------
    wards : list of shapely.geometry.Polygon
        Input polygons (modified in place).

    Note
    ----
    This function modifies the ward list in place by finding
    intersection points and adding them to both touching polygons.
    """
    for g0, g1 in itertools.combinations(wards, 2):
        if g0.touches(g1):
            g0_ring = LineString(list(g0.exterior.coords))
            g1_ring = LineString(list(g1.exterior.coords))

            # Union creates a linestring with all intersection points
            union = g0_ring.union(g1_ring)

            # Polygonize reconstructs polygons with shared vertices
            result = [g for g in polygonize(union)]
            if len(result) >= 2:
                g0 = result[0]
                g1 = result[1]


def voronoi_finite_polygons_2d(
    vor: Voronoi,
    radius: float | None = None
) -> tuple[list[list[int]], NDArray[np.floating]]:
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

