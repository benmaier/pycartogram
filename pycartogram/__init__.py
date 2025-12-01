"""
pycartogram - Tools for generating cartograms using the Gastner-Newman algorithm.

This package provides classes for creating density-equalizing cartograms
from geographic data using diffusion-based methods.
"""

from .ward_cartogram import WardCartogram
from .point_cartogram import PointCartogram
from .voronoi_cartogram import VoronoiCartogram
from .google_loader import GoogleShapeProject
from .tools import (
    polygon_patch,
    enrich_polygon_with_points,
    enrich_polygon_to_n_points,
    match_vertex_count,
    coarse_grain_wards,
    get_json,
    get_polygon_network,
    voronoi_finite_polygons_2d,
)

__version__ = "0.1.0"

__all__ = [
    # Main cartogram classes
    "WardCartogram",
    "PointCartogram",
    "VoronoiCartogram",
    # Data loaders
    "GoogleShapeProject",
    # Utility functions
    "polygon_patch",
    "enrich_polygon_with_points",
    "enrich_polygon_to_n_points",
    "match_vertex_count",
    "coarse_grain_wards",
    "get_json",
    "get_polygon_network",
    "voronoi_finite_polygons_2d",
]
