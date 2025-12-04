"""
Google location history integration for cartogram generation.

This module provides the GoogleShapeProject class for loading Google
location history data and projecting it onto geographic regions (wards)
for use with pycartogram's cartogram generation.
"""

from __future__ import annotations

from typing import Any
import cartopy.crs as ccrs
import numpy as np
from cartopy.io.shapereader import Reader
import json
import shapely.geometry as sgeom
from shapely.geometry import Polygon, Point
from matplotlib.axes import Axes
from pycartogram.tools import (
    polygon_patch,
    enrich_polygon_with_points,
    enrich_polygon_to_n_points,
    get_json,
)
from shapely.ops import unary_union
from datetime import date
from datetime import datetime


class GoogleShapeProject():
    """
    Load and project Google location history onto geographic regions.

    Combines shapefile-defined regions (wards) with Google location history
    or arbitrary lon/lat point data, projecting everything to a common
    coordinate system for cartogram generation.

    Parameters
    ----------
    shape_file : str
        Path to shapefile containing ward/region boundaries.
    google_file : str or None
        Path to Google location history JSON file. If None, lon_lat_list
        must be provided.
    lon_lat_list : list of tuple, optional
        Alternative to google_file: list of (longitude, latitude) pairs.
    shape_source_proj : cartopy.crs.CRS, optional
        Source projection of shapefile (default: PlateCarree/WGS84).
    google_source_proj : cartopy.crs.CRS, optional
        Source projection of Google data (default: PlateCarree/WGS84).
    target_proj : cartopy.crs.CRS, optional
        Target projection for all data (default: UTM zone 33N).
    minimum_time_as_unix_seconds : float, optional
        Filter Google data to only include points after this timestamp
        (default: 0, includes all data).

    Attributes
    ----------
    wards : list of shapely.geometry.Polygon
        Ward boundaries in target projection.
    records : list
        Shapefile record metadata for each ward.
    whole_shape : shapely.geometry.Polygon
        Union of all ward boundaries.
    orig_bbox : shapely.geometry.Polygon
        Bounding box of all wards.
    relevant_points : list
        Points from Google data within the bounding box.
    relevant_points_by_hour : list of list
        Points grouped by hour of day (0-23).
    all_relevant_points : list
        Flattened list of all relevant points.
    """

    def __init__(
        self,
        shape_file: str,
        google_file: str | None,
        lon_lat_list: list[tuple[float, float]] | None = None,
        shape_source_proj: ccrs.CRS = ccrs.PlateCarree(),
        google_source_proj: ccrs.CRS = ccrs.PlateCarree(),
        target_proj: ccrs.CRS = ccrs.UTM(zone=33, southern_hemisphere=False),
        minimum_time_as_unix_seconds: float = 0,
    ) -> None:

        # load the single (the first) polygon from
        # each multipolygon represnting an electoral district
        # so the end is a list of polygons
        #self.wards = [\
        #    geom[0] \
        #    for geom in Reader(shape_file).geometries()\
        #  ]
        shape_reader = Reader(shape_file)
        shape_data = []
        shape_records = []
        for rec in shape_reader.records():
            shape_data.append(rec.geometry)
            shape_records.append(rec)

        new_shape_data = []
        for geom in shape_data:
            poly_coords = []
            for x,y in geom.exterior.coords:
                poly_coords.append(
                          target_proj.transform_point(
                                   x,
                                   y,
                                   shape_source_proj
                                )
                )
            new_shape_data.append(Polygon(poly_coords))
        shape_data = new_shape_data

        self.wards = shape_data
        self.records = shape_records

        # cascade each polygon s.t. we get a big polygon representing berlin
        self.whole_shape = unary_union(self.wards)
        # get Berlin bounding box as polygon
        x_ = (self.whole_shape.bounds[0],self.whole_shape.bounds[2])
        y_ = (self.whole_shape.bounds[1],self.whole_shape.bounds[3])
        self.orig_bbox = Polygon([
                                   (x_[0],y_[0]),
                                   (x_[1],y_[0]),
                                   (x_[1],y_[1]),
                                   (x_[0],y_[1]),
                                 ])

        # load the google data
        if google_file is not None:
            with open(google_file) as f:
                google_data = json.load(f)
                google_data = google_data['locations']


            # extract the relevant google data
            ggl_lon_lat_timestamp = [ ( dat['longitudeE7'] / 1e7,
                                        dat['latitudeE7']  / 1e7,
                                        float(dat['timestampMs']) / 1e3,
                                      ) for dat in google_data if float(dat['timestampMs']) / 1e3 > minimum_time_as_unix_seconds ]

            # convert google data to our polygon metric
            ggl_utm33n = [ (
                            sgeom.Point(
                                   *target_proj.transform_point(
                                           d[0],
                                           d[1],
                                           google_source_proj
                                        )
                                       ),
                            d[2]
                           ) for d in ggl_lon_lat_timestamp ]
            # filter all points for relevance (throw away those wich are not
            # in the berlin bounding box)
            self.relevant_points = [ d for d in ggl_utm33n if self.orig_bbox.contains(d[0]) ]
            self.relevant_points_by_hour = [ [] for i in range(24) ]

            # sort by hour of timestamp but also cumulate
            self.all_relevant_points = []
            for point, timestamp in self.relevant_points:
                hour = datetime.fromtimestamp(timestamp).hour
                self.relevant_points_by_hour[hour].append(point)
                self.all_relevant_points.append(point)
        elif lon_lat_list is not None:
            ggl_utm33n = [
                            sgeom.Point(
                                   *target_proj.transform_point(
                                           d[0],
                                           d[1],
                                           google_source_proj
                                        )
                                       )\
                            for d in lon_lat_list ]
            # filter all points for relevance (throw away those wich are not
            # in the berlin bounding box)
            self.relevant_points = [ d for d in ggl_utm33n if self.orig_bbox.contains(d) ]
            self.all_relevant_points = self.relevant_points


    def label_this(
        self,
        ax: Axes,
        wards: list[Polygon],
        labels: dict[str, str],
        key: str = 'PLZ99',
    ) -> None:
        """
        Add text labels to wards on a matplotlib axes.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            Axes to add labels to.
        wards : list of shapely.geometry.Polygon
            Ward polygons to label at centroids.
        labels : dict
            Mapping from ward key (e.g., postal code) to label text.
        key : str, optional
            Attribute name to use for ward identification (default: 'PLZ99').
        """
        for iward, ward in enumerate(wards):
            PLZ = self.records[iward].attributes['PLZ99']
            ha = 'center'
            va = 'center'
            if PLZ in labels:
                ax.text(ward.centroid.x,
                        ward.centroid.y,
                        labels[PLZ],
                        ha=ha,
                        va=va,
                        fontsize='small'
                       )

    def export_json(
        self,
        carto: Any,
        additional_data: dict[str, Any] | None = None,
        export_rec_attributes: list[str] | None = None,
    ) -> dict[str, Any]:
        """
        Export cartogram data to JSON format.

        Parameters
        ----------
        carto : WardCartogram
            Computed cartogram object with old_ward_coords and new_ward_coords.
        additional_data : dict, optional
            Extra data to include in export.
        export_rec_attributes : list of str, optional
            Shapefile attribute names to include for each ward.

        Returns
        -------
        dict
            JSON-serializable dictionary with original and transformed
            ward coordinates, domain bounds, and ward attributes.
        """
        if export_rec_attributes is None:
            export_rec_attributes = []
        export = {
                   'attributes': [],
                   'domain':{ 'x': (self.orig_bbox.bounds[0],self.orig_bbox.bounds[2]),
                              'y': (self.orig_bbox.bounds[1],self.orig_bbox.bounds[3
]),
                            },
                 }

        export['original'] = carto.old_ward_coords
        export['new'] = carto.new_ward_coords

        for iward in range(len(self.wards)):
            export['attributes'].append({})
            for attr in export_rec_attributes:
                export['attributes'][-1][attr] = self.records[iward].attributes[attr]

        if additional_data is not None:
            export.update(additional_data)

        return export


    def _export_json(
        self,
        new_wards: list[Polygon],
        delta: float,
        old_wards: list[Polygon] | None = None,
        additional_data: dict[str, Any] | None = None,
        export_rec_attributes: list[str] | None = None,
    ) -> dict[str, Any]:
        """
        Export ward data to JSON with polygon enrichment.

        Internal method that enriches polygons with additional vertices
        before export to enable smooth animations.

        Parameters
        ----------
        new_wards : list of shapely.geometry.Polygon
            Transformed ward polygons.
        delta : float
            Maximum distance between vertices when enriching polygons.
        old_wards : list of shapely.geometry.Polygon, optional
            Original ward polygons. If None, uses self.wards.
        additional_data : dict, optional
            Extra data to include in export.
        export_rec_attributes : list of str, optional
            Shapefile attribute names to include for each ward.

        Returns
        -------
        dict
            JSON-serializable dictionary with enriched original and
            transformed ward coordinates.
        """
        if export_rec_attributes is None:
            export_rec_attributes = []
        export = {
                   'attributes': [],
                   'domain':{ 'x': (self.orig_bbox.bounds[0],self.orig_bbox.bounds[2]),
                              'y': (self.orig_bbox.bounds[1],self.orig_bbox.bounds[3
]),
                            },
                 }

        no_old_wards_given = old_wards is None

        if no_old_wards_given:
            old_wards = self.wards

        enriched_wards = []
        for iward, ward in enumerate(old_wards):
            new_ward = enrich_polygon_with_points(ward,delta)
            n_points = len(new_ward.exterior.coords.xy[0])
            enriched_wards.append(enrich_polygon_to_n_points(new_ward,n_points))

        export['original'] = get_json(enriched_wards)
        export['new'] = get_json(new_wards)

        for iward in range(len(self.wards)):
            export['attributes'].append({})
            for attr in export_rec_attributes:
                export['attributes'][-1][attr] = self.records[iward].attributes[attr]

        if additional_data is not None:
            export.update(additional_data)

        return export






