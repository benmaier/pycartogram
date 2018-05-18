import cartopy.crs as ccrs
import numpy as np
from cartopy.io.shapereader import Reader
import json
import shapely.geometry as sgeom
from shapely.geometry import Polygon
from descartes.patch import PolygonPatch
from shapely.ops import cascaded_union
from datetime import date
from datetime import datetime
from pycartogram.tools import *

class GoogleShapeProject():

    def __init__(self,
                 shape_file,
                 google_file,
                 shape_source_proj = ccrs.PlateCarree(),
                 google_source_proj = ccrs.PlateCarree(),
                 target_proj = ccrs.UTM('33N'),
                 ):

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
            shape_data.append(rec.geometry[0])
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
        self.whole_shape = cascaded_union(self.wards) 
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
        with open(google_file) as f:
            google_data = json.load(f)
            google_data = google_data['locations']


        # extract the relevant google data
        ggl_lon_lat_timestamp = [ ( dat['longitudeE7'] / 1e7,
                                    dat['latitudeE7']  / 1e7,
                                    float(dat['timestampMs']) / 1e3,
                                  ) for dat in google_data ]

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

    def label_this(ax,wards,labels,key='PLZ99'):
        for iward, ward in enumerate(wards):
            PLZ = self.records[iward].attributes['PLZ99']
            ha = 'center'
            va = 'center'
            if PLZ in labels:
            #if iward in indices:
                ax.text(ward.centroid.x,
                        ward.centroid.y,
                        labels[PLZ],
                        ha = ha,
                        va = va,
                        fontsize='small'
                       )       

    def export_json(self,carto,additional_data=None,export_rec_attributes=[]):
        
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
        

    def _export_json(self,new_wards,delta,old_wards=None,additional_data=None,export_rec_attributes=[]):
        
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
        


        

        
