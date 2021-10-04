#GDAL_LIBRARY_PATH = "/usr/local/lib/libgdal.dylib"
#import ctypes
#ctypes.CDLL(GDAL_LIBRARY_PATH)
import numpy as np
import geopandas as gpd
import pandas as pd
from pycartogram.ward_cartogram import WardCartogram
from shapely.geometry import Polygon, MultiPolygon
from easing_functions import QuadEaseInOut, CubicEaseInOut

from pycartogram.tools import enrich_polygon_to_n_points, match_vertex_count

from copy import deepcopy

class GeoDataFrameWardCartogram(WardCartogram):

    def __init__(self,
                 geo_df,
                 ward_density_column,
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

        # loop through 

        self.gdf = geo_df.copy()
        polygons = geo_df.geometry.to_list()
        orig_ward_density = geo_df[ward_density_column].to_list()
        self.ward_indices = [ [] for i in range(len(geo_df)) ]

        i = 0
        wards = []
        ward_density = []
        for ipoly, (poly, dens) in enumerate(zip(polygons, orig_ward_density)):
            if isinstance(poly, Polygon):
                wards.append(poly)
                ward_density.append(dens)
                self.ward_indices[ipoly].append(i)
                i += 1
            elif isinstance(poly, MultiPolygon):
                for _poly in list(poly):
                    wards.append(_poly)
                    ward_density.append(dens)
                    self.ward_indices[ipoly].append(i)
                    i += 1

        self._equal_length_enriched_wards = None
        self._equal_length_new_wards = None

        WardCartogram.__init__(self,wards,
                                    ward_density,
                                    False,
                                    margin_ratio,
                                    map_orientation,
                                    x_raster_size,
                                    y_raster_size,
                                )

    def get_enriched_original_geo_df(self):
        """
        Get a copy of the original geo dataframe where each polygon
        has more vertices (such that transitions to the cartogram
        can be smoother)
        """
        return self._get_geo_df(self.new_old_wards)

    def get_cartogram_geo_df(self):
        """
        Get the computed cartogram as a geo data frame.
        """
        return self._get_geo_df(self.new_wards)

    def get_interpolated_geo_df(self,t,ease='QuadEaseInOut'):

        if self._equal_length_enriched_wards is None:
            self._equal_length_enriched_wards = list(self.new_old_wards)
            self._equal_length_new_wards = list(self.new_wards)
            def shape_length(shape):
                return len(shape.exterior.xy[0])

            for i, (old, new) in enumerate(zip(self._equal_length_enriched_wards, self._equal_length_new_wards)):
                _old, _new = match_vertex_count(old, new)
                self._equal_length_enriched_wards[i] = _old
                self._equal_length_new_wards[i] = _new

        if ease == 'QuadEaseInOut':
            t = QuadEaseInOut()(t)
        elif ease == 'CubicEaseInOut':
            t = CubicEaseInOut()(t)

        assert(0<=t<=1)

        if t == 0:
            return self.get_enriched_original_geo_df()
        elif t == 1:
            return self.get_cartogram_geo_df()

        wards = []
        _a = np.array
        for old, new in zip(self._equal_length_enriched_wards, self._equal_length_new_wards):
            xo, yo = old.exterior.xy
            xn, yn = new.exterior.xy
            x = (1-t) * _a(xo) + t * _a(xn)
            y = (1-t) * _a(yo) + t * _a(yn)
            shape = Polygon(list(zip(x,y)))
            wards.append(shape)

        return self._get_geo_df(wards)


    def _get_geo_df(self,wards):
        """
        Get the computed cartogram as a geo data frame.
        """
        gdf = self.gdf.copy()
        new_geometry = []
        for indices in self.ward_indices:
            if len(indices) == 1:
                new_geometry.append(wards[indices[0]])
            else:
                these_polygons = [ wards[ndx] for ndx in indices ]
                this_shape = MultiPolygon(these_polygons)
                new_geometry.append(this_shape)
        gdf['geometry'] = new_geometry
        return gdf

    def _get_custom_json(self,wards,label):
        """
        .. code:: python

            {
                'type': 'single_map',
                'map_label': 'cartogram',
                'xlim': [0,4],
                'ylim': [0,4],
                'polygons' : [
                    {
                        'name': 'foo',
                        'attr1': 'value1',
                        'attr2': 2.8,
                        'polygon': {
                            'x': [0,1,2,3,4,0],
                            'y': [2,3,4,0,0,1],
                        }
                        'affiliated_polygon_indices': [],
                    },
                    {
                        'name': 'bar',
                        'attr1': 'value3',
                        'attr2': 3.2,
                        'polygon': {
                            'x': [1.5,1.5,0.5,0.5,1.5],
                            'y': [1.5,0.5,0.5,1.5,1.5],
                        }
                        'affiliated_polygon_indices': [2],
                    },
                    {
                        'name': 'bar',
                        'attr1': 'value3',
                        'attr2': 3.2,
                        'polygon': {
                            'x': [2.5,2.5,1.5,1.5,2.5],
                            'y': [2.5,1.5,1.5,2.5,2.5],
                        }
                        'affiliated_polygon_indices': [1],
                    },
                ]
            }
        """

        gdf = self.gdf.drop('geometry')

        cols = gdf.columns.to_list()

        polygons = [ None for i in range(len(wards)) ]

        xlim, ylim = [1e300,-1e300], [1e300, -1e300]

        for irow, row in enumerate(gdf.to_records(index=False)):

            this_dict = dict(zip(cols, row))
            indices = self.ward_indices[irow]

            for i in indices:
                this_polygon = deepcopy(this_dict)
                _indices = deepcopy(indices)
                _indices.pop(_indices.index(i))
                x, y = zip(*list(wards[i].exterior.coords))

                xlim[0] = min(xlim[0], min(x))
                xlim[1] = max(xlim[1], max(x))
                ylim[0] = min(ylim[0], min(y))
                ylim[1] = max(ylim[1], max(y))

                this_polygon['polygon'] = {
                            'x': x,
                            'y': y,
                        }

                this_polygon['affiliated_polygon_indices'] = _indices
                polygons[i] = this_polygon

        this_data = {
                    'type': 'single_map',
                    'map_label': label,
                    'xlim': xlim,
                    'ylim': ylim,
                    'polygons': polygons
                }

        return this_data

    def get_enriched_original_as_custom_json(self,df,label='original'):
        """
        """
        return self._get_custom_json(self.new_old_wards)

    def get_cartogram_as_custom_json(self,df,label='cartogram'):
        """
        """
        return self._get_custom_json(self.new_wards,label)


def merge_custom_jsons(*list_of_custom_jsons):

    """
    .. code:: python

        {
            'type': 'multiple_maps',
            'map_labels': [ 'original', 'cartogram' ],
            'xlim': [0,4],
            'ylim': [0,4],
            'polygons' : [
                {
                    'name': 'foo',
                    'attr1': 'value1',
                    'attr2': 3.3,
                    'polygons': [
                        {
                            'x': [0,1,2,3,4,0],
                            'y': [2,3,4,0,0,1],
                        },
                        {
                            'x': [0,0.5,1,1.5,2,0],
                            'y': [0.5,1.5,2,0,0,0.5],
                        }
                    ]
                },
                {
                    'name': 'bar',
                    'attr1': 'value3',
                    'attr2': 4.47,
                    'polygons': [
                        {
                            'x': [1.5,1.5,0.5,0.5,1.5],
                            'y': [1.5,0.5,0.5,1.5,1.5],
                        },
                        {
                            'x': [0.75,0.75,0.25,0.25,0.75],
                            'y': [0.75,0.25,0.25,0.75,0.75],
                        }
                    ]
                }
            ]
        }
    """


    data = list_of_custom_jsons
    xmin = min(filter(lambda d: d['xlim'][0],data))
    xmax = max(filter(lambda d: d['xlim'][1],data))
    ymin = min(filter(lambda d: d['ylim'][0],data))
    ymax = max(filter(lambda d: d['ylim'][1],data))
    xlim = [xmin, xmax]
    ylim = [ymin, ymax]
    labels = list(filter(lambda d: d['label'],data))

    polygons = deepcopy(data[0]['polygons'])
    for i, entry in enumerate(polygons):
        polygons[i]['polygons'] = [ polygons[i].pop('polygon') ]

    for other in data[1:]:
        for i, entry in enumerate(other):
            polygons[i]['polygons'].append(deepcopy(entry['polygon']))

    this_data = {
                'type' : 'multiple_maps',
                'map_labels' : labels,
                'xlim' : xlim,
                'ylim' : ylim,
                'polygons' : polygons
            }

    return this_data




if __name__ == "__main__":

    np.random.seed(30)
    gdf = gpd.read_file("/Users/bfmaier/Seafile/german_geo_data/construct_nuts3_plus_berliner_bezirke/targetdata/nuts3_and_berliner_bezirke_simplified.json")
    print(gdf.head())

    from rocsDB import rocsDB
    import matplotlib as mpl
    import matplotlib.pyplot as pl


    db = rocsDB()
    pop = dict(db.submit_query("""
        select germanid, population from censusdata.german_counties_info where tl_merge_id is not null
    """))
    db.close()
    _pop = []
    for _id, pol in zip(gdf.countyid.to_list(), gdf.geometry.to_list()):
        _pop.append(pop[int(_id)] / pol.area)
    gdf['population'] = _pop

    carto = GeoDataFrameWardCartogram(gdf,'population',y_raster_size=2048,x_raster_size=768*2,map_orientation='portrait')


    carto.cast_density_to_matrix(verbose=True,set_boundary_to='400percent_min')
    pl.figure()
    carto.plot(show_new_wards=False,show_density_matrix=False)
    pl.show()


    carto.compute(verbose=True)

    
    new_gdf = carto.get_cartogram_geo_df()
    #new_gdf.plot()

    pl.figure()
    carto.plot()

    pl.show()
