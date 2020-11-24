GDAL_LIBRARY_PATH = "/usr/local/lib/libgdal.dylib"
import ctypes
ctypes.CDLL(GDAL_LIBRARY_PATH)
import numpy as np
import geopandas as gpd
import pandas as pd
from pycartogram.ward_cartogram import WardCartogram
from shapely.geometry import Polygon, MultiPolygon

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

if __name__ == "__main__":

    np.random.seed(30)
    gdf = gpd.read_file("/Users/bfmaier/Seafile/german_geo_data/construct_nuts3_plus_berliner_bezirke/targetdata/nuts3_and_berliner_bezirke_simplified.json")
    print(gdf.head())

    from rocsDB import rocsDB

    db = rocsDB()
    pop = dict(db.submit_query("""
        select germanid, population from censusdata.german_counties_info where tl_merge_id is not null
    """))
    db.close()
    _pop = []
    for _id, pol in zip(gdf.countyid.to_list(), gdf.geometry.to_list()):
        _pop.append(pop[int(_id)] / pol.area)
    gdf['population'] = _pop

    carto = GeoDataFrameWardCartogram(gdf,'population',y_raster_size=2048,x_raster_size=2048,map_orientation='portrait')
    carto.compute(verbose=True)

    import matplotlib as mpl
    import matplotlib.pyplot as pl

    
    new_gdf = carto.get_cartogram_geo_df()
    new_gdf.plot()

    #carto.plot(show_new_wards=False,show_density_matrix=True)
    pl.show()
