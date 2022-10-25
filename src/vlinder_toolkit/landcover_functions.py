#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:28:36 2022

@author: thoverga
"""
import pandas as pd
import geopandas as gpd
# from shapely.geometry import Point
import rasterstats

# =============================================================================
# Geometry functions
# =============================================================================


def coordinate_to_point_geometry(lat, lon, crs):
    """ This function returns a shapely point object in the given coordinate referece frame. """

    point_geo = gpd.GeoDataFrame(pd.DataFrame(data={'lat': [lat], 'lon': [lon]}), geometry=gpd.points_from_xy([lon], [lat])) #to geopandas df
    point_geo = point_geo.set_crs(epsg = 4326) #inpunt are gps coordinates
    point_geo = point_geo.to_crs(crs) #coordinate transform
    point_geo = point_geo.iloc[0]['geometry']

    return point_geo


# =============================================================================
# geodataset functions
# =============================================================================




def geotiff_point_extraction(lat, lon, geotiff_location, geotiff_crs, class_to_human_mapper,
                             band=1, interpolate='nearest'):
    

    #transform that lat-lon coordinates to the coordinatesystem of the rasterfile, and make geometry point object.
    point = coordinate_to_point_geometry(lat=lat,
                                         lon=lon,
                                         crs=geotiff_crs)

    #extract raster value for the point object
    raster_value = rasterstats.point_query(vectors=point,
                                           raster=geotiff_location,
                                           band=band,
                                           nodata=None,
                                           affine=None,
                                           interpolate=interpolate, #nearest or bilinear (beliniear --> no direct mapping to a class possible)
                                           geojson_out=False)
    if isinstance(raster_value[0], type(None)):
        print("Location ", lat, '; ', lon, " not found on map --> return 'unknown'.")
        cover = 'Unknown'
    else:
        cover = class_to_human_mapper[raster_value[0]]
    return cover