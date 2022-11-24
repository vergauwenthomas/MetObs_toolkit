#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:28:36 2022

@author: thoverga
"""
import pandas as pd
# from shapely.geometry import Point
import rasterstats


    
# =============================================================================
# geodataset functions
# =============================================================================




def geotiff_point_extraction(geodf, geotiff_location, geotiff_crs, class_to_human_mapper,
                             band=1, interpolate='nearest'):

    coords_geodf = geodf[~geodf['geometry'].isnull()]
    missing_coords_geodf = geodf[geodf['geometry'].isnull()]
    


    #extract raster value for the point object
    coords_geodf['_numeric_label'] = rasterstats.point_query(vectors=coords_geodf['geometry'],
                                           raster=geotiff_location,
                                           band=band,
                                           nodata=None,
                                           affine=None,
                                           interpolate=interpolate, #nearest or bilinear (beliniear --> no direct mapping to a class possible)
                                           geojson_out=False)
    coords_geodf['_human_label'] = coords_geodf['_numeric_label'].map(class_to_human_mapper)
    missing_coords_geodf['_human_label'] = 'Unknown'
    
    comb_geodf = pd.concat([coords_geodf, missing_coords_geodf])
    comb_geodf = comb_geodf.sort_index()
    
    
    return comb_geodf['_human_label']
  