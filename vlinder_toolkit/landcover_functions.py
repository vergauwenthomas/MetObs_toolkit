#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:28:36 2022

@author: thoverga
"""


import sys
import pandas as pd
import ee


# =============================================================================
#  Connection functions
# =============================================================================

def connect_to_gee():
    if not ee.data._credentials: #check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return



# =============================================================================
# Object convertors
# =============================================================================



def get_ee_obj(mapinfo):
    if mapinfo['is_image']:
        obj = ee.Image(mapinfo['location'])
    elif mapinfo['is_imagecollection']:
        if isinstance(mapinfo['band_of_use'], type(None)):
            obj = ee.ImageCollection(mapinfo['location'])
        else:
            obj = ee.ImageCollection(mapinfo['location']).select(mapinfo['band_of_use'])
        
    else:
        sys.exit('Map type is not an Image or Imagecollection.')
    return obj




def coords_to_geometry(lat=[], lon=[]):
    if len(lat) == 1:
        return ee.Geometry.Point(lon, lat)
    else:
        return ee.Geometry.MultiPoint(list(zip(lon, lat)))



# =============================================================================
# Data extractors
# =============================================================================



def extract_pointvalues(mapinfo, lat=[], lon=[]):
    
    # make google earth objects 
    raster = get_ee_obj(mapinfo) #dataset
    point = coords_to_geometry(lat, lon) #location
    
    #extract point out of dataset
    nested_values = raster.getRegion(geometry=point, scale=mapinfo['scale']).getInfo()
   
    #convert to dataframe
    val_df = pd.DataFrame(nested_values[1:], columns=nested_values[0])
    
    
    #map to human space if categorical
    
    if mapinfo['value_type'] == 'categorical':
        val_df[mapinfo['band_of_use']] = val_df[mapinfo['band_of_use']].map(mapinfo['categorical_mapper'])
    
    
    return val_df[mapinfo['band_of_use']].to_list(), val_df