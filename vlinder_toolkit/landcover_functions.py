#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 11:28:36 2022

@author: thoverga
"""


import sys
import pandas as pd
import ee

    

# gee_datasets = {
#     'global_lcz_map':{'location': "RUB/RUBCLIM/LCZ/global_lcz_map/v1", #GEE location
#                       'usage': 'LCZ', #Human readable application domain
#                       'band_of_use': 'LCZ_Filter', #band to use for imagecollections (or None if no band available)
#                       'value_type': 'categorical', #categorical or numeric
#                       'dynamical': False, #time evolution? To be used for timeseries
#                       'scale': 1000, 
#                       'is_image': False,
#                       'is_imagecollection': True,
#                       'categorical_mapper': {
#                          1: 'Compact highrise', #mapvalue: (color, human class)
#                          2: 'Compact midrise',
#                          3:'Compact lowrise',
#                          4:	'Open highrise',
#                          5:	'Open midrise',
#                          6:'Open lowrise',
#                          7:	'Lightweight lowrise',
#                          8:	'Large lowrise',
#                          9:	'Sparsely built',
#                          10: 'Heavy industry',
#                          11: 'Dense Trees (LCZ A)',
#                          12: 'Scattered Trees (LCZ B)',
#                          13:'Bush, scrub (LCZ C)',
#                          14:'Low plants (LCZ D)',
#                          15:'Bare rock or paved (LCZ E)',
#                          16:'Bare soil or sand (LCZ F)',
#                          17:'Water (LCZ G)',                          
#                           }
#                       },
    
#     }



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
# geodataset functions
# =============================================================================



def geotiff_point_extraction():
    print("not yet implemented")
  



# =============================================================================
#  Connection functions
# =============================================================================

def connect_to_gee():
    if not ee.data._credentials: #check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return





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