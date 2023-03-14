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



def get_ee_obj(mapinfo, band=None):
    if mapinfo['is_image']:
        obj = ee.Image(mapinfo['location'])
    elif mapinfo['is_imagecollection']:
        if isinstance(band, type(None)):
            obj = ee.ImageCollection(mapinfo['location'])
        else:
            obj = ee.ImageCollection(mapinfo['location']).select(band)
        
    else:
        sys.exit('Map type is not an Image or Imagecollection.')
    return obj




def coords_to_geometry(lat=[], lon=[], proj='EPSG:4326'):
    if len(lat) == 1:
        return ee.Geometry.Point(coords=[lon[0], lat[0]], proj=proj)
    else:
        return ee.Geometry.MultiPoint(list(zip(lon, lat)), proj=proj)

# =============================================================================
# Helpers
# =============================================================================

def _addDate(image):
    """ add the image datetime as a band """
    img_date = ee.Date(image.date())
    img_date = ee.Number.parse(img_date.format('YYYYMMddHHmmss'))
    return image.addBands(ee.Image(img_date).rename('datetime'))



def _df_to_featurescollection(df, loncolname, latcolname):
    """ Convert a dataframe to a featurecollections row-wise"""
    features=[]
    for index, row in df.iterrows():
    #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row[loncolname], row[latcolname]])
    #     construct the attributes (properties) for each point 
        poi_properties = dict(row)
    #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)
    
    return ee.FeatureCollection(features)







# =============================================================================
# Data extractors
# =============================================================================



def extract_pointvalues(metadf, mapinfo, output_column_name, latcolname='lat', loncolname='lon'):
    """
    Extract values for point locations from a GEE dataset. 
    The pointlocations are defined in a dataframe by EPSG:4326 lat lon coordinates. 
    

    A dataframe with the extracted values is returned. 
    The values are mapped to human classes if the dataset value type is labeld as categorical.


    Parameters
    ----------
    metadf : pd.DataFrame
        dataframe containing coordinates and a column "name", representing the name for each location.
    mapinfo : Dict
        The information about the GEE dataset.
    output_column_name : String
        Column name for the extracted values.
    latcolname : String, optional
        Columnname of latitude values. The default is 'lat'.
    loncolname : String, optional
        Columnname of longitude values. The default is 'lon'.

    Returns
    -------
    pd.DataFrame
        A dataframe with name as index, all columns from the metadf + extracted extracted values column.

    """
    scale=mapinfo['scale']
    
    # =============================================================================
    # df to featurecollection
    # =============================================================================
    
    ee_fc = _df_to_featurescollection(metadf, loncolname, latcolname)


        
    
    # =============================================================================
    # extract raster values
    # =============================================================================

    def rasterExtraction(image):
        feature = image.sampleRegions(
            collection = ee_fc, # feature collection here
            scale = scale # Cell size of raster
        )
        return feature

    
    raster = get_ee_obj(mapinfo, mapinfo['band_of_use']) #dataset
    results = raster.map(rasterExtraction) \
                    .flatten() \
                    .getInfo()         
                    
                    
    # =============================================================================
    # to dataframe
    # =============================================================================
    
    # extract properties
    properties = [x['properties'] for x in results['features']]
    df = pd.DataFrame(properties)
    
    #map to human space if categorical
    if mapinfo['value_type'] == 'categorical':
        df[mapinfo['band_of_use']] = df[mapinfo['band_of_use']].map(mapinfo['categorical_mapper'])
        
    #rename to values to toolkit space
    df = df.rename(columns={mapinfo['band_of_use']: output_column_name})
    
    # #format index
    df = df.set_index(['name'])
   
   
    return df




def gee_extract_timeseries(metadf, mapinfo, startdt, enddt, obstype='temp', latcolname='lat', loncolname='lon'):
    """
    Extract a timeseries, for a given obstype, for point locations from a GEE dataset. The pointlocations are defined in a dataframe by EPSG:4326 lat lon coordinates. 
    
    The startdate is included, the enddate is excluded. 
    
    A multi-index dataframe with the timeseries is returned

    Parameters
    ----------
    metadf : pd.DataFrame
        dataframe containing coordinates and a column "name", representing the name for each location.
    mapinfo : Dict
        The information about the GEE dataset.
    startdt : datetime obj
        Start datetime for timeseries (included). 
    enddt : datetime obj
        End datetime for timeseries (excluded).
    obstype : String, optional
        toolkit observation type. The default is 'temp'.
    latcolname : String, optional
        Columnname of latitude values. The default is 'lat'.
    loncolname : String, optional
        Columnname of longitude values. The default is 'lon'.

    Returns
    -------
    pd.DataFrame
        A dataframe with name - datetime multiindex, all columns from the metadf + extracted timeseries
        column with the same name as the obstype.

    """
    
    
    
    scale=mapinfo['scale']
    bandname=mapinfo['band_of_use'][obstype]['name']
    

    # =============================================================================
    # df to featurecollection
    # =============================================================================
    
    ee_fc = _df_to_featurescollection(metadf, loncolname, latcolname)


        
    
    # =============================================================================
    # extract raster values
    # =============================================================================

    def rasterExtraction(image):
        feature = image.sampleRegions(
            collection = ee_fc, # feature collection here
            scale = scale # Cell size of raster
        )
        return feature

    
    raster = get_ee_obj(mapinfo, bandname) #dataset
    results = raster.filter(ee.Filter.date(startdt, enddt)) \
                    .map(_addDate) \
                    .map(rasterExtraction) \
                    .flatten() \
                    .getInfo()         

    # =============================================================================
    # to dataframe
    # =============================================================================
    
    # extract properties
    properties = [x['properties'] for x in results['features']]
    df = pd.DataFrame(properties)
    
    #format datetime
    df['datetime'] = pd.to_datetime(df['datetime'], format='%Y%m%d%H%M%S')
    
    #format index
    df = df.set_index(['name', 'datetime'])
    df = df.sort_index()
    
    #rename to values to toolkit space
    df = df.rename(columns={bandname: obstype})

    return df
