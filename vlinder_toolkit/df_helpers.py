#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of functions on dataframe that are often used.

Created on Thu Mar  2 16:00:59 2023

@author: thoverga
"""
import pandas as pd
import geopandas as gpd





def expand_gabs_to_multiind(df, gapsdf, obstype='temp'):
    """
    Expand the gaps to timestamps, that are missing in the observations. 
    
    Parameters
    ----------
    df : dataset.df
        The observations where gaps are present as Nan values.
    gapsdf : dataset.gapsdf
        The dataframe with detaild information on the start and end of a gap.
    obstype : String, optional
        Observation type. The default is 'temp'.
    Returns
    -------
    expanded_gabsidx : pd.MultiIndex
        A Station-datetime multiindex of missing gap-records.
    """
    
    expanded_gabsidx = pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])
    
    
    for sta, row in gapsdf.iterrows():
    
        sta_df = df.xs(sta, level='name') #filter by name 
        posible_gaps_dt = sta_df[sta_df[obstype].isnull()].index #filter by missing observations
        gaps_dt = posible_gaps_dt[(posible_gaps_dt >= row['start_gap']) & #filter if the observations are within a gap
                                  (posible_gaps_dt <= row['end_gap'])]
        
        gaps_multiidx = pd.MultiIndex.from_arrays(arrays=[[sta]*len(gaps_dt),
                                                          gaps_dt],
                                                  names=[u'name', u'datetime'])
        
        expanded_gabsidx = expanded_gabsidx.append(gaps_multiidx)
    
    return expanded_gabsidx


    
    





def metadf_to_gdf(df, crs=4326):
    """
    Function to convert a dataframe with 'lat' en 'lon' columnst to a geopandas 
    dataframe with a geometry column containing points.
    
    Special care for stations with missing coordinates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with a 'lat' en 'lon' column.
    crs : Integer, optional
        The epsg number of the coordinates. The default is 4326.

    Returns
    -------
    geodf : geopandas.GeaDataFrame
        The geodataframe equivalent of the df.

    """
    
    # only conver to points if coordinates are present
    coordsdf = df[(~df['lat'].isnull()) & (~df['lon'].isnull())]
    missing_coords_df =  df[(df['lat'].isnull()) | (df['lon'].isnull())]
    
    geodf = gpd.GeoDataFrame(coordsdf,
                              geometry=gpd.points_from_xy(coordsdf.lon,
                                                          coordsdf.lat)) 
    geodf = geodf.set_crs(epsg = crs)
    geodf = pd.concat([geodf, missing_coords_df])
    
    geodf = geodf.sort_index()
    return geodf





def datetime_subsetting(df, starttime, endtime):
    """
    Wrapper function for subsetting a dataframe with datetimeindex with a start- and 
    endtime. 

    Parameters
    ----------
    df : pandas.DataFrame with datetimeindex
        The dataframe to apply the subsetting to.
    starttime : datetime.Datetime
        Starttime for the subsetting period (included).
    endtime : datetime.Datetime
        Endtime for the subsetting period (included).

    Returns
    -------
    pandas.DataFrame
        Subset of the df.

    """
    
    stand_format = '%Y-%m-%d %H:%M:%S'
    
    if isinstance(starttime, type(None)):
        startstring = None #will select from the beginning of the df
    else:
        startstring = starttime.strftime(stand_format)
    if isinstance(endtime, type(None)):
        endstring = None
    else: 
        endstring = endtime.strftime(stand_format)

    return df[startstring: endstring]


