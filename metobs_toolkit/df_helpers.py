#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of functions on dataframe that are often used.

Created on Thu Mar  2 16:00:59 2023

@author: thoverga
"""

import pandas as pd
import numpy as np
import geopandas as gpd
import itertools






def init_multiindex():
     return pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])
def init_multiindexdf():
    return pd.DataFrame(index = init_multiindex())

def init_triple_multiindex():
    my_index = pd.MultiIndex(levels=[['name'],['datetime'],['obstype']],
                             codes=[[],[],[]],
                             names=[u'name', u'datetime', u'obstype'])
    return my_index
def init_triple_multiindexdf():
    return pd.DataFrame(index=init_triple_multiindex())


def format_outliersdf_to_doubleidx(outliersdf):
    """
    Convert outliersdf to multiindex dataframe if needed.

    This is applied when the obstype level in the index is not relevant.


    Parameters
    ----------
    ouliersdf : Dataset.outliersdf
        The outliers dataframe to format to name - datetime index.

    Returns
    -------
    pandas.DataFrame
        The outliersdfdataframe where the 'obstype' level is dropped, if it was present.

    """


    if 'obstype' in outliersdf.index.names:
        return outliersdf.droplevel('obstype')
    else:
        return outliersdf





def remove_outliers_from_obs(obsdf, outliersdf):
    #TODO this function can only be used with care!!!
    # because all timestamps will be removed that have an oulier in one specific obstype !!!!
    return obsdf.loc[~obsdf.index.isin(outliersdf.index)]


def conv_tz_multiidxdf(df, timezone):

    df.index = df.index.set_levels(df.index.levels[1].tz_convert(timezone), level=1)
    return df



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





def multiindexdf_datetime_subsetting(df, starttime, endtime):
    " The multiindex equivalent of datetime_subsetting"
    dt_df = df.reset_index().set_index('datetime')
    subset_dt_df = datetime_subsetting(dt_df, starttime, endtime)

    # back to multiindex name-datetime
    subset_dt_df = subset_dt_df.reset_index()
    idx = pd.MultiIndex.from_frame(subset_dt_df[['name', 'datetime']])
    returndf = subset_dt_df.set_index(idx).drop(columns=['name', 'datetime'], errors='ignore')

    if returndf.empty:
        print(f'Warning: No observations left after subsetting datetime {starttime} -- {endtime} ')

    return returndf

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


def conv_applied_qc_to_df(obstypes, ordered_checknames):
    if isinstance(obstypes, str):
        obstypes = [obstypes]
    if isinstance(ordered_checknames, str):
        ordered_checknames = [ordered_checknames]

    obslist = list(itertools.chain.from_iterable(itertools.repeat(item, len(ordered_checknames)) for item in obstypes))

    checknamelist = list(itertools.chain.from_iterable(itertools.repeat(ordered_checknames, len(obstypes))))

    df = pd.DataFrame({'obstype': obslist,
                       'checkname': checknamelist})
    return df

