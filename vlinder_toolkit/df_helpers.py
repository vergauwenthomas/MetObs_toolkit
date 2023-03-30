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




def init_multiindexdf():
    my_index = pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])


    df = pd.DataFrame(index=my_index)
    return df


def init_multiindex():
     return pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])



def add_final_label_to_outliersdf(outliersdf, data_res_series, observation_types, checks_info):
    """
    V3
        This function creates a final label based on de individual qc labels. The final label will be that of the individual qc-label
        which rejected the obseration.

        This functions converts labels to numeric values, algebra to get final label, and inversly
        convert to labels. This is faster than looping over the rows.

        Parameters
        ----------
        outliersdf : pandas.DataFrame
            The dataset outliers dataframe containing the observations and QC labels.

        data_res_series : Pandas.Series
            The series that contain the dataset resolution (values) per station (index). This
            is stored in the dataset.metadf as column 'dataset_resolution'. These are used to explode the gaps.

        Returns
        -------
        outliersdf : pd.DataFrame
            The outliersdf with extra columns indicated by example 'temp_final_label' and 'humid_final_label'.

        """




    # order columns
    labels_columns = [column for column in outliersdf.columns if not column in observation_types]
    #drop final columns if they are in the outliersdf
    labels_columns = [column for column in labels_columns if not column.endswith('_final_label')]

    checked_obstypes = [obstype for obstype in observation_types if any([qc_column.startswith(obstype+'_') for qc_column in labels_columns])]
    columns_on_record_lvl = [info['label_columnname'] for checkname, info in checks_info.items() if info['apply_on'] == 'record']


    # Construct numeric mapper
    labels_to_numeric_mapper = {info['outlier_flag']:info['numeric_flag'] for info in checks_info.values()}

    # add 'ok' and 'not checked' labels
    labels_to_numeric_mapper['ok'] = 0
    labels_to_numeric_mapper['not checked'] = np.nan
    #invert numeric mapper
    inv_label_to_num = {v: k for k, v in labels_to_numeric_mapper.items()}


    #generete final label per obstype
    for obstype in checked_obstypes:
        # logger.debug(f'Generating final QC labels for {obstype}.')
        #Get qc column namse specific for this obstype
        specific_columns = [col for col in labels_columns if col.startswith(obstype+'_')]
        #add qc labels that are applicable on all obstypes
        specific_columns.extend(columns_on_record_lvl)

        #Drop columns that are not present
        specific_columns = [colmn for colmn in specific_columns if colmn in outliersdf.columns]


        #get labels dataframe
        qc_df = outliersdf[specific_columns]
        num_qc_df = pd.DataFrame()
        num_qc_df = qc_df.applymap(labels_to_numeric_mapper.get )

        outliersdf[obstype+'_final_label'] = num_qc_df.sum(axis=1, skipna=True).map(inv_label_to_num)


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
