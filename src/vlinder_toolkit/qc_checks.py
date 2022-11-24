#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""
import pandas as pd
import numpy as np
import logging

# from .stations import Station, Dataset
from .settings import Settings
from .settings_files.qc_settings import check_settings, outlier_values, observation_labels

logger = logging.getLogger(__name__)

def flag_series_when_no_settings_found(inputseries, obstype, checkname):
    """
    Apply this function when the check could not be executed. The check-labels series is returned.

    Parameters
    ----------
    inputseries : pd.Series
        The observations to be checked
    obstype : String
        Name of the observationtype that is checked. This should be the name of the inputseries
    checkname : String
        Default name of the check

    Returns
    -------
    flag_series : pd.Series
        The flags for the observations with the same index as the inputseries.

    """
    flag_series = pd.Series(data='not checked',
                            index=inputseries.index
                            )
    flag_series.name = obstype + '_' + checkname + '_' + 'label'
    return flag_series


def split_to_check_to_ignore(input_series, ignore_val=np.nan):
    if np.isnan(ignore_val):
         #better handling for nan values
         ignore_input_series = input_series[input_series.isna()]
         to_check_input_series = input_series[input_series.notna()]
    else:  
         ignore_input_series = input_series[input_series == ignore_val]
         to_check_input_series = input_series[input_series != ignore_val]
    
    return to_check_input_series, ignore_input_series


def make_checked_obs_and_labels_series(checked_obs_series,
                                          ignored_obs_series, 
                                          outlier_obs,
                                          checkname,
                                          outlier_label,
                                          outlier_value,
                                          obstype,
                                          not_checked_label='not checked',
                                          ok_label='ok'):
    
    flag_column_name = obstype + '_' + checkname + '_' + 'label'
    
    #create labels for the checked observations
    df_checked = checked_obs_series.to_frame()
    df_checked[flag_column_name] = ok_label #Set all checked labels as ok to start
    
    
    #update flags of outlier observations
    df_checked.at[df_checked.index.isin(outlier_obs), flag_column_name] = outlier_label
    #convert observations of outliers
    df_checked.at[df_checked.index.isin(outlier_obs), obstype] = outlier_value
    
    
    #create labels for ignored observations
    df_ignored = ignored_obs_series.to_frame()
    df_ignored[flag_column_name] = not_checked_label #Set all checked labels as ok to start
    
    
    #concat checked and ignored and sort
    df_tot = pd.concat([df_checked, df_ignored])
    df_tot = df_tot.sort_index()    
    
    
    return df_tot[obstype], df_tot[flag_column_name]
  


    
    
    
# =============================================================================
# Quality assesment checks on data import
# =============================================================================
def missing_timestamp_check(df):
    """
    V2
    Looking for missing timestaps by assuming an observation frequency. The assumed frequency is the most occuring frequency.
    If missing observations are detected, the observations dataframe is extended by these missing timestamps with Nan's as filling values.

    Parameters
    ----------
    station : Station object
        The station you whant to apply this check on.

    Returns
    -------
    df : pandas.DataFrame()
        The observations dataframe (same as Station.df()).
    missing_datetimes : list of datetimes
        The list of the missing timestamps.

    """     
    
    checkname = 'missing_timestamp'
    
   
    flag_column = checkname + '_' + 'label'
    df[flag_column] = observation_labels['ok']
    
    
    #missing timestamp per station (because some stations can have other frequencies!)
    stationnames = df.index.get_level_values(level='name').unique()
    for station in stationnames:
        timestamps = df.xs(station, level='name').index
        likely_freq =timestamps.to_series().diff().value_counts().idxmax()
        
        
        missing_datetimeindices = pd.date_range(start = timestamps.min(),
                                                end = timestamps.max(),
                                                freq=likely_freq).difference(timestamps)
        
        if not missing_datetimeindices.empty:
            logging.warning(f'{len(missing_datetimeindices)} missing records ({missing_datetimeindices[:10]} ...) found for {station}. These will be filled with Nans.')
            missing_records_df = pd.DataFrame(columns=df.columns,
                                              index=pd.MultiIndex.from_arrays([[station]*len(missing_datetimeindices),
                                                                              missing_datetimeindices.to_list()]))
            missing_records_df[flag_column] = observation_labels[checkname]
            
            
            df = pd.concat([df, missing_records_df])
    
            
    df = df.sort_index()
    return df
    

        
def duplicate_timestamp_check(df):

    checkname = 'duplicate_timestamp'
    
   
    duplicates = pd.Series(data=df.index.duplicated(keep=check_settings[checkname]['keep']),
                           index=df.index)
    
    if not df.loc[duplicates].empty:
        logging.warning(f' Following records are labeld as duplicates: {df.loc[duplicates]}, and are removed')
    
    df = df[~df.index.duplicated(keep=check_settings[checkname]['keep'])]

    return df
# =============================================================================
# Quality assesment checks on dataset
# =============================================================================


def gross_value_check(input_series, obstype, ignore_val):
    checkname = 'gross_value'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        flag_series = flag_series_when_no_settings_found(input_series, obstype, checkname)
        return input_series, flag_series
   
    
    #Split into two sets
    to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
                                                                 ignore_val)
    
    
    
    
    #find outlier observations as a list of tuples [(name, datetime), (name, datetime)]
    outl_obs = to_check_series.loc[(to_check_series <= specific_settings['min_value']) | 
                                          (to_check_series >= specific_settings['max_value'])
                                          ].index.to_list()
    
    
         
    # #Update observations and quality flags
    updated_obs_series, qc_flags_series = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                                ignored_obs_series=to_ignore_series, 
                                                                outlier_obs=outl_obs,
                                                                checkname=checkname,
                                                                outlier_label=observation_labels[checkname],
                                                                outlier_value=outlier_values[checkname],
                                                                obstype=obstype)
    return updated_obs_series, qc_flags_series






def persistance_check(input_series, obstype, ignore_val=np.nan):
    checkname = 'persistance'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        flag_series = flag_series_when_no_settings_found(input_series, obstype, checkname)
        return input_series, flag_series
    
    
    
    
    #Split into two sets
    # to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
    #                                                              ignore_val)
    #  Because np.nan != np.nan, the nan's have always an persistance counter of 1
    #  To make shure that the shift() is one hour late, it is thus easier to apply this
    #  check on all input observations, and split the dataset later to be consistent
    #  for the labels.
    
    
    #find outlier datetimes
    
    #make consec groups
    grouped = input_series.groupby(['name', (input_series.shift() != input_series).cumsum()]) 
    #the above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = grouped.size()
    outlier_groups = group_sizes[group_sizes > specific_settings['max_valid_repetitions']]

    
    #add to outl_obs.
    outl_obs = []
    for group_idx in outlier_groups.index:
        outl_obs.extend(grouped.get_group(group_idx).index.to_list())
    
    
    #only used for labeling consistency
    to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
                                                                  ignore_val)
    
    
    
    
    # #Update observations and quality flags
    updated_obs_series, qc_flags_series = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                                ignored_obs_series=to_ignore_series, 
                                                                outlier_obs=outl_obs,
                                                                checkname=checkname,
                                                                outlier_label=observation_labels[checkname],
                                                                outlier_value=outlier_values[checkname],
                                                                obstype=obstype)
    
    return updated_obs_series, qc_flags_series
