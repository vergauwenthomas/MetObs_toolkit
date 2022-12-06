#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""
import pandas as pd
import numpy as np

from datetime import timedelta

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
    """
    This functions splits an input_series in two sets based on if the value is
    equal or not equal to the ignore_val.

    Parameters
    ----------
    input_series : pd.Series
        The series to be splits in two sets.
    ignore_val : Value, optional
        The value used to split the input series on. The default is np.nan.

    Returns
    -------
    to_check_input_series : pandas.Series
        Subset of the input_series for which the value != ignore_val.
    ignore_input_series : pandas.Series
        Subset of the input_series for which the value != ignore_val.

    """
    if np.isnan(ignore_val):
         #better handling for nan values
         ignore_input_series = input_series[input_series.isna()]
         to_check_input_series = input_series[input_series.notna()]
    else:  
         ignore_input_series = input_series[input_series == ignore_val]
         to_check_input_series = input_series[input_series != ignore_val]
    
    return to_check_input_series, ignore_input_series


def make_checked_obs_and_labels_series(checked_obs_series,ignored_obs_series, 
                                       outlier_obs, checkname, outlier_label,
                                       outlier_value, obstype,
                                       not_checked_label='not checked', ok_label='ok'):
    """
    This function combines all the input series to :
        * updated observation series
        * series containing the quality flags

    Parameters
    ----------
    checked_obs_series : pandas.Series
        Series with observations that are checked.
    ignored_obs_series : pandas.Series
        Series with observations that are not checked.
    outlier_obs : list
        List of indices of the observation series that are flagged as outliers.
    checkname : str
        Default name of the check.
    outlier_label : str
        Default name of the label series.
    outlier_value : str
        Observation value of an outlier.
    obstype : str
        Default name of the observationtype.
    not_checked_label : str, optional
        Label for observations that ar not checked. The default is 'not checked'.
    ok_label : str, optional
        Label for the observations that are not flagged as outliers. The default is 'ok'.

    Returns
    -------
    pandas.Series
        The updated observation values.
    pandas.Series
        The quality labels for the observations.

    """
    flag_column_name = obstype + '_' + checkname + '_' + 'label'
    
    #create labels for the checked observations
    df_checked = checked_obs_series.to_frame()
    df_checked[flag_column_name] = ok_label #Set all checked labels as ok to start
    
    
    #The at gives problems when running in python >= 3.8  
    # #update flags of outlier observations
    # df_checked.at[df_checked.index.isin(outlier_obs), flag_column_name] = outlier_label
    # #convert observations of outliers
    # df_checked.at[df_checked.index.isin(outlier_obs), obstype] = outlier_value
    
    #update flags of outlier observations
    df_checked.loc[df_checked.index.isin(outlier_obs), flag_column_name] = outlier_label
    #convert observations of outliers
    df_checked.loc[df_checked.index.isin(outlier_obs), obstype] = outlier_value
    
    
    
    
    
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
    Looking for missing timestaps by assuming an observation frequency. The assumed frequency is the most occuring frequency PER STATION.
    If missing observations are detected, the observations dataframe is extended by these missing timestamps with Nan's as filling values.

    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)

    Returns
    -------
    df : pandas.DataFrame()
        The observations dataframe updated for missing timestamps (values updated + quality flag column added).

    """     
    
    checkname = 'missing_timestamp'
    
   
    flag_column = checkname + '_' + 'label'
    df[flag_column] = observation_labels['ok']
    
    stationnames = df.index.get_level_values(level='name').unique()
    for station in stationnames:
        timestamps = df.xs(station, level='name').index
        likely_freq = min(abs(timestamps.to_series().diff().value_counts().index))

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
    """
    Looking for duplcate timestaps per station. Duplicated records are removed by the method specified in the qc_settings. 

    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)

    Returns
    -------
    df : pandas.DataFrame()
        The observations dataframe updated for missing timestamps (values updated + quality flag column added).

    """     
    
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



def step_check(input_series, obstype='temp', ignore_val=np.nan):   
    """
    @MICHIEL

    Parameters
    ----------
    input_series : TYPE
        DESCRIPTION.
    obstype : TYPE, optional
        DESCRIPTION. The default is 'temp'.
    ignore_val : TYPE, optional
        DESCRIPTION. The default is np.nan.

    Returns
    -------
    TYPE
        DESCRIPTION.
    TYPE
        DESCRIPTION.

    """
    checkname='step'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        # logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        flag_series = flag_series_when_no_settings_found(input_series, obstype, checkname)
        return input_series, flag_series
    
    
    
    
    #Do not split in advace, so no gabs appear in the dataset!! Split afterwards for labeling consistency
    
    outl_obs = []
    for _station, obs in input_series.groupby(level='name'):
        increaments = input_series.shift(1) - input_series
        outl_obs.extend(increaments[abs(increaments) > specific_settings['max_value']].index.to_list())
        
    # Split into two sets
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




def compute_dew_point(df, spec_settings):
    dew_temp = spec_settings['c']*np.log(df['humidity']/100 *np.exp((spec_settings['b'] - df['temp']/spec_settings['d']) * (df['temp']/(spec_settings['c'] + df['temp']))))/(spec_settings['b'] - np.log(df['humidity']/100 *np.exp((spec_settings['b'] - df['temp']/spec_settings['d']) * (df['temp']/(spec_settings['c'] + df['temp'])))))
    return dew_temp
  




def internal_consistency_check(input_series, humidity_series, obstype='temp', ignore_val=np.nan):
    checkname = 'internal_consistency'
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        flag_series = flag_series_when_no_settings_found(input_series, obstype, checkname)
        return input_series, flag_series
    
    #Check if humidity is available
    if humidity_series.isnull().all():
        print(f'{checkname} could not be performed, no humidity observations found.')
        logger.warning(f'{checkname} could not be performed, no humidity observations found.')
        flag_series = flag_series_when_no_settings_found(input_series, obstype, checkname)
        return input_series, flag_series
    
    #Split series when humidity is not zero and when observatiosn are not Nan
    comb_df = input_series.to_frame()
    comb_df['humidity'] = humidity_series
    
    if np.isnan(ignore_val):
         #better handling for nan values
         ignore_mask = comb_df[(comb_df[obstype].isna()) | (comb_df['humidity']==0)].index
         check_mask =  comb_df[(~comb_df[obstype].isna()) & (comb_df['humidity']!=0)].index
    else:  
         ignore_mask = comb_df[(comb_df[obstype]==ignore_val) | (comb_df['humidity']==0)].index
         check_mask = comb_df[(comb_df[obstype]!=ignore_val) & (comb_df['humidity']!=0)].index
    
    
    to_check_df= comb_df.loc[check_mask]
    to_ignore_df = comb_df.loc[ignore_mask]
    
    #compute dewpoint temperature
    to_check_df['dew_temp'] = compute_dew_point(to_check_df, specific_settings)
    
    #compute hourly rolling max and min temp for groupby station (to avoid overlap between stations)
    rolling_agg= to_check_df.reset_index(level=0).groupby('name')[obstype].rolling('H', center=True).agg([np.max, np.min])
    to_check_df['rolling_max'], to_check_df['rolling_min'] = rolling_agg['amax'], rolling_agg['amin']
    
    #get outliers
    #TODO @Michiel : Zeker dat deze defenitie klopt? 
    outl_obs = to_check_df[(to_check_df['rolling_min'] > to_check_df[obstype]) |
                                 (to_check_df['rolling_max'] < to_check_df[obstype]) |
                                 (to_check_df[obstype] < to_check_df['dew_temp'])].index.to_list()
    
         
    updated_obs_series, qc_flags_series = make_checked_obs_and_labels_series(checked_obs_series=to_check_df[obstype],
                                                                    ignored_obs_series=to_ignore_df[obstype], 
                                                                    outlier_obs=outl_obs,
                                                                    checkname=checkname,
                                                                    outlier_label=observation_labels[checkname],
                                                                    outlier_value=outlier_values[checkname],
                                                                    obstype=obstype)
    
    return updated_obs_series, qc_flags_series

# def dew_point(temp_hum_dataframe, spec_settings):
#     if not (temp_hum_dataframe['humidity'] == 0):
#         dew_temp = spec_settings['c']*np.log(temp_hum_dataframe['humidity']/100 *np.exp((spec_settings['b'] - temp_hum_dataframe['temp']/spec_settings['d']) * (temp_hum_dataframe['temp']/(spec_settings['c'] + temp_hum_dataframe['temp']))))/(spec_settings['b'] - np.log(temp_hum_dataframe['humidity']/100 *np.exp((spec_settings['b'] - temp_hum_dataframe['temp']/spec_settings['d']) * (temp_hum_dataframe['temp']/(spec_settings['c'] + temp_hum_dataframe['temp'])))))
#         return dew_temp
#     else:
#         return -999

# def internal_consistency(input_series, humidity_series, obstype='temp', ignore_val=np.nan):
    
#     try:
#          specific_settings = check_settings['internal_consistency'][obstype]
#     except:
#         print('No internal_consistency settings found for obstype=', obstype, '. Check is skipped!') 
#          # return station
#         qc_flags = pd.Series('not checked', index=input_series.index)
#         qc_flags.name = 'internal_consistency'
#         return input_series, qc_flags.name
    
#     #Split into two sets
#     to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
#                                                                  ignore_val)

#     if not to_check_series.empty:
#         max_value = to_check_series.rolling("H", center=True).max().rename('max')
#         min_value = to_check_series.rolling("H", center=True).min().rename('min')
        
#         temp_hum_data = pd.concat([to_check_series, humidity_series, max_value, min_value], axis=1)
#         temp_hum_data['dew_temp'] = temp_hum_data.apply(dew_point, spec_settings=specific_settings, axis=1)
#         datetimes_no_hum = temp_hum_data[temp_hum_data['humidity'] == 0].index.to_list()
#         outl_datetimes = temp_hum_data[(temp_hum_data['min'] > temp_hum_data['temp']) | (temp_hum_data['max'] < temp_hum_data['temp']) | (temp_hum_data['temp'] < temp_hum_data['dew_temp'])].index.to_list()

#     if to_check_series.empty:
#          datetimes_no_hum = []
#          outl_datetimes = []
         
                
                
#     updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
#                                                                    ignored_obs_series=to_ignore_series, 
#                                                                    outlier_dt_list=outl_datetimes,
#                                                                    checkname='internal_consistency',
#                                                                    outlier_label=observation_labels['internal_consistency'],
#                                                                    outlier_value=outlier_values['internal_consistency'],
#                                                                    obstype=obstype)
    
#     updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
#                                                                    ignored_obs_series=to_ignore_series, 
#                                                                    outlier_dt_list=datetimes_no_hum,
#                                                                    checkname='internal_consistency',
#                                                                    outlier_label='humidity 0',
#                                                                    outlier_value=outlier_values['internal_consistency'],
#                                                                    obstype=obstype)
        
#     return updated_obs, qc_flags


def qc_info(qc_dataframe):
    
    number_gross_error_outliers = qc_dataframe.persistance.str.contains('gross value outlier').sum()
    number_persistance_outliers = qc_dataframe.persistance.str.contains('persistance outlier').sum()
    number_step_outliers = qc_dataframe.persistance.str.contains('step outlier').sum()
    number_int_consistency_outliers = qc_dataframe.persistance.str.contains('internal consistency outlier').sum()
    
    fraction_gross_outliers = number_gross_error_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
    fraction_persistance_outliers = number_persistance_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
    fraction_step_outliers = number_step_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
    fraction_int_consistency_outliers = number_int_consistency_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
    
    percentage_of_data_gross_outlier = number_gross_error_outliers/len(qc_dataframe)
    percentage_of_data_persistance_outlier = number_persistance_outliers/len(qc_dataframe)
    percentage_of_data_step_outlier = number_step_outliers/len(qc_dataframe)
    percentage_of_data_consistency_outlier = number_int_consistency_outliers/len(qc_dataframe)
    
    print("Number of gross error outliers: ", number_gross_error_outliers)
    print("Number of persistance outliers: ", number_persistance_outliers)
    print("Number of step outliers: ", number_step_outliers)
    print("Number of internal consistency outliers: ", number_int_consistency_outliers)
    
    print("Fraction of gross error outliers: ", fraction_gross_outliers)
    print("Fraction of persistance outliers: ", fraction_persistance_outliers)
    print("Fraction of step outliers: ", fraction_step_outliers)
    print("Fraction of internal consistency outliers: ", fraction_int_consistency_outliers)
    
    print("Percentage of data gross outlier: ", percentage_of_data_gross_outlier)
    print("Percentage of data persistance outlier: ", percentage_of_data_persistance_outlier)
    print("Percentage of data step outlier: ", percentage_of_data_step_outlier)
    print("Percentage of data internal consistency outlier: ", percentage_of_data_consistency_outlier)
    
    print("Gross error outliers at: ", qc_dataframe[qc_dataframe['gross_value'] == 'gross value outlier'].index)
    print("Persistance outliers at: ", qc_dataframe[qc_dataframe['persistance'] == 'persistance outlier'].index)
    print("Step outliers at: ", qc_dataframe[qc_dataframe['step'] == 'step outlier'].index)
    print("Internal consistency outliers at: ", qc_dataframe[qc_dataframe['internal_consistency'] == 'internal consistency outlier'].index)
    print("Humidity 0 at: ", qc_dataframe[qc_dataframe['internal_consistency'] == 'humidity 0'].index)
    
    


