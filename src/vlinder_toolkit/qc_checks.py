#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""
import pandas as pd
import numpy as np

# from .stations import Station, Dataset
from .settings import Settings
from .settings_files.qc_settings import check_settings, outlier_values, observation_labels


def make_checked_obs_and_labels_series(checked_obs_series,
                      ignored_obs_series, 
                      outlier_dt_list,
                      checkname,
                      outlier_label,
                      outlier_value,
                      obstype,
                      not_checked_label='not checked',
                      ok_label='ok'):
    
    
    
    #create labels for the checked observations
    df_checked = checked_obs_series.to_frame()
    df_checked['label'] = [outlier_label if dt in outlier_dt_list else 
                           ok_label for dt in df_checked.index]
    
    #create labels for ignored observations
    df_ignored = ignored_obs_series.to_frame()
    df_ignored['label'] = not_checked_label
    
    #merge labels df together
    df = pd.concat([df_checked, df_ignored])
    df = df.sort_index()
    
    #convert outl observations to outlier values
    df.loc[df['label'] == outlier_label, obstype] = outlier_value
    

    return df[obstype], df['label']

def split_to_check_to_ignore(input_series, ignore_val=np.nan):
    if np.isnan(ignore_val):
         #better handling for nan values
         ignore_input_series = input_series[input_series.isna()]
         to_check_input_series = input_series[input_series.notna()]
    else:  
         ignore_input_series = input_series[input_series == ignore_val]
         to_check_input_series = input_series[input_series != ignore_val]
    
    return to_check_input_series, ignore_input_series
    

  


    
    
    
# =============================================================================
# Quality assesment checks
# =============================================================================

    

def gross_value(input_series, obstype, ignore_val=np.nan):
   
    try:
        specific_settings = check_settings['gross_value'][obstype]
    except:
        print('No gross_value settings found for obstype=', obstype, '. Check is skipped!') 
        qc_flags = pd.Series('not checked', index=input_series.index)
        qc_flags.name = 'gross_value'
        return input_series, qc_flags.name
    
    #Split into two sets
    to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
                                                                 ignore_val)
    



    #find outlier datetimes
    outl_datetimes = to_check_series.loc[(to_check_series <= specific_settings['min_value']) | 
                                         (to_check_series >= specific_settings['max_value'])
                                         ].index.to_list()
         
    #Update observations and quality flags
    updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                               ignored_obs_series=to_ignore_series, 
                                                               outlier_dt_list=outl_datetimes,
                                                               checkname='gross_value',
                                                               outlier_label=observation_labels['gross_value'],
                                                               outlier_value=outlier_values['gross_value'],
                                                               obstype=obstype)

    return updated_obs, qc_flags





def persistance(input_series, obstype='temp', ignore_val=np.nan):
   
    
    try:
         specific_settings = check_settings['persistance'][obstype]
    except:
        print('No persistance settings found for obstype=', obstype, '. Check is skipped!') 
         # return station
        qc_flags = pd.Series('not checked', index=input_series.index)
        qc_flags.name = 'persistance'
        return input_series, qc_flags.name
    
    

    
    #Split into two sets
    to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
                                                                 ignore_val)


    #find outlier datetimes
    
    #make consec groups
    grouped = to_check_series.groupby((to_check_series.shift() != to_check_series).cumsum()) 
    #the above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = grouped.size()
    outlier_groups = group_sizes[group_sizes > specific_settings['max_valid_repetitions']]
    
    #add datetimes of outlier groups to outl_datetime
    outl_datetimes = []
    for group_idx in outlier_groups.index:
        outl_datetimes.extend(grouped.get_group(group_idx).index.to_list())
    
    
    #Update observations and quality flags
    updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                               ignored_obs_series=to_ignore_series, 
                                                               outlier_dt_list=outl_datetimes,
                                                               checkname='persistance',
                                                               outlier_label=observation_labels['persistance'],
                                                               outlier_value=outlier_values['persistance'],
                                                               obstype=obstype)
    
    return updated_obs, qc_flags
