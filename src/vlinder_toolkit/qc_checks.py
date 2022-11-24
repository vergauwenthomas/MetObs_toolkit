#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""
import pandas as pd
import numpy as np
from datetime import timedelta
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



def step(input_series, obstype='temp', ignore_val=np.nan):
    
    try:
         specific_settings = check_settings['step'][obstype]
    except:
        print('No step settings found for obstype=', obstype, '. Check is skipped!') 
         # return station
        qc_flags = pd.Series('not checked', index=input_series.index)
        qc_flags.name = 'step'
        return input_series, qc_flags.name
    
    

    
    #Split into two sets
    to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
                                                                 ignore_val)
    
    
    outl_datetimes = []
    increaments = to_check_series.shift(1) - to_check_series
    outl_datetimes = increaments[abs(increaments) > specific_settings['max_value']].index.to_list()
   
    
    updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                               ignored_obs_series=to_ignore_series, 
                                                               outlier_dt_list=outl_datetimes,
                                                               checkname='step',
                                                               outlier_label=observation_labels['step'],
                                                               outlier_value=outlier_values['step'],
                                                               obstype=obstype)
    if not to_check_series.empty:
        updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                                   ignored_obs_series=to_ignore_series, 
                                                                   outlier_dt_list=[to_check_series.index[0]],
                                                                   checkname='step',
                                                                   outlier_label='not checked',
                                                                   outlier_value=outlier_values['step'],
                                                                   obstype=obstype)
    
    return updated_obs, qc_flags



def dew_point(temp_hum_dataframe, spec_settings):
    if not (temp_hum_dataframe['humidity'] == 0):
        dew_temp = spec_settings['c']*np.log(temp_hum_dataframe['humidity']/100 *np.exp((spec_settings['b'] - temp_hum_dataframe['temp']/spec_settings['d']) * (temp_hum_dataframe['temp']/(spec_settings['c'] + temp_hum_dataframe['temp']))))/(spec_settings['b'] - np.log(temp_hum_dataframe['humidity']/100 *np.exp((spec_settings['b'] - temp_hum_dataframe['temp']/spec_settings['d']) * (temp_hum_dataframe['temp']/(spec_settings['c'] + temp_hum_dataframe['temp'])))))
        return dew_temp
    else:
        return -999

def internal_consistency(input_series, humidity_series, obstype='temp', ignore_val=np.nan):
    
    try:
         specific_settings = check_settings['internal_consistency'][obstype]
    except:
        print('No internal_consistency settings found for obstype=', obstype, '. Check is skipped!') 
         # return station
        qc_flags = pd.Series('not checked', index=input_series.index)
        qc_flags.name = 'internal_consistency'
        return input_series, qc_flags.name
    
    #Split into two sets
    to_check_series, to_ignore_series = split_to_check_to_ignore(input_series,
                                                                 ignore_val)

    if not to_check_series.empty:
        max_value = to_check_series.rolling("H", center=True).max().rename('max')
        min_value = to_check_series.rolling("H", center=True).min().rename('min')
        
        temp_hum_data = pd.concat([to_check_series, humidity_series, max_value, min_value], axis=1)
        temp_hum_data['dew_temp'] = temp_hum_data.apply(dew_point, spec_settings=specific_settings, axis=1)
        datetimes_no_hum = temp_hum_data[temp_hum_data['humidity'] == 0].index.to_list()
        outl_datetimes = temp_hum_data[(temp_hum_data['min'] > temp_hum_data['temp']) | (temp_hum_data['max'] < temp_hum_data['temp']) | (temp_hum_data['temp'] < temp_hum_data['dew_temp'])].index.to_list()

    if to_check_series.empty:
         datetimes_no_hum = []
         outl_datetimes = []
         
                
                
    updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                                   ignored_obs_series=to_ignore_series, 
                                                                   outlier_dt_list=outl_datetimes,
                                                                   checkname='internal_consistency',
                                                                   outlier_label=observation_labels['internal_consistency'],
                                                                   outlier_value=outlier_values['internal_consistency'],
                                                                   obstype=obstype)
    
    updated_obs, qc_flags = make_checked_obs_and_labels_series(checked_obs_series=to_check_series,
                                                                   ignored_obs_series=to_ignore_series, 
                                                                   outlier_dt_list=datetimes_no_hum,
                                                                   checkname='internal_consistency',
                                                                   outlier_label='humidity 0',
                                                                   outlier_value=outlier_values['internal_consistency'],
                                                                   obstype=obstype)
        
    return updated_obs, qc_flags


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
    
    
    