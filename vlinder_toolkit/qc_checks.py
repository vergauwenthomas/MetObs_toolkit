#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""

import sys
import pandas as pd
import numpy as np
import logging


from vlinder_toolkit.settings_files.qc_settings import (check_settings,
                                                        checks_info)





logger = logging.getLogger(__name__)



# =============================================================================
# Helper functions
# =============================================================================
def init_outlier_multiindexdf():
    my_index = pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])

   
    df = pd.DataFrame(index=my_index)
    return df

def make_outlier_df_for_check(station_dt_list, values_in_dict, flagcolumnname, flag, stationname=None, datetimelist=None):
    """
    V3
    Helper function to create an outlier dataframe for the given station(s) and datetimes. This will be returned by 
    a quality control check and later added to the dastes.outlierdf. 
    
    Multiple commum inputstructures can be handles
    
    A multiindex dataframe with the relevant observationtypes i.e. the values_in_dict and a specific quality flag column (i.g. the labels) is returned.
    Parameters
    ----------
    station_dt_list : MultiIndex or list of tuples: (name, datetime)
        The stations with corresponding datetimes that are labeled as outliers.
    values_in_dict : Dictionary
        A Dictionary {columnname: pd.Series()} with the values on which this check is applied.
    flagcolumnname : String
        Name of the labels column
    flag : String
        The label for all the outliers.
    stationname : String, optional
        It is possible to give the name of one station. The default is None.
    datetimelist : DatetimeIndex or List, optional
        The outlier timestamps for the stationname. The default is None.
    Returns
    -------
    check_outliers : pd.Dataframe
        A multiindex (name -- datetime) datatframe with one column (columname) and all
        values are the flag.
    """

        
    columnorder =list(values_in_dict.keys())
    columnorder.append(flagcolumnname)
        
    
    if isinstance(station_dt_list, pd.MultiIndex):
        check_outliers = pd.DataFrame(data=flag, index=station_dt_list, columns=[flagcolumnname])
        for obstype, obsvalues in values_in_dict.items():    
            check_outliers[obstype] = obsvalues
        check_outliers = check_outliers[columnorder]
        return check_outliers
    elif isinstance(station_dt_list, list): #list of tuples: (name, datetime)
        multi_idx = pd.MultiIndex.from_tuples(station_dt_list, names=['name', 'datetime'])
        check_outliers = pd.DataFrame(data=flag, index=multi_idx, columns=[flagcolumnname])
        for obstype, obsvalues in values_in_dict.items():    
            check_outliers[obstype] = obsvalues
        check_outliers = check_outliers[columnorder]
        return check_outliers
    elif not isinstance(stationname, type(None)):
        if isinstance(datetimelist, pd.DatetimeIndex):
            datetimelist = datetimelist.to_list()
        if isinstance(datetimelist, list):
            indexarrays = list(zip([stationname]*len(datetimelist), datetimelist))
            multi_idx = pd.MultiIndex.from_tuples(indexarrays, names=['name', 'datetime'])
            check_outliers = pd.DataFrame(data=flag, index=multi_idx, columns=[flagcolumnname])
            for obstype, obsvalues in values_in_dict.items():    
                check_outliers[obstype] = obsvalues
            check_outliers = check_outliers[columnorder]
            return check_outliers
        else:
            sys.exit(f'Type of datetimelist: {type(datetimelist)} is not implemented.')




    
# =============================================================================
# Quality assesment checks on data import
# =============================================================================


def invalid_input_check(df, obstype):
    
    checkname = 'invalid_input'
    
    nan_df = df[df[obstype].isnull()]
    nan_indices = nan_df.index

    outlierdf = make_outlier_df_for_check(station_dt_list = nan_indices,
                                          values_in_dict = nan_df.to_dict(orient='series'),
                                          flagcolumnname=obstype+'_'+checks_info[checkname]['label_columnname'],
                                          flag=checks_info[checkname]['outlier_flag'])
    
    df = df.drop(nan_indices)
    
    return df, outlierdf
  
def duplicate_timestamp_check(df):
    """
    V3
    Looking for duplcate timestaps per station. Duplicated records are removed by the method specified in the qc_settings. 

    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)

    Returns
    -------
    df : pandas.DataFrame()
        The observations dataframe updated for duplicate timestamps. Duplicated timestamps are removed.

    """     
    
    checkname = 'duplicated_timestamp'

    
    
    duplicates = pd.Series(data=df.index.duplicated(keep=check_settings[checkname]['keep']),
                           index=df.index)
    
    if not df.loc[duplicates].empty:
        logging.warning(f' Following records are labeld as duplicates: {df.loc[duplicates]}, and are removed')
    

    #Fill the outlierdf with the duplicates
    outliers = df[df.index.duplicated(keep=check_settings[checkname]['keep'])]
   
    outlierdf = make_outlier_df_for_check(station_dt_list = outliers.index,
                                          values_in_dict = outliers.to_dict(orient='series'),
                                          flagcolumnname=checks_info[checkname]['label_columnname'],
                                          flag=checks_info[checkname]['outlier_flag'])
    

    #Remove duplicates from the observations
    df = df[~df.index.duplicated(keep=check_settings[checkname]['keep'])]
    
    
    
    return df, outlierdf

# =============================================================================
# Quality assesment checks on dataset
# =============================================================================



def gross_value_check(input_series, obstype):
    """
    Looking for values of an observation type that are not physical. These values are labeled and the physical limits are specified in the qc_settings. 

    Parameters
    ----------
    input_series : pandas.Series
        The observations series of the dataset object
        
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.
        

    Returns
    -------
    updated_obs_series : pandas.Series
        The observations series updated for this check. Observations that didn't pass are removed.
        
    outlier_df : pandas.DataFrame
        The collection of records flagged as outliers by this check.

    """  
    
    checkname = 'gross_value'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        
        return input_series, init_outlier_multiindexdf()
   

    
    
    #find outlier observations as a list of tuples [(name, datetime), (name, datetime)]
    outl_obs = input_series.loc[(input_series <= specific_settings['min_value']) | 
                                (input_series >= specific_settings['max_value'])
                                ].index.to_list()
    
    #make outlierdf
    outlier_df = make_outlier_df_for_check(station_dt_list=outl_obs,
                                           values_in_dict={obstype:input_series.loc[outl_obs]},
                                           flagcolumnname=obstype+'_'+ checks_info[checkname]['label_columnname'],
                                           flag=checks_info[checkname]['outlier_flag'])
    
    
    #drop outliers from input series
    input_series = input_series.drop(outl_obs)
    
    return input_series, outlier_df



def persistance_check(station_frequencies, input_series, obstype):

    """
    V3
    Looking for values of an observation type that do not change during a timewindow. These are flagged as outliers.
    
    In order to perform this check, at least N observations chould be in that time window.
    

    Parameters
    ----------
    input_series : pandas.Series
        The observations series of the dataset object
        
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.
        

        Returns
        -------
        updated_obs_series : pandas.Series
            The observations series updated for this check. Observations that didn't pass are removed.
            
        outlier_df : pandas.DataFrame
            The collection of records flagged as outliers by this check.

    """  
    
    checkname = 'persistance'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        
        return input_series, init_outlier_multiindexdf()
    
    #apply persistance
    def is_unique(window):   #comp order of N (while using the 'unique' function is Nlog(N))
        a = window.values
        a = a[~np.isnan(a)]
        return (a[0] == a).all()
    
    #TODO: Tis is very expensive if no coarsening is applied !!!! Can we speed this up? 
    window_output = input_series.reset_index(level=0).groupby('name').rolling(window= specific_settings['time_window_to_check'],
                                                                            closed='both',
                                                                            center=True,
                                                                            min_periods=specific_settings['min_num_obs']).apply(is_unique)
    
    
    list_of_outliers = []
    outl_obs = window_output.loc[window_output[obstype] == True].index
    for outlier in outl_obs:
        outliers_list = get_outliers_in_daterange(input_series, outlier[1], outlier[0], specific_settings['time_window_to_check'], station_frequencies)
      
        list_of_outliers.extend(outliers_list)
        
    list_of_outliers = list(set(list_of_outliers))
    
    #Create outlier df
    outlier_df = make_outlier_df_for_check(station_dt_list=list_of_outliers,
                                           values_in_dict={obstype:input_series.loc[list_of_outliers]},
                                           flagcolumnname=obstype+'_'+ checks_info[checkname]['label_columnname'],
                                           flag=checks_info[checkname]['outlier_flag'])
    
  
    #drop outliers from input series
    input_series = input_series.drop(list_of_outliers)
    return input_series, outlier_df
      


def repetitions_check(input_series, obstype):
    """
    Looking for values of an observation type that are repeated at least with the frequency specified in the qc_settings. These values are labeled.

    Parameters
    ----------
    input_series : pandas.Series
        The observations series of the dataset object
        
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.
    

    Returns
    -------
    updated_obs_series : pandas.Series
        The observations series updated for this check. Observations that didn't pass are removed.
        
    outlier_df : pandas.DataFrame
        The collection of records flagged as outliers by this check.


    """  
    
    checkname = 'repetitions'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        return input_series, init_outlier_multiindexdf()

    
    #find outlier datetimes
    
    #add time interval between two consecutive records, group by consecutive records without missing records
   
    time_diff = input_series.index.get_level_values('datetime').to_series().diff()
    time_diff.index = input_series.index #back to multiindex
    
    persistance_filter = ((input_series.shift() != input_series)).cumsum()
    

    grouped = input_series.groupby(['name', persistance_filter]) 
    #the above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = grouped.size()
    outlier_groups = group_sizes[group_sizes > specific_settings['max_valid_repetitions']]
    

    
    #add to outl_obs.
    outl_obs = []
    for group_idx in outlier_groups.index:
        groupseries = grouped.get_group(group_idx)
        if len(set(groupseries)) == 1: #Check if all observations are equal in group
            outl_obs.extend(groupseries.index.to_list())
    
    #Create outlier df
    outlier_df = make_outlier_df_for_check(station_dt_list=outl_obs,
                                           values_in_dict={obstype:input_series.loc[outl_obs]},
                                           flagcolumnname=obstype+'_'+ checks_info[checkname]['label_columnname'],
                                           flag=checks_info[checkname]['outlier_flag'])
    
   
    #drop outliers from input series
    input_series = input_series.drop(outl_obs)
    
    return input_series, outlier_df


def step_check(input_series, obstype):
    """

    V3
    Looking for jumps of the values of an observation type that are larger than the limit specified in the qc_settings. These values are removed from 
    the input series and combined in the outlier df.
    
    The purpose of this check is to flag observations with a value that is too much different compared to the previous (not flagged) recorded value.

    Parameters
    ----------
    input_series : pandas.Series
        The observations series of the dataset object
        
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.
        
  

     Returns
     -------
     updated_obs_series : pandas.Series
         The observations series updated for this check. Observations that didn't pass are removed.
         
     outlier_df : pandas.DataFrame
         The collection of records flagged as outliers by this check.

    """ 
    
    checkname='step'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        return input_series, init_outlier_multiindexdf()
    
    list_of_outliers = []
    
    for name in input_series.index.droplevel('datetime').unique():
        subdata = input_series.xs(name, level='name', drop_level=False)
        
        time_diff = subdata.index.get_level_values('datetime').to_series().diff()
        time_diff.index = subdata.index #back to multiindex
        #define filter
        step_filter = (((subdata - subdata.shift(1)) > (specific_settings['max_increase_per_second']*time_diff.dt.total_seconds())) | ((subdata - subdata.shift(1)) < (specific_settings['max_decrease_per_second']*time_diff.dt.total_seconds()))) #& 
                       #(time_diff == station_frequencies[name]))
        outl_obs = step_filter[step_filter==True].index
        
     
        list_of_outliers.extend(outl_obs)
        
    outlier_df = make_outlier_df_for_check(station_dt_list=list_of_outliers,
                                           values_in_dict={obstype:input_series.loc[list_of_outliers]},
                                           flagcolumnname=obstype+'_'+ checks_info[checkname]['label_columnname'],
                                           flag=checks_info[checkname]['outlier_flag'])
    
    input_series = input_series.drop(list_of_outliers)
    
    return input_series, outlier_df

def window_variation_check(station_frequencies, input_series, obstype):   
    """

    V3
    Looking for jumps of the values of an observation type that are larger than the limit specified in the qc_settings. These values are removed from 
    the input series and combined in the outlier df.
    
    There is a increament threshold (that is if there is a max value difference and the maximum value occured later than the minimum value occured.) 
    And vice versa is there a decreament threshold.
    
    The check is only applied if there are at leas N observations in the time window.


    Parameters
    ----------
    station_frequencies: pandas.Series
        The series containing the dataset time-resolution for all stations (as index).
    
    input_series : pandas.Series
        The observations series of the dataset object
        
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.
        
  

     Returns
     -------
     updated_obs_series : pandas.Series
         The observations series updated for this check. Observations that didn't pass are removed.
         
     outlier_df : pandas.DataFrame
         The collection of records flagged as outliers by this check.

    """ 
    checkname='window_variation'
    
    try:
        specific_settings = check_settings[checkname][obstype]
    except:
        print(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        return input_series, init_outlier_multiindexdf()

    # Calculate window thresholds (by linear extarpolation)
    windowsize_seconds = pd.Timedelta(specific_settings['time_window_to_check']).total_seconds()
    max_window_increase = specific_settings['max_increase_per_second'] * windowsize_seconds
    max_window_decrease = specific_settings['max_decrease_per_second'] * windowsize_seconds
    

    #apply steptest
    def variation_test(window):
        if ((max(window) - min(window) > max_window_increase) & 
            (window.idxmax() > window.idxmin())):
            return 1

        if ((max(window) - min(window) > max_window_decrease) & 
            (window.idxmax() < window.idxmin())):
            return 1
        else:
            return 0
        
    window_output = input_series.reset_index(level=0).groupby('name').rolling(window=specific_settings['time_window_to_check'],
                                                                            closed='both',
                                                                            center=True,
                                                                            min_periods=specific_settings['min_window_members']).apply(variation_test)

    list_of_outliers = []
    outl_obs = window_output.loc[window_output[obstype] == 1].index

    for outlier in outl_obs:
        outliers_list = get_outliers_in_daterange(input_series, outlier[1], outlier[0], specific_settings['time_window_to_check'], station_frequencies)
      
        list_of_outliers.extend(outliers_list)
        
    list_of_outliers = list(set(list_of_outliers))
    
    #Create outlier df
    outlier_df = make_outlier_df_for_check(station_dt_list=list_of_outliers,
                                           values_in_dict={obstype:input_series.loc[list_of_outliers]},
                                           flagcolumnname=obstype+'_'+ checks_info[checkname]['label_columnname'],
                                           flag=checks_info[checkname]['outlier_flag'])
   
   
    #drop outliers from input series
    input_series = input_series.drop(list_of_outliers)
    
    return input_series, outlier_df

def get_outliers_in_daterange(input_data, date, name, time_window, station_freq):
    end_date = date + (pd.Timedelta(time_window)/2).floor(station_freq[name])
    start_date = date - (pd.Timedelta(time_window)/2).floor(station_freq[name])
    
    daterange = pd.date_range(start=start_date, end = end_date, freq=station_freq[name])
    
    multi_idx = pd.MultiIndex.from_arrays(arrays=[[name]*len(daterange), daterange.to_list()],
                                          sortorder=1,
                                          names=['name', 'datetime'])
    outlier_sub_df = pd.DataFrame(data=None,
                                         index=multi_idx, 
                                         columns=None)
    
    
    intersection = outlier_sub_df.index.intersection(input_data.dropna().index).values
    
    return intersection




# def qc_info(qc_dataframe):
    
#     number_gross_error_outliers = qc_dataframe.persistance.str.contains('gross value outlier').sum()
#     number_persistance_outliers = qc_dataframe.persistance.str.contains('persistance outlier').sum()
#     number_step_outliers = qc_dataframe.persistance.str.contains('step outlier').sum()
#     number_int_consistency_outliers = qc_dataframe.persistance.str.contains('internal consistency outlier').sum()
    
#     fraction_gross_outliers = number_gross_error_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
#     fraction_persistance_outliers = number_persistance_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
#     fraction_step_outliers = number_step_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
#     fraction_int_consistency_outliers = number_int_consistency_outliers/(number_gross_error_outliers+number_persistance_outliers+number_step_outliers+number_int_consistency_outliers)
    
#     percentage_of_data_gross_outlier = number_gross_error_outliers/len(qc_dataframe)
#     percentage_of_data_persistance_outlier = number_persistance_outliers/len(qc_dataframe)
#     percentage_of_data_step_outlier = number_step_outliers/len(qc_dataframe)
#     percentage_of_data_consistency_outlier = number_int_consistency_outliers/len(qc_dataframe)
    
#     print("Number of gross error outliers: ", number_gross_error_outliers)
#     print("Number of persistance outliers: ", number_persistance_outliers)
#     print("Number of step outliers: ", number_step_outliers)
#     print("Number of internal consistency outliers: ", number_int_consistency_outliers)
    
#     print("Fraction of gross error outliers: ", fraction_gross_outliers)
#     print("Fraction of persistance outliers: ", fraction_persistance_outliers)
#     print("Fraction of step outliers: ", fraction_step_outliers)
#     print("Fraction of internal consistency outliers: ", fraction_int_consistency_outliers)
    
#     print("Percentage of data gross outlier: ", percentage_of_data_gross_outlier)
#     print("Percentage of data persistance outlier: ", percentage_of_data_persistance_outlier)
#     print("Percentage of data step outlier: ", percentage_of_data_step_outlier)
#     print("Percentage of data internal consistency outlier: ", percentage_of_data_consistency_outlier)
    
#     print("Gross error outliers at: ", qc_dataframe[qc_dataframe['gross_value'] == 'gross value outlier'].index)
#     print("Persistance outliers at: ", qc_dataframe[qc_dataframe['persistance'] == 'persistance outlier'].index)
#     print("Step outliers at: ", qc_dataframe[qc_dataframe['step'] == 'step outlier'].index)
#     print("Internal consistency outliers at: ", qc_dataframe[qc_dataframe['internal_consistency'] == 'internal consistency outlier'].index)
#     print("Humidity 0 at: ", qc_dataframe[qc_dataframe['internal_consistency'] == 'humidity 0'].index)
    
    

