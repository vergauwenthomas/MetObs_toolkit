#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022

@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
import numpy as np
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

# % Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

# lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')


#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                          # geotiff_lcz_file=lcz_map
                          output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                          )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)



dataset.apply_quality_control(persistance=False, repetitions=False, step=True)


# dataset.write_to_csv(filename='remove_me', add_final_labels=True)




#%%

check_settings = {
    
    "gaps_finder": {'gapsize_n': 40}, #gaps defined as n times the highest frequency on IO. 
    
    #checks on all observation types
    "duplicated_timestamp": {'keep': False}, #No numeric settings (False: drop all duplicates)

    "missing_timestamp": {},
    
    "persistance": {'temp': {'time_window_of_assumed_change': 5400,#in seconds
                             'minimum_numer': 5}}, #Minimum numer of records in window to perform check
    
    "repetitions": {'temp': {'max_valid_repetitions': 5}},
    
    #checks on specific observation types
    "gross_value": {'temp': {'min_value': -15.0,
                             'max_value': 39.0},
                    },
   
    
    "step": {'temp': {'max_increase_per_second': 8.0/3600.0, #== max 8° change in one hour
                      'max_decrease_per_second': 10.0/3600.0}
             },
    
    "internal_consistency": {'temp': {'b': 18.678,
                             'c': 257.14, 'd': 234.5}} 
    }


checks_info={
    # 1. --> on data import
    'duplicated_timestamp':{'label_columnname': 'duplicated_timestamp_label',
                            'outlier_flag': 'duplicated timestamp outlier',
                            'numeric_flag': 1,
                            'apply_on': 'record'
                            },
    # 2(A). --> on data import
    'gaps_finder':{'label_columnname': 'gap_timestamp_label',
                            'outlier_flag': 'missing timestamp (gap)',
                            'numeric_flag': 2,
                            'apply_on': 'record'
                            },
    # 2(B). --> on data import
    'missing_timestamp':{'label_columnname': 'missing_timestamp_label',
                            'outlier_flag': 'missing timestamp',
                            'numeric_flag': 3,
                            'apply_on': 'record'
                            },
    # 3. --> on observed values
    'gross_value':{'label_columnname': 'gross_value_label', #Obstype_ is prefix
                            'outlier_flag': 'gross value outlier',
                            'numeric_flag': 4,
                            'apply_on': 'obstype'
                            },
    # 4(A). --> on observed values
    'persistance':{'label_columnname': 'persistance_label', #Obstype_ is prefix
                            'outlier_flag': 'persistance outlier',
                            'numeric_flag': 5,
                            'apply_on': 'obstype'
                            },
    # 4(B). --> on observed values
    'repetitions':{'label_columnname': 'repetitions_label', #Obstype_ is prefix
                            'outlier_flag': 'repetitions outlier',
                            'numeric_flag': 6,
                            'apply_on': 'obstype'
                            },
    # 5. --> on observed values
    'step':{'label_columnname': 'step_label', #Obstype_ is prefix
                            'outlier_flag': 'step outlier',
                            'numeric_flag': 7,
                            'apply_on': 'obstype'
                            },
    
    }

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



#%%


def step_check(input_series, dataset_resolution, obstype='temp'):   
    """

    V3
    Looking for jumps of the values of an observation type that are larger than the limit specified in the qc_settings. These values are removed from 
    the input series and combined in the outlier df.


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
        # logger.warning(f'No {checkname} settings found for obstype={obstype}. Check is skipped!')
        return input_series, init_outlier_multiindexdf()

    
    
    # Define window settings
    highest_res = min(dataset_resolution)
    windowsize = specific_settings['max_window_members'] * highest_res
    
    max_window_increase = specific_settings['max_increase_per_second'] * windowsize.total_seconds()
    max_window_decrease = specific_settings['max_decrease_per_second'] * windowsize.total_seconds()

    def steptest(window):
        if ((max(window) - min(window) > max_window_increase) & 
            (window.idxmax() > window.idxmin())):
            return 1

        if ((max(window) - min(window) > max_window_decrease) & 
            (window.idxmax() < window.idxmin())):
            return 1
        else:
            return 0

    
    #apply steptest
    step_output = input_series.reset_index(level=0).groupby('name').rolling(windowsize,
                                                                            closed='both',
                                                                            center=True,
                                                                            min_periods=specific_settings['min_window_members']).apply(steptest)

    outl_obs = step_output.loc[step_output == 1].index

    
    #Create outlier df
    outlier_df = make_outlier_df_for_check(station_dt_list=outl_obs,
                                           values_in_dict={obstype:input_series.loc[outl_obs]},
                                           flagcolumnname=obstype+'_'+ checks_info[checkname]['label_columnname'],
                                           flag=checks_info[checkname]['outlier_flag'])
    
    #drop outliers from input series
    input_series = input_series.drop(outl_obs)
    return input_series, outlier_df




#%%

resolution = dataset.metadf['dataset_resolution']

highest_res = min(resolution)


min_obs_per_window = 2
input_series = dataset.df['temp']

input_series.iloc[12] = 22.3


min_window_members = 3
max_window_members = 5


windowsize = max_window_members * highest_res
max_increase_per_second = 8.0/3600.0


max_increase_per_window = max_increase_per_second * windowsize.total_seconds()

decrease_dif_window = 3.
increase_dif_window = 2.





def steptest(window):
    if ((max(window) - min(window) > increase_dif_window) & 
        (window.idxmax() > window.idxmin())):
            return 1

    if ((max(window) - min(window) > decrease_dif_window) & (window.idxmax() < window.idxmin())):
        return 2
    else:
        return 0
# df = input_series.to_frame().reset_index()

# test = input_series.groupby(level='name').droplevel(0).rolling('1H', min_periods=min_obs_per_window).apply(steptest)
# test = input_series.reset_index(level=0).groupby('name').rolling('1H', center=True).apply(steptest)
test = input_series.reset_index(level=0).groupby('name').rolling(windowsize,closed='both', center=True).apply(steptest)

testa = input_series.to_frame()
testa['label'] = test


#%%
testa = input_series.groupby(level='name')

# test = testa.get_group('vlinder01').droplevel(0)
test = input_series.reset_index(level=0).groupby('name').rolling('H', center=True).apply(steptest)