#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022

@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

#% Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')


#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                          input_metadata_file=static_data,
                          geotiff_lcz_file=lcz_map
                         )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)




dataset.apply_quality_control(obstype='temp')


sta1 = dataset.get_station('vlinder05')
sta1.make_plot()




#%%



sta1.make_plot()

dataset.update_dataset_df_with_stations()

test=dataset.df


#%%
dataset.make_plot(['vlinder02', 'vlinder05'])

#%%



#%% 
import numpy as np


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


dat = [3,5,66,4,2,np.nan, 17, np.nan]
# dat = [3,5,66,4,2,123, 17, 77]
ind = [1,2,3,4,5,6,7,8]

data = pd.Series(data=dat, index=ind)
data.name= 'temp'
check_series, ignore_series = split_to_check_to_ignore(data)


outl_datetimes = check_series.loc[(check_series <= 2) | (check_series >= 50)].index.to_list()


df = make_checked_obs_and_labels_series(checked_obs_series=check_series,
                      ignored_obs_series = ignore_series, 
                      outlier_dt_list=outl_datetimes,
                      checkname='checkname',
                      outlier_label='fout',
                      outlier_value = -999,
                      obstype='temp',
                      not_checked_label='not checked',
                      ok_label='ok')

print(df)

#%%
import numpy as np



check_settings = {
    
    #checks on all observation types
    "duplicate_timestamp": {}, #No numeric settings
    
    
    #checks on specific observation types
    "gross_value": {'temp': {'min_value': 8.0,
                             'max_value': 18.0},
                    },
    "persistance": {'temp': {'max_valid_repetitions': 5}}
    
    }


outlier_values = {
    "duplicate_timestamp": np.nan, 
    "gross_value": np.nan,
    "persistance": np.nan    
    }



observation_labels={
    'ok': 'ok',
    'duplicated_timestamp': 'duplicated timestamp outlier',
    'gross_value': 'gross value outlier',
    'persistance': 'percistance outlier'
    }




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



def persistance(input_series, obstype='temp', ignore_val=np.nan):
   
    
    try:
         specific_settings = check_settings['persistance'][obstype]
    except:
        print('No persistance settings found for obstype=', obstype, '. Check is skipped!') 
         # return station
        qc_flags = pd.Series('not checked', index=input_series.index)
        qc_flags.name = 'gross_value'
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



test = persistance(input_series=sta1.temp,
                   obstype='temp',
                   ignore_val=np.nan)


#%%

dataset.make_geo_plot(variable='temp')

#%%

dataset.make_geo_plot(variable='temp', vmin=10, vmax=13)

#%%
# =============================================================================
# checks
# =============================================================================

sta = dataset.get_station('vlinder02')

df_init = sta.df()
sta.make_plot(title='init temp')


sta = vlinder_toolkit.qc_checks.duplicate_timestamp(sta)
sta.make_plot(title='after timstamp dub qc')
sta = vlinder_toolkit.qc_checks.gross_value_check(sta)
sta.make_plot(title='after gross value qc')
sta = vlinder_toolkit.qc_checks.persistance(sta)
sta.make_plot(title='after persistance qc')

# df = sta.df()
# sta.make_plot()


# sta.make_plot()



