#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""
import pandas as pd
import numpy as np

from .stations import Station, Dataset
from .settings import Settings
from .settings_files.qc_settings import check_settings, outlier_values, observation_labels




def duplicate_timestamp(station, obstype='temp'):
    """
    This check looks for duplicate timestamps. If a duplicate timestamp is found, it will be handled as specified in the outlier_values.
    Either is will be removed or replaced by a value.

    Parameters
    ----------
    station : Station object
        The Station on which to apply this check.
    obstype : string, optional
        The observationtype on which to apply this check. The default is 'temp'.
    
    Returns
    -------
    station : Station object
        The observations of the station object are filtered by this check. The duplicate timestamps are added in the qc_info attribute of the station.

    """
    
    try:
        specific_settings = check_settings['duplicate_timestamp'][obstype]
    except:
       print('No duplicate_timestamp settings found for obstype=', obstype, '. Check is skipped!') 
       return station
    
    
    #extract observations
    obs = getattr(station, obstype)
    #sort observations
    obs = obs.sort_index()
    
   #Locate duplicates
    dub_mask = obs.index.duplicated(keep='first')
    dub_timestamps = obs.index[dub_mask]
    
    
    #update qc_info
    station.qc_info['duplicate_timestamp'] = {'Original duplicate timestamps': dub_timestamps.to_list()}
    
    
    
    
    #Update observations
    if outlier_values['duplicate_timestamp'] == 'drop':
        obs = obs.loc[~dub_mask] #drop duplicated timestamps
    else:
        obs.loc[dub_mask] = outlier_values['duplicate_timestamp']
    
    #Update labels
    
    
    #Update station object
    setattr(station, obstype, obs)
    return station
    
    


def gross_value_check(station, obstype='temp', ignore_val=np.nan):
   try:
       specific_settings = check_settings['gross_value'][obstype]
   except:
      print('No gross_value settings found for obstype=', obstype, '. Check is skipped!') 
      # return station

   #extract observations into dataframe
   obs = getattr(station, obstype)
   obs = pd.DataFrame({'obs': obs})



   #get observations to ignore
   ignore_obs = obs[obs['obs'] == ignore_val]
   to_check_obs = obs[obs['obs'] != ignore_val]


                            
   #create labels
   to_check_obs['labels'] = [observation_labels['ok'] if (specific_settings['min_value'] <= obs_value <=specific_settings['max_value'])
                                 else observation_labels['gross_value'] for obs_value in to_check_obs['obs']]

   ignore_obs['labels'] = 'not_checked'



   #Merge together
   checked_obs = to_check_obs.append(ignore_obs).sort_index()


   #Update station qc labels
   station.qc_labels_df[obstype]['gross_value'] = checked_obs['labels']


   #update observations
   def new_obs(row):
       if row['labels'] != observation_labels['gross_value']:
           return row['obs']
       else:
           return outlier_values['gross_value']

   checked_obs['new_obs'] = checked_obs.apply(lambda row: new_obs(row), axis=1)
   setattr(station, obstype, checked_obs['new_obs'])

   return station