#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:01:35 2022

@author: thoverga
"""
import os
from pathlib import Path
main_folder = Path(__file__).resolve().parents[1]
testdata_file = os.path.join(str(main_folder), 'tests', 'test_data',  'vlinderdata.csv' )


import vlinder_toolkit



#%%

# =============================================================================
# Settings
# =============================================================================


# 1. Initiate settings object. This object contains all settings needed for furthur analysis
settings = vlinder_toolkit.Settings()




# 3. If the output data folder and input file are not exported as system variables, you need to update them:
settings.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                         output_data_folder='/home/$USER/output/')





# 4. Check the setting by using the .show() or .check_settings() function 
settings.show()





# =============================================================================
#  Import data
# =============================================================================



#1. Importing a dataset containing mulitple different stations is a function in the Dataset class. First we need to initiate a Dataset with a name of choise.

sept_2022_all_vlinders = vlinder_toolkit.Dataset()


# ---------------- Importing from CSV file -------------------------------


#The dataset is initiated but still empty. Filling it with the data from a csv file is simply done by:
    
sept_2022_all_vlinders.import_data_from_file(settings) #Rember that you added the input file in the settings object, this file will be used.


# =============================================================================
# Applying quality control
# =============================================================================

#  a number of quality control methods are available in the toolkit. We can classify them in two groups:
    # 1) Quality control for missing/duplicated timestamps 
    
        #These are automatically performed when the dataset is created
        
    # 2) Quality control for bad observations
    
        #These are not automatically executed. These checks are performd in a sequence 
        # of specific checks, that are looking for signatures of typically bad observations.
        # The following checks are available:
            # Gross value check: a threshold check that observations should be between the thresholds
            # Persistance check: a check that looks for repetitive observation values (indicating a connection error.)
            
        #All settings, labels, replacement values are defind in the default settings 
        # in /settings_files/qc_settings.py
        
        # Here a copy:
                    
        # check_settings = {
            
        #     #checks on all observation types
        #     "duplicate_timestamp": {}, #No numeric settings
            
            
        #     #checks on specific observation types
        #     "gross_value": {'temp': {'min_value': 8.0,
        #                              'max_value': 18.0},
        #                     },
        #     "persistance": {'temp': {'max_valid_repetitions': 5}}
            
        #     }
        
        
        # outlier_values = {
        #     "duplicate_timestamp": nan, 
        #     "gross_value": nan,
        #     "persistance": nan    
        #     }
        
        # observation_labels={
        #     'ok': 'ok',
        #     'duplicated_timestamp': 'duplicated timestamp outlier',
        #     'gross_value': 'gross value outlier',
        #     'persistance': 'percistance outlier'
        #     }
        
        
# performing a quality control check on a station

station = sept_2022_all_vlinders.get_station('vlinder05')

station.apply_gross_value_check(obstype='temp')
station.apply_persistance_check(obstype='temp')
#Now the temperature observations are checked by the gross_value_check and persistance check.
 # You can maybe see this in station.temp or by plotting (gaps in the observations)


# You can also see an overvieuw of all performed checks on a station

qc_overview = station.qc_labels_df

#This is a dictionary where each observation type is a key and an overview dataframe as values

print(qc_overview['temp'])



#To apply quality control checks on the whole dataset you can do this:
    
sept_2022_all_vlinders.apply_quality_control(obstype='temp',
                                            gross_value=True, #apply this check 
                                            persistance=True, #apply this check
                                            )




# =============================================================================
# Timeseries plots
# =============================================================================

#NOTE: All styling settings are defined in the src/vlinder_toolkit/settings_files/default_formats_settings.py


#To plot timeseries for one station you can use the make_plot function on the station object:
favorite_station = sept_2022_all_vlinders.get_station(stationname='vlinder02')


#Possible obstypes to plot:
    # 'temp','radiation_temp','humidity','precip','precip_sum','wind_speed',
    # 'wind_gust','wind_direction','pressure','pressure_at_sea_level'
    
favorite_station.make_plot(variable='temp', 
                           title=None) #if title=None, a title will be generated


#Nan values will be not shown in the plot




#If you want to compair multiple stations by a timeseries you can use the make_plot function on the dataset:

from datetime import datetime    

sept_2022_all_vlinders.make_plot(stationnames=['vlinder02', 'vlinder05', 'vlinder07'],
                                variable='humidity',
                                starttime=datetime(2022, 9,4), # 2022/09/04 00:00:00
                                endtime=datetime(2022,9,12,12,45), #2022/09/12 12:45:00
                                title=None,
                                legend=True
                                )
#If you do not specify a start and endtime, all available timestamps are used.





# =============================================================================
# Geospatial plots
# =============================================================================

# geospatial plots can be made for a given moment (datetimeinstance) by 
# applying the make_geo_plot() on a dataset object:
    
sept_2022_all_vlinders.make_geo_plot(variable='temp',
                                    timeinstance=datetime(2022, 9,4), # 2020/09/04 00:00:00
                                    title=None,
                                    legend=True,
                                    vmin=None, #value corresponding to the minimum of the colorscheme
                                    vmax=None) #value corresponding to the maximum of the colorscheme

# Notes:
    # * If no timeinstance is given, the first timestamp of the network is used
    # * vmin and vmax are default (None) set to the minimum and maximum of the observations. 
    #   You can change this if you whant to compair geoplots using the same color scheme
    












