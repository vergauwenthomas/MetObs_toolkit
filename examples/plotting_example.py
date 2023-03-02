#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 12:01:35 2022

@author: thoverga
"""
import os
from pathlib import Path
main_folder = Path(__file__).resolve().parents[1]
testdata_file = os.path.join(str(main_folder), 'tests', 'test_data',  'vlinderdata_small.csv' )
vlinders_metadatafile = os.path.join(str(main_folder), 'static_data', 'vlinder_metadata.csv' )


#import sys
#sys.path.append(str(main_folder))


import vlinder_toolkit



#%%

# =============================================================================
# Settings
# =============================================================================


# 1. Initiate settings object. This object contains all settings needed for furthur analysis
settings = vlinder_toolkit.Settings()




# 3. If the output data folder and input file are not exported as system variables, you need to update them:
settings.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                         output_folder='/home/$USER/output/',
                         input_metadata_file=vlinders_metadatafile)





# 4. Check the setting by using the .show() or .check_settings() function 
settings.show()





# =============================================================================
#  Import data
# =============================================================================



#1. Importing a dataset containing mulitple different stations is a function in the Dataset class. First we need to initiate a Dataset with a name of choise.

aug_2020_all_vlinders = vlinder_toolkit.Dataset()


# ---------------- Importing from CSV file -------------------------------


#The dataset is initiated but still empty. Filling it with the data from a csv file is simply done by:
    
aug_2020_all_vlinders.import_data_from_file(settings) #Rember that you added the input file in the settings object, this file will be used.




# =============================================================================
# Timeseries plots
# =============================================================================

#NOTE: All styling settings are defined in the src/vlinder_toolkit/settings_files/default_formats_settings.py


#To plot timeseries for one station you can use the make_plot function on the station object:
favorite_station = aug_2020_all_vlinders.get_station(stationname='vlinder02')


#Possible obstypes to plot:
    # 'temp','radiation_temp','humidity','precip','precip_sum','wind_speed',
    # 'wind_gust','wind_direction','pressure','pressure_at_sea_level'
    
favorite_station.make_plot(variable='temp', 
                           title=None) #if title=None, a title will be generated


#Nan values will be not shown in the plot




#If you want to compair multiple stations by a timeseries you can use the make_plot function on the dataset:

from datetime import datetime    

aug_2020_all_vlinders.make_plot(stationnames=['vlinder02', 'vlinder05', 'vlinder07'],
                                variable='humidity',
                                title=None,
                                legend=True
                                )
#If you do not specify a start and endtime, all available timestamps are used.





# =============================================================================
# Geospatial plots
# =============================================================================

# geospatial plots can be made for a given moment (datetimeinstance) by 
# applying the make_geo_plot() on a dataset object:
    
aug_2020_all_vlinders.make_geo_plot(variable='temp',
                                    timeinstance=datetime(2022, 9,6), # 2022/09/06 00:00:00
                                    title=None,
                                    legend=True,
                                    vmin=None, #value corresponding to the minimum of the colorscheme
                                    vmax=None) #value corresponding to the maximum of the colorscheme

# Notes:
    # * If no timeinstance is given, the first timestamp of the network is used
    # * vmin and vmax are default (None) set to the minimum and maximum of the observations. 
    #   You can change this if you whant to compair geoplots using the same color scheme
    












