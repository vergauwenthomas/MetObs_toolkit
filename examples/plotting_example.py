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


# import sys
# sys.path.append(str(main_folder))


import vlinder_toolkit




# =============================================================================
#  Import data
# =============================================================================



#1. Importing a dataset containing mulitple different stations is a function in the Dataset class. First we need to initiate a Dataset with a name of choise.

aug_2020_all_vlinders = vlinder_toolkit.Dataset()
aug_2020_all_vlinders.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                         output_folder='/home/$USER/output/',
                         input_metadata_file=vlinders_metadatafile)


# ---------------- Importing from CSV file -------------------------------


#The dataset is initiated but still empty. Filling it with the data from a csv file is simply done by:

aug_2020_all_vlinders.import_data_from_file() #Rember that you added the input file in the settings object, this file will be used.



# =============================================================================
# Timeseries plots
# =============================================================================

#NOTE: All styling settings are defined in the src/vlinder_toolkit/settings_files/default_formats_settings.py

# plotting the timeseries of all stations can simply be done by
from datetime import datetime

aug_2020_all_vlinders.make_plot(
                                #specify the names of the stations in a list, or use None to plot them all.
                                stationnames = None,
                                #what variable to plot (default is 'temp')
                                variable='temp',
                                #choose how to color the timeseries:
                                    #'name' : a specific color per station
                                    #'label': a specific color per quality control label
                                colorby='name',
                                #choose a start en endtime for the series (datetime).
                                #Default is None, wich uses all available data
                                starttime=None,
                                endtime = datetime(2022, 9, 5),
                                #specify title, by default a title is created
                                title=None,
                                # Add legend to plot?, by default true
                                legend=True,
                                # Plot observations that are labeld as outliers?
                                show_outliers=True
                                )


#To plot timeseries for one station you can use the make_plot function on the station object:
favorite_station = aug_2020_all_vlinders.get_station(stationname='vlinder02')
favorite_station.apply_quality_control()

#Possible obstypes to plot:
    # 'temp','radiation_temp','humidity','precip','precip_sum','wind_speed',
    # 'wind_gust','wind_direction','pressure','pressure_at_sea_level'


favorite_station.make_plot(variable='temp',
                           colorby='label',
                           title=None) #if title=None, a title will be generated

#Nan values will be not shown in the plot

#If you want to compair multiple stations by a timeseries you can use the make_plot function on the dataset:

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













