#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

import metobs_toolkit


#%%
import pandas as pd
import datetime

#%%

# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template # Contains also the metadata mapping
                        )

# Load the data from the demo data files
dataset.import_data_from_file()

#%%

# The outliers are stored in the outliersdf object of the dataset:
outliers = dataset.outliersdf

# Print this object to see what is stored in this data frame:
print(outliers)

#%%

# All settings, labels, replacement values are defined in the default settings in /settings_files/qc_settings.py
# To inspect (and change) these quality control settings, you can extract them:
qc_settings = dataset.settings.qc["qc_check_settings"]

# These settings are in a dictionary which contains multiple levels.
# The first level concerns the specific quality control check which the settings are for.
# You can print the keys of the dictionary to get an idea of the different available checks:
print(qc_settings.keys())


#%%
dataset.apply_quality_control(
    obstype="temp",         # choose which observations you want to check
    gross_value=True,       # set True if you want to perform the gross value check
    persistance=True,       # set True if you want to perform the persistence check
    step=True,              # set True if you want to perform the spike check
    window_variation=True,  # set True if you want to perform the window variation check
)


#%%

# Print the outliers dataframe. Are there more outliers than before?
dataset.outliersdf.xs('temp', level='obstype') #select only the temperature outliers
#%%

qc_statistics = dataset.get_qc_stats(
    obstype="temp",     # Specify which observation variable you want to get the statistics for; here we choose temperature
    stationname=None,  # None means all stations are plotted. You can also plot a specific station by saying 'station_A'
    make_plot=True,     # Set True to make a plot
)


#%%

specific_station = 'vlinder01' #the name of the station

station = dataset.get_station(specific_station)

station.apply_quality_control(
    obstype="temp",         # choose which observations you want to check
    gross_value=True,       # set True if you want to perform the gross value check
    persistance=True,       # set True if you want to perform the persistence check
    step=True,              # set True if you want to perform the spike check
    window_variation=True,  # set True if you want to perform the window variation check
)

qc_statistics = station.get_qc_stats(
    obstype="temp",     # Specify which observation variable you want to get the statistics for; here we choose temperature
    make_plot=True,     # Set True to make a plot
)


#%%

dataset.make_plot(colorby = 'label')


#%%


# dataset = metobs_toolkit.Dataset()

# # Add the demo data files to the dataset settings
# dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
#                         input_metadata_file = metobs_toolkit.demo_metadatafile,
#                         data_template_file = metobs_toolkit.demo_template,
#                         metadata_template_file = metobs_toolkit.demo_template # Contains also the metadata mapping
#                         )

# # Load the data from the demo data files
# dataset.import_data_from_file()

#update the settings
dataset.update_qc_settings(obstype='temp',
                           gross_value_max_value=27.2,
                           win_var_time_win_to_check='3H', #3 hour
                           step_max_decrease_per_sec=7.2/3600)



print(dataset._applied_qc)

#coarsen timeresolution

dataset.coarsen_time_resolution(freq='1H')

#apply QC
dataset.apply_quality_control(
    obstype="temp",
)

#Visualise the effect
dataset.make_plot(obstype='temp', colorby='label')




#%%
# TESTING






# #update the settings
# dataset.update_qc_settings(obstype='temp',
#                            gross_value_max_value=27.2,
#                            win_var_time_win_to_check='3H', #3 hour
#                            step_max_decrease_per_sec=7.2/3600)
# #coarsen timeresolution
# dataset.coarsen_time_resolution(freq='1H')

# #apply QC
# dataset.apply_quality_control(
#     obstype="temp",
# )

# #Visualise the effect
# dataset.make_plot(obstype='temp', colorby='label')



