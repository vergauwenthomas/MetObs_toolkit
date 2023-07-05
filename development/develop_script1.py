#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit
#
import os
import sys
from pathlib import Path
import pandas as pd
import time
import math


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%
import pandas as pd


# # Specify the path to your Drive folder
# # BASE_DIR = '/content/drive/MyDrive/FAIRNESS_summerschool_23/'
# # print('BASE_DIR: ',BASE_DIR)


# kwargs_data = {'sep':';'}
# kwargs_metadata = {'sep':';'}



# # Construct the path, by adding the name of the folder and file
data_path = '/home/thoverga/Desktop/testdata/fNCEI_BCG_2010.csv'
metadata = "/home/thoverga/Desktop/testdata/Stations_coordinates.csv"



df = pd.read_csv(data_path, sep=',')

df['TMP'] = df['TMP'].astype(str)

df['TMP'] = [dat if len(dat) == 8 else '0'+dat for dat in df['TMP'] ]

newdata = '/home/thoverga/Desktop/testdata/fNCEI_BCG_2010_cleanup.csv'
df.to_csv(newdata)


#%%
dataf = '/home/thoverga/Desktop/testdata/fNCEI_BCG_2010_cleanup.csv'
metaf ="/home/thoverga/Desktop/testdata/Stations_coordinates.csv"

temp ="/home/thoverga/Desktop/testdata/default_template.csv"


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
                        input_data_file= dataf,
                        input_metadata_file = metaf,
                        data_template_file = temp,
                        metadata_template_file = temp
                        )

dataset.import_data_from_file()

dataset.sync_observations(tollerance='20T', _drop_target_nan_dt=True, _force_resolution_minutes=30)


after_sync = dataset.combine_all_to_obsspace()
after_sync = after_sync.xs('temp', level='obstype')
pure_input = dataset.input_df
#%%


from datetime import timedelta

test = timedelta(minutes=30.)

stations = dataset.metadf.index


freq_series = pd.Series(index=dataset.metadf.index,
                        data = [test]*len(stations))

# idee is om de freqencie series mee te geven aan de syncronize functie, om
# te forcen naar een specifieke resolutie







#%%
# READ IN INCOMPLETE DATA SET

# Make an empty data set
demo_incomplete = metobs_toolkit.Dataset()

# Update settings with path to data files
demo_incomplete.update_settings(
                        input_data_file= data_path,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template
                        )

# Update the gap definition if required
demo_incomplete.update_qc_settings(
                                  gapsize_in_records = 12,
                                 )

# Import the data set
demo_incomplete.import_data_from_file()
#%%

import datetime as dt
# # Print the first five lines
# print("First five lines of data set:")
# print(demo_incomplete.df.head(5))

# # Print the column 'assumed_import_frequency' of the metadf dataframe, to look at the assumed time resolution
# print("Metadf dataframe:")
# print(demo_incomplete.metadf.assumed_import_frequency)

# # Print the outliersdf
# print(demo_incomplete.outliersdf)

# # Update the missing observations and gaps
# demo_incomplete.update_gaps_and_missing_from_outliers(n_gapsize=12)

# # Get general information about all the gaps of the data set
# print("Information about all gaps in the data set:")
# print(demo_incomplete.get_gaps_df())
# print(demo_incomplete.get_gaps_info())

# # Get detailed information about one of the gaps of the data set (! Start counting from 0 to get the first gap !)
# print("Information on the first gap:")
# print(demo_incomplete.gaps[0].to_df())
# print(demo_incomplete.gaps[0].get_info())

# # Get info about missing observations
# print("Information about the missing observations")
# print(demo_incomplete.get_missing_obs_info())

# # Get detailed information of the gaps of one station
# your_chosen_station = demo_incomplete.get_station('vlinder01')
# print(your_chosen_station.get_gaps_df())

# # Make timeplot of data set

# # Select a station and time period
# stationname = "vlinder28"
# begin_plot = dt.datetime.strptime("2022-09-04 20:00:00", '%Y-%m-%d %H:%M:%S')
# end_plot = dt.datetime.strptime("2022-09-07 15:00:00", '%Y-%m-%d %H:%M:%S')

# # Make timeplot
# demo_incomplete.make_plot(
#                           stationnames=[stationname],
#                           obstype='temp',
#                           colorby='name',
#                           starttime=begin_plot,
#                           endtime=end_plot,
#                           legend=False,
#                           show_outliers=False
#                           )
# Change the time resolution
demo_incomplete.coarsen_time_resolution(freq = "30T")

# Select one station
Station = demo_incomplete.get_station('vlinder04')
# #%%
# # Print the first 5 lines of the result
# print(Station.df.head(5))

# # Fill the missing observations through linear interpolation
# # (you don't have to run this code cel, since vlinder04 does not have any missing obsevations)
# Station.fill_missing_obs_linear()

# # Perform linear interpolation
# Station.fill_gaps_linear(obstype='temp')

# # Select a station and time period
# stationname = "vlinder04"
# begin_plot = dt.datetime.strptime("2022-09-03 00:00:00", '%Y-%m-%d %H:%M:%S')
# end_plot = dt.datetime.strptime("2022-09-03 23:00:00", '%Y-%m-%d %H:%M:%S')

# # Make timeplot
# Station.make_plot(
#                   stationnames=[stationname],
#                   obstype='temp',
#                   colorby='label',
#                   starttime=begin_plot,
#                   endtime=end_plot,
#                   legend=False,
#                   show_outliers=False
#                   )
#%%
# Get the ERA5 data
era_model = demo_incomplete.get_modeldata(
                    modelname='ERA5_hourly',  # Name of the model data to get from GEE
                    stations='vlinder04',     # If only for one station: specify its name, if None the data is downloaded for all stations
                    startdt=None,             # None means the starttime of the observations is taken
                    enddt=None)               # None means the endtime of the observations is taken

# Print the first five lines
# print(era_model.df.head())
#%%

Station.fill_gaps_era5(modeldata=era_model,   # The name of the Modeldata object
                       method='debias',       # Before filling with ERA5, debias first
                       obstype='temp',        # Perform gap-filling on temperature
                       overwrite_fill=True)   # Overwrite the previous results

# Select a station and time period
stationname = "vlinder04"
begin_plot = dt.datetime.strptime("2022-09-03 00:00:00", '%Y-%m-%d %H:%M:%S')
end_plot = dt.datetime.strptime("2022-09-03 23:00:00", '%Y-%m-%d %H:%M:%S')

# Make timeplot
Station.make_plot(
                  stationnames=[stationname],
                  obstype='temp',
                  colorby='label',
                  starttime=begin_plot,
                  endtime=end_plot,
                  legend=False,
                  show_outliers=False
                  )


#%%

# Fill gaps with the hybrid method
Station.fill_gaps_automatic(modeldata=era_model,
                            obstype="temp",
                            max_interpolate_duration_str='6H',
                            overwrite_fill=True)

# Select a station and time period
stationname = "vlinder04"
begin_plot = dt.datetime.strptime("2022-09-03 00:00:00", '%Y-%m-%d %H:%M:%S')
end_plot = dt.datetime.strptime("2022-09-03 23:00:00", '%Y-%m-%d %H:%M:%S')

# Make timeplot
Station.make_plot(
                  stationnames=[stationname],
                  obstype='temp',
                  colorby='label',
                  starttime=begin_plot,
                  endtime=end_plot,
                  legend=False,
                  show_outliers=False
                  )

#%%
Station.fill_gaps_automatic(modeldata=era_model,
                            obstype="temp",
                            max_interpolate_duration_str='6H', # Change this for your choice of parameter!
                            overwrite_fill=True)

# Select a station and time period
stationname = "vlinder04"
begin_plot = dt.datetime.strptime("2022-09-03 00:00:00", '%Y-%m-%d %H:%M:%S') # Change this for your chosen gap!
end_plot = dt.datetime.strptime("2022-09-03 23:00:00", '%Y-%m-%d %H:%M:%S') # Change this for you chosen gap!

# Make timeplot
Station.make_plot(
                  stationnames=[stationname],
                  obstype='temp',
                  colorby='label',
                  starttime=begin_plot,
                  endtime=end_plot,
                  legend=False,
                  show_outliers=False
                  )


Station.fill_gaps_automatic(modeldata=era_model,
                            obstype="temp",
                            max_interpolate_duration_str='6H', # Change this for your choice of parameter!
                            overwrite_fill=True)

# Select a station and time period
stationname = "vlinder04"
begin_plot = dt.datetime.strptime("2022-09-03 00:00:00", '%Y-%m-%d %H:%M:%S') # Change this for your chosen gap!
end_plot = dt.datetime.strptime("2022-09-03 23:00:00", '%Y-%m-%d %H:%M:%S') # Change this for you chosen gap!

# Make timeplot
Station.make_plot(
                  stationnames=[stationname],
                  obstype='temp',
                  colorby='label',
                  starttime=begin_plot,
                  endtime=end_plot,
                  legend=False,
                  show_outliers=False
                  )
