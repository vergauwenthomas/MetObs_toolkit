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
metadata = os.path.join(str(main_folder), 'static_data', 'vlinder_metadata.csv')
templatefile = os.path.join(str(main_folder), 'vlinder_toolkit', 'data_templates',
                            'template_defaults', 'default_template.csv')

import sys
sys.path.append(str(main_folder))


import vlinder_toolkit


# =============================================================================
# Dataset
# =============================================================================

# A dataset is a collection of all observational data. Most of the methods are
# applied directly on a dataset. Always start by initializing a dataset object:

aug_2020_all_vlinders = vlinder_toolkit.Dataset()


# we created an dataset and stored in under the variable 'aug_2020_all_vlinder'.
# The show function prints out an overview of data in the dataset:
aug_2020_all_vlinders.show()

# =============================================================================
# Importing data
# =============================================================================

# To import data, you must specify where your data (and additional templates) are.
# All this information is stored in the Settings of a dataset. To get an overview
# of all the settings:
aug_2020_all_vlinders.show_settings()

# you can see that the settings contain default values. But there are no default
# values for where your data is. So you have to update the settings.


# Start by updating the settings of your dataset so it knows here the data files are:


aug_2020_all_vlinders.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                         input_metadata_file=metadata,
                         output_folder='/home/$USER/output/',
                         data_template_file=templatefile,
                         metadata_template_file=templatefile)



# ---------------- Importing from CSV file -------------------------------


#The dataset is initiated but still empty. Filling it with the data from a csv file is simply done by:

aug_2020_all_vlinders.import_data_from_file() #Rember that you added the input file in the settings object, this file will be used.

# Now the dataset is filled with your data.
aug_2020_all_vlinders.show()





# ---------------- Importing from Database -------------------------------

#First of all, make shure that your VNP connection (UGent) is on, or you ar working from within the Ugent !!!
# 4. Updating DB settings
# To extract data directly from the database you need an active UGent VNP connection (or from inside the UGent)
# To connect to the database you need a specific login and password. (Contact Thomas (thomas.vergauwen@meteo.be) for this account.)

#The user and password are extracted from the environment-variables: VLINDER_DB_USER_NAME and VLINDER_DB_USER_PASW

#If this is done correctly than the user and password should appear when running settings.show()

# We have to start by first making an (empty) dataset
# sept_2022_all_vlinders = vlinder_toolkit.Dataset()


#We need to specify the start-moment and end-moment for this period. To do this we need to import the datetime module (base python module)
# from datetime import datetime

#To get the data it is a simple as:
# sept_2022_all_vlinders.import_data_from_database(start_datetime=datetime(2022, 9,1), # 2022/09/01 00:00:00
                                                #  end_datetime=datetime(2022,9,2,12,45)) #2022/09/25 12:45:00




# =============================================================================
# Analysing one station
# =============================================================================

#if you are interested in one station, you can extract all the info for that one station from the dataset by:

favorite_station = aug_2020_all_vlinders.get_station(stationname='vlinder02')

#the variable favorite_station now contains all info about that station. If you would like to have the observation in a tabular form (pandas dataframe):

observations = favorite_station.df

print(observations.head())


#Or you can make a timeseries plot for a field of chose:

favorite_station.make_plot(variable='temp')

#You can extract attributes from the station:

temperature = favorite_station.df['temp']




# =============================================================================
# Writing to a file
# =============================================================================

# to write the dataset to a file, specify an outputfolder first. The data will be written there.
aug_2020_all_vlinders.update_settings(output_folder=os.getcwd())


# To write the output to a file the following function can be used:
aug_2020_all_vlinders.write_to_csv(filename='testdata.csv',
                                   include_outliers=True,
                                   add_final_labels=True,
                                   use_tlk_obsnames=True,
                                   )