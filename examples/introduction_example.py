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

import vlinder_toolkit



#%%

# =============================================================================
# Settings
# =============================================================================





# Settings contains all information about all the paths that will be used. Without proper initiating the settings, a script will not run.

# 1. Initiate settings object. This object contains all settings needed for furthur analysis
settings = vlinder_toolkit.Settings()

#At any time you can see all the used settings:
settings.show()

# 2. If the output data folder and input file are not exported as system variables, you need to update them:
settings.update_settings(input_data_file=testdata_file, #A demo data file, downloaded with brian tool: https://vlinder.ugent.be/vlinderdata/multiple_vlinders.php
                         input_metadata_file=metadata,
                         output_folder='/home/$USER/output/')



# 3. Check the setting again by using the .show() or .check_settings() function 
settings.show()


# 4. Updating DB settings
# To extract data directly from the database you need an active UGent VNP connection (or from inside the UGent)
# To connect to the database you need a specific login and password. (Contact Thomas (thomas.vergauwen@meteo.be) for this account.)

#The user and password are extracted from the environment-variables: VLINDER_DB_USER_NAME and VLINDER_DB_USER_PASW

#If this is done correctly than the user and password should appear when running settings.show()

 



# =============================================================================
#  Import data
# =============================================================================





#The toolkit is build around two classes: Station and Dataset.
    
    # The station class contains all observations and meta data of a station, and functions that apply to it.

    # The Dataset class is a collection of all station objects. All function that apply to multiple stations are defined in the Dataset class.
    

#1. Importing a dataset containing mulitple different stations is a function in the Dataset class. First we need to initiate a Dataset with a name of choise.

aug_2020_all_vlinders = vlinder_toolkit.Dataset()


# ---------------- Importing from CSV file -------------------------------


#The dataset is initiated but still empty. Filling it with the data from a csv file is simply done by:
    
aug_2020_all_vlinders.import_data_from_file(settings) #Rember that you added the input file in the settings object, this file will be used.


#You can see at any time what is in the dataset by:
aug_2020_all_vlinders.show()



# ---------------- Importing from Database -------------------------------

#First of all, make shure that your VNP connection (UGent) is on, or you ar working from within the Ugent !!!

# We have to start by first making an (empty) dataset
sept_2022_all_vlinders = vlinder_toolkit.Dataset()


#We need to specify the start-moment and end-moment for this period. To do this we need to import the datetime module (base python module)
from datetime import datetime

#To get the data it is a simple as:
sept_2022_all_vlinders.import_data_from_database(start_datetime=datetime(2022, 9,1), # 2022/09/01 00:00:00
                                                 end_datetime=datetime(2022,9,25,12,45)) #2022/09/25 12:45:00



# =============================================================================
# Analysing one station
# =============================================================================







#if you are interested in one station, you can extract all the info for that one station from the dataset by:

favorite_station = aug_2020_all_vlinders.get_station(stationname='vlinder02')

#the variable favorite_station now contains all info about that station. If you would like to have the observation in a tabular form (pandas dataframe):
    
observations = favorite_station.df()

print(observations.head())


#Or you can make a timeseries plot for a field of chose:

favorite_station.plot(variable='temp')

#You can extract attributes from the station:

temperature = favorite_station.temp








