#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:28:15 2022

@author: thoverga
"""

from pathlib import Path
import os, sys

lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
from src import vlinder_toolkit

#%% define inputfiles
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata.csv')



#%% Define setting 

settings = vlinder_toolkit.Settings()
settings.show()


settings.update_settings(input_file=testdatafile)
settings.check_settings()

#%% import data from file

# dataset = vlinder_toolkit.Dataset()

# dataset.import_data_from_file()


# station = dataset.get_station('vlinder02')
# stationdf = station.df()
# print(stationdf.head())

#%% import data from DB
from datetime import datetime

dataset2 = vlinder_toolkit.Dataset()
dataset2.import_data_from_database(start_datetime=datetime(2022, 6,12),
                                    end_datetime=datetime(2022,6,19,12,45)) #2022/7/19 12:45:00


station = dataset2.get_station('vlinder02')
stationdf = station.df()
# print(stationdf.head())

print('Description: ',station.obs_description)
print('units: ',station.units)


#%%

# sta = dataset.get_station('vlinder02')

# ax =sta.make_plot()






