#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:28:15 2022

@author: thoverga
"""

import sys,os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit


#%% define inputfiles
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata_small.csv')



#%% Define setting 

settings = vlinder_toolkit.Settings()
settings.show()


settings.update_settings(input_data_file=testdatafile)
settings.check_settings()

#%% import data from file

dataset = vlinder_toolkit.Dataset()

dataset.import_data_from_file()


station = dataset.get_station('vlinder02')
stationdf = station.df()
print(stationdf.head())

ax =station.make_plot()








