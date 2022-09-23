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
settings.update_settings(input_file=testdatafile)


#%% import data from file

dataset = vlinder_toolkit.Dataset()

dataset.import_data_from_file(Settings=settings)



#%%
station = dataset.get_station('vlinder02').df()
print(station.head())




#%%

sta = dataset.get_station('vlinder02')

ax =sta.plot()






