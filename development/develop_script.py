#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022

@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

#% Import

# testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt.csv')

# static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

# lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')


# #% Setup dataset
# settings = vlinder_toolkit.Settings()
# settings.update_settings(input_data_file=testdatafile,
#                           input_metadata_file=static_data,
#                           geotiff_lcz_file=lcz_map
#                          )


# dataset = vlinder_toolkit.Dataset()
# dataset.import_data_from_file(coarsen_timeres=True)




# dataset.apply_quality_control(obstype='temp')


# sta1 = dataset.get_station('vlinder05')
# sta1.make_plot()

# dataset.make_geo_plot()


#% Mocca data for Amber
inputdatafile = '/home/thoverga/Documents/Thesis studenten/Amber/data/merged_data.csv'

df = pd.read_csv(inputdatafile, sep=';')

settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=inputdatafile,
                         )




#%

#% Setup dataset




dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)
#%% apply QC
sta = dataset.get_station('SNZ')
sta.make_plot()


dataset.apply_quality_control()

sta = dataset.get_station('SNZ')
sta.make_plot()


#%%

export_csv_file = '/home/thoverga/Documents/Thesis studenten/Amber/data/checked_coarsed_merged_data.csv'


dataset.df.to_csv(export_csv_file)
