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

# %Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt.csv')

static_data =  os.path.join(str(lib_folder), 'static_data', 'vlinder_metadata.csv')

lcz_map = os.path.join(str(lib_folder), 'physiograpy', 'lcz_filter_v1.tif')

 
#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile,
                           # input_metadata_file=static_data,
                          geotiff_lcz_file=lcz_map,
                          output_folder=os.path.join(str(lib_folder), 'temp_output')
                          )


dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)

dataset.apply_quality_control(obstype='temp')


dataset.make_geo_plot()


#%%

vlinderlist = []
for i in range(28):
    # print(i+1)
    vlindername = 'vlinder' + str(i+1).zfill(2)
    
    vlinderlist.append(vlindername)

print(vlinderlist)
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata.csv')

df= pd.read_csv(testdatafile, sep=';')

df =df[df['Vlinder'].isin(vlinderlist)]

outputfile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata_small.csv')

df.to_csv(outputfile, sep=';', index=False)
