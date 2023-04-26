#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%
import metobs_toolkit
import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))




#%% % Import


# testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')


template = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')



# #% Setup dataset



dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        # input_metadata_file=static_data,
                        data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )

# dataset.apply_quality_control()


dataset.import_data_from_file()
#%%
dataset.apply_quality_control()

test = dataset.combine_all_to_obsspace()


#%%
# dataset.make_plot(colorby='label')

sta = dataset.get_station('1')

