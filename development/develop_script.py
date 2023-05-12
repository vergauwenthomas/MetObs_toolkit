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


testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/vlinder_metadata.csv'
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')




# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        input_metadata_file=static_data,
                        # data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )

dataset.import_data_from_file()
dataset.coarsen_time_resolution()
dataset.apply_quality_control()
#%%
import copy
testdataset =copy.deepcopy(dataset)

# era = testdataset.get_modeldata()

modelfile = '/home/thoverga/Downloads/era5_data_(3).csv'

era = metobs_toolkit.Modeldata('era5')
era.set_model_from_csv(csvpath=modelfile,
                       )


#%%
dataset.fill_gaps_era5(modeldata=era,
                       )
#%%




