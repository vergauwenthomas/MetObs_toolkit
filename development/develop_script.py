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

dataset.update_gaps_and_missing_from_outliers()


print(dataset.gaps)
sta = dataset.get_station('vlinder01')


dataset.fill_gaps_linear('temp')

for gap in dataset.gaps:
    gap.get_info()







#%%
# import pandas as pd
# from metobs_toolkit.gap import _gap_collection_from_list_of_gaps


# import copy
# testdataset =copy.deepcopy(dataset)

# # era = testdataset.get_modeldata()

# modelfile = '/home/thoverga/Downloads/era5_data__(3).csv'

# era = metobs_toolkit.Modeldata('era5')
# era.set_model_from_csv(csvpath=modelfile,
#                         )

# #%%
# dataset.update_gap_and_missing_fill_settings(gap_debias_minimum_leading_period_hours=10,
#                                               gap_debias_minimum_trailing_period_hours=0,
#                                               gap_debias_prefered_leading_period_hours=30,
#                                               gap_debias_prefered_trailing_period_hours=30)



# test = dataset.fill_gaps_automatic(modeldata=era,
#                             obstype='temp',
#                             max_interpolate_duration_str='4H')






#%%
# dataset.fill_gaps_era5(modeldata=era,
#                        )
#%%






