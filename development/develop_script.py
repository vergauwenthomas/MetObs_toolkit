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

dataset.fill_gaps_linear()
dataset.update_gaps_and_missing_from_outliers()
dataset.fill_missing_obs_linear()





#%%

test = dataset.write_to_csv(obstype = 'temp', use_tlk_obsnames=False)

print(test)




#%%

# comb_df = dataset.combine_all_to_obsspace()
# comb_df = comb_df[~comb_df.index.duplicated(keep='first')]
# comb_df.unstack('obstype')

# # to one level for the columns
# comb_df.columns = [' : '.join(col).strip() for col in comb_df.columns.values]






#%%
# nobs_orig = len(dataset.missing_obs.idx)
# ngaps_orig = len(dataset.gaps.list)

# comb_init = dataset.combine_all_to_obsspace()

# dataset.update_gaps_and_missing_from_outliers(obstype='temp', n_gapsize = 10)

# # nobs = len(dataset.missing_obs.idx)
# # ngaps = len(dataset.gaps.list)

# # assert (nobs == 43 & nobs_orig == 26), 'Something wrong with the update gaps and missing from outliers'
# # assert (ngaps == 5 & nobs_orig == 2), 'Something wrong with the update gaps and missing from outliers'


# # check if the mergedf does not contain them as duplicates
# comb2 = dataset.combine_all_to_obsspace()

# # assert (comb2.index.duplicated().any()), 'duplicated indexes in comb df after the outliers updated to gaps/missing'


# print(f'init shape {comb_init.shape}, after update: {comb2.shape}')

#%%
# outliers = dataset.outliersdf.xs('temp', level='obstype')




















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






