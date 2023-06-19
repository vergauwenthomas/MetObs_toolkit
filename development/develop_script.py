#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit



# # data
# era5_congo_file = '/home/thoverga/Downloads/era5_data_kongo.csv'
data_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/ATHTS01_all.csv'
# metadata_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/CONGO_meta.csv'
template_file ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/template.csv'
#%%

dataset = metobs_toolkit.Dataset()

dataset.update_settings(output_folder=None,
                        # input_data_file=metobs_toolkit.demo_datafile,
                        input_data_file = data_file,
                        # input_metadata_file=metobs_toolkit.demo_metadatafile,
                        # data_template_file=metobs_toolkit.demo_template,
                        data_template_file = template_file,
                        # metadata_template_file=metobs_toolkit.demo_template,
                        )


dataset.import_data_from_file()
# dataset.coarsen_time_resolution()

dataset.apply_quality_control()

#%%
# dataset.make_plot(colorby='label')
# dataset.make_plot(colorby='label')

dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)
dataset.make_plot(colorby='label')


dataset.fill_gaps_linear()
dataset.make_plot(colorby='label')


# from metobs_toolkit.df_helpers import xs_save




# mergedf = dataset.combine_all_to_obsspace()
# mergedf = mergedf.xs('temp', level='obstype')

# #%%

# import pandas as pd
# import numpy as np



# # line_labels = ['ok', 'gap_interpolation', 'gap_debiased_era5', 'missing_obs_interpolation']

# sta = 'vlinder05'
# # for sta in mergedf.index.get_level_values('name').unique():
# stadf = xs_save(mergedf, sta, 'name')
# stadf_idx = stadf.index

# linedf = stadf[stadf['label'].isin(line_labels)]
# stadf.loc[~stadf.index.isin(linedf.index), 'value'] = np.nan



#%%
# dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)
# dataset.fill_missing_obs_linear()




#%%

staname = 'n13'
# print(dataset.outliersdf.shape)
dataset.make_plot(colorby='label', title='after qc')
dataset.get_station(staname).make_plot(colorby='label', title='station: after qc')
dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)
dataset.make_plot(colorby='label', title='after qc updating to gaps and missing')
dataset.get_station(staname).make_plot(colorby='label', title='station: after qc updateing to gaps and missing')
dataset.fill_missing_obs_linear()
dataset.make_plot(colorby='label', title='after fill missing')
dataset.get_station(staname).make_plot(colorby='label', title='station: after fill missing')
dataset.fill_gaps_linear()
dataset.make_plot(colorby='label', title='after fill missing and gapsfill')
dataset.get_station(staname).make_plot(colorby='label', title='station: after fill missing and gapsfill')




# dataset.update_qc_settings(obstype='temp',
#                            step_max_decrease_per_sec=0.5,
#                            step_max_increase_per_sec=0.5)


# dataset.apply_quality_control(
#     obstype="temp",         # choose which observations you want to check
#     gross_value=False,       # set True if you want to perform the gross value check
#     persistance=False,       # set True if you want to perform the persistence check
#     step=True,              # set True if you want to perform the spike check
#     window_variation=False,  # set True if you want to perform the window variation check
# )



#%%
# dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)
# dataset.fill_missing_obs_linear()

#%%
# dataset.make_plot(colorby='label', title='after fill missing')


#%%







# combdf = dataset.combine_all_to_obsspace()
# #%%
# line_labels = ['ok', 'gap_interpolation', 'gap_debiased_era5', 'missing_obs_interpolation']

# test = combdf[combdf['label'].isin(line_labels)]
# print(test['label'].value_counts())









