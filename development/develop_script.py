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


#%%

# # data
# era5_congo_file = '/home/thoverga/Downloads/era5_data_kongo.csv'
data_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/ATHTS01_all.csv'
# metadata_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/CONGO_meta.csv'
template_file ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/template.csv'
#%%

dataset = metobs_toolkit.Dataset()

dataset.update_settings(output_folder=None,
                        input_data_file=data_file,
                        # input_metadata_file=metobs_toolkit.demo_metadatafile,
                        data_template_file=template_file,
                        # metadata_template_file=metobs_toolkit.demo_template,
                        )


dataset.import_data_from_file()


#%%




dataset.update_qc_settings(obstype='temp',
                           step_max_decrease_per_sec=0.5,
                           step_max_increase_per_sec=0.5)


dataset.apply_quality_control(
    obstype="temp",         # choose which observations you want to check
    gross_value=False,       # set True if you want to perform the gross value check
    persistance=False,       # set True if you want to perform the persistence check
    step=True,              # set True if you want to perform the spike check
    window_variation=False,  # set True if you want to perform the window variation check
)



#%%

# dataset.make_plot(colorby='label')
comb1 = dataset.combine_all_to_obsspace()
comb1 = comb1.xs('temp', level='obstype')
print(f' bfore: {comb1["label"].value_counts()}')

#%%

from metobs_toolkit.gap import remove_gaps_from_outliers


# print(dataset.outliersdf.shape)

dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)


# print(dataset.outliersdf.shape)


#%%

# #%%
dataset.fill_missing_obs_linear()



dataset.make_plot(colorby='label')

#%%
