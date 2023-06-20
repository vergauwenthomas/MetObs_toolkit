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
                        input_data_file=metobs_toolkit.demo_datafile,
                        # input_data_file = data_file,
                        input_metadata_file=metobs_toolkit.demo_metadatafile,
                        data_template_file=metobs_toolkit.demo_template,
                        # data_template_file = template_file,
                        metadata_template_file=metobs_toolkit.demo_template,
                        )


dataset.import_data_from_file()
dataset.coarsen_time_resolution()

# dataset.apply_quality_control()

#%%


# print(dataset.settings.qc['titan_check_settings']['titan_buddy_check'])


# dataset.update_titan_qc_settings(obstype='temp', buddy_radius = 1567.2)


# print(dataset.settings.qc['titan_check_settings']['titan_buddy_check'])




#%%
# dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)
# dataset.fill_missing_obs_linear()


# dataset.get_altitude()


dataset.update_titan_qc_settings(obstype='temp', sct_tneg = 12.3,
                                 buddy_radius=2345.6, sct_outer_radius=30200)



dataset.apply_titan_buddy_check(use_constant_altitude=True)
# dataset.apply_titan_sct_resistant_check(use_constant_altitude=True)
# dataset.make_plot(colorby='label')





#%%

# import numpy as np
# import titanlib
# # Set up some fake data
# lats = [60,60.1,60.2]
# lons = [10,10,10]
# elevs = [0,0,0]
# temp_obs = [0, 1, -111]
# points = titanlib.Points(lats, lons, elevs)

# radius = np.full(points.size(), 5000)
# num_min = np.full(points.size(), 5)
# num_min = 5
# threshold = 2
# max_elev_diff = 200
# elev_gradient = -0.0065
# min_std = 1
# num_iterations = 5






