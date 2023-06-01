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

import metobs_toolkit






#%%
# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
                        )

dataset.import_data_from_file()

dataset.coarsen_time_resolution()

dataset.get_altitude()
obs, outl = dataset.apply_titan_buddy_check()

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

# flags = titanlib.buddy_check(
#     points,
#     temp_obs,
#     radius,
#     num_min,
#     threshold,
#     max_elev_diff,
#     elev_gradient,
#     min_std,
#     num_iterations,
# )



# def create_titanlib_points_series(dataset, obstype, datetime):
#     import titanlib

#     # Set up some fake data
#     lats = [60,60.1,60.2]
#     lons = [10,10,10]
#     elevs = [0,0,0]
#     obs = [0, 1, -111]
#     points = titanlib.Points(lats, lons, elevs)




# dataset.apply_quality_control()

#%%

#%%
# dataset.make_plot(stationnames=['vlinder01', 'vlinder05'], colorby='label')

