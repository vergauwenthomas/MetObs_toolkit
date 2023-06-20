#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:29:37 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))


import metobs_toolkit

# %% IO testdata



dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=metobs_toolkit.demo_datafile,
                        input_metadata_file=metobs_toolkit.demo_metadatafile)
dataset.import_data_from_file()
dataset.coarsen_time_resolution()


# %% Apply Qc on dataset level

dataset.apply_quality_control(obstype="temp")


outliersdf = dataset.combine_all_to_obsspace()
dataset.get_qc_stats(make_plot=False)
dataset.get_qc_stats(obstype="humidity", make_plot=False)


# %% Apply Qc on obstype not specified in settings


dataset.apply_quality_control(obstype="humidity")
dataset.get_qc_stats(obstype="humidity", make_plot=False)
# %% Apply QC on station level

sta = dataset.get_station("vlinder05")

sta.apply_quality_control()

test = sta.get_qc_stats(make_plot=True)

# sta.get_qc_stats(make_plot=True)


#%% Apply titan checks

#  ------ Buddy check --------------
dataset.update_titan_qc_settings(obstype='temp',
                                 buddy_radius=50000,
                                 buddy_num_min=3,
                                 buddy_max_elev_diff=200,
                                 buddy_threshold=3)



dataset.apply_titan_buddy_check(use_constant_altitude=True)


# count test
assert dataset.outliersdf['label'].value_counts()['buddy check outlier'] == 57, 'The buddy check did not perfom good.'

# test if a check does not overwrite itself
dataset.update_titan_qc_settings(obstype='temp',
                                 buddy_radius=80000,
                                 buddy_num_min=3,
                                 buddy_max_elev_diff=200,
                                 buddy_threshold=0.5)
dataset.apply_titan_buddy_check(use_constant_altitude=True)

assert dataset.outliersdf['label'].value_counts()['buddy check outlier'] == 57, 'The buddy check did overwrite itself!'





#  ------------- SCT check ---------------


dataset.apply_titan_sct_resistant_check(use_constant_altitude=True)

