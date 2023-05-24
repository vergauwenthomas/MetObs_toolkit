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
# data Ian
# datafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/WBGTdata_totaal.csv'

# data Wout
datafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/Outdoor_module_Netatmo_Sara_new.csv'
metafile ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/metadata_Outdoor_module_Netatmo_Sara_new.csv'



# metobs_toolkit.build_template_prompt()



#%%

data_file = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/Outdoor_module_Netatmo_Sara_new.csv"
meta_data_file = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/metadata_Outdoor_module_Netatmo_Sara_new.csv"
data_template = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/template.csv"
meta_data_template = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/template.csv"

your_dataset = metobs_toolkit.Dataset()


your_dataset.update_settings(
    input_data_file = data_file,
    data_template_file = data_template,
    input_metadata_file = meta_data_file,
    metadata_template_file = meta_data_template,
    )

your_dataset.update_timezone(timezonestr = 'Europe/Brussels')
your_dataset.update_qc_settings(gapsize_in_records = 40)


your_dataset.import_data_from_file(
    long_format=True,
    freq_estimation_method='median', #None for default(="highest"), "highest" or "median"
    freq_estimation_simplify=None, #None for default(=True), True or False
    freq_estimation_simplify_error='5T', #None for default(="2T"), or timedelta string.
    )

# your_dataset.coarsen_time_resolution(freq='1H')

#%%
import copy

dataset = copy.deepcopy(your_dataset)

dataset.make_plot()

