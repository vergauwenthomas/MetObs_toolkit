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

file = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide.csv"
file_template ="/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_template.csv"
file_metadata = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_metadata.csv"

# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = file,
                        input_metadata_file = file_metadata,
                        data_template_file = file_template,
                        metadata_template_file = file_template, # Contains also the metadata mapping
                        )

# Load the data from the demo data files
dataset.import_data_from_file(long_format=False,
                              obstype='temp',
                              obstype_dtype='float')

# dataset.coarsen_time_resolution()
#%%


dataset.make_plot()
