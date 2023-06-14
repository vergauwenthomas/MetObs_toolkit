#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

# %%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(lib_folder))

tmp_pickle = os.path.join(lib_folder, "development", "tmp", "dev_pickle.pkl")

import metobs_toolkit


# %%
import pandas as pd
import datetime


# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    data_template_file=metobs_toolkit.demo_template,
    metadata_template_file=metobs_toolkit.demo_template,  # Contains also the metadata mapping
    output_folder="/home/thoverga/Documents/VLINDER_github/MetObs_toolkitss",
)


# Load the data from the demo data files
dataset.import_data_from_file()

dataset.coarsen_time_resolution()
dataset.apply_quality_control()
dataset.update_gaps_and_missing_from_outliers(n_gapsize=4)

# %%
modeldata = metobs_toolkit.Modeldata("era5")

modeldata.set_model_from_csv(
    csvpath="/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/era5_data.csv"
)


# %%


gapfilldf1 = dataset.fill_gaps_linear()

# %%
gapfilldf2 = dataset.fill_gaps_automatic(modeldata, overwrite_fill=True)
