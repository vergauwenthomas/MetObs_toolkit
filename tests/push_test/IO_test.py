#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:28:15 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))


import metobs_toolkit


# %% define inputfiles


# %% import data from file (long standard format)

testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "vlinderdata_small.csv"
)


dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile)
dataset.show_settings()

dataset.import_data_from_file()

dataset.show()


station = dataset.get_station("vlinder02")


ax = station.make_plot()

ax = dataset.make_plot()


# %% import default dataset.


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    data_template_file=metobs_toolkit.demo_template,
)


dataset.show_settings()

dataset.import_data_from_file()

assert dataset.df.shape == (120957, 10), "Shape of demo data is not correct."
# dataset.show()


# %% Import wide dataset + syncronize

widedatafile = os.path.join(str(lib_folder), "tests", "test_data", "wide_test_data.csv")
widetemplate = os.path.join(
    str(lib_folder), "tests", "test_data", "wide_test_template.csv"
)


# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=widedatafile,
    # input_metadata_file=static_data,
    data_template_file=widetemplate,
)


dataset.import_data_from_file(long_format=False, obstype="temp")

assert dataset.df.shape == (597, 1), "Shape of unsynced widedata is not correct."

# Sycnronize data
test = dataset.sync_observations(tollerance="5T", verbose=True)


assert dataset.df.shape == (180, 1), "Shape after syncronizing widedata is not correct."

assert dataset.missing_obs.series.shape == (
    18,
), "Number of missing obs after sync wide data not correct"
