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
testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "vlinderdata_small.csv"
)


# %% import data from file

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile)
dataset.show_settings()

dataset.import_data_from_file()

dataset.show()


station = dataset.get_station("vlinder02")


ax = station.make_plot()

ax = dataset.make_plot()
