#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:32:45 2023

@author: thoverga
"""


import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]

import metobs_toolkit

# print(metobs_toolkit.__version__)

# %% Import dataset

dataset = metobs_toolkit.Dataset()

dataset.update_settings(
    output_folder=None,
    input_data_file=metobs_toolkit.demo_datafile,
    data_template_file=metobs_toolkit.demo_template,
)


dataset.import_data_from_file()
dataset.coarsen_time_resolution()


# %% Import modeldata
model_data = metobs_toolkit.Modeldata("ERA5")

csv_file = os.path.join(lib_folder, "tests", "test_data", "era5_data.csv")

model_data.set_model_from_csv(csv_file)

# %% Test repr

print(model_data)


# %% Test plotting

a = model_data.df.shape

model_data.make_plot(stationnames=["vlinder01", "vlinder02"])


assert model_data.df.shape == (5348, 1), "Shape of modeldata df changed after plotting."


model_data.make_plot(dataset=dataset, show_outliers=False)
