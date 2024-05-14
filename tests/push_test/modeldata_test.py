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
    template_file=metobs_toolkit.demo_template,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
)


dataset.import_data_from_file()
dataset.coarsen_time_resolution()

# dataset.get_modeldata()

# %% test adding gee information
model_data = metobs_toolkit.Modeldata("ERA5_hourly")

# Define a regular obstype
new_obstype = metobs_toolkit.Obstype(
    obsname="special_pressure",
    std_unit="pa",
    description="just for testing",
    unit_aliases={
        "pa": ["Pascal", "Pa", "N/m²"],
    },
    unit_conversions={"hpa": ["x * 100"]},
)

# add new obstype to model_data
model_data.add_obstype(
    Obstype=new_obstype,
    bandname="surface_pressure",
    band_units="hpa",
)


model_data.get_info()
from datetime import datetime

tstart = datetime(2022, 10, 3, 23)
tend = datetime(2022, 10, 4, 4)
model_data = dataset.get_modeldata(
    modeldata=model_data, obstype="special_pressure", startdt=tstart, enddt=tend
)

assert (
    model_data.df.shape[0] == 168
), "No modeldata extracted from gee for new unit and obstype!"
assert model_data.df.columns.to_list() == [
    "special_pressure"
], "Something is wrong with column names"

model_data.make_plot(obstype_model="special_pressure")
# %% Test 2D vector fields

model_data = dataset.get_modeldata(
    modeldata=model_data, obstype="wind_speed", startdt=tstart, enddt=tend
)

print(model_data)

assert model_data.df.columns.to_list() == [
    "wind_speed_amplitude",
    "wind_speed_direction",
], "Something is wrong with column names"


# %% Testing multiple field extraction
model_data.get_gee_dataset_data(
    mapname=model_data.modelname,
    metadf=dataset.metadf,
    obstypes=["temp", "wind_speed"],
    startdt_utc=tstart,
    enddt_utc=tend,
)

assert model_data.df.columns.to_list() == [
    "temp",
    "wind_speed_amplitude",
    "wind_speed_direction",
], "Something is wrong with column names"


# %% Import modeldata
model_data = metobs_toolkit.Modeldata("ERA5_hourly")
# mutliple observations and vector components
csv_file = os.path.join(lib_folder, "tests", "test_data", "era5_modeldata_test.csv")

model_data.set_model_from_csv(csv_file)

assert model_data.df.columns.to_list() == [
    "temp",
    "wind_speed_amplitude",
    "wind_speed_direction",
], "something wrong with reading modeldata from csv (drive)."
model_data.make_plot(obstype_model="wind_speed_amplitude")
# %% Test repr

print(model_data)

# %% test saving and importing
outfolder = os.path.join(lib_folder, "tests", "test_data")
pkl_file = "delete_me_if_you_see_me"
# save
model_data.save_modeldata(outputfolder=outfolder, filename=pkl_file)

# read it again
newmod = metobs_toolkit.Modeldata("ERA5_hourly")
newmod2 = newmod.import_modeldata(folder_path=outfolder, filename=pkl_file + ".pkl")

# delete file
fullpath = os.path.join(outfolder, pkl_file + ".pkl")
if os.path.exists(fullpath):
    os.remove(fullpath)


# %% test interpolation

dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

dataset.import_data_from_file()

interpdf = model_data.interpolate_modeldata(dataset.df.index)


# test that there are no nan values
if not interpdf[interpdf["temp"].isnull()].empty:
    sys.exit("Error in modeldata interpolation")


assert interpdf.shape == (120957, 3), "Error in modeldata interpolation"

# check if other obstypes are interpolated as well
assert (
    interpdf["wind_speed_amplitude"].shape[0] == 120957
), "Error in modeldata interpolation"
assert (
    interpdf["wind_speed_amplitude"][interpdf["wind_speed_amplitude"].isnull()].shape[0]
    == 0
), "Error in modeldata interpolation"


# %% Test plotting

a = model_data.df.shape

model_data.make_plot(stationnames=["vlinder01", "vlinder02"])


assert model_data.df.shape == (
    10108,
    3,
), "Shape of modeldata df changed after plotting."


model_data.make_plot(dataset=dataset, show_outliers=False)
