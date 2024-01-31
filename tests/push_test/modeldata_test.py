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

#%% Import dataset

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

#%% test adding gee information
model_data = metobs_toolkit.Modeldata(
    metadf=dataset.metadf, extractor=metobs_toolkit.GeeExtractor()
)


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

# Convert to modelobstype
new_obstype_modelform = metobs_toolkit.ModelObstype(
    obstype=new_obstype,
    band_name="surface_pressure",
    band_unit="hpa",
    band_description="oijmjimo",
)


# add new obstype to model_data
model_data.add_obstype(new_obstype_modelform)

assert set(model_data.obstypes.keys()) == set(
    ["temp", "pressure", "wind", "special_pressure"]
), "New obstype not added to the modeldata"


model_data.get_info()
from datetime import datetime

tstart = datetime(2022, 9, 3, 23)
tend = datetime(2022, 9, 4, 4)

model_data.extractor.activate_ERA5()


model_data.import_from_gee(
    target_obstypes=model_data.obstypes["special_pressure"],
    start_utc=tstart,
    end_utc=tend,
    gdrive_filename="era5_data",
)


assert (
    model_data.df.shape[0] == 168
), "No modeldata extracted from gee for new unit and obstype!"
assert model_data.df.columns.to_list() == [
    "special_pressure"
], "Something is wrong with column names"

model_data.make_plot(obstype_model="special_pressure")
#%% Test 2D vector fields

model_data.import_from_gee(
    target_obstypes=model_data.obstypes["wind"],
    start_utc=tstart,
    end_utc=tend,
)

print(model_data)

assert set(model_data.df.columns) == set(
    ["u_comp_wind", "v_comp_wind", "wind_amplitude", "wind_direction"]
), "something wrong with exploiting vector fields"
assert model_data.df.sum().to_dict() == {
    "u_comp_wind": 75.64337158203125,
    "v_comp_wind": 359.2424736022949,
    "wind_amplitude": 378.537614761677,
    "wind_direction": 9559.611737920999,
}, "something wrong with exploiting vector fields."

#%% Testing multiple field extraction
model_data.import_from_gee(
    target_obstypes=[model_data.obstypes["wind"], model_data.obstypes["temp"]],
    start_utc=tstart,
    end_utc=tend,
)

assert set(model_data.df.columns) == set(
    ["u_comp_wind", "v_comp_wind", "wind_amplitude", "wind_direction", "temp"]
), "Something is wrong with column names"


#%% Import modeldata
dummy_model_data = metobs_toolkit.Modeldata(
    metadf=dataset.metadf, extractor=metobs_toolkit.GeeExtractor()
)
# mutliple observations and vector components
csv_file = os.path.join(lib_folder, "tests", "test_data", "era5_modeldata_test.csv")

dummy_model_data.import_gee_data_from_csv(csv_file)

assert set(dummy_model_data.df.columns) == set(
    ["temp", "u_comp_wind", "v_comp_wind", "wind_amplitude", "wind_direction"]
), "something wrong with reading modeldata from csv (drive)."
dummy_model_data.make_plot(obstype_model="wind_amplitude")
#%% Test repr

print(model_data)

#%% test saving and importing
outfolder = os.path.join(lib_folder, "tests", "test_data")
pkl_file = "delete_me_if_you_see_me"
# save
model_data.save_modeldata(outputfolder=outfolder, filename=pkl_file, overwrite=True)

# read it again
# newmod = metobs_toolkit.Modeldata(metadf=dataset.)
newmod2 = metobs_toolkit.import_modeldata(
    target_pkl_file=os.path.join(outfolder, pkl_file + ".pkl")
)

# delete file
fullpath = os.path.join(outfolder, pkl_file + ".pkl")
if os.path.exists(fullpath):
    os.remove(fullpath)


#%% test interpolation

dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

dataset.import_data_from_file()

interpdf = model_data.sample_data_as(dataset.df)


assert (
    interpdf.columns == model_data.df.columns
).all(), "Error in modeldata interpolation"
assert interpdf.shape[0] == dataset.df.shape[0], "Error in modeldata interpolation"
assert (
    interpdf[~interpdf["temp"].isnull()].shape[0] == 1708
), "Error in modeldata interpolation"


#%% Test plotting


model_data.make_plot(stationnames=["vlinder01", "vlinder02"])


assert model_data.df.shape == (
    168,
    5,
), "Shape of modeldata df changed after plotting."

dataset.coarsen_time_resolution(freq="1H")
model_data.make_plot(dataset=dataset, show_outliers=False)
