#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:28:15 2022

Test importing a variaty of data/metadata combinations and formats


@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]

import metobs_toolkit

# print(metobs_toolkit.__version__)

# %%
# %% import data from file (long standard format)

testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "vlinderdata_small.csv"
)


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=testdatafile, template_file=metobs_toolkit.demo_template
)
dataset.show_settings()

dataset.import_data_from_file()

dataset.show()


station = dataset.get_station("vlinder02")


# %% import default dataset.


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)


dataset.show_settings()

dataset.import_data_from_file()

assert dataset.df.shape == (120957, 4), "Shape of demo data is not correct."


# %% Import wide dataset (Multiple stations) + syncronize

widedatafile = os.path.join(str(lib_folder), "tests", "test_data", "wide_test_data.csv")
widetemplate = os.path.join(
    str(lib_folder), "tests", "test_data", "wide_test_template.json"
)


# #% Setup dataset
dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=widedatafile,
    # input_metadata_file=static_data,
    template_file=widetemplate,
)


dataset.import_data_from_file(
    # long_format=False,
    # obstype="temp",
    # obstype_description="2mT",
    # obstype_unit="Celcius"
)

assert dataset.df.shape == (597, 1), "Shape of unsynced widedata is not correct."


# %% Test syncronizing wide

# Sycnronize data
test = dataset.sync_observations(tolerance="5T", verbose=True)


assert dataset.df.shape == (182, 1), "Shape after syncronizing widedata is not correct."

assert dataset.missing_obs.series.shape == (
    15,
), "Number of missing obs after sync wide data not correct"


# %% import wide dataset (One station)

singlestationdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station.csv"
)
singlestationtemplate = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_template.json"
)
singlestationmetadata = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_metadata.csv"
)


# #% Setup dataset

dataset_single = metobs_toolkit.Dataset()
dataset_single.update_settings(
    input_data_file=singlestationdatafile,
    input_metadata_file=singlestationmetadata,
    template_file=singlestationtemplate,
)


dataset_single.import_data_from_file()

assert dataset_single.df.shape == (13, 2), "Shape singlestation dataset is not correct."

assert (
    dataset_single.df.index.get_level_values("name")[0] == "whats_the_name"
), "The single station name in the metadata is not set for the data."

assert (
    dataset_single.metadf.shape[0] == 1
), "Shape metadf for single station is not correct"


assert (
    dataset_single.metadf["lat"].iloc[0] == 51.558
), "Metadf latitde is not merged correct."

assert (
    dataset_single.df.index.get_level_values("name").unique()[0] == "whats_the_name"
), "single station name not represented correctly."


# %%

# helper


def del_file(file_path):
    if os.path.isfile(file_path):
        os.remove(file_path)
        print(f"{file_path} deleted.")
    else:
        print(f"{file_path} not found.")


# %% Pickle save and read dataset
outfolder = os.path.join(str(lib_folder), "tests", "test_data")
file = "dataset_IO_test"


del_file(os.path.join(outfolder, file + ".pkl"))


# save dataset as pickle


dataset.update_default_name("this_is_a_test_name")

dataset.save_dataset(outputfolder=outfolder, filename=file)


del dataset  # remove from kernel


# read dataset
new_dataset = metobs_toolkit.Dataset()
new_dataset = new_dataset.import_dataset(folder_path=outfolder, filename=file + ".pkl")

del_file(os.path.join(outfolder, file + ".pkl"))

assert (
    new_dataset.settings.app["default_name"] == "this_is_a_test_name"
), "some attributes are not correctly saved when pickled."


# =============================================================================
# Testing the IO properties for new observation types and units
# =============================================================================

dataset = metobs_toolkit.Dataset()

n_obstypes = len(dataset.obstypes)
# add unit to unexisting obstype
dataset.add_new_unit(
    obstype="wetbulptem", new_unit="fake_wbtemp", conversion_expression=["x+100"]
)

new_n_obstypes = len(dataset.obstypes)

assert (
    n_obstypes == new_n_obstypes
), "Adding a new unit to unexisting obstype creates and obstype!"


# test addition of obstype and unit
dataset.add_new_unit(
    obstype="temp", new_unit="fake_temp", conversion_expression=["x+100"]
)

wetbulp_obstype = metobs_toolkit.Obstype(
    obsname="wetbulptemp",
    std_unit="Celcius",
    description="THe wet bulb temperature",
    unit_aliases={"Celcius": ["°C"], "Kelvin": ["K"]},
    unit_conversions={"Kelvin": ["x-273"]},
)
dataset.add_new_observationtype(wetbulp_obstype)
new_n_obstypes = len(dataset.obstypes)

assert n_obstypes == new_n_obstypes - 1, "Adding a new obstype not stored in dataset!"


# test if data can be imported with the new obstype and the new unit

testdata = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_new_obstypes.csv"
)
testmetadata = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_metadata.csv"
)
testtemplate = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_new_obstype_template.json"
)

dataset.update_settings(
    input_data_file=testdata,
    input_metadata_file=testmetadata,
    template_file=testtemplate,
)

dataset.import_data_from_file()


# test if all obstypes are present in the dataset
assert list(dataset.df.columns) == [
    "temp",
    "wetbulptemp",
], "New obstype not use when importing data"

# check if the unist of the obstypes are correct (the default)
assert (
    dataset.obstypes["temp"].get_standard_unit() == "Celsius"
), "Standard unit not correct"
assert (
    dataset.obstypes["wetbulptemp"].get_standard_unit() == "Celcius"
), "Standard unit not correct"

# Check if unit conversion is done

assert (
    dataset.df["temp"].mean() > 100.0
), "THe units of the temperature observations are not converted to std units"
