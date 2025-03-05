#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:28:15 2022

Test importing a variaty of data/metadata combinations and formats


@author: thoverga
"""

import sys, os
from pathlib import Path
import pandas as pd

repodir = str(Path(__file__).resolve().parents[2])

# add the solutions
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
print(sys.path)
import solutions.solutions_creator as solution

# point to current version of the toolkit
lib_folder = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(lib_folder))
import metobs_toolkit

# %%

# point to current version of the toolkit
lib_folder = repodir


# %% import data from file (long standard format)

testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "vlinderdata_small.csv"
)


dataset = metobs_toolkit.Dataset()

dataset.import_data_from_file(
    input_data_file=testdatafile, template_file=metobs_toolkit.demo_template
)
dataset.get_info()
dataset.template.get_info()
station = dataset.get_station("vlinder02")

# %% Import dataset from online location
dataloc = "https://raw.githubusercontent.com/vergauwenthomas/MetObs_toolkit/master/tests/test_data/vlinderdata_small.csv"
template_loc = "https://raw.githubusercontent.com/vergauwenthomas/MetObs_toolkit/master/metobs_toolkit/datafiles/demo_template.json"


# dataset = metobs_toolkit.Dataset()
# dataset.import_data_from_file(input_data_file=dataloc,
#                                template_file=template_loc,
#                               templatefile_is_url=True)
# assert not dataset.df.empty, "something wrong with importing from onlin locations"


# %% import default dataset.


dataset = metobs_toolkit.Dataset()

dataset.import_data_from_file(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

assert dataset.df.shape == (483840, 2), "Shape of demo data is not correct."

# %% Test template class methods

dataset.template.get_info()


# %% Import wide dataset (Multiple stations) + syncronize

widedatafile = os.path.join(str(lib_folder), "tests", "test_data", "wide_test_data.csv")
widetemplate = os.path.join(
    str(lib_folder), "tests", "test_data", "wide_test_template.json"
)


# #% Setup dataset
dataset = metobs_toolkit.Dataset()

dataset.import_data_from_file(
    input_data_file=widedatafile,
    template_file=widetemplate,
    freq_estimation_method="median",
    freq_estimation_simplify_tolerance="2min",
    origin_simplify_tolerance="5min",
    timestamp_tolerance="4min",
)


assert int(dataset.df["value"].min()) == -268, "Unit is not converted properly "

assert dataset.df.shape == (197, 2), "Shape of synced widedata is not correct."

for sta in dataset.stations:
    assert sta.obsdata["temp"].freq == pd.Timedelta(
        "1h"
    ), "frequency of records not correct synchronized"


# %% Test syncronizing wide


syncedwidedf_file = "synced_wide_test_combdf.pkl"


def _create_widesync_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    # Sycnronize data
    dataset.sync_records(
        timestamp_shift_tolerance="5min2s", freq_shift_tolerance="2min"
    )
    dataset.df.to_pickle(os.path.join(solution.solutions_dir, syncedwidedf_file))


# _create_widesync_solutions()


def get_synced_solutions():
    sol_combdf = pd.read_pickle(os.path.join(solution.solutions_dir, syncedwidedf_file))
    return sol_combdf


# Sycnronize data
dataset.sync_records(timestamp_shift_tolerance="5min2s", freq_shift_tolerance="2min")

assert dataset.df.shape == (196, 2), "something wrong with syncronization"

diff_df = solution.test_df_are_equal(
    testdf=dataset.df, solutiondf=get_synced_solutions()
)
assert diff_df is None

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


# %% Setup dataset

dataset_single = metobs_toolkit.Dataset()
dataset_single.import_data_from_file(
    input_data_file=singlestationdatafile,
    input_metadata_file=singlestationmetadata,
    template_file=singlestationtemplate,
)

assert dataset_single.df.shape == (26, 2), "Shape singlestation dataset is not correct."

assert (
    dataset_single.df.index.get_level_values("name")[0]
    == "dummy2_station_name"  # Must be the name defined in the template/data rather than in the matedata
), "The single station name in the data/template is not set for the metadata."

assert (
    dataset_single.metadf.shape[0] == 1
), "Shape metadf for single station is not correct"


assert (
    dataset_single.metadf["lat"].iloc[0] == 51.558
), "Metadf latitde is not merged correct."

assert (
    dataset_single.df.index.get_level_values("name").unique()[0]
    == "dummy2_station_name"
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


dataset.rename_stations({"MolenhofWeer - Temperatuur (°C)": "this_is_a_test_name"})

dataset.save_dataset_to_pkl(target_folder=outfolder, filename=file)


del dataset  # remove from kernel


# read dataset
new_dataset = metobs_toolkit.import_dataset_from_pkl(
    target_path=os.path.join(outfolder, file + ".pkl")
)

del_file(os.path.join(outfolder, file + ".pkl"))


# %%
# =============================================================================
# Testing the IO properties for new observation types and units
# =============================================================================

dataset = metobs_toolkit.Dataset()

n_obstypes = len(dataset.obstypes)

# test addition of obstype

wetbulp_obstype = metobs_toolkit.Obstype(
    obsname="wetbulptemp",
    std_unit="°F",
    description="THe wet bulb temperature",
)
wetbulp_obstype.get_info()


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


dataset.import_data_from_file(
    input_data_file=testdata,
    input_metadata_file=testmetadata,
    template_file=testtemplate,
)
dataset.template.get_info()


# test if all obstypes are present in the dataset
assert dataset.df.index.get_level_values("obstype").unique().to_list() == [
    "temp",
    "wetbulptemp",
], "New obstype not use when importing data"

# check if the unist of the obstypes are correct (the default)
assert (
    dataset.obstypes["temp"].std_unit == "degree_Celsius"
), "Standard unit not correct"
assert (
    dataset.obstypes["wetbulptemp"].std_unit == "degree_Fahrenheit"
), "Standard unit not correct"

# Check if unit conversion is done

assert (
    dataset.df.xs("temp", level="obstype")["value"].mean() < -200.0
), "THe units of the temperature observations are not converted to std units"


# %% Testing importing a "metadata-only" dataset with

testmetadata = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_metadata.csv"
)

testtemplate = os.path.join(
    str(lib_folder), "tests", "test_data", "single_station_new_obstype_template.json"
)

dataset = metobs_toolkit.Dataset()

dataset.import_data_from_file(
    input_metadata_file=testmetadata, template_file=testtemplate
)

assert dataset.df.empty, "metadata only dataset has data"
assert dataset.metadf.shape == (1, 3), "metadata not correct in metadata only dataset"


dataset = metobs_toolkit.Dataset()

dataset.import_data_from_file(
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

dataset.get_info()


assert dataset.df.empty, "metadata only dataset has data"
assert dataset.metadf.shape == (28, 4), "metadata not correct in metadata only dataset"

# %%
