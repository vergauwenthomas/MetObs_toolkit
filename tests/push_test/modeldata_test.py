#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:32:45 2023

@author: thoverga
"""


import sys, os
import pandas as pd
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
import solutions.solutions_creator as solution

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

surfpres_sol_file = "surfpress_era5.pkl"


def _create_surfpres_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    # % test adding gee information
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

    model_data.df.to_pickle(os.path.join(solution.solutions_dir, surfpres_sol_file))


# _create_surfpres_solutions()


def get_unfilled_solutions():
    surf_pres_model = pd.read_pickle(
        os.path.join(solution.solutions_dir, surfpres_sol_file)
    )
    return surf_pres_model


diff_df = solution.test_df_are_equal(
    testdf=model_data.df, solutiondf=get_unfilled_solutions()
)
assert diff_df is None

model_data.make_plot(obstype_model="special_pressure")
# %% Test 2D vector fields


d2_sol_model = "d2_sol_model.pkl"


def _create_2d_sol(model_data):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    model_data = dataset.get_modeldata(
        modeldata=model_data, obstype="wind_speed", startdt=tstart, enddt=tend
    )

    model_data.df.to_pickle(os.path.join(solution.solutions_dir, d2_sol_model))
    return


# _create_2d_sol(model_data)


def get_2d_solutions():
    d2_model = pd.read_pickle(os.path.join(solution.solutions_dir, d2_sol_model))
    return d2_model


model_data = dataset.get_modeldata(
    modeldata=model_data, obstype="wind_speed", startdt=tstart, enddt=tend
)


diff_df = solution.test_df_are_equal(
    testdf=model_data.df, solutiondf=get_2d_solutions()
)
assert diff_df is None

# %% Testing multiple field extraction


multifield_modelfile = "multifield_sol_model.pkl"


def _create_multifield_sol(model_data):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    model_data.get_gee_dataset_data(
        mapname=model_data.modelname,
        metadf=dataset.metadf,
        obstypes=["temp", "wind_speed"],
        startdt_utc=tstart,
        enddt_utc=tend,
    )

    model_data.df.to_pickle(os.path.join(solution.solutions_dir, multifield_modelfile))
    return


# _create_multifield_sol(model_data)


def get_multifield_solutions():
    multifield_model = pd.read_pickle(
        os.path.join(solution.solutions_dir, multifield_modelfile)
    )
    return multifield_model


model_data.get_gee_dataset_data(
    mapname=model_data.modelname,
    metadf=dataset.metadf,
    obstypes=["temp", "wind_speed"],
    startdt_utc=tstart,
    enddt_utc=tend,
)


diff_df = solution.test_df_are_equal(
    testdf=model_data.df, solutiondf=get_multifield_solutions()
)
assert diff_df is None


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

mod_interp_pkl = "modl_interp_df.pkl"


def _create_interp_model_sol():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset = metobs_toolkit.Dataset()
    dataset.update_settings(
        input_data_file=metobs_toolkit.demo_datafile,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        template_file=metobs_toolkit.demo_template,
    )

    dataset.import_data_from_file()

    target = dataset.df.xs("temp", level="obstype").index
    interpdf = model_data._interpolate_modeldata(target)
    interpdf.to_pickle(os.path.join(solution.solutions_dir, mod_interp_pkl))


# _create_interp_model_sol()


def get_interp_mmodel_sol():
    return pd.read_pickle(os.path.join(solution.solutions_dir, mod_interp_pkl))


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

dataset.import_data_from_file()

target = dataset.df.xs("temp", level="obstype").index
interpdf = model_data._interpolate_modeldata(target)


diff_df = solution.test_df_are_equal(
    testdf=interpdf, solutiondf=get_interp_mmodel_sol()
)
assert diff_df is None


# %% Test plotting

a = model_data.df.shape

model_data.make_plot(stationnames=["vlinder01", "vlinder02"])


assert model_data.df.shape == (
    10108,
    3,
), "Shape of modeldata df changed after plotting."


model_data.make_plot(dataset=dataset, show_outliers=False)
