#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 20 20:32:45 2023

@author: thoverga
"""


import sys, os
import pandas as pd
from pathlib import Path

# add the solutions
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
print(sys.path)
import solutions.solutions_creator as solution

# point to current version of the toolkit
lib_folder = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(lib_folder))
import metobs_toolkit


# print(metobs_toolkit.__version__)

# %% Import dataset

dataset = metobs_toolkit.Dataset()

dataset.update_file_paths(
    input_data_file=metobs_toolkit.demo_datafile,
    template_file=metobs_toolkit.demo_template,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
)


dataset.import_data_from_file()
dataset.coarsen_time_resolution()

# dataset.get_modeldata()

# %% test adding gee information
model_data = dataset.gee_datasets["ERA5-land"]

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
new_modobstype = metobs_toolkit.ModelObstype(
    obstype=new_obstype, model_unit="hpa", model_band="surface_pressure"
)
model_data.add_modelobstype(new_modobstype)

model_data.get_info()
from datetime import datetime

tstart = datetime(2022, 10, 3, 23)
tend = datetime(2022, 10, 4, 4)
model_data = dataset.get_modeldata(
    Model=model_data,
    obstypes=["special_pressure"],
    startdt=tstart,
    enddt=tend,
    force_direct_transfer=True,
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

    model_data.modeldf.to_pickle(
        os.path.join(solution.solutions_dir, surfpres_sol_file)
    )


# _create_surfpres_solutions()


def get_unfilled_solutions():
    surf_pres_model = pd.read_pickle(
        os.path.join(solution.solutions_dir, surfpres_sol_file)
    )
    return surf_pres_model


diff_df = solution.test_df_are_equal(
    testdf=model_data.modeldf, solutiondf=get_unfilled_solutions()
)
assert diff_df is None

model_data.make_plot(obstype_model="special_pressure")
# %% Test 2D vector fields


d2_sol_model = "d2_sol_model.pkl"


def _create_2d_sol(model_data):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    model_data = dataset.get_modeldata(
        Model=model_data,
        obstypes=["wind"],
        startdt=tstart,
        enddt=tend,
        force_direct_transfer=True,
    )

    model_data.modeldf.to_pickle(os.path.join(solution.solutions_dir, d2_sol_model))
    return


# _create_2d_sol(model_data)


def get_2d_solutions():
    d2_model = pd.read_pickle(os.path.join(solution.solutions_dir, d2_sol_model))
    return d2_model


model_data = dataset.get_modeldata(
    Model=model_data,
    obstypes=["wind"],
    startdt=tstart,
    enddt=tend,
    force_direct_transfer=True,
)


diff_df = solution.test_df_are_equal(
    testdf=model_data.modeldf, solutiondf=get_2d_solutions()
)
assert diff_df is None

# %% Testing multiple field extraction


multifield_modelfile = "multifield_sol_model.pkl"


def _create_multifield_sol(model_data):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    dataset.get_modeldata(
        Model=model_data,
        obstypes=["temp", "wind"],
        startdt=tstart,
        enddt=tend,
        force_direct_transfer=True,
    )

    model_data.modeldf.to_pickle(
        os.path.join(solution.solutions_dir, multifield_modelfile)
    )
    return


# _create_multifield_sol(model_data)


def get_multifield_solutions():
    multifield_model = pd.read_pickle(
        os.path.join(solution.solutions_dir, multifield_modelfile)
    )
    return multifield_model


dataset.get_modeldata(
    Model=model_data,
    obstypes=["temp", "wind"],
    startdt=tstart,
    enddt=tend,
    force_direct_transfer=True,
)


diff_df = solution.test_df_are_equal(
    testdf=model_data.modeldf, solutiondf=get_multifield_solutions()
)
assert diff_df is None
# %% Test get all bands


allfield_modelfile = "allfield_sol_model.pkl"


def _create_allfield_sol(model_data):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    dataset.get_modeldata(
        Model=model_data,
        obstypes=["temp", "wind"],
        startdt=tstart,
        enddt=tend,
        get_all_bands=True,
        force_direct_transfer=True,
    )

    model_data.modeldf.to_pickle(
        os.path.join(solution.solutions_dir, allfield_modelfile)
    )
    return


# _create_allfield_sol(model_data)


def get_allfield_solutions():
    allfield_model = pd.read_pickle(
        os.path.join(solution.solutions_dir, allfield_modelfile)
    )
    return allfield_model


dataset.get_modeldata(
    Model=model_data,
    obstypes=["temp", "wind"],
    startdt=tstart,
    get_all_bands=True,
    enddt=tend,
    force_direct_transfer=True,
)


diff_df = solution.test_df_are_equal(
    testdf=model_data.modeldf, solutiondf=get_allfield_solutions()
)
assert diff_df is None


# %% Import modeldata
model_data = metobs_toolkit.Dataset().gee_datasets["ERA5-land"]
# mutliple observations and vector components
csv_file = os.path.join(
    lib_folder, "tests", "test_data", "ERA5-land_timeseries_data.csv"
)

model_data.set_modeldata_from_csv(csv_file)

assert model_data.modeldf.shape == (
    960,
    71,
), "something wrong with reading modeldata from csv (drive)."
model_data.make_plot(obstype_model="wind_speed")
# %% Test repr

print(model_data)
model_data.get_info()

# %% test saving and importing
outfolder = os.path.join(lib_folder, "tests", "test_data")
pkl_file = "delete_me_if_you_see_me"
# save
model_data.set_metadf(dataset.metadf)
model_data.save_modeldata(outputfolder=outfolder, filename=pkl_file)

# read it again
newmod2 = metobs_toolkit.import_modeldata_from_pkl(
    folder_path=outfolder, filename=pkl_file + ".pkl"
)


diff_df = solution.test_df_are_equal(
    testdf=newmod2.modeldf, solutiondf=model_data.modeldf
)
assert diff_df is None


diff_df = solution.test_df_are_equal(
    testdf=newmod2.metadf, solutiondf=model_data.metadf
)
assert diff_df is None

# delete file
fullpath = os.path.join(outfolder, pkl_file + ".pkl")
if os.path.exists(fullpath):
    os.remove(fullpath)


# %%
model_data._subset_to_obstypes(
    [model_data.modelobstypes["temp"], model_data.modelobstypes["wind"]]
)

assert model_data.modeldf.shape[1] == 3, "something wrong with subset to obstype"

# %% test interpolation

mod_interp_pkl = "modl_interp_df.pkl"


def _create_interp_model_sol():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    # Create target
    dataset = metobs_toolkit.Dataset()
    dataset.update_file_paths(
        input_data_file=metobs_toolkit.demo_datafile,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        template_file=metobs_toolkit.demo_template,
    )

    dataset.import_data_from_file()
    target = dataset.df.xs("temp", level="obstype").index

    # Create modeltimeseries
    dataset.get_modeldata(
        Model=model_data,
        obstypes=["temp", "wind"],
        startdt=datetime(2022, 9, 4),
        enddt=datetime(2022, 9, 5, 16, 24),
        force_direct_transfer=True,
    )

    interpdf = model_data._interpolate_modeldata(target)
    interpdf.to_pickle(os.path.join(solution.solutions_dir, mod_interp_pkl))


# _create_interp_model_sol()


def get_interp_mmodel_sol():
    return pd.read_pickle(os.path.join(solution.solutions_dir, mod_interp_pkl))


# Create target
dataset = metobs_toolkit.Dataset()
dataset.update_file_paths(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

dataset.import_data_from_file()
target = dataset.df.xs("temp", level="obstype").index

# Create modeltimeseries
dataset.get_modeldata(
    Model=model_data,
    obstypes=["temp", "wind"],
    startdt=datetime(2022, 9, 4),
    enddt=datetime(2022, 9, 5, 16, 24),
    force_direct_transfer=True,
)

interpdf = model_data._interpolate_modeldata(target)


diff_df = solution.test_df_are_equal(
    testdf=interpdf, solutiondf=get_interp_mmodel_sol()
)
assert diff_df is None

# %% Test station Geemodel

sta = dataset.get_station("vlinder02")
assert sta.gee_datasets[
    "ERA5-land"
].metadf.empty, f"non empty metadf for gee model passing to Station"
assert sta.gee_datasets[
    "ERA5-land"
].modeldf.empty, f"non empty modeldf for gee model passing to Station"

# %% Test temeseries plotting

model_data.make_plot(
    obstype_model="wind_speed",
    Dataset=dataset,
    obstype_dataset="temp",
    show_outliers=False,
)
