#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:18:09 2022

@author: thoverga
"""
import copy
import sys, os
import pandas as pd
from datetime import datetime
import numpy as np

from metobs_toolkit.df_helpers import init_multiindexdf, conv_tz_multiidxdf


from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]

# sys.path.append(str(lib_folder))


import metobs_toolkit

# %% Import data

testdatafile = os.path.join(
    str(lib_folder), "tests", "test_data", "testdata_okt_small.csv"
)

static_data = os.path.join(str(lib_folder), "static_data", "vlinder_metadata.csv")


# #% import data


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=testdatafile,
    input_metadata_file=static_data,
    output_folder="/home/thoverga/Documents/VLINDER_github/metobs_toolkit",
)
dataset.import_data_from_file()
dataset.coarsen_time_resolution()
# %% Basic tests on the gaps
gapsdf = dataset.get_gaps_df()
# check if two gaps are found
assert (
    gapsdf.shape[0] == 2
), f"There are assumed 2 gaps, but the tlkit found {gapsdf.shape[0]}"

assert list(gapsdf.index.unique()) == [
    "vlinder01"
], f"Only gaps assumed in vlinder01. Tlkit found gaps for these {list(gapsdf.index.unique())}"


# %% Basic tests on missing obs
missingobs = dataset.missing_obs.series
# check if two gaps are found
assert (
    missingobs.shape[0] == 26
), f"There are assumed 26 missing obs, but the tlkit found {missingobs.shape[0]}"

assert list(missingobs.index.unique()) == [
    "vlinder01",
    "vlinder02",
    "vlinder03",
], f"Only missing obs assumed in vl01, vl02, vl03. Tlkit found missing obs for these {list(missingobs.index.unique())}"

#%% Test linear interpolation on missing obs
dataset.fill_missing_obs_linear()

solution = {
    "temp": {
        (
            "vlinder03",
            pd.Timestamp("2022-10-07 11:00:00+0000", tz="UTC"),
        ): 15.333333333333334,
        (
            "vlinder03",
            pd.Timestamp("2022-10-07 12:00:00+0000", tz="UTC"),
        ): 16.566666666666666,
    },
    "temp_final_label": {
        (
            "vlinder03",
            pd.Timestamp("2022-10-07 11:00:00+0000", tz="UTC"),
        ): "missing_obs_interpolation",
        (
            "vlinder03",
            pd.Timestamp("2022-10-07 12:00:00+0000", tz="UTC"),
        ): "missing_obs_interpolation",
    },
}


assert dataset.missing_fill_df.equals(
    pd.DataFrame(solution)
), "something wrong with the missing obs fill!"

# %% Test functions on gaps
from metobs_toolkit.gap import (
    get_station_gaps,
    get_gaps_indx_in_obs_space,
    remove_gaps_from_obs,
)

get_station_gaps(dataset.gaps, "vlinder01")


dataset.gaps[0].get_info()

remove_gaps_from_obs(dataset.gaps, dataset.df)

get_gaps_indx_in_obs_space(
    dataset.gaps, dataset.df, dataset.outliersdf, dataset.metadf["dataset_resolution"]
)


dataset.gaps[0].update_leading_trailing_obs(dataset.df, dataset.outliersdf)

# %% Test gapfilling

# ------------------ linear interpolation ----------------------------

dataset.fill_gaps_linear(obstype="temp")

gapsfilldf = dataset.gapfilldf

gapsfilldf["manual"] = [
    9.487500,
    9.875000,
    10.262500,
    10.650000,
    11.037500,
    11.425000,
    11.812500,
    15.833333,
    14.766667,
    13.700000,
    12.633333,
    11.566667,
]

gapsfilldf["diff"] = gapsfilldf["temp"] - gapsfilldf["manual"]


assert (
    not gapsfilldf["temp"].isnull().any()
), f"np.nan value found in tlk interpolation ! \n {gapsfilldf}"
assert (
    gapsfilldf["diff"]
).sum() < 1e-5, f"Tlk interpolation differs from manual: \n {gapsfilldf}"


#%% Test if filled values are present in the combined df

comb_df = dataset.combine_all_to_obsspace()
comb_df = comb_df.xs("temp", level="obstype")


comb_gaps = comb_df.loc[dataset.gapfilldf.index]
comb_missing = comb_df.loc[dataset.missing_fill_df.index]

assert (
    comb_gaps["value"].eq(dataset.gapfilldf["temp"]).all()
), "Something wrong with the filled gaps in the combined df"
assert (
    comb_missing["value"].eq(dataset.missing_fill_df["temp"]).all()
), "Something wrong with the filled missing in the combined df"

#%% Test the update of outliers to gaps
nobs_orig = len(dataset.missing_obs.idx)
ngaps_orig = len(dataset.gaps)

dataset2 = copy.deepcopy(dataset)
dataset2.apply_quality_control()
outliersbefore = copy.deepcopy(dataset2.outliersdf.xs("temp", level="obstype"))
missingbefore = copy.deepcopy(dataset2.missing_obs.series)
dataset2.update_gaps_and_missing_from_outliers(obstype="temp", n_gapsize=10)
missingafter = copy.deepcopy(dataset2.missing_obs.series)

nobs = len(dataset2.missing_obs.idx)
ngaps = len(dataset2.gaps)

assert (nobs == 28) & (
    nobs_orig == 26
), "Something wrong with the update gaps and missing from outliers"
assert (ngaps == 5) & (
    ngaps_orig == 2
), "Something wrong with the update gaps and missing from outliers"


# check if the mergedf does not contain them as duplicates
comb2 = dataset2.combine_all_to_obsspace()

assert (
    comb2[comb2.index.duplicated()].shape[0] == 0
), "duplicated indexes in comb df after the outliers updated to gaps/missing"


# %%
# ------------------ ERA debias fill
obstype = "temp"


dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=testdatafile,
    input_metadata_file=static_data,
    output_folder="/home/thoverga/Documents/VLINDER_github/metobs_toolkit",
)


dataset.import_data_from_file()
dataset.coarsen_time_resolution(freq="30T")


# offline mode
era = metobs_toolkit.Modeldata(
    metadf=dataset.metadf, extractor=metobs_toolkit.GeeExtractor()
)

era_datafile = os.path.join(str(lib_folder), "tests", "test_data", "era5_data.csv")


era.import_gee_data_from_csv(era_datafile)

assert era.df.shape[0] == 5348, "Something wrong with importing era data from csv."

#%%
output = dataset.fill_gaps_automatic(
    era, max_interpolate_duration_str="5H", overwrite_fill=True
)


checked = {
    "temp": {
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:00:00+0000", tz="UTC"),
        ): 15.760000000000002,
        ("vlinder01", pd.Timestamp("2022-10-06 17:30:00+0000", tz="UTC")): 15.22,
        ("vlinder01", pd.Timestamp("2022-10-06 18:00:00+0000", tz="UTC")): 14.68,
        ("vlinder01", pd.Timestamp("2022-10-06 18:30:00+0000", tz="UTC")): 14.14,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:00:00+0000", tz="UTC"),
        ): 13.600000000000001,
        ("vlinder01", pd.Timestamp("2022-10-06 19:30:00+0000", tz="UTC")): 13.06,
        ("vlinder01", pd.Timestamp("2022-10-06 20:00:00+0000", tz="UTC")): 12.52,
        ("vlinder01", pd.Timestamp("2022-10-06 20:30:00+0000", tz="UTC")): 11.98,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 21:00:00+0000", tz="UTC"),
        ): 11.440000000000001,
        ("vlinder01", pd.Timestamp("2022-10-04 02:00:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 02:30:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 03:00:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 03:30:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 04:00:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 04:30:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 05:00:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 05:30:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 06:00:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 06:30:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 07:00:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 07:30:00+0000", tz="UTC")): np.nan,
        ("vlinder01", pd.Timestamp("2022-10-04 08:00:00+0000", tz="UTC")): np.nan,
    },
    "temp_final_label": {
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:00:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:30:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 18:00:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 18:30:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:00:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:30:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 20:00:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 20:30:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 21:00:00+0000", tz="UTC"),
        ): "gap_interpolation",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 02:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 02:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 03:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 03:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 04:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 04:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 05:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 05:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 06:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 06:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 07:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 07:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-04 08:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
    },
}

checkeddf = pd.DataFrame(checked)


assert checkeddf.equals(output), "something wrong with the automatic gapfill"


#%%
# # Fill gaps using era5 data:
dataset.fill_gaps_using_debiased_modeldata(
    modeldata=era, obstype="temp", overwrite_fill=True
)


# validate

checked = {
    "temp": {
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:00:00+0000", tz="UTC"),
        ): 14.719558715820341,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:30:00+0000", tz="UTC"),
        ): 14.105651664733898,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 18:00:00+0000", tz="UTC"),
        ): 13.523644638061523,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 18:30:00+0000", tz="UTC"),
        ): 13.272471523284917,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:00:00+0000", tz="UTC"),
        ): 12.417074584960952,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:30:00+0000", tz="UTC"),
        ): 12.083017063140879,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 20:00:00+0000", tz="UTC"),
        ): 12.194173049926757,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 20:30:00+0000", tz="UTC"),
        ): 11.708401107788086,
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 21:00:00+0000", tz="UTC"),
        ): 11.562461853027344,
    },
    "temp_final_label": {
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 17:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 18:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 18:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 19:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 20:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 20:30:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
        (
            "vlinder01",
            pd.Timestamp("2022-10-06 21:00:00+0000", tz="UTC"),
        ): "gap_debiased_era5",
    },
}


checkeddf = pd.DataFrame(checked)
checkeddf = checkeddf.reset_index()
checkeddf = checkeddf.rename(columns={"level_0": "name", "level_1": "datetime"})

# checkeddf["datetime"] = pd.to_datetime(checkeddf["datetime"]).dt.tz_localize(
#     tz='UTC'
# )

checkeddf = checkeddf.set_index(["name", "datetime"])


test = dataset.gapfilldf[obstype].astype(float).eq(checkeddf[obstype])
if not test.all():
    print("Gapfill for era debias not equal to manual labels! Here is the difference")
    print(f"manual: {checkeddf}, tlk: {dataset.gapfilldf}")
    print(f"equal at: {test}")
    sys.exit("Error in debias gapfilling")

print("Done! ")
