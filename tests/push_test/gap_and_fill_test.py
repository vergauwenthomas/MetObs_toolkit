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

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
import solutions.solutions_creator as solution

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
# %% Basic tests on the gaps location and gaps functions


combdf_gaps_file = "gaps_test_combdf.pkl"
gapsdf_file = "gaps_gapdf.pkl"


def _create_gap_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.combine_all_to_obsspace().to_pickle(
        os.path.join(solution.solutions_dir, combdf_gaps_file)
    )
    dataset.get_gaps_fill_df().to_pickle(
        os.path.join(solution.solutions_dir, gapsdf_file)
    )


# _create_gap_solutions()


def get_unfilled_solutions():
    sol_combdf = pd.read_pickle(os.path.join(solution.solutions_dir, combdf_gaps_file))
    sol_gapsdf = pd.read_pickle(os.path.join(solution.solutions_dir, gapsdf_file))
    return sol_combdf, sol_gapsdf


sol_combdf, sol_gapsdf = get_unfilled_solutions()
# Check if gaps are well found and converted to dataframes


solution.test_df_are_equal(
    testdf=dataset.combine_all_to_obsspace(), solutiondf=sol_combdf
)

solution.test_df_are_equal(testdf=dataset.get_gaps_fill_df(), solutiondf=sol_gapsdf)


assert len(dataset.gaps) == 14, f"the number of gaps is not correct"

# %% Real vs unreal gaps


upsample_gapdf = "upsample_gapsdf.pkl"
downsample_gapdf = "downsample_gapsdf.pkl"


def _create_samplegap_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.coarsen_time_resolution(freq="30T")
    dataset.get_gaps_fill_df().to_pickle(
        os.path.join(solution.solutions_dir, upsample_gapdf)
    )
    dataset.coarsen_time_resolution(freq="4H")
    dataset.get_gaps_fill_df().to_pickle(
        os.path.join(solution.solutions_dir, downsample_gapdf)
    )


# _create_samplegap_solutions()


def get_sample_solution():
    upsampledf = pd.read_pickle(os.path.join(solution.solutions_dir, upsample_gapdf))
    downsampledf = pd.read_pickle(
        os.path.join(solution.solutions_dir, downsample_gapdf)
    )
    return upsampledf, downsampledf


upsample_sol, downsample_sol = get_sample_solution()


# to higher resolution
dataset.coarsen_time_resolution(
    freq="30T"
)  # since some gaps are smaller than 30 minuts this must create unreal gaps

solution.test_df_are_equal(testdf=dataset.get_gaps_fill_df(), solutiondf=upsample_sol)

assert len(dataset.gaps) == 14, f"the number of gaps is not correct after upsampling"

# to lower resolution
dataset.coarsen_time_resolution(
    freq="4H"
)  # since some gaps are smaller than 30 minuts this must create unreal gaps

solution.test_df_are_equal(testdf=dataset.get_gaps_fill_df(), solutiondf=downsample_sol)

assert len(dataset.gaps) == 14, f"the number of gaps is not correct after downsampling"

# %%
# =============================================================================
# Test interpolation on demo data
# =============================================================================

dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)
dataset.import_data_from_file()
dataset.coarsen_time_resolution(freq="1H")
dataset.apply_quality_control()
dataset.convert_outliers_to_gaps()

combdf_demo_startpoin_file = "gapfill_startpoint_combdf.pkkl"


def _create_startpoint_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.combine_all_to_obsspace().to_pickle(
        os.path.join(solution.solutions_dir, combdf_demo_startpoin_file)
    )


# _create_startpoint_solutions()
def get_startpoint_solution():
    return pd.read_pickle(
        os.path.join(solution.solutions_dir, combdf_demo_startpoin_file)
    )


solution.test_df_are_equal(
    testdf=dataset.combine_all_to_obsspace(), solutiondf=get_startpoint_solution()
)


# %% Test linear interpolation on missing obs
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


# %% Test if filled values are present in the combined df

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

# %% Test the update of outliers to gaps
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
era = metobs_toolkit.Modeldata("ERA5_hourly")

era_datafile = os.path.join(str(lib_folder), "tests", "test_data", "era5_data.csv")


era.set_model_from_csv(era_datafile)

assert era.df.shape[0] == 5348, "Something wrong with importing era data from csv."

# %%
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


# %%
# # Fill gaps using era5 data:
dataset.fill_gaps_era5(
    modeldata=era, method="debias", obstype="temp", overwrite_fill=True
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
