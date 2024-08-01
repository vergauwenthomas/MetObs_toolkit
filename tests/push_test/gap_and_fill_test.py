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
    str(lib_folder), "tests", "test_data", "testdata_with_gaps.csv"
)

static_data = os.path.join(str(lib_folder), "static_data", "vlinder_metadata.csv")

templatefile = os.path.join(
    str(lib_folder), "tests", "test_data", "template_testdata_with_gaps.json"
)


# #% import data
dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=testdatafile,
    input_metadata_file=static_data,
    template_file=templatefile,
)
dataset.import_data_from_file()
dataset.coarsen_time_resolution()


# %% Basic tests on the gaps location and gaps functions


combdf_gaps_file = "gaps_test_combdf.pkl"
gapsdf_file = "gaps_gapdf.pkl"


def _create_gap_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.get_full_status_df().to_pickle(
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

_ = dataset.get_full_status_df()


diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=sol_combdf
)
assert diff_df is None

diff_df = solution.test_df_are_equal(
    testdf=dataset.get_gaps_fill_df(), solutiondf=sol_gapsdf
)
assert diff_df is None


# %% Test gap relate methods (execution tests)

_ = dataset.get_gaps_fill_df()

_ = dataset.gaps[1].get_info()


# %%
# =============================================================================
# Gaps on station subsetting
# =============================================================================

station_combdf_gaps_file = "station_gapsdf.pkl"


def _create_station_gap_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.get_station("vlinder01").get_gaps_fill_df().to_pickle(
        os.path.join(solution.solutions_dir, station_combdf_gaps_file)
    )


# _create_station_gap_solutions()


def get_unfilled_station_solutions():
    sol_sta_gapsdf = pd.read_pickle(
        os.path.join(solution.solutions_dir, station_combdf_gaps_file)
    )
    return sol_sta_gapsdf


sol_sta_gapsdf = get_unfilled_station_solutions()
# Check if gaps are well found and converted to dataframes
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_station("vlinder01").get_gaps_fill_df(),
    solutiondf=sol_sta_gapsdf,
)


# %%
# =============================================================================
# Test interpolation
# =============================================================================

dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=testdatafile,
    input_metadata_file=static_data,
    template_file=templatefile,
)

dataset.import_data_from_file()
dataset.coarsen_time_resolution(freq="1H")
dataset.apply_quality_control()
dataset.convert_outliers_to_gaps()


dataset.get_qc_stats()
combdf_demo_startpoin_file = "gapfill_startpoint_combdf.pkl"


def _create_startpoint_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.get_full_status_df().to_pickle(
        os.path.join(solution.solutions_dir, combdf_demo_startpoin_file)
    )


# _create_startpoint_solutions()
def get_startpoint_solution():
    return pd.read_pickle(
        os.path.join(solution.solutions_dir, combdf_demo_startpoin_file)
    )


diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_startpoint_solution()
)
assert diff_df is None


inter_sol_df_file = "interp_sol_df_file.pkl"
station_inter_sol_df = "interp_sol_sta_df_file.pkl"


def _create_interpolation_solution():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.interpolate_gaps(
        obstype="temp",
        overwrite_fill=True,
        method="time",
        max_consec_fill=10,
        max_lead_to_gap_distance="2H",
        max_trail_to_gap_distance=None,
    )

    trg_file = os.path.join(solution.solutions_dir, inter_sol_df_file)
    dataset.get_full_status_df().to_pickle(trg_file)

    sta = dataset.get_station("vlinder04")
    sta.interpolate_gaps(
        obstype="temp",
        overwrite_fill=True,
        method="nearest",
        max_consec_fill=4,
        max_lead_to_gap_distance="2H",
        max_trail_to_gap_distance=None,
    )
    trg_sta_file = os.path.join(solution.solutions_dir, station_inter_sol_df)
    sta.get_full_status_df().to_pickle(trg_sta_file)


# _create_interpolation_solution()


def get_interp_solution():
    dataset_comb_solution = pd.read_pickle(
        os.path.join(solution.solutions_dir, inter_sol_df_file)
    )
    sta_comb_solution = pd.read_pickle(
        os.path.join(solution.solutions_dir, station_inter_sol_df)
    )
    return dataset_comb_solution, sta_comb_solution


dataset_sol, sta_solution = get_interp_solution()

# interpolate the gaps
dataset.interpolate_gaps(
    obstype="temp",
    overwrite_fill=True,
    method="time",
    max_consec_fill=10,
    max_lead_to_gap_distance="2H",
    max_trail_to_gap_distance=None,
)


# test interpolation on dataset level
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=dataset_sol
)
assert diff_df is None


sta = dataset.get_station("vlinder04")
sta.interpolate_gaps(
    obstype="temp",
    overwrite_fill=True,
    method="linear",
    max_consec_fill=4,
    max_lead_to_gap_distance="2H",
    max_trail_to_gap_distance=None,
)
# check if this methods overwrites the linear fille
sta.interpolate_gaps(
    obstype="temp",
    overwrite_fill=True,
    method="nearest",
    max_consec_fill=4,
    max_lead_to_gap_distance="2H",
    max_trail_to_gap_distance=None,
)


# test interpolation on station level
diff_df = solution.test_df_are_equal(
    testdf=sta.get_full_status_df(), solutiondf=sta_solution
)
assert diff_df is None

# %% Test higher order interpolation


highorder_inter_sol_df_file = "highorderinterp_sol_df_file.pkl"


def _create_polynom_interpolation_solution():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    sta = dataset.get_station("vlinder04")
    sta.interpolate_gaps(
        obstype="temp",
        overwrite_fill=True,
        method="polynomial",
        max_consec_fill=5,
        n_trailing_anchors=4,
        n_leading_anchors=3,
        max_lead_to_gap_distance="4h",
        max_trail_to_gap_distance=None,
        method_kwargs={"order": 3},
    )

    trg_sta_file = os.path.join(solution.solutions_dir, highorder_inter_sol_df_file)
    sta.get_gaps_fill_df().to_pickle(trg_sta_file)


# _create_polynom_interpolation_solution()


def get_polynom_interp_solution():
    sol = pd.read_pickle(
        os.path.join(solution.solutions_dir, highorder_inter_sol_df_file)
    )

    return sol


sta = dataset.get_station("vlinder04")
sta.interpolate_gaps(
    obstype="temp",
    overwrite_fill=True,
    method="polynomial",
    max_consec_fill=5,
    n_trailing_anchors=4,
    n_leading_anchors=3,
    max_lead_to_gap_distance="4h",
    max_trail_to_gap_distance=None,
    method_kwargs={"order": 3},
)


# test interpolation on station level
diff_df = solution.test_df_are_equal(
    testdf=sta.get_gaps_fill_df(), solutiondf=get_polynom_interp_solution()
)
assert diff_df is None


sta.make_plot(colorby="label")


# %% Get modeldata and test model gapfill methods


modeldatafile = "era_modeldata_for_gapfill_dataset.pkl"


def _create_modeldata():

    era = dataset.get_modeldata()
    era.save_modeldata(outputfolder=solution.solutions_dir, filename=modeldatafile)


# _create_modeldata()


def get_modeldata():
    era = metobs_toolkit.Modeldata("ERA5_hourly")
    era = era.import_modeldata(
        folder_path=solution.solutions_dir, filename=modeldatafile
    )
    return era


era = get_modeldata()


# %% Test regular model gapfill


raw_gapfill_sol_file = "raw_gapfill_solution.pkl"


def _create_raw_gapfil_solution():
    print("SOLUTION WILL BE OVERWRITTEN !!!!!")
    trg_path = os.path.join(solution.solutions_dir, raw_gapfill_sol_file)
    dataset.fill_gaps_with_raw_modeldata(
        Modeldata=era, obstype="temp", overwrite_fill=True
    )
    dataset.get_full_status_df().to_pickle(trg_path)


# _create_raw_gapfil_solution()


def get_raw_gapfill_sol():
    trg_path = os.path.join(solution.solutions_dir, raw_gapfill_sol_file)
    return pd.read_pickle(trg_path)


dataset.fill_gaps_with_raw_modeldata(Modeldata=era, obstype="temp", overwrite_fill=True)
# test raw gapfill on dataset
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_raw_gapfill_sol()
)
assert diff_df is None


# %% Test debias gapfill

debias_gapfill_sol_file = "debias_gapfill_solution.pkl"


def _create_debias_gapfil_solution():
    print("SOLUTION WILL BE OVERWRITTEN !!!!!")
    trg_path = os.path.join(solution.solutions_dir, debias_gapfill_sol_file)
    dataset.fill_gaps_with_debiased_modeldata(
        Modeldata=era,
        obstype="temp",
        overwrite_fill=True,
        leading_period_duration="15h",
        min_leading_records_total=12,
        trailing_period_duration="3h",
        min_trailing_records_total=3,
    )
    dataset.get_full_status_df().to_pickle(trg_path)


# _create_debias_gapfil_solution()


def get_debias_gapfill_sol():
    trg_path = os.path.join(solution.solutions_dir, debias_gapfill_sol_file)
    return pd.read_pickle(trg_path)


dataset.fill_gaps_with_debiased_modeldata(
    Modeldata=era,
    obstype="temp",
    overwrite_fill=True,
    leading_period_duration="15h",
    min_leading_records_total=12,
    trailing_period_duration="3h",
    min_trailing_records_total=3,
)


# test raw gapfill on dataset
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_debias_gapfill_sol()
)
assert diff_df is None


# %% Test diurnal debias gapfill

diur_debias_gapfill_sol_file = "diurnal_debias_gapfill_solution.pkl"


def _create_diur_debias_gapfil_solution():
    print("SOLUTION WILL BE OVERWRITTEN !!!!!")
    trg_path = os.path.join(solution.solutions_dir, diur_debias_gapfill_sol_file)
    dataset.fill_gaps_with_diurnal_debiased_modeldata(
        Modeldata=era,
        obstype="temp",
        overwrite_fill=True,
        leading_period_duration="48h",
        min_debias_sample_size=2,
        trailing_period_duration="12h",
    )
    dataset.get_full_status_df().to_pickle(trg_path)


# _create_diur_debias_gapfil_solution()


def get_diur_debias_gapfill_sol():
    trg_path = os.path.join(solution.solutions_dir, diur_debias_gapfill_sol_file)
    return pd.read_pickle(trg_path)


dataset.fill_gaps_with_diurnal_debiased_modeldata(
    Modeldata=era,
    obstype="temp",
    overwrite_fill=True,
    leading_period_duration="48h",
    min_debias_sample_size=2,
    trailing_period_duration="12h",
)


# test raw gapfill on dataset
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_diur_debias_gapfill_sol()
)
assert diff_df is None

# %% Test weighted diurnal debias gapfill

weight_diur_debias_gapfill_sol_file = "weight_diurnal_debias_gapfill_solution.pkl"
weight_diur_debias_gapfilldf_sol = "weight_diurnal_debias_gapfilldf_sol.pkl"


def _create_weight_diur_debias_gapfil_solution():
    print("SOLUTION WILL BE OVERWRITTEN !!!!!")
    trg_path = os.path.join(solution.solutions_dir, weight_diur_debias_gapfill_sol_file)
    dataset.fill_gaps_with_weighted_diurnal_debias_modeldata(
        Modeldata=era,
        obstype="temp",
        overwrite_fill=True,
        leading_period_duration="58h",
        min_lead_debias_sample_size=1,
        trailing_period_duration="24h",
        min_trail_debias_sample_size=1,
    )
    dataset.get_full_status_df().to_pickle(trg_path)
    dataset.get_gaps_fill_df().to_pickle(
        os.path.join(solution.solutions_dir, weight_diur_debias_gapfilldf_sol)
    )


# _create_weight_diur_debias_gapfil_solution()


def get_weight_diur_debias_gapfill_sol():
    fulldf = pd.read_pickle(
        os.path.join(solution.solutions_dir, weight_diur_debias_gapfill_sol_file)
    )
    gapfilldf = pd.read_pickle(
        os.path.join(solution.solutions_dir, weight_diur_debias_gapfilldf_sol)
    )
    return fulldf, gapfilldf


dataset.fill_gaps_with_weighted_diurnal_debias_modeldata(
    Modeldata=era,
    obstype="temp",
    overwrite_fill=True,
    leading_period_duration="58H",
    min_lead_debias_sample_size=1,
    trailing_period_duration="24H",
    min_trail_debias_sample_size=1,
)


fulldf_sol, gapfilldf_sol = get_weight_diur_debias_gapfill_sol()

# test raw gapfill on dataset
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(),
    solutiondf=fulldf_sol,
)
assert diff_df is None


# test raw gapfill on dataset
diff_df = solution.test_df_are_equal(
    testdf=dataset.get_gaps_fill_df(),
    solutiondf=gapfilldf_sol,
)
assert diff_df is None


# %% Find gap

import datetime

filledgap = dataset.find_gap(
    stationname="vlinder01",
    obstype="temp",
    in_gap_timestamp=datetime.datetime(2022, 10, 6, 18, 53, 2),
)

assert filledgap is not None, "find gap does not work accordingly"

filledgap.get_info()


nogap = dataset.find_gap(
    stationname="vlinder01",
    obstype="temp",
    in_gap_timestamp=datetime.datetime(2022, 10, 6, 13, 53, 2),
)


assert nogap is None, "find gap does not work accordingly"
