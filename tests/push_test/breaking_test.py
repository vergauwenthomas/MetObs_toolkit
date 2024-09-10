#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# PROGRAM FOR TESTING THE BREAKING DATAFILE
"""
Created on Tue Nov 29 12:19:03 2022

@author: mivieijra
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

# %%


lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
print(str(lib_folder))
testdata = os.path.join(str(lib_folder), "tests", "test_data", "testdata_breaking.csv")


####### Create dataset ######
# % add template
template_file = os.path.join(
    str(lib_folder), "tests", "test_data", "template_breaking.json"
)

dataset_coarsened = metobs_toolkit.Dataset()
dataset_coarsened.update_file_paths(
    input_data_file=testdata, template_file=template_file
)


#####################################################################
# Set settings for QC
dupl_dropping = False  # method used to drop duplicated timestamps

persistence_time_window_to_check = "1h"  # Use this format as example: "1h20min50s"
min_num_obs = 3  # Minimum number of records in window to perform persistence check

max_valid_repetitions = 5  # Maximal number of repetitions that is allowed

min_value = -15.0  # Minimal allowed value
max_value = 29.0  # Maximal allowed value

max_increase_per_second = (
    8.0 / 3600.0
)  # Maximal allowed increase per second (for window variation check)
max_decrease_per_second = (
    10.0 / 3600.0
)  # Maximal allowed decrease per second (for window variation check)
time_window_to_check = "1h"  # Use this format as example: "1h20min50s"
min_window_members = 3  # Minimal number of records in window to perform check

max_increase_per_second_step = (
    8.0 / 3600.0
)  # Maximal allowed increase per second (for step check)
max_decrease_per_second_step = (
    -10.0 / 3600.0
)  # Maximal allowed increase per second (for step check)


# %%
dataset = metobs_toolkit.Dataset()
dataset.update_file_paths(input_data_file=testdata, template_file=template_file)
dataset.update_qc_settings(
    obstype="temp",
    dupl_timestamp_keep=dupl_dropping,
    persis_time_win_to_check=persistence_time_window_to_check,
    persis_min_num_obs=min_num_obs,
    rep_max_valid_repetitions=max_valid_repetitions,
    gross_value_min_value=min_value,
    win_var_max_increase_per_sec=max_increase_per_second,
    win_var_max_decrease_per_sec=max_decrease_per_second,
    win_var_time_win_to_check=time_window_to_check,
    win_var_min_num_obs=min_window_members,
    step_max_increase_per_sec=max_increase_per_second_step,
    step_max_decrease_per_sec=max_decrease_per_second_step,
)

dataset.import_data_from_file(origin_simplify_tolerance="1min")
dataset.apply_quality_control()


_ = dataset.get_qc_stats()

dataset.make_plot(stationnames=["Fictional"], colorby="label", show_outliers=True)

# %% Debug

combdf = dataset.get_full_status_df()
tlk_temp_df = dataset.get_full_status_df()["temp"]

# %% Copare with manual labels

man_df = pd.read_csv(testdata)
# create datetime coumn
man_df["datetime"] = man_df["date"] + " " + man_df["time"]
man_df["datetime"] = pd.to_datetime(man_df["datetime"])
man_df["datetime"] = man_df["datetime"].dt.tz_localize(dataset._get_tz())

# create double index
man_df["name"] = man_df["station"]
man_df = man_df.set_index(["name", "datetime"]).sort_index()

# format the label column
man_df = man_df.rename(columns={"qc_flags": "label_manual"})

# check if labels are the same

compare = man_df[["label_manual"]]

compare = man_df.merge(tlk_temp_df, how="left", left_index=True, right_index=True)

compare["correct_labeled"] = compare["label_manual"].eq(compare["label"])


diff = compare.loc[compare["correct_labeled"] == False, :]
assert diff.empty, f"Incorrect labels compared to the manual labels: \n {diff}"


# %% test if dataframes are equal to the solution


breaking_df_csv = "breaking_df.pkl"
breaking_outliers_csv = "breaking_outliersdf.pkl"
breaking_combined_df_csv = "breaking_combined_df.pkl"
breaking_gapsfill_df = "breaking_gapsfill_df.pkl"


def _create_qc_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.df.to_pickle(os.path.join(solution.solutions_dir, breaking_df_csv))
    dataset.outliersdf.to_pickle(
        os.path.join(solution.solutions_dir, breaking_outliers_csv)
    )
    dataset.get_full_status_df().to_pickle(
        os.path.join(solution.solutions_dir, breaking_combined_df_csv)
    )

    dataset.get_gaps_fill_df().to_pickle(
        os.path.join(solution.solutions_dir, breaking_gapsfill_df)
    )


# _create_qc_solutions()


def get_solutions():
    sol_df = pd.read_pickle(os.path.join(solution.solutions_dir, breaking_df_csv))
    sol_outliersdf = pd.read_pickle(
        os.path.join(solution.solutions_dir, breaking_outliers_csv)
    )
    sol_combdf = pd.read_pickle(
        os.path.join(solution.solutions_dir, breaking_combined_df_csv)
    )
    sol_gapsfilldf = pd.read_pickle(
        os.path.join(solution.solutions_dir, breaking_gapsfill_df)
    )
    return sol_df, sol_outliersdf, sol_combdf, sol_gapsfilldf


sol_df, sol_outliersdf, sol_combdf, sol_gapsfilldf = get_solutions()


# 1 Test the df attribute

solution.test_df_are_equal(testdf=dataset.df, solutiondf=sol_df)


# 2 Test the outliersdf attribute

solution.test_df_are_equal(testdf=dataset.outliersdf, solutiondf=sol_outliersdf)


# 3 Test the combined df
solution.test_df_are_equal(testdf=dataset.get_full_status_df(), solutiondf=sol_combdf)


# 4 Test the gapsfilldf
solution.test_df_are_equal(testdf=dataset.get_gaps_fill_df(), solutiondf=sol_gapsfilldf)


print("Done")
