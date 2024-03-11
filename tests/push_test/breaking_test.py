#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# PROGRAM FOR TESTING THE BREAKING DATAFILE
"""
Created on Tue Nov 29 12:19:03 2022

@author: mivieijra
"""

import sys, os

from pathlib import Path

import metobs_toolkit

lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
print(str(lib_folder))


testdata = os.path.join(str(lib_folder), "tests", "test_data", "testdata_breaking.csv")


####### Create dataset ######
# % add template
template_file = os.path.join(
    str(lib_folder), "tests", "test_data", "template_breaking.csv"
)

dataset_coarsened = metobs_toolkit.Dataset()
dataset_coarsened.update_settings(input_data_file=testdata, template_file=template_file)


#####################################################################
# Set settings for QC
dupl_dropping = False  # method used to drop duplicated timestamps

persistance_time_window_to_check = "1h"  # Use this format as example: "1h20min50s"
min_num_obs = 3  # Minimum number of records in window to perform persistance check

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

dataset_coarsened.update_qc_settings(
    obstype="temp",
    dupl_timestamp_keep=dupl_dropping,
    persis_time_win_to_check=persistance_time_window_to_check,
    persis_min_num_obs=min_num_obs,
    rep_max_valid_repetitions=max_valid_repetitions,
    gross_value_min_value=min_value,
    gross_value_max_value=max_value,
    win_var_max_increase_per_sec=max_increase_per_second,
    win_var_max_decrease_per_sec=max_decrease_per_second,
    win_var_time_win_to_check=time_window_to_check,
    win_var_min_num_obs=min_window_members,
    step_max_increase_per_sec=max_increase_per_second_step,
    step_max_decrease_per_sec=max_decrease_per_second_step,
)


dataset_coarsened.import_data_from_file()


# %%
dataset_coarsened.coarsen_time_resolution()
dataset_coarsened.apply_quality_control()

# _ = dataset_coarsened.get_qc_stats()
# %%
dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdata, template_file=template_file)
dataset.update_qc_settings(
    obstype="temp",
    dupl_timestamp_keep=dupl_dropping,
    persis_time_win_to_check=persistance_time_window_to_check,
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

dataset.import_data_from_file()
dataset.apply_quality_control()


_ = dataset.get_qc_stats()

dataset.make_plot(stationnames=["Fictional"], colorby="label", show_outliers=True)

# %% Debug

combdf = dataset.combine_all_to_obsspace()


# %% Compare manual and toolkit labeling


man_df = dataset.input_df  # manual label

tlk_df = dataset.combine_all_to_obsspace()


tlk_temp_df = tlk_df.xs("temp", level="obstype").sort_index()

# %%


print(tlk_temp_df.iloc[35:39])


# %%
all_manual_labels = list(man_df["flags"].unique())
manual_to_tlkit_label_map = {
    "ok": "ok",
    "in step outlier group": dataset_coarsened.settings.qc["qc_checks_info"]["step"][
        "outlier_flag"
    ],
    "repetitions outlier": dataset_coarsened.settings.qc["qc_checks_info"][
        "repetitions"
    ]["outlier_flag"],
    "duplicated timestamp outlier": dataset_coarsened.settings.qc["qc_checks_info"][
        "duplicated_timestamp"
    ]["outlier_flag"],
    "gross value outlier": dataset_coarsened.settings.qc["qc_checks_info"][
        "gross_value"
    ]["outlier_flag"],
    "in window variation outlier group": dataset_coarsened.settings.qc[
        "qc_checks_info"
    ]["window_variation"]["outlier_flag"],
    "persistance outlier": dataset_coarsened.settings.qc["qc_checks_info"][
        "persistance"
    ]["outlier_flag"],
}

# check if the mapper is still up to date
assert all(
    [True for label in all_manual_labels if label in manual_to_tlkit_label_map.keys()]
), "Update the manual to toolkit mapper"


# =============================================================================
# iterate over all labels and validate if the indices are equal between manual and toolkit
# =============================================================================

to_check = [
    "ok",
    "in step outlier group",
    "repetitions outlier",
    # 'duplicated timestamp outlier',
    "gross value outlier",
    "in window variation outlier group",
    "persistance outlier",
]
for man_label, tlk_label in manual_to_tlkit_label_map.items():
    if not man_label in to_check:
        continue

    print(
        f" Testing equality of the {tlk_label} with the manual labeling ({man_label})."
    )

    man_idx = man_df[man_df["flags"] == man_label].index.sort_values()
    tlk_idx = (
        tlk_df[tlk_df["label"] == tlk_label]
        .xs("temp", level="obstype")
        .index.sort_values()
    )

    if not tlk_idx.equals(man_idx):
        print(f"ERROR: wrong labels for {tlk_label}")

        print(f"differences tlkit --> manual: { tlk_idx.difference(man_idx)}")
        print(f"differences manual --> tlkit: {man_idx.difference(tlk_idx)}")
        sys.exit(1)

    else:
        print("OK!")


# =============================================================================
# test duplicates
# =============================================================================
# tested seperatly because duplicates are in tlk stored as one record, to avoid
# duplicate index errors. So we have to do the same for the manual labeling
man_label = "duplicated timestamp outlier"
tlk_label = manual_to_tlkit_label_map[man_label]

print(f" Testing equality of the {tlk_label} with the manual labeling ({man_label}).")

man_df_no_duplic = man_df[~man_df.index.duplicated(keep="first")]

man_idx = man_df_no_duplic[man_df_no_duplic["flags"] == man_label].index.sort_values()
tlk_idx = (
    tlk_df[tlk_df["label"] == tlk_label].xs("temp", level="obstype").index.sort_values()
)

if not tlk_idx.equals(man_idx):
    print(f"ERROR: wrong labels for {tlk_label}")

    print(f"differences tlkit --> manual: { tlk_idx.difference(man_idx)}")
    print(f"differences manual --> tlkit: {man_idx.difference(tlk_idx)}")
    sys.exit(1)

else:
    print("OK!")
# =============================================================================
# test missing Gaps
# =============================================================================


from datetime import datetime
import pandas as pd


manual_missing_gaps = [
    {
        "name": "Fictional",
        "start_gap": datetime(2020, 9, 14, 22, 30),
        "end_gap": datetime(2020, 9, 14, 23, 55),
    }
]  # UPDATE MANUALLY !!!!!!!!!!

print("Testing the gaps")

man_gapsdf = pd.DataFrame().from_records(manual_missing_gaps)
man_gapsdf = man_gapsdf.set_index("name")

# Set timezone
man_gapsdf["start_gap"] = pd.to_datetime(man_gapsdf["start_gap"]).dt.tz_localize(
    tz="UTC"
)
man_gapsdf["end_gap"] = pd.to_datetime(man_gapsdf["end_gap"]).dt.tz_localize(tz="UTC")


# %%

tlk_gapsdf = dataset.get_gaps_df()
tlk_gapsdf = tlk_gapsdf[list(man_gapsdf.columns)]


if not tlk_gapsdf.equals(man_gapsdf):
    print(f"ERROR: wrong gaps detection")

    print(
        f"differences tlkit --> manual: {tlk_gapsdf[~tlk_gapsdf.apply(tuple,1).isin(man_gapsdf.apply(tuple,1))]}"
    )
    print(
        f"differences manual --> tlkit: {man_gapsdf[~man_gapsdf.apply(tuple,1).isin(tlk_gapsdf.apply(tuple,1))]}"
    )
    sys.exit(1)

else:
    print("OK!")


# =============================================================================
# test missing missing timestamps
# =============================================================================


# This has to be done properly, there are too much missing timestamps for manual labelling.
# as a dirty fix, now only count the number of missing timestamps and see it that is oke

number_missing_timestamps = {"1": 1, "Fictional": 307}


print("Testing the missing obs")

man_missing_timestamps_df = pd.DataFrame(
    [
        ("Fictional", datetime(2020, 9, 15, 2, 50)),
        ("Fictional", datetime(2020, 9, 15, 3, 00)),
        ("Fictional", datetime(2020, 9, 15, 3, 5)),
        ("Fictional", datetime(2020, 9, 15, 3, 10)),
        ("Fictional", datetime(2020, 9, 15, 3, 15)),
        ("Fictional", datetime(2020, 9, 15, 3, 20)),
        ("Fictional", datetime(2020, 9, 15, 5, 10)),
        ("Fictional", datetime(2020, 9, 15, 6, 45)),
        ("Fictional", datetime(2020, 9, 15, 7, 40)),
        ("Fictional", datetime(2020, 9, 15, 12, 40)),
        ("Fictional", datetime(2020, 9, 15, 17, 35)),
        ("Fictional", datetime(2020, 9, 15, 20, 40)),
        ("Fictional", datetime(2020, 9, 15, 23, 50)),
        ("1", datetime(2020, 9, 16, 21, 30)),
    ],
    columns=["name", "datetime"],
)

# set timezone
man_missing_timestamps_df["datetime"] = pd.to_datetime(
    man_missing_timestamps_df["datetime"]
).dt.tz_localize(tz="UTC")

# get index
man_missing_timestamps_idx = pd.MultiIndex.from_frame(
    man_missing_timestamps_df
).sort_values()

tlk_missing_series = dataset.missing_obs.series
tlk_missing_series = tlk_missing_series.reset_index().rename(
    columns={"index": "name", 0: "datetime"}
)
tlk_missing_timestamps_idx = pd.MultiIndex.from_frame(tlk_missing_series).sort_values()


if not tlk_missing_timestamps_idx.equals(man_missing_timestamps_idx):
    print(f"ERROR: wrong missing timestamps detection")

    print(
        f"differences tlkit --> manual: { tlk_missing_timestamps_idx.difference(man_missing_timestamps_idx)}"
    )
    print(
        f"differences manual --> tlkit: {man_missing_timestamps_idx.difference(tlk_missing_timestamps_idx)}"
    )
    sys.exit(1)

else:
    print("OK!")
