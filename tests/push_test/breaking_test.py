#!/usr/bin/env python3
# -*- coding: utf-8 -*-


# PROGRAM FOR TESTING THE BREAKING DATAFILE
"""
Created on Tue Nov 29 12:19:03 2022

@author: mivieijra
"""

import sys, os

from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
import solutions.solutions_creator as solution
import pandas as pd
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


# %% test if dataframes are equal to the solution


breaking_df_csv = "breaking_df.pkl"
breaking_outliers_csv = "breaking_outliersdf.pkl"
breaking_combined_df_csv = "breaking_combined_df.pkl"


def _create_qc_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.df.to_pickle(os.path.join(solution.solutions_dir, breaking_df_csv))
    dataset.outliersdf.to_pickle(
        os.path.join(solution.solutions_dir, breaking_outliers_csv)
    )
    dataset.combine_all_to_obsspace().to_pickle(
        os.path.join(solution.solutions_dir, breaking_combined_df_csv)
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
    return sol_df, sol_outliersdf, sol_combdf


sol_df, sol_outliersdf, sol_combdf = get_solutions()


# 1 Test the df attribute

solution.test_df_are_equal(testdf=dataset.df, solutiondf=sol_df)


# 2 Test the outliersdf attribute

solution.test_df_are_equal(testdf=dataset.outliersdf, solutiondf=sol_outliersdf)


# 3 Test the combined df
solution.test_df_are_equal(
    testdf=dataset.combine_all_to_obsspace(), solutiondf=sol_combdf
)
# %%


# all_manual_labels = list(man_df["flags"].unique())
# manual_to_tlkit_label_map = {
#     "ok": "ok",
#     "in step outlier group": dataset_coarsened.settings.qc["qc_checks_info"]["step"][
#         "outlier_flag"
#     ],
#     "repetitions outlier": dataset_coarsened.settings.qc["qc_checks_info"][
#         "repetitions"
#     ]["outlier_flag"],
#     "duplicated timestamp outlier": dataset_coarsened.settings.qc["qc_checks_info"][
#         "duplicated_timestamp"
#     ]["outlier_flag"],
#     "gross value outlier": dataset_coarsened.settings.qc["qc_checks_info"][
#         "gross_value"
#     ]["outlier_flag"],
#     "in window variation outlier group": dataset_coarsened.settings.qc[
#         "qc_checks_info"
#     ]["window_variation"]["outlier_flag"],
#     "persistance outlier": dataset_coarsened.settings.qc["qc_checks_info"][
#         "persistance"
#     ]["outlier_flag"],
# }

# # check if the mapper is still up to date
# assert all(
#     [True for label in all_manual_labels if label in manual_to_tlkit_label_map.keys()]
# ), "Update the manual to toolkit mapper"


# # =============================================================================
# # iterate over all labels and validate if the indices are equal between manual and toolkit
# # =============================================================================

# to_check = [
#     "ok",
#     "in step outlier group",
#     "repetitions outlier",
#     # 'duplicated timestamp outlier',
#     "gross value outlier",
#     "in window variation outlier group",
#     "persistance outlier",
# ]
# for man_label, tlk_label in manual_to_tlkit_label_map.items():
#     if not man_label in to_check:
#         continue

#     print(
#         f" Testing equality of the {tlk_label} with the manual labeling ({man_label})."
#     )

#     man_idx = man_df[man_df["flags"] == man_label].index.sort_values()
#     tlk_idx = (
#         tlk_df[tlk_df["label"] == tlk_label]
#         .xs("temp", level="obstype")
#         .index.sort_values()
#     )

#     if not tlk_idx.equals(man_idx):
#         print(f"ERROR: wrong labels for {tlk_label}")

#         print(f"differences tlkit --> manual: { tlk_idx.difference(man_idx)}")
#         print(f"differences manual --> tlkit: {man_idx.difference(tlk_idx)}")
#         sys.exit(1)

#     else:
#         print("OK!")


# # =============================================================================
# # test duplicates
# # =============================================================================
# # tested seperatly because duplicates are in tlk stored as one record, to avoid
# # duplicate index errors. So we have to do the same for the manual labeling
# man_label = "duplicated timestamp outlier"
# tlk_label = manual_to_tlkit_label_map[man_label]

# print(f" Testing equality of the {tlk_label} with the manual labeling ({man_label}).")

# man_df_no_duplic = man_df[~man_df.index.duplicated(keep="first")]

# man_idx = man_df_no_duplic[man_df_no_duplic["flags"] == man_label].index.sort_values()
# tlk_idx = (
#     tlk_df[tlk_df["label"] == tlk_label].xs("temp", level="obstype").index.sort_values()
# )

# if not tlk_idx.equals(man_idx):
#     print(f"ERROR: wrong labels for {tlk_label}")

#     print(f"differences tlkit --> manual: { tlk_idx.difference(man_idx)}")
#     print(f"differences manual --> tlkit: {man_idx.difference(tlk_idx)}")
#     sys.exit(1)

# else:
#     print("OK!")
# =============================================================================
# test missing Gaps
# =============================================================================


from datetime import datetime
import pandas as pd

# all_gaps = dataset.get_gaps_fill_df()
# temp_gaps = all_gaps.xs('temp', level='obstype').sort_index()


# # manual_gaps_creator
# man_gaps = []
# for gap in dataset.gaps:
#     if gap.obstype.name == 'temp':
#         man_gaps.append(tuple([gap.name, gap.startdt, gap.enddt]))


from pandas import Timestamp

man_gaps = [
    (
        "1",
        Timestamp("2020-09-16 21:30:00+0000", tz="UTC"),
        Timestamp("2020-09-16 21:30:00+0000", tz="UTC"),
    ),  # ok
    (
        "1",
        Timestamp("2020-09-16 22:30:00+0000", tz="UTC"),
        Timestamp("2020-09-16 22:30:00+0000", tz="UTC"),
    ),  # ok
    (
        "1",
        Timestamp("2020-09-16 23:00:00+0000", tz="UTC"),
        Timestamp("2020-09-16 23:45:00+0000", tz="UTC"),
    ),  # ok
    (
        "Fictional",
        Timestamp("2020-09-14 22:30:00+0000", tz="UTC"),
        Timestamp("2020-09-14 23:55:00+0000", tz="UTC"),
    ),  # ok
    (
        "Fictional",
        Timestamp("2020-09-15 01:15:00+0000", tz="UTC"),
        Timestamp("2020-09-15 01:15:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 02:50:00+0000", tz="UTC"),
        Timestamp("2020-09-15 02:50:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 03:00:00+0000", tz="UTC"),
        Timestamp("2020-09-15 03:20:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 04:15:00+0000", tz="UTC"),
        Timestamp("2020-09-15 04:15:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 05:10:00+0000", tz="UTC"),
        Timestamp("2020-09-15 05:15:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 06:00:00+0000", tz="UTC"),
        Timestamp("2020-09-15 06:00:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 06:45:00+0000", tz="UTC"),
        Timestamp("2020-09-15 06:45:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 07:35:00+0000", tz="UTC"),
        Timestamp("2020-09-15 07:40:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 08:35:00+0000", tz="UTC"),
        Timestamp("2020-09-15 08:35:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 09:40:00+0000", tz="UTC"),
        Timestamp("2020-09-15 09:40:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 10:55:00+0000", tz="UTC"),
        Timestamp("2020-09-15 10:55:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 12:40:00+0000", tz="UTC"),
        Timestamp("2020-09-15 12:40:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 14:55:00+0000", tz="UTC"),
        Timestamp("2020-09-15 14:55:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 17:25:00+0000", tz="UTC"),
        Timestamp("2020-09-15 17:25:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 17:35:00+0000", tz="UTC"),
        Timestamp("2020-09-15 17:35:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 19:35:00+0000", tz="UTC"),
        Timestamp("2020-09-15 19:35:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 20:40:00+0000", tz="UTC"),
        Timestamp("2020-09-15 20:40:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 21:50:00+0000", tz="UTC"),
        Timestamp("2020-09-15 21:50:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 22:50:00+0000", tz="UTC"),
        Timestamp("2020-09-15 22:50:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 23:25:00+0000", tz="UTC"),
        Timestamp("2020-09-15 23:25:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-15 23:50:00+0000", tz="UTC"),
        Timestamp("2020-09-16 00:00:00+0000", tz="UTC"),
    ),
    (
        "Fictional",
        Timestamp("2020-09-16 00:45:00+0000", tz="UTC"),
        Timestamp("2020-09-16 00:45:00+0000", tz="UTC"),
    ),
]


# create a list of found temperature gaps
found_gaps = [
    tuple([gap.name, gap.startdt, gap.enddt])
    for gap in dataset.gaps
    if gap.obstype.name == "temp"
]
# Check if all found gaps are also manual gap
gaps_are_ok = True
for gap in found_gaps:
    if gap not in man_gaps:
        print(f"{gap} Not found in the manual gaps!")
        gaps_are_ok = False

for gap in man_gaps:
    if gap not in found_gaps:
        print(f"{gap} is a manual gap but not found by the toolkit")
        gaps_are_ok = False


assert gaps_are_ok, "problem with the gaps"
print("Gaps are checked")


# manual_missing_gaps = [
#     {
#         "name": "Fictional",
#         "start_gap": datetime(2020, 9, 14, 22, 30),
#         "end_gap": datetime(2020, 9, 14, 23, 55),
#     }
# ]  # UPDATE MANUALLY !!!!!!!!!!

# print("Testing the gaps")

# man_gapsdf = pd.DataFrame().from_records(manual_missing_gaps)
# man_gapsdf = man_gapsdf.set_index("name")

# # Set timezone
# man_gapsdf["start_gap"] = pd.to_datetime(man_gapsdf["start_gap"]).dt.tz_localize(
#     tz="UTC"
# )
# man_gapsdf["end_gap"] = pd.to_datetime(man_gapsdf["end_gap"]).dt.tz_localize(tz="UTC")


# %%

# tlk_gapsdf = dataset.get_gaps_df()
# tlk_gapsdf = tlk_gapsdf[list(man_gapsdf.columns)]


# if not tlk_gapsdf.equals(man_gapsdf):
#     print(f"ERROR: wrong gaps detection")

#     print(
#         f"differences tlkit --> manual: {tlk_gapsdf[~tlk_gapsdf.apply(tuple,1).isin(man_gapsdf.apply(tuple,1))]}"
#     )
#     print(
#         f"differences manual --> tlkit: {man_gapsdf[~man_gapsdf.apply(tuple,1).isin(tlk_gapsdf.apply(tuple,1))]}"
#     )
#     sys.exit(1)

# else:
#     print("OK!")

print("Done")
