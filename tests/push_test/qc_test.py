#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:29:37 2022

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


# %% IO testdata


dataset = metobs_toolkit.Dataset()

dataset.import_data_from_file(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)

dataset.resample(target_freq="1h")


# %% Apply Qc on dataset level
# general_qc_sol = "general_qc.pkl"


# def _create_general_qc_solution():
#     print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
#     dataset.apply_quality_control(obstype="temp")

#     dataset.get_full_status_df().to_pickle(
#         os.path.join(solution.solutions_dir, general_qc_sol)
#     )


# # _create_general_qc_solution()


# def get_general_qc_sol():
#     return pd.read_pickle(os.path.join(solution.solutions_dir, general_qc_sol))


# dataset.apply_quality_control(obstype="temp")

# diff_df = solution.test_df_are_equal(
#     testdf=dataset.get_full_status_df(), solutiondf=get_general_qc_sol()
# )
# assert diff_df is None


# dataset.get_qc_stats(make_plot=False)
# dataset.get_qc_stats(obstype="humidity", make_plot=False)


# %% Apply buddy check


buddy_qc_sol = "buddy_qc.pkl"


def _create_buddy_qc_solution(dataset):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.update_qc_settings(
        buddy_radius=17000,
        buddy_min_sample_size=3,
        buddy_max_elev_diff=150,
        buddy_min_std=1.2,
        buddy_threshold=2.4,
        buddy_elev_gradient=None,
    )

    dataset.apply_buddy_check(use_constant_altitude=True)

    dataset.get_full_status_df().to_pickle(
        os.path.join(solution.solutions_dir, buddy_qc_sol)
    )


# _create_buddy_qc_solution(dataset)


def get_buddy_qc_sol():
    return pd.read_pickle(os.path.join(solution.solutions_dir, buddy_qc_sol))


import copy

datset1 = copy.deepcopy(dataset)
datset2 = copy.deepcopy(dataset)

datset1.buddy_check(
    target_obstype="temp",
    buddy_radius=25000,
    min_sample_size=3,
    max_alt_diff=None,
    min_std=1.0,
    std_threshold=2.1,
    N_iter=1,
    instantanious_tolerance=pd.Timedelta("4min"),
    lapserate=None,  # -0.0065
    use_mp=False,
)


datset2.buddy_check(
    target_obstype="temp",
    buddy_radius=25000,
    min_sample_size=3,
    max_alt_diff=None,
    min_std=1.0,
    std_threshold=2.1,
    N_iter=2,
    instantanious_tolerance=pd.Timedelta("4min"),
    lapserate=None,  # -0.0065
    use_mp=False,
)


diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_buddy_qc_sol()
)
# assert diff_df is None


# %% Apply Qc on obstype not specified in settings


dataset.apply_quality_control(obstype="humidity")
dataset.get_qc_stats(obstype="humidity", make_plot=False)


# %% Apply QC on station level

sta = dataset.get_station("vlinder05")

sta.apply_quality_control()
sta.apply_buddy_check(use_constant_altitude=True)
test = sta.get_qc_stats(make_plot=True)


# %% Apply titan checks


titan_buddy_qc_sol = "titan_buddy_qc.pkl"


def _create_titan_buddy_qc_solution(dataset):
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")
    dataset.update_titan_qc_settings(
        obstype="temp",
        buddy_radius=50000,
        buddy_num_min=3,
        buddy_max_elev_diff=200,
        buddy_threshold=2,
    )

    dataset.apply_titan_buddy_check(use_constant_altitude=True)
    dataset.get_full_status_df().to_pickle(
        os.path.join(solution.solutions_dir, titan_buddy_qc_sol)
    )


# _create_titan_buddy_qc_solution(dataset)


def get_titan_buddy_qc_sol():
    return pd.read_pickle(os.path.join(solution.solutions_dir, titan_buddy_qc_sol))


#  ------ Buddy check --------------
dataset.update_titan_qc_settings(
    obstype="temp",
    buddy_radius=50000,
    buddy_num_min=3,
    buddy_max_elev_diff=200,
    buddy_threshold=2,
)


dataset.apply_titan_buddy_check(use_constant_altitude=True)

diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_titan_buddy_qc_sol()
)
assert diff_df is None


# count test
assert (
    dataset.outliersdf["label"].value_counts()["titan buddy check outlier"] == 277
), "The TITAN buddy check did not perfom good."

# test if a check does not overwrite itself
dataset.update_titan_qc_settings(
    obstype="temp",
    buddy_radius=80000,
    buddy_num_min=3,
    buddy_max_elev_diff=200,
    buddy_threshold=0.5,
)
dataset.apply_titan_buddy_check(use_constant_altitude=True)

diff_df = solution.test_df_are_equal(
    testdf=dataset.get_full_status_df(), solutiondf=get_titan_buddy_qc_sol()
)
assert diff_df is None
# %%


# import numpy as np


# #  ------------- SCT check ---------------
# dataset.update_titan_qc_settings(obstype='temp',
#                                 sct_basic=True,
#                                 sct_eps2=0.5,
#                                 sct_inner_radius=20000,
#                                 sct_kth_closest_obs_horizontal_scale=2,
#                                 sct_max_horizontal_scale=100000,
#                                 sct_maxa_deviation=7,
#                                 sct_maxv_deviation=0.5,
#                                 sct_min_elev_diff=100,
#                                 sct_min_horizontal_scale=250,
#                                 sct_mina_deviation=7,
#                                 sct_minv_deviation=0.7,
#                                 sct_num_iterations=2,
#                                 sct_num_max_outer=10,
#                                 sct_num_min_outer=3,
#                                 sct_num_min_prof=1,
#                                 sct_outer_radius=50000,
#                                 sct_tneg=16,
#                                 sct_tpos=16,
#                                 sct_vertical_scale=200,
#                                 sct_debug=False,
#                                 )


# dataset.metadf['altitude'] = np.random.randint(1, 250, dataset.metadf.shape[0])


# dataset.apply_titan_sct_resistant_check()
