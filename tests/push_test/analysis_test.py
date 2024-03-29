#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  8 15:55:13 2023

@author: thoverga
"""


import sys, os

from pathlib import Path
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
print(sys.path)
import solutions.solutions_creator as solution

import metobs_toolkit


lib_folder = Path(__file__).resolve().parents[2]
# sys.path.append(str(lib_folder))
# print(str(lib_folder))


# %% Create startpoint
dataset = metobs_toolkit.Dataset()
dataset.update_settings(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)
dataset.import_data_from_file()
dataset.coarsen_time_resolution(freq="1H")
dataset.apply_quality_control()
# get LCZ
# dataset.get_lcz()
# Outliers as gaps
dataset.convert_outliers_to_gaps()
# interpolate the gaps
dataset.interpolate_gaps(max_consec_fill=20)

dataset.metadf = solution.get_all_metadata_solution()


# %% Creation Analysis

combinedf = dataset.combine_all_to_obsspace()

num_of_records = (
    combinedf["label"].value_counts()["ok"]
    + combinedf["label"].value_counts()["interpolation"]
)

an = metobs_toolkit.Analysis(dataset)
an.convert_filled_to_obs()

an_df_solution = "analysis_solution.pkl"


def create_df_attr_solution():
    print("SOLUTION WILL BE OVERWRITTEN!! ")
    an.df.to_pickle(os.path.join(solution.solutions_dir, an_df_solution))
    return None


# create_df_attr_solution()

# Test equality
solution.test_df_are_equal(
    testdf=an.df,
    solutiondf=pd.read_pickle(os.path.join(solution.solutions_dir, an_df_solution)),
)


# %%
# =============================================================================
# test diurnal methods
# =============================================================================

teststa = ["vlinder01", "vlinder02", "vlinder03"]

from datetime import datetime

startdt = datetime(2022, 9, 4)

# Test plotting and functions
temp_diurnal = an.get_diurnal_statistics(
    colorby="lcz", stations=teststa, startdt=startdt
)


test2 = an.get_diurnal_statistics_with_reference(
    refstation="vlinder08", colorby="name", errorbands=True
)

test3 = an.get_aggregated_cycle_statistics(aggregation=["lcz"])

# =============================================================================
# test anual cycle
# =============================================================================

# test4 = an.get_anual_statistics(groupby=['name'])


diurnal_solution_df_file = "diurnal_demo_df_solution.pkl"


def create_diurnal_cycle_solution():
    print("SOLUTION WILL BE OVERWRITTEN!! ")
    temp_diurnal.to_pickle(
        os.path.join(solution.solutions_dir, diurnal_solution_df_file)
    )
    return None


# create_diurnal_cycle_solution() #RUN ONLY IF MANUALLY CHECKED !!!!


def get_diurnal_cycle_solution():
    return pd.read_pickle(
        os.path.join(solution.solutions_dir, diurnal_solution_df_file)
    )


# Test equality
solution.test_df_are_equal(testdf=temp_diurnal, solutiondf=get_diurnal_cycle_solution())


# %%


# =============================================================================
# aggregate method
# =============================================================================

agg_df = an.aggregate_df(agg=["lcz", "hour"])
assert agg_df.shape == (216, 10), "aggregate on analysis problem"


# =============================================================================
# Correlation check
# =============================================================================
import numpy as np

test = an.get_lc_correlation_matrices(
    obstype=["temp", "humidity"], groupby_labels=["lcz", "season"]
)


# plot test
an.plot_correlation_heatmap(groupby_value=("Open lowrise", "autumn"))


cor_solution_file = "cor_mat_solution.pkl"


def _create_cor_mat_solution():
    print("SOLUTION WILL BE OVERWRITTEN")
    trg = os.path.join(solution.solutions_dir, cor_solution_file)
    an.lc_cor_dict[("Open lowrise", "autumn")]["cor matrix"].to_pickle(trg)


# _create_cor_mat_solution()
def get_cor_mat_solution():
    return pd.read_pickle(os.path.join(solution.solutions_dir, cor_solution_file))


# Test equality
solution.test_df_are_equal(
    testdf=an.lc_cor_dict[("Open lowrise", "autumn")]["cor matrix"],
    solutiondf=get_cor_mat_solution(),
)


# scatter plot test

an.plot_correlation_variation()


# =============================================================================
# Test represetation
# =============================================================================
print(an)

# =============================================================================
# Test filter method
# =============================================================================

alldf = an.get_full_dataframe()

assert alldf.shape == (91080, 32), "returned df not correct"

# filter the data
filterdf = alldf[
    (alldf["obstype"] == "temp")
    & (alldf["value"] > 12.2)
    & (alldf["hour"] > 16)
    & (alldf["lcz"] == "Open midrise")
]

# set data
an.set_data(df=filterdf)


assert an.df.shape == (313, 1), "set data not correct"
assert an.metadf.shape[0] == 3, "metadata not filterd after data set"
# %%
test = alldf.set_index(["name", "datetime", "obstype"]).unstack(level="obstype")
