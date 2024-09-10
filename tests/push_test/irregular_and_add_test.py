#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 15:59:45 2024

@author: thoverga
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 26 11:30:06 2024

@author: thoverga
"""


import sys, os

from pathlib import Path
import pandas as pd

# add the solutions
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
print(sys.path)
import solutions.solutions_creator as solution

# point to current version of the toolkit
lib_folder = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(lib_folder))
import metobs_toolkit


# %% Importing tricky datast (no metadata + single station + irregular timestaps)
data2 = os.path.join(
    lib_folder,
    "tests",
    "test_data",
    "Outdoor_module_Netatmo_Sara_small.csv",
)
template2 = os.path.join(lib_folder, "tests", "test_data", "sara_template.json")

irr_combdf_file = "irr_combdf.pkl"
irr_metadf_file = "irr_metadf.pkl"


def _create_irr_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    dataset2 = metobs_toolkit.Dataset()
    dataset2.update_file_paths(
        input_data_file=data2,
        input_metadata_file=None,
        template_file=template2,
    )
    dataset2.import_data_from_file(
        freq_estimation_method="median",
        freq_estimation_simplify_tolerance="1min",
        origin_simplify_tolerance="3min",
        timestamp_tolerance="4min",
    )
    trgfile = os.path.join(solution.solutions_dir, irr_combdf_file)
    dataset2.get_full_status_df().to_pickle(trgfile)
    dataset2.metadf.to_pickle(os.path.join(solution.solutions_dir, irr_metadf_file))

    return


# _create_irr_solutions()


def get_irr_solutions():
    irr_combdf = pd.read_pickle(os.path.join(solution.solutions_dir, irr_combdf_file))
    irr_metadf = pd.read_pickle(os.path.join(solution.solutions_dir, irr_metadf_file))
    return irr_combdf, irr_metadf


irr_sol_combdf, irr_sol_metadf = get_irr_solutions()


dataset2 = metobs_toolkit.Dataset()
dataset2.update_file_paths(
    input_data_file=data2,
    input_metadata_file=None,
    template_file=template2,
)
dataset2.import_data_from_file(
    freq_estimation_method="median",
    freq_estimation_simplify_tolerance="1min",
    origin_simplify_tolerance="3min",
    timestamp_tolerance="4min",
)


# test the observation records
diff_df = solution.test_df_are_equal(
    testdf=dataset2.get_full_status_df(), solutiondf=irr_sol_combdf
)
assert diff_df is None

# test the metadata
diff_df = solution.test_df_are_equal(testdf=dataset2.metadf, solutiondf=irr_sol_metadf)
assert diff_df is None


# %% Combining two dataset without station overalp


add_combdf_file = "add_combdf.pkl"
add_metadf_file = "add_metadf.pkl"


def _create_add_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    dataset1 = metobs_toolkit.Dataset()
    dataset1.update_file_paths(
        input_data_file=metobs_toolkit.demo_datafile,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        template_file=metobs_toolkit.demo_template,
    )
    dataset1.import_data_from_file()
    dataset1.apply_quality_control()  # to see if outliers are combined aswell
    comb = dataset1 + dataset2

    trgfile = os.path.join(solution.solutions_dir, add_combdf_file)
    comb.get_full_status_df().to_pickle(trgfile)
    comb.metadf.to_pickle(os.path.join(solution.solutions_dir, add_metadf_file))

    return


# _create_add_solutions()


def read_add_solution():
    add_combdf = pd.read_pickle(os.path.join(solution.solutions_dir, add_combdf_file))
    add_metadf = pd.read_pickle(os.path.join(solution.solutions_dir, add_metadf_file))
    return add_combdf, add_metadf


add_sol_comb, add_sol_meta = read_add_solution()

dataset1 = metobs_toolkit.Dataset()
dataset1.update_file_paths(
    input_data_file=metobs_toolkit.demo_datafile,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)
dataset1.import_data_from_file()
dataset1.apply_quality_control()  # to see if outliers are combined aswell
comb = dataset1 + dataset2


# test the observation records
diff_df = solution.test_df_are_equal(
    testdf=comb.get_full_status_df(), solutiondf=add_sol_comb
)
assert diff_df is None

# test the metadata
diff_df = solution.test_df_are_equal(testdf=comb.metadf, solutiondf=add_sol_meta)
assert diff_df is None
# %%Combining two dataset with station overlap (no time overlap)

data3 = os.path.join(lib_folder, "tests", "test_data", "testdata_okt_small.csv")
add_sta_overlap_combdf_file = "add_vlinders_combdf.pkl"
add_sta_overlap_metadf_file = "add_vlinders_metadf.pkl"


def _create_add_sta_overlap_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    dataset3 = metobs_toolkit.Dataset()
    dataset3.update_file_paths(
        input_data_file=data3,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        template_file=metobs_toolkit.demo_template,
    )
    dataset3.import_data_from_file()

    comb = dataset1 + dataset3
    trgfile = os.path.join(solution.solutions_dir, add_sta_overlap_combdf_file)
    comb.get_full_status_df().to_pickle(trgfile)
    comb.metadf.to_pickle(
        os.path.join(solution.solutions_dir, add_sta_overlap_metadf_file)
    )
    return


# _create_add_sta_overlap_solutions()


def read_add_sta_overlap_solution():
    add_combdf = pd.read_pickle(
        os.path.join(solution.solutions_dir, add_sta_overlap_combdf_file)
    )
    add_metadf = pd.read_pickle(
        os.path.join(solution.solutions_dir, add_sta_overlap_metadf_file)
    )
    return add_combdf, add_metadf


add_sol_comb, add_sol_meta = read_add_sta_overlap_solution()

dataset3 = metobs_toolkit.Dataset()
dataset3.update_file_paths(
    input_data_file=data3,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)
dataset3.import_data_from_file()

comb = dataset1 + dataset3


# test the observation records
diff_df = solution.test_df_are_equal(
    testdf=comb.get_full_status_df(), solutiondf=add_sol_comb
)
assert diff_df is None

# test the metadata
diff_df = solution.test_df_are_equal(testdf=comb.metadf, solutiondf=add_sol_meta)
assert diff_df is None


# %%Combining two dataset with station overlap WITH time overlap

data4 = os.path.join(
    lib_folder, "tests", "test_data", "vlinders_okt_small_part_two.csv"
)
add_sta_period_overlap_combdf_file = "add_vlinders_period_combdf.pkl"
add_sta_period_overlap_metadf_file = "add_vlinders_period_metadf.pkl"


def _create_add_sta_period_overlap_solutions():
    print("WARNING!!! THE SOLUTION WILL BE OVERWRITTEN!")

    dataset4 = metobs_toolkit.Dataset()
    dataset4.update_file_paths(
        input_data_file=data4,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        template_file=metobs_toolkit.demo_template,
    )
    dataset4.import_data_from_file()

    comb = dataset3 + dataset4
    trgfile = os.path.join(solution.solutions_dir, add_sta_period_overlap_combdf_file)
    comb.get_full_status_df().to_pickle(trgfile)
    comb.metadf.to_pickle(
        os.path.join(solution.solutions_dir, add_sta_period_overlap_metadf_file)
    )
    return


# _create_add_sta_period_overlap_solutions()


def read_add_sta_overlap_solution():
    add_combdf = pd.read_pickle(
        os.path.join(solution.solutions_dir, add_sta_period_overlap_combdf_file)
    )
    add_metadf = pd.read_pickle(
        os.path.join(solution.solutions_dir, add_sta_period_overlap_metadf_file)
    )
    return add_combdf, add_metadf


add_sol_comb, add_sol_meta = read_add_sta_overlap_solution()

dataset4 = metobs_toolkit.Dataset()
dataset4.update_file_paths(
    input_data_file=data4,
    input_metadata_file=metobs_toolkit.demo_metadatafile,
    template_file=metobs_toolkit.demo_template,
)
dataset4.import_data_from_file()

comb = dataset3 + dataset4


# test the observation records
diff_df = solution.test_df_are_equal(
    testdf=comb.get_full_status_df(), solutiondf=add_sol_comb
)
assert diff_df is None

# test the metadata
diff_df = solution.test_df_are_equal(testdf=comb.metadf, solutiondf=add_sol_meta)
assert diff_df is None
