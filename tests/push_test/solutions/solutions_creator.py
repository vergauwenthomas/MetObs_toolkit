#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 16:36:38 2024

@author: thoverga
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 16:55:54 2024

@author: thoverga
"""

import sys, os
from pathlib import Path
import pandas as pd
import numpy as np

solutions_dir = Path(__file__).resolve().parents[0]

lib_dir = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(lib_dir))
import metobs_toolkit


def test_df_are_equal(testdf, solutiondf):
    # test is columns are the same
    if not set(testdf.columns) == set(solutiondf.columns):
        print("The columns of the test and solutiondf are not the same")
        print(f"Columns in the testdf : {testdf.columns}")
        print(f"Columns in the solutionsdf : {solutiondf.columns}")
        raise SystemError("The columns of the test and solutiondf are not the same")

    # order columns alphabetically
    testdf = testdf.reindex(sorted(testdf.columns), axis=1)
    solutiondf = solutiondf.reindex(sorted(solutiondf.columns), axis=1)
    # common error/bug signal checks

    # check for duplicates in index
    if (not solutiondf.index.duplicated().any()) & (testdf.index.duplicated().any()):
        print("The tested dataframe contains duplicated indexes!: ")
        print(testdf.loc[testdf.index.duplicated()])
        raise SystemError("The tested dataframe contains duplicated indexes")

    # overlap_data = testdf[testdf.index.isin(solutiondf.index)]
    added_data = testdf[~testdf.index.isin(solutiondf.index)]
    missing_data = solutiondf[~solutiondf.index.isin(testdf.index)]
    # situation 1: testdf has additive rows wrt solution
    if (not added_data.empty) & (missing_data.empty):
        print("These rows are found IN the test but NOT IN the solution: \n")
        print(added_data)
        raise SystemError("Rows are found IN the test but NOT IN the solution")

    # situation 2 tes # situation 1: testdf has additive rows wrt solution
    if (added_data.empty) & (not missing_data.empty):
        print("These rows are missing in the test, that are in the solution : \n")
        print(missing_data)
        raise SystemError("Rows are missing in the test")

    # situation 3 some rows are lacking, and some rows are added
    if (not added_data.empty) & (not missing_data.empty):
        print(
            "there is a difference between the test and the solution : \n  --------------------"
        )
        print(
            f"The following rows are IN the test but NOT in the solution: \n {added_data}"
        )
        print(" .................................")
        print(f"The following rows missing in  the test: \n {missing_data}")
        raise SystemError("There is a difference between the test and the solution")

    # test is columns are the same
    if not (testdf.index == solutiondf.index).all():
        print("The index of the test and solutiondf are not the same!")
        # check the names and hierarcy of index levels
        if list(testdf.index.names) != list(solutiondf.index.names):
            print(
                "The index structure is not the same between the test and the solution"
            )
            print(f"Index structure of the test: {testdf.index.names}")
            print(f"Index structure of the solution: {solutiondf.index.names}")
            raise SystemError("The index of the test and solutiondf are not the same!")

    are_equal = testdf.equals(solutiondf)
    if not are_equal:
        # apply a floating point mapping
        # for col in testdf.columns:
        #     if np.issubdtype(testdf[col].dtype, np.number):
        #         testdf[col] = testdf[col].astype('float16')
        #         solutiondf[col] = solutiondf[col].astype('float16')

        diffdf = testdf.compare(
            solutiondf, keep_shape=False, result_names=("TEST", "SOLUTION")
        )

        print("The stucture of datasets is equal, but the content differs!! ")
        print(f"See the diff_df for more details: {diffdf}")
        return diffdf
        # raise SystemError("There is a difference between the test and the solution"

    print("ok")
    return None


# Fix filenames
qc_applied_demo_pkl = "qc_applied_demo_dataset_5min_res.pkl"
gaps_from_qc_demo_pkl = "gaps_from_qc_applied_demo_dataset_5min_res.pkl"
itnerp_gaps_from_qc_demo_pkl = "interp_gaps_from_qc_applied_demo_dataset_5min_res.pkl"

qc_applied_demo_comb_df_pkl = "comb_df_qc_applied_demo_5min_res.pkl"


hourly_qc_applied_demo_pkl = "hourly_qc_applied_demo_dataset.pkl"
hourly_gaps_from_qc_demo_pkl = "hourly_gaps_from_qc_applied_demo_dataset.pkl"
hourly_itnerp_gaps_from_qc_demo_pkl = (
    "hourly_interp_gaps_from_qc_applied_demo_dataset.pkl"
)

hourly_qc_applied_demo_comb_df_pkl = "hourly_comb_df_qc_applied_demo.pkl"


# =============================================================================
# Create general solutions
# =============================================================================

# dataset pickles 5min resolution


def get_demo_qc_applied_5_min_dataset():
    dataset = metobs_toolkit.Dataset()
    new_dataset = dataset.import_dataset(
        folder_path=solutions_dir, filename=qc_applied_demo_pkl
    )
    return new_dataset


def get_demo_outliers_as_gaps_5_min_dataset():
    dataset = metobs_toolkit.Dataset()
    new_dataset = dataset.import_dataset(
        folder_path=solutions_dir, filename=gaps_from_qc_demo_pkl
    )
    return new_dataset


def get_demo_outliers_as_gaps_and_interp_5_min_dataset():
    dataset = metobs_toolkit.Dataset()
    new_dataset = dataset.import_dataset(
        folder_path=solutions_dir, filename=itnerp_gaps_from_qc_demo_pkl
    )
    return new_dataset


# dataset pickles hourly resolution
def get_demo_qc_applied_hourly_dataset():
    dataset = metobs_toolkit.Dataset()
    new_dataset = dataset.import_dataset(
        folder_path=solutions_dir, filename=hourly_qc_applied_demo_pkl
    )
    return new_dataset


def get_demo_outliers_as_gaps_hourly_dataset():
    dataset = metobs_toolkit.Dataset()
    new_dataset = dataset.import_dataset(
        folder_path=solutions_dir, filename=hourly_gaps_from_qc_demo_pkl
    )
    return new_dataset


def get_demo_outliers_as_gaps_and_interp_hourly_dataset():
    dataset = metobs_toolkit.Dataset()
    new_dataset = dataset.import_dataset(
        folder_path=solutions_dir, filename=hourly_itnerp_gaps_from_qc_demo_pkl
    )
    return new_dataset


# Dataframes


def get_demo_qc_applied_5_min_comb_df():
    targetfile = os.path.join(solutions_dir, qc_applied_demo_comb_df_pkl)
    df = pd.read_pickle(targetfile)
    return df


def get_demo_qc_applied_hourly_comb_df():
    targetfile = os.path.join(solutions_dir, hourly_qc_applied_demo_comb_df_pkl)
    df = pd.read_pickle(targetfile)
    return df


# =============================================================================
# solutions of metadata
# =============================================================================

demo_lcz_series_file = "demo_lcz_series.pkl"
demo_altitude_series_file = "demo_altiude_series.pkl"
demo_all_metadata = "demo_all_metadata_df.pkl"


def create_metadata_solution():
    print("THIS WILL OVERWRITE THE SOLUTIONS")
    dataset = metobs_toolkit.Dataset()
    dataset.update_file_paths(
        input_data_file=metobs_toolkit.demo_datafile,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        template_file=metobs_toolkit.demo_template,
    )
    dataset.import_data_from_file()
    dataset.get_lcz()
    dataset.get_altitude()
    dataset.get_landcover(buffers=[50, 100, 250], aggregate=True)
    # write series to pickles

    dataset.metadf["lcz"].to_pickle(os.path.join(solutions_dir, demo_lcz_series_file))
    dataset.metadf["altitude"].to_pickle(
        os.path.join(solutions_dir, demo_altitude_series_file)
    )
    dataset.metadf.to_pickle(os.path.join(solutions_dir, demo_all_metadata))


# create_metadata_solution()


def get_all_metadata_solution():
    return pd.read_pickle(os.path.join(solutions_dir, demo_all_metadata))


def get_lcz_series_solution():
    return pd.read_pickle(os.path.join(solutions_dir, demo_lcz_series_file))


def get_altitude_series_solution():
    return pd.read_pickle(os.path.join(solutions_dir, demo_altitude_series_file))


# # # ----- QC demo dataset ------------
# print('BE AWARE THAT THE SOLUTIONS WILL BE OVERWRITTEN !!! ')
# dataset = metobs_toolkit.Dataset()
# # ------- 5 min demo dataset ---------------
# dataset.update_file_paths(
#                           input_data_file=metobs_toolkit.demo_datafile,
#                           input_metadata_file=metobs_toolkit.demo_metadatafile,
#                           template_file=metobs_toolkit.demo_template,
#                           )
# dataset.import_data_from_file()
# dataset.apply_quality_control()

# # Dataset as pickle
# dataset.save_dataset(outputfolder=solutions_dir,
#                       filename=qc_applied_demo_pkl,
#                       overwrite=True)

# # Outliers as gaps
# dataset.convert_outliers_to_gaps()
# # Dataset as pickle
# dataset.save_dataset(outputfolder=solutions_dir,
#                      filename=gaps_from_qc_demo_pkl,
#                      overwrite=True)

# # interpolate the gaps
# dataset.interpolate_gaps(max_consec_fill=20)
# # Dataset as pickle

# dataset.save_dataset(outputfolder=solutions_dir,
#                      filename=itnerp_gaps_from_qc_demo_pkl,
#                      overwrite=True)

# # Dataframe as pickle
# combdf = dataset.combine_all_to_obsspace()
# targetfile = os.path.join(solutions_dir, qc_applied_demo_comb_df_pkl)
# combdf.to_pickle(targetfile)


# # -------  hourly demo dataset ---------------
# dataset = metobs_toolkit.Dataset()
# dataset.update_file_paths(
#                           input_data_file=metobs_toolkit.demo_datafile,
#                           input_metadata_file=metobs_toolkit.demo_metadatafile,
#                           template_file=metobs_toolkit.demo_template,
#                           )
# dataset.import_data_from_file()
# dataset.coarsen_time_resolution(freq='1H')
# dataset.apply_quality_control()

# # Dataset as pickle
# dataset.save_dataset(outputfolder=solutions_dir,
#                       filename=hourly_qc_applied_demo_pkl,
#                       overwrite=True)

# # Outliers as gaps
# dataset.convert_outliers_to_gaps()
# # Dataset as pickle
# dataset.save_dataset(outputfolder=solutions_dir,
#                       filename=hourly_gaps_from_qc_demo_pkl,
#                       overwrite=True)

# # interpolate the gaps
# dataset.interpolate_gaps(max_consec_fill=20)
# # Dataset as pickle
# dataset.save_dataset(outputfolder=solutions_dir,
#                       filename=hourly_itnerp_gaps_from_qc_demo_pkl,
#                       overwrite=True)


# # Dataframe as pickle
# combdf = dataset.combine_all_to_obsspace()
# targetfile = os.path.join(solutions_dir, hourly_qc_applied_demo_comb_df_pkl)
# combdf.to_pickle(targetfile)


# =============================================================================
# Analysis test
# =============================================================================


# -----------   Create solutions  -------------------
def get_analysis_dataset_solution():
    return get_demo_qc_applied_5_min_dataset()
