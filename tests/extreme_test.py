#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Idea is to test all testfiles with all kind of datasets

@author: thoverga
"""

import copy
import sys, os

from pathlib import Path

lib_folder = Path(__file__).resolve().parents[1]

# print(str(lib_folder))

sys.path.insert(0, str(lib_folder))
import metobs_toolkit

# %%
testfolder = os.path.join(str(lib_folder), "tests", "push_test")

from tests.push_test.test_data_paths import testdata


# %%
def read_in_the_dataset(dataname, testdatadict):
    print(f"\n ------ read dataset ({dataname}) ---------\n")
    datafile = testdatadict[dataname]["datafile"]
    metafile = testdatadict[dataname]["metadatafile"]
    template = testdatadict[dataname]["template"]
    # kwargsdict = testdatadict[dataname]["kwargs"]

    dataset = metobs_toolkit.Dataset()
    dataset.update_settings(
        input_data_file=datafile,
        input_metadata_file=metafile,
        template_file=template,
    )

    dataset.import_data_from_file()
    return dataset


def IO_test(dataset, name):
    print(f"\n ------ IO tests ({name}) ---------\n")

    def del_file(file_path):
        if os.path.isfile(file_path):
            os.remove(file_path)
            print(f"{file_path} deleted.")
        else:
            print(f"{file_path} not found.")

    # Sycnronize data
    test = dataset.sync_observations(tolerance="5T")

    # pickel test
    outfolder = os.path.join(str(lib_folder), "tests", "test_data")
    file = "dataset_IO_test"

    del_file(os.path.join(outfolder, file + ".pkl"))

    # save dataset as pickle
    dataset.update_default_name("this_is_a_test_name")

    dataset.save_dataset(outputfolder=outfolder, filename=file)

    del dataset  # remove from kernel

    # read dataset
    new_dataset = metobs_toolkit.Dataset()
    new_dataset = new_dataset.import_dataset(
        folder_path=outfolder, filename=file + ".pkl"
    )

    del_file(os.path.join(outfolder, file + ".pkl"))

    assert (
        new_dataset.settings.app["default_name"] == "this_is_a_test_name"
    ), "some attributes are not correctly saved when pickled."


def qc_testing(dataset, name):
    print(f"\n ------ QC tests ({name}) ---------\n")

    # on get station
    stationname = dataset.metadf.index[0]
    station = dataset.get_station(stationname)

    station.apply_quality_control(obstype="temp")
    station.get_qc_stats(make_plot=False)

    # on dataset
    dataset.get_qc_stats(make_plot=False)
    dataset.apply_quality_control(obstype="temp")
    dataset.get_qc_stats(make_plot=True)

    # titan test
    dataset.update_titan_qc_settings(
        obstype="temp",
        buddy_radius=50000,
        buddy_num_min=3,
        buddy_max_elev_diff=200,
        buddy_threshold=3,
    )


def gapfill_testing(dataset, name):
    print(f"\n ------ gaps missing tests ({name})---------\n")

    # testing conversion to df + update from outliers
    _ = dataset.get_gaps_df()

    init_outl_shape = dataset.outliersdf.shape
    dataset.update_gaps_and_missing_from_outliers(n_gapsize=3)
    if init_outl_shape[0] > 0:

        assert (
            init_outl_shape != dataset.outliersdf.shape
        ), "outliers still the same as before updateing to gaps"

    _ = dataset.get_gaps_df()

    dataset.get_gaps_info()
    dataset.get_missing_obs_info()
    # filling
    dataset.fill_missing_obs_linear()
    dataset.fill_gaps_linear()

    dataset.get_gaps_info()
    dataset.get_missing_obs_info()

    # plot test
    # dataset.make_plot(colorby='label', title='AFTER GAP AND MISSING FILL')


def plot_testing(dataset, name):
    print(f"\n ------ plot tests ({name})---------\n")

    dataset.make_plot(colorby="name", title=name)
    dataset.make_plot(colorby="label", title=name)

    if not dataset.metadf.empty:
        dataset.make_geo_plot(variable="temp", title=name)


def analysis_test(dataset, name):
    print(f"\n ------ Analysis testing({name})---------\n")

    an = dataset.get_analysis()

    # Test plotting and functions
    temp_diurnal = an.get_diurnal_statistics(colorby="lcz", title=name)
    an.get_anual_statistics(agg_method="median", plot=False)
    test3 = an.get_aggregated_cycle_statistics(aggregation=["lcz"], title=name)

    print(an)

    filter_an = an.apply_filter("temp < 15.5 &  hour <= 19")

    agg_df = an.aggregate_df(agg=["lcz", "hour"])

    if "humidity" in dataset.df.columns:
        an.get_lc_correlation_matrices(
            obstype=["temp", "humidity"], groupby_labels=["lcz", "season"]
        )
    else:
        an.get_lc_correlation_matrices(
            obstype=["temp"], groupby_labels=["lcz", "season"]
        )


def get_lcz_and_lc(name, dataset):
    print(f"\n ------ gee lcz and lc extraction ({name})---------\n")

    dataset.get_lcz()
    dataset.get_landcover(buffers=[50, 100])

    # save for later use if needed
    metadf = dataset.metadf.copy()

    # relevant columns:
    rel_columns = [
        col for col in metadf.columns if (col.endswith("50m") | col.endswith("100m"))
    ]
    rel_columns.append("lcz")

    metadf = metadf[rel_columns]
    metadf = metadf.reset_index()
    filename = meta_path_generator(name)

    metadf.to_csv(filename)


def meta_path_generator(name):
    metafolder = os.path.join(
        str(lib_folder), "tests", "test_data/meta_data_extreme_test"
    )

    filename = name.replace(" ", "_") + "_lc_info.csv"

    return os.path.join(metafolder, filename)


# %%

for name in testdata:
    print(f"\n ************ {name} *************\n")
    dataset = read_in_the_dataset(name, testdata)
    print(f"Initial df shape: {dataset.df.shape}")
    dataset.coarsen_time_resolution(freq=testdata[name]["coarsen"])
    # print(f"after coarsening df shape: {dataset.df.shape}")
    # qc_testing(dataset, name)
    # plot_testing(dataset, name)
    # gapfill_testing(dataset, name)
    # get_lcz_and_lc(name, dataset)
    # analysis_test(dataset, name)
    # IO_test(dataset, name)

# dataset = read_in_the_dataset(name, testdata)
# dataset.coarsen_time_resolution(freq=testdata[name]['coarsen'])
# get_lcz_and_lc(name, dataset)
# an = dataset.get_analysis()
# #%%
# an.get_lc_correlation_matrices(obstype=['temp', 'humidity'], groupby_labels=['lcz', 'hour'])


# an.plot_correlation_variation()
