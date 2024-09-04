#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 09:53:55 2023

@author: thoverga
"""

import sys, os

from pathlib import Path

from os.path import join

# add the solutions
sys.path.insert(0, str(Path(__file__).resolve().parents[0]))
print(sys.path)
import solutions.solutions_creator as solution

# point to current version of the toolkit
lib_folder = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(lib_folder))
import metobs_toolkit


test_data_dir = os.path.join(str(lib_folder), "tests", "test_data")


testdata = {
    # demo
    "demo": {
        "datafile": metobs_toolkit.demo_datafile,
        "metadatafile": metobs_toolkit.demo_metadatafile,
        "template": metobs_toolkit.demo_template,
        # "kwargs": {},
        "coarsen": "20min",
    },
    # paper dataset (based on the demo dataset)
    "paper_dataset": {
        "datafile": join(test_data_dir, "paper_dataset", "datafile.csv"),
        "metadatafile": join(test_data_dir, "paper_dataset", "metadatafile.csv"),
        "template": join(test_data_dir, "paper_dataset", "paper_dataset_template.json"),
        # "kwargs": {},
        "coarsen": "20min",
    },
    # wide test dataset
    "debug_wide": {
        "datafile": join(test_data_dir, "debug_wide.csv"),
        "metadatafile": join(test_data_dir, "debug_wide_metadata.csv"),
        "template": join(test_data_dir, "debug_wide_template.json"),
        # "kwargs": {"long_format": False, "obstype": "temp"},
        "coarsen": "20min",
    },
    # Single station dataset
    "single_station": {
        "datafile": join(test_data_dir, "single_station.csv"),
        "metadatafile": join(test_data_dir, "single_station_metadata.csv"),
        "template": join(test_data_dir, "single_station_template.json"),
        # "kwargs": {"long_format": False, "obstype": "temp"},
        "coarsen": "20min",
    },
    # breaking
    "breaking data": {
        "datafile": join(test_data_dir, "testdata_breaking.csv"),
        "metadatafile": None,
        "template": join(test_data_dir, "template_breaking.json"),
        # "kwargs": {},
        "coarsen": "60min",
    },
    # Kobe congo (single station)
    "Congo_single_station": {
        "datafile": join(
            test_data_dir, "testdata_testday", "Kobe", "meteo_soil_clean_2023-01-19.csv"
        ),
        "metadatafile": join(
            test_data_dir, "testdata_testday", "Kobe", "CONGO_meta.csv"
        ),
        "template": join(
            test_data_dir, "testdata_testday", "Kobe", "kongo_template.json"
        ),
        # "kwargs": {},
        "coarsen": "60min",
    },
    # Single Netatmo station Sara
    "single_netatmo_sara_station": {
        "datafile": join(
            test_data_dir,
            "testdata_testday",
            "Sara",
            "Outdoor_module_Netatmo_Sara_small.csv",
        ),
        "metadatafile": join(
            test_data_dir,
            "testdata_testday",
            "Sara",
            "metadata_Outdoor_module_Netatmo_Sara_new.csv",
        ),
        "template": join(
            test_data_dir, "testdata_testday", "Sara", "sara_template.json"
        ),
        # "kwargs": {"freq_estimation_method": "median"},
        "coarsen": "60min",
    },
    # Vlinders 2022
    "vlindergent2022": {
        "datafile": join(
            test_data_dir, "testdata_testday", "Sara", "Vlinder_gent_2022.csv"
        ),
        "metadatafile": join(
            test_data_dir, "testdata_testday", "Sara", "all_vlinders_metadata.csv"
        ),
        "template": join(
            test_data_dir, "testdata_testday", "Sara", "vlinders22_template.json"
        ),
        # "kwargs": {"freq_estimation_method": "median"},
        "coarsen": "60min",
    },
    # amsterdam
    "amsterdam": {
        "datafile": join(
            test_data_dir,
            "testdata_testday",
            "amsterdam",
            "Amsterdam_D2222z6together.csv",
        ),
        "metadatafile": join(
            test_data_dir,
            "testdata_testday",
            "amsterdam",
            "Latlon_stations_Amsterdam.csv",
        ),
        "template": join(
            test_data_dir, "testdata_testday", "amsterdam", "amsterdam_template.json"
        ),
        # "kwargs": {},
        "coarsen": "60min",
    },
}
