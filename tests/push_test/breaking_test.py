#!/usr/bin/env python3
# -*- coding: utf-8 -*-



# PROGRAM FOR TESTING THE BREAKING DATAFILE
"""
Created on Tue Nov 29 12:19:03 2022

@author: mivieijra
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

testdata = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')

settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdata)
settings.check_settings()
dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=False)
dataset.apply_quality_control(persistance = False, repetitions=True)
station = dataset.get_station('Fictional')
