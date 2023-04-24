#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:29:37 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
#sys.path.append(str(lib_folder))


import metobs_toolkit

#%% IO testdata

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata_small.csv')


dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile)
dataset.import_data_from_file(coarsen_timeres=True)


#%% Apply Qc on dataset level

dataset.apply_quality_control(obstype='temp')


outliersdf = dataset.combine_all_to_obsspace()
dataset.get_qc_stats(make_plot = False)


#%% Apply Qc on obstype not specified in settings



dataset.apply_quality_control(obstype='humidity')

#%% Apply QC on station level

sta = dataset.get_station('vlinder05')

sta.apply_quality_control()
sta.get_qc_stats(make_plot=True)
