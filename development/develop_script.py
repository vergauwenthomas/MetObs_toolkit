#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022

@author: thoverga
"""

# import vlinder_toolkit
import os, sys
import pandas as pd
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
print(str(lib_folder))

from src import vlinder_toolkit

#%% Import

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata.csv')

#% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_file=testdatafile)

dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)


# sta = dataset.get_station('vlinder05')

#%%
# =============================================================================
# checks
# =============================================================================

sta = dataset.get_station('vlinder05')

df_init = sta.df()
sta.make_plot(title='init temp')


sta = vlinder_toolkit.qc_checks.duplicate_timestamp(sta)
sta.make_plot(title='after timstamp dub qc')
sta = vlinder_toolkit.qc_checks.gross_value_check(sta)
sta.make_plot(title='after gross value qc')
sta = vlinder_toolkit.qc_checks.persistance(sta)
sta.make_plot(title='after persistance qc')

# df = sta.df()
# sta.make_plot()


# sta.make_plot()



