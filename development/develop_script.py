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

#%% Setup dataset
settings = vlinder_toolkit.Settings()
settings.update_settings(input_file=testdatafile)

dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file()

#%%
# =============================================================================
# developments
# =============================================================================

sta = dataset.get_station('vlinder05')

stachecked = vlinder_toolkit.qc_checks.duplicate_timestamp(sta)


temp = sta.temp.append(pd.Series(index=[sta.temp.index[2]], data=[sta.temp.iloc[2]]))






sta.temp = temp
sta= vlinder_toolkit.qc_checks.duplicate_timestamp(sta)
sta= vlinder_toolkit.qc_checks.gross_value_check(sta)



print(sta.qc_labels_df['temp'])


