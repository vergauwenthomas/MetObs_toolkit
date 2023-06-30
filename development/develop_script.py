#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit
#
import os
import sys
from pathlib import Path
import pandas as pd
import time
import math


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%
# # use_dataset = 'debug_wide'
# use_dataset = 'single_netatmo_sara_station'
use_dataset = 'demo'
dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
dataset.apply_quality_control()
dataset.update_gaps_and_missing_from_outliers()
#%%


station = dataset.get_station('vlinder05')
# station.fill_gaps_linear()
# station.fill_missing_obs_linear()

station.make_plot(colorby='label')



# mergedf = station.combine_all_to_obsspace()
# mergedf = mergedf.xs('temp', level='obstype')

#%%
station.fill_gaps_linear()
#%%
station.make_plot(colorby='label', title='after fix')

# station.get_gaps_info()
mergedf = station.combine_all_to_obsspace()
mergedf = mergedf.xs('temp', level='obstype')


#%%

# missing = station.missing_obs

# misrec = missing.idx
# mis_series = missing.series
# misfil = missing.fill_df

# # cheat method
# unfilled = misfil[misfil['temp'].isnull()]


# misfil = misfil.dropna(subset='temp')

# test_unfilled = misrec[~misrec.isin(misfil.index)]




