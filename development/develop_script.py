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


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit




#%%

# # use_dataset = 'debug_wide'
use_dataset = 'single_netatmo_sara_station'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )
tz ='Europe/Brussels'

dataset.update_timezone(timezonestr=tz)


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# dataset.apply_quality_control()

#%%
dataset.sync_observations(tollerance='3T')
dataset.coarsen_time_resolution(freq='30T')


#%%
from datetime import datetime
import pytz



tstart = datetime(2021, 2, 28,9, tzinfo=pytz.timezone('UTC'))

dataset.make_plot(starttime=tstart)
#%%

from metobs_toolkit.df_helpers import multiindexdf_datetime_subsetting
import pytz





