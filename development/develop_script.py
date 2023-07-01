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
use_dataset = 'vlindergent2022'
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
# dataset.update_gaps_and_missing_from_outliers()
#%%

ann = dataset.get_analysis()

#%%
# Compute mean annual cycle for each station + plot
ann_mean_df = ann.get_anual_statistics(groupby=['name'], #each line represents a station
                                                obstype='temp', # on temperatures
                                                agg_method='mean', #value of the line is the means of the aggregation
                                                stations=None, #use all stations
                                                startdt=None, #use the full analysis
                                                enddt=None, #use the full analysis
                                                plot=True,
                                                errorbands=False, #Display the std as a band
                                                title = None,
                                                )

