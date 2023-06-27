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


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit



#%%

import pandas as pd
import pytz


mmm = metobs_toolkit.Modeldata('cordex')

csvpath= '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/Alaro_2.5_melle.csv'

mmm.set_alaro_25_model_from_csv(csvpath)




#%%

# # use_dataset = 'debug_wide'
use_dataset = 'demo'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )


dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq ='60T')
# dataset.apply_quality_control()
#%%


mmm.df = mmm.df.reset_index()
mmm.df['name']= 'vlinder01'
mmm.df = mmm.df.set_index(['name', 'datetime'])

mmm.make_plot(obstype_model='temp_TEB', dataset=dataset,
              obstype_dataset = 'temp')


