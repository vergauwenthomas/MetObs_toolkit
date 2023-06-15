#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%

# single station
file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/single_station.csv'
metadfile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/single_station_metadata.csv'


# wide
file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide.csv'
metadfile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_metadata.csv'





metobs_toolkit.build_template_prompt(debug=True)
