#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:28:15 2022

@author: thoverga
"""

from pathlib import Path
import os, sys

lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))
from src import vlinder_toolkit

#%% 
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata.csv')

#%% 
settings = vlinder_toolkit.Settings()
settings.update_settings(input_file=testdatafile)


df = vlinder_toolkit.import_data(settings)



