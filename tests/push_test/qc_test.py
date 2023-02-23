#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:29:37 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
sys.path.append(str(lib_folder))


import vlinder_toolkit

#%% IO testdata

testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'vlinderdata_small.csv')



settings = vlinder_toolkit.Settings()
settings.update_settings(input_data_file=testdatafile)
settings.check_settings()

dataset = vlinder_toolkit.Dataset()
dataset.import_data_from_file(coarsen_timeres=True)


#%% Apply Qc on dataset level

dataset.apply_quality_control(obstype='temp',
                                            gross_value=True, #apply this check 
                                            persistance=True, #apply this check
                                            )


dataset.get_final_qc_labels()



#%% Apply Qc on obstype not specified in settings



dataset.apply_quality_control(obstype='humidity',
                                            gross_value=True, #apply this check 
                                            persistance=True, #apply this check
                                            )


