#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.append(str(lib_folder))




#%% % Import


testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/vlinder_metadata.csv'
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')




# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        input_metadata_file=static_data,
                        # data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )

obstype = 'temp'



print('before : ', dataset.settings.qc['qc_check_settings']["gross_value"][obstype]['min_value'])

dataset.update_qc_settings(obstype=obstype, gross_value_min_value=-100)


print('after: ', dataset.settings.qc['qc_check_settings']["gross_value"][obstype]['min_value'])



#%% Datset2


dataset2 = metobs_toolkit.Dataset()
dataset2.update_settings(input_data_file=testdatafile,
                        output_folder='/home/thoverga/Documents'
                        )

print('before2 : ', dataset2.settings.qc['qc_check_settings']["gross_value"][obstype]['min_value'])

dataset2.update_qc_settings(obstype=obstype, gross_value_min_value=-100)


print('after2: ', dataset2.settings.qc['qc_check_settings']["gross_value"][obstype]['min_value'])




#%%

from metobs_toolkit.settings import Settings

testset = Settings()
print('after settings creation: ', testset.qc['qc_check_settings']["gross_value"][obstype]['min_value'])

