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


# testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_data.csv')
testdatafile = '/home/thoverga/Downloads/testdataset.csv'



template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')




# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        # input_metadata_file=static_data,
                        data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )



dataset.import_data_from_file(long_format=False, obstype='temp')

dataset.show()

dfinit = dataset.df.copy()
metadfinit = dataset.metadf.copy()
outliersdfinit = dataset.outliersdf.copy()
#%%
test = dataset.sync_observations(tollerance='5T', verbose=True)


dataset.show()



#%%
import pandas as pd


# df = pd.DataFrame({'A': [1, 2, 3, 4, 5],
#                    'B': [1, 1, 1, pd.np.nan, pd.np.nan],
#                    'C': [pd.np.nan, 2, 3, pd.np.nan, pd.np.nan],
#                    'D': [pd.np.nan, pd.np.nan, pd.np.nan, pd.np.nan, pd.np.nan]})

# column_names = df.nunique().loc[lambda x: x <= 1].index.tolist()