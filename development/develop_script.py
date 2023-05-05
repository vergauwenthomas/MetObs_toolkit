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


# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')
# template = os.path.join(str(lib_folder), 'tests', 'test_data',  'wide_test_template.csv')




# #% Setup dataset

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                        # input_metadata_file=static_data,
                        # data_template_file= template,
                        output_folder='/home/thoverga/Documents'
                        )



dataset.import_data_from_file()


#%%
from metobs_toolkit.analysis import Analysis
an = Analysis(obsdf = dataset.df,
              metadf = dataset.metadf,
              settings = dataset.settings)


df = an.get_diurnal_statistics(refstation='vlinder08')



# #%% Make plot
# import matplotlib.pyplot as plt

# refstation = 'vlinder08'


# # Make title
# if title is None:
#     if start


#     if refstation is None:
#         title = "Diurnal cycle "




# fig, ax = plt.subplots
