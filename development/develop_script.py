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

import metobs_toolkit


#%%
datafile = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_okt.csv"
# metafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/vlinder_metadata.csv'



# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        # data_template_file = metobs_toolkit.demo_template,
                        # metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
                        )

dataset.import_data_from_file()

#%%
import copy
dataset.coarsen_time_resolution(freq='30T')

print(f' gaps before qc: {len(dataset.gaps)}')
print(f' missing before qc: {len(dataset.missing_obs)}')

dataset.apply_quality_control()

print(f' gaps after qc: {len(dataset.gaps)}')
print(f' missing after qc: {len(dataset.missing_obs)}')
dataset.update_gaps_and_missing_from_outliers(n_gapsize=4)


print(f' gap  after update: {len(dataset.gaps)}')
print(f' missing after update: {len(dataset.missing_obs)}')




teststation = 'vlinder09'

sta = dataset.get_station(teststation)
sta = copy.deepcopy(sta)
print(f' station gap  after update: {len(sta.gaps)}')
print(f' station missing after update: {len(sta.missing_obs)}')
print(sta.gaps[0].get_info())


dataset.fill_gaps_linear()


sta2 = dataset.get_station(teststation)

print('after fill', sta2.gaps[0].get_info())



#%%

# sta = dataset.get_station('vlinder05')
# sta.make_plot(title = 'before qc')

#%%


# #%%
# dataset.apply_quality_control()
# sta.make_plot(colorby='label', title='after qc')

#%%
# dataset.update_gaps_and_missing_from_outliers()
# sta.make_plot(colorby='label', title='after updating gaps')


#%%
# sta.fill_gaps_linear()
# sta.make_plot(colorby='label', title='after gaps fill')

# print(sta.gaps[0].get_info())

