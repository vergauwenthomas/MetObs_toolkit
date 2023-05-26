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
metafile = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/data.csv'



# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = datafile,
                        input_metadata_file = metafile,
                        # data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template, # Contains also the metadata mapping
                        output_folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/development'
                        )

dataset.import_data_from_file()

#%%

dataset.coarsen_time_resolution(freq='30T')
dataset.apply_quality_control()
dataset.update_gaps_and_missing_from_outliers(n_gapsize=4)


# era = dataset.get_modeldata()
era = metobs_toolkit.Modeldata('era5')
era.set_model_from_csv(csvpath='/home/thoverga/Downloads/era_5_data.csv')

dataset.fill_gaps_linear()
# dataset.fill_gaps_era5(modeldata=era)

#%%
sta = dataset.get_station('vlinder36')
sta.gaps[4].get_info()


# sta.gaps[4].gapfill_info


# i = 0
# while i < len(sta.gaps):
#     print(f' ITEM {i}')
#     sta.gaps[i].get_info()
#     i+=1

# def format_gapfill_info(gap, obstype):
#     df = gap.gapfill_info
#     if gap.gapfill_technique == 'gap_debiased_era5': #extract from settings
#         if not df.empty:
#             # make time column
#             df['time'] = (df['hours'].astype(str).str.zfill(2) + ':' +
#                           df['minutes'].astype(str).str.zfill(2) + ':' +
#                           df['seconds'].astype(str).str.zfill(2))

#             # rename
#             df = df.rename(columns={obsytpe:f'{obstype}_model',
#                                     f'{obstype}_fill'})
#         else:







# test = format_gapfill_info(sta.gaps[4])









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

