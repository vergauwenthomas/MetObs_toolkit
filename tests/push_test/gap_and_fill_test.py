#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 08:18:09 2022

@author: thoverga
"""

import sys, os

from pathlib import Path


lib_folder = Path(__file__).resolve().parents[2]
#sys.path.append(str(lib_folder))


import vlinder_toolkit

#%% Import data

testdatafile = os.path.join(
    str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')


# #% import data



dataset = vlinder_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                         input_metadata_file=static_data,
                         output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                         )
dataset.import_data_from_file(coarsen_timeres=True)

#%% Basic tests on the gaps
gapsdf = dataset.gaps.to_df()
# check if two gaps are found
assert gapsdf.shape[0] == 2, f'There are assumed 2 gaps, but the tlkit found {gapsdf.shape[0]}'

assert list(gapsdf.index.unique()) == ['vlinder01'], f'Only gaps assumed in vlinder01. Tlkit found gaps for these {list(gapsdf.index.unique())}'


#%% Basic tests on missing obs
missingobs = dataset.missing_obs.series
# check if two gaps are found
assert missingobs.shape[0] == 26, f'There are assumed 26 missing obs, but the tlkit found {missingobs.shape[0]}'

assert list(missingobs.index.unique()) == ['vlinder01', 'vlinder02', 'vlinder03'], \
             f'Only missing obs assumed in vl01, vl02, vl03. Tlkit found missing obs for these {list(missingobs.index.unique())}'


#%% Test functions on gaps

stagap = dataset.gaps.get_station_gaps('vlinder01')

dataset.gaps.remove_gaps_from_obs(dataset.df)

dataset.gaps.get_gaps_indx_in_obs_space(dataset.df, dataset.outliersdf, dataset.metadf['dataset_resolution'])


dataset.gaps.list[0].update_leading_trailing_obs(dataset.df, dataset.outliersdf)

#%% Test gapfilling


interp_series = dataset.gaps.apply_interpolate_gaps(dataset.df, dataset.outliersdf, dataset.metadf['dataset_resolution'])
interp_series.name = 'tlk'

interpdf = interp_series.to_frame()

interpdf['manual']=[     
                9.487500,
                9.875000,
                10.262500,
                10.650000,
                11.037500,
                11.425000,
                11.812500,
                15.833333,
                14.766667,
                13.700000,
                12.633333,
                11.566667]

interpdf['diff'] = interpdf['tlk'] - interpdf['manual']


assert not interpdf['tlk'].isnull().any(), f'Nan value found in tlk interpolation ! \n {interpdf}'
assert (interpdf['diff']).sum() < 1e-5, f'Tlk interpolation differs from manual: \n {interpdf}'



print('Done! ')

