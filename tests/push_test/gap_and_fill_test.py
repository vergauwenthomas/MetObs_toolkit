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


import metobs_toolkit

#%% Import data

testdatafile = os.path.join(
    str(lib_folder), 'tests', 'test_data',  'testdata_okt_small.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')

# DEBUG
print('sys path: ', sys.path)

# #% import data

dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                         input_metadata_file=static_data,
                         output_folder='/home/thoverga/Documents/VLINDER_github/metobs_toolkit'
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

# ------------------ linear interpolation ----------------------------

dataset.fill_gaps_linear(obstype='temp')

gapsfilldf = dataset.gapfilldf

gapsfilldf['manual']=[
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

gapsfilldf['diff'] = gapsfilldf['temp'] - gapsfilldf['manual']


assert not gapsfilldf['temp'].isnull().any(), f'Nan value found in tlk interpolation ! \n {gapsfilldf}'
assert (gapsfilldf['diff']).sum() < 1e-5, f'Tlk interpolation differs from manual: \n {gapsfilldf}'

#%%
# ------------------ ERA debias fill
obstype='temp'


dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=testdatafile,
                         input_metadata_file=static_data,
                         output_folder='/home/thoverga/Documents/VLINDER_github/metobs_toolkit'
                         )

dataset.settings.time_settings['target_time_res'] = '30T'

dataset.import_data_from_file(coarsen_timeres=True)





#offline mode
era = metobs_toolkit.Modeldata('ERA5_hourly')

era_datafile = os.path.join(
    str(lib_folder), 'tests', 'test_data',  'era5_data.csv')


era.set_model_from_csv(era_datafile, 'ERA5_hourly', obstype)

assert era.df.shape[0] == 5348, 'Something wrong with importing era data from csv.'



# # Fill gaps using era5 data:
dataset.fill_gaps_era5(modeldata=era,
                        method='debias',
                        obstype='temp',
                        overwrite=True)


#%%
from datetime import datetime
import numpy as np
import pandas as pd
from metobs_toolkit.df_helpers import init_multiindexdf, conv_tz_multiidxdf


# validate

checked = {'temp': {('vlinder01',
      datetime.strptime('2022-10-04 02:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 02:30:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 03:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 03:30:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 04:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 04:30:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 05:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 05:30:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 06:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 06:30:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 07:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 07:30:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-04 08:00:00', '%Y-%m-%d %H:%M:%S')): np.nan,
     ('vlinder01',
      datetime.strptime('2022-10-06 17:00:00', '%Y-%m-%d %H:%M:%S')): 15.046638488769531,
     ('vlinder01',
      datetime.strptime('2022-10-06 17:30:00', '%Y-%m-%d %H:%M:%S')): 14.839948940277099,
     ('vlinder01',
      datetime.strptime('2022-10-06 18:00:00', '%Y-%m-%d %H:%M:%S')): 14.560950469970702,
     ('vlinder01',
      datetime.strptime('2022-10-06 18:30:00', '%Y-%m-%d %H:%M:%S')): 13.707459545135507,
     ('vlinder01',
      datetime.strptime('2022-10-06 19:00:00', '%Y-%m-%d %H:%M:%S')): 12.349718475341811,
     ('vlinder01',
      datetime.strptime('2022-10-06 19:30:00', '%Y-%m-%d %H:%M:%S')): 11.568361473083502,
     ('vlinder01',
      datetime.strptime('2022-10-06 20:00:00', '%Y-%m-%d %H:%M:%S')): 11.356404495239257,
     ('vlinder01',
      datetime.strptime('2022-10-06 20:30:00', '%Y-%m-%d %H:%M:%S')): 10.896783733367933,
     ('vlinder01',
      datetime.strptime('2022-10-06 21:00:00', '%Y-%m-%d %H:%M:%S')): 10.707939147949247}}


checkeddf = pd.DataFrame(checked)
checkeddf = checkeddf.reset_index()
checkeddf = checkeddf.rename(columns={'level_0':'name', 'level_1': 'datetime'})

checkeddf['datetime'] = pd.to_datetime(checkeddf['datetime']).dt.tz_localize(tz='Europe/Brussels')

checkeddf = checkeddf.set_index(['name', 'datetime'])

# fill na by -1 (no equality operation on nans)
checkeddf = checkeddf.fillna(-1)
dataset.gapfilldf = dataset.gapfilldf.fillna(-1)



test = dataset.gapfilldf[obstype].astype(float).eq(checkeddf[obstype] )
if not test.all():
    print('Gapfill for era debias not equal to manual labels! Here is the difference')
    print(f'manual: {checkeddf}, tlk: {dataset.gapfilldf}')
    print(f'equal at: {test}')
    sys.exit('Error in debias gapfilling')

print('Done! ')

