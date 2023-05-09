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



dataset.import_data_from_file()

dataset.coarsen_time_resolution()

#%%
test = dataset.missing_obs.get_station_missingobs('vlinder03')
print(test)
#%%
# model = dataset.get_modeldata()
# dataset.fill_gaps_linear()
dataset.fill_missing_obs_linear()



# test = dataset.combine_all_to_obsspace(repr_outl_as_nan=True)


# print(dataset.gapfilldf)
#%%
# dataset.fill_missing_obs_linear()
# print(dataset.missing_fill_df)


#%%
# missing_obs = dataset.missing_obs
# missing_obs.interpolate_missing(dataset.df, dataset.metadf['dataset_resolution'])



#%%
# gaps = dataset.gaps

# print(gaps)

# #%%
# method='time'
# obsdf = dataset.df
# obstype = 'temp'
# missing_obs = dataset.missing_obs

# missing_obsspace = missing_obs.get_missing_indx_in_obs_space(dataset.df, dataset.metadf['dataset_resolution'])


# from metobs_toolkit.gap import _find_closes_occuring_date

# import numpy as np

# for staname, missingdt in missing_obsspace:
#     staobs = obsdf.xs(staname, level='name')[obstype]
#     # exclude nan values because they are no good leading/trailing
#     staobs = staobs[~staobs.isnull()]

#     # find leading and trailing datetimes
#     leading_dt = missingdt - timedelta(
#         seconds=_find_closes_occuring_date(
#             refdt = missingdt,
#             series_of_dt = staobs.index,
#             where='before')
#         )
#     trailing_dt = missingdt + timedelta(
#         seconds=_find_closes_occuring_date(
#             refdt = missingdt,
#             series_of_dt = staobs.index,
#             where='after')
#         )

#     # extract the values and combine them in a dataframe
#     leading_val = staobs.loc[leading_dt]
#     trailing_val = staobs.loc[trailing_dt]

#     stadf = pd.DataFrame(
#         index=[leading_dt, missingdt, trailing_dt],
#         data={obstype: [leading_val, np.nan, trailing_val]}
#     )


#     # interpolate the missing obs
#     stadf['interp'] = stadf[obstype].interpolate(
#                                         method=method,
#         )






