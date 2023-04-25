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
testdatafile = os.path.join(str(lib_folder), 'tests', 'test_data',  'testdata_breaking.csv')


template = os.path.join(str(lib_folder), 'tests', 'test_data',  'template_breaking.csv')

static_data = os.path.join(
    str(lib_folder), 'static_data', 'vlinder_metadata.csv')



# #% Setup dataset



raw_dataset = metobs_toolkit.Dataset()
raw_dataset.update_settings(input_data_file=testdatafile,
                        # input_metadata_file=static_data,
                        data_template_file= template,
                        output_folder='/home/thoverga/Documents/VLINDER_github/vlinder_toolkit'
                        )

# dataset.apply_quality_control()


raw_dataset.import_data_from_file()

raw_dataset.apply_quality_control()

test = raw_dataset.combine_all_to_obsspace()


#%%
import pandas as pd
def init_triple_multiindex():
    my_index = pd.MultiIndex(levels=[['name'],['datetime'],['obstype']],
                             codes=[[],[],[]],
                             names=[u'name', u'datetime', u'obstype'])
    return my_index
def init_triple_multiindexdf():
    return pd.DataFrame(index=init_triple_multiindex())


def init_multiindex():
     return pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])
def init_multiindexdf():
    return pd.DataFrame(index = init_multiindex())


#%%

dataset = raw_dataset



#%%
outliersdf = dataset.outliersdf
#remove duplicate indixes (needed for update)
outliersdf = outliersdf[~outliersdf.index.duplicated(keep='first')]

outliersdf_values = outliersdf['value'].unstack()






#%%
import numpy as np
repr_outl_as_nan = False
dataset._qc_checked_obstypes = []


outliersdf = dataset.outliersdf


#remove duplicate indixes (needed for update)
outliersdf = outliersdf[~outliersdf.index.duplicated(keep='first')]


if not outliersdf.empty:
    outliersdf_values = outliersdf['value'].unstack() # for later use
    # convert to wide df with labels
    outliersdf = outliersdf['label'].unstack()

    # convert to final label names for columns
    outliersdf = outliersdf.rename(columns={col: col+'_final_label' for col in outliersdf.columns})


else:
    outliersdf = init_multiindexdf()
    outliersdf_values = init_multiindexdf() # for later use
    outliercolumns = [col+'_final_label' for col in dataset.df if col in dataset.settings.app['observation_types']]
    for column in outliercolumns:
        outliersdf[column] = 'not checked'



# =============================================================================
# Combine observations and outliers
# =============================================================================
# get observations
df = dataset.df


# 1. Merge the label columns
df_and_outl = df.merge(outliersdf, how='outer', left_index=True, right_index=True)

# 2. fill the missing labels

# split between obstype that are checked by qc and obstypes that are not checked

checked_cols = [col+'_final_label' for col in dataset._qc_checked_obstypes]
not_checked_cols = [col for col in df_and_outl.columns if ((col.endswith('_final_label')) and (not col in checked_cols))]

# if obstype checked, and value is nan --> label ok
df_and_outl[checked_cols] = df_and_outl[checked_cols].fillna('ok')
# if obstype is not checked and label is missing --> label 'not checked'
df_and_outl[not_checked_cols] = df_and_outl[not_checked_cols].fillna('not checked')


# 3. Update the values if needed
if not repr_outl_as_nan:

# Merge obs and outliers, where obs values will be updated by outliers
    df_and_outl.update(other=outliersdf_values,
                 join='left',
                 overwrite=True,
                 errors='ignore')


# =============================================================================
# Make gaps, gapsfill and missing dataframes
# =============================================================================

# add gaps observations and fill with default values
gapsidx = dataset.gaps.get_gaps_indx_in_obs_space(
       dataset.df, dataset.outliersdf, dataset.metadf['dataset_resolution'])
gapsdf = gapsidx.to_frame()

# add missing observations if they occure in observation space
missingidx = dataset.missing_obs.get_missing_indx_in_obs_space(
        dataset.df, dataset.metadf['dataset_resolution'])
missingdf = missingidx.to_frame()

# add gapfill and remove the filled records from gaps
gapsfilldf = dataset.gapfilldf.copy()

gapsdf = gapsdf.drop(gapsfilldf.index, errors='ignore')




# initiate default values
for col in df_and_outl.columns:
    if col in dataset.settings.app['observation_types']:
        default_value_gap = np.nan  # nan for observations
        default_value_missing = np.nan

    elif col.endswith('_final_label'):
        # 'gap' for final label
        default_value_gap = dataset.settings.gap['gaps_info']['gap']['outlier_flag']
        # 'is_missing_timestamp' for final label
        default_value_missing = dataset.settings.gap['gaps_info']['missing_timestamp']['outlier_flag']

    else:
        default_value_gap = 'not checked'
        default_value_missing = 'not checked'
        gapsfilldf[col] = 'not checked'

    gapsdf[col] = default_value_gap
    missingdf[col] = default_value_missing

# sort columns
gapsdf = gapsdf[list(df_and_outl.columns)]
missingdf = missingdf[list(df_and_outl.columns)]


# Merge all together
comb_df = pd.concat([df_and_outl, gapsdf, missingdf, gapsfilldf]).sort_index()

test = comb_df



