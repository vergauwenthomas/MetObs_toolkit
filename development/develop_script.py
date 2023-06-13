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
import pandas as pd
import datetime


# Make an empty dataset
dataset = metobs_toolkit.Dataset()

# Add the demo data files to the dataset settings
dataset.update_settings(input_data_file = metobs_toolkit.demo_datafile,
                        input_metadata_file = metobs_toolkit.demo_metadatafile,
                        data_template_file = metobs_toolkit.demo_template,
                        metadata_template_file = metobs_toolkit.demo_template # Contains also the metadata mapping
                        )

# Load the data from the demo data files
dataset.import_data_from_file()

dataset.coarsen_time_resolution()

dataset.get_landcover()

#%%

anal = dataset.get_analysis()


new_anal = anal.apply_filter('temp < 15.2')

import copy

backup = copy.deepcopy(anal)

#%%

import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np
from metobs_toolkit.analysis import _make_time_derivatives



from scipy.stats import pearsonr
import numpy as np



def get_lc_correlation_matrices(anal, obstype=['temp'], groupby_labels=['hour']):



    # TODO: visualisation ??

    if not isinstance(obstype, list):
        obstype = [obstype]

    # get data
    df = anal.df[obstype].reset_index()
    df = _make_time_derivatives(df, groupby_labels)

    # subset columns
    relev_columns = [label for label in groupby_labels] #to avoid deep copy import
    relev_columns.append('name')
    relev_columns.extend(obstype)
    df = df[relev_columns]

    # find landcover columnnames in the metadf
    lc_columns = [col for col in anal.metadf.columns if (('_' in col ) & (col.endswith('m')))]

    # get landcover data
    lc_df = anal.metadf[lc_columns]

    if lc_df.empty:
        print('WARNING: No landcover columns found in the metadf. Landcover correlations cannot be computed.')
        return None


    # merge together
    df = df.merge(lc_df, how='left', left_on='name', right_index = True)

    # remove name column if it is not explicit in the groupby labels
    if 'name' not in groupby_labels:
        df = df.drop(columns=['name'])

    # create return
    cor_dict = {}

    #Iterate over all groups

    # avoid futur pandas warning for groupby labels of len==1
    if len(groupby_labels) == 1:
        groups = df.groupby(groupby_labels[0])
    else:
        groups = df.groupby(groupby_labels)


    for group_lab, groupdf in groups:

        # drop groupby labels
        groupdf = groupdf.drop(columns=groupby_labels, errors='ignore')

        rho = groupdf.corr(method='pearson')
        pval = groupdf.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(*rho.shape)
        # represent p values by stars
        p_stars = pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x<=t]))

        cor_dict[group_lab] = {'cor matrix': rho,
                               'significance matrix': pval,
                               'combined matrix': rho.astype(str) +' ' +  p_stars}



    return cor_dict

test = get_lc_correlation_matrices(anal,obstype=['temp', 'humidity'], groupby_labels=['hour'])


#%%

import matplotlib.pyplot as plt

cor =test[5]


# make heatmap of cor

fig, ax = plt.subplots(figsize=(8,8))
im = ax.imshow(cor['cor matrix'], interpolation='nearest')
fig.colorbar(im, orientation='vertical', fraction = 0.05)

# Show all ticks and label them with the dataframe column name
ax.set_xticklabels(cor['cor matrix'].columns, rotation=65, fontsize=15)
ax.set_yticklabels(cor['cor matrix'].index, rotation=0, fontsize=15)

# Loop over data dimensions and create text annotations
for i in range(len(cor['cor matrix'].columns)-1):
    for j in range(len(cor['cor matrix'].index)-1):
        text = ax.text(j, i,
                       f'{cor["combined matrix"].to_numpy()[i, j]:.2f}',
                       # round(pear_corr.to_numpy()[i, j], 2),
                       # ha="center", va="center", color="black",
                       )

plt.show()





