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


#%%

import matplotlib.pyplot as plt
# import seaborn as sns
import numpy as np


hour = 5
obstype = 'temp'
buffer_rad = 100


def filter_data(df, metadf, quarry_str):
    """
    Function to filter a dataframe by a user definde string expression. This
    can be used to filter the observation to specific meteorological conditions
    (i.e. low windspeeds, high humidity, cold temperatures, ...)

    The filter expression contains only columns present in the df and/or the
    metadf.

    The filtered df and metadf are returned

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing all the observations to be filterd.
    metadf : pandas.DataFrame
        The dataframe containig all the metadata per station.
    quarry_str : str
        A filter expression using columnnames present in either df or metadf,
        number and expressions like <, >, ==, >=, *, +, .... Multiple filters
        can be combine to one expression by using & (AND) and | (OR).

    Returns
    -------
    filter_df : pandas.DataFrame
        The filtered df.
    filter_metadf : pandas.DataFrame
        The filtered metadf.

    """


    # save index order and names for reconstruction
    df_init_idx = list(df.index.names)
    metadf_init_idx = list(metadf.index.names)

    # easyer for sperationg them
    df = df.reset_index()
    metadf = metadf.reset_index()


    # save columns orders
    df_init_cols = df.columns
    metadf_init_cols = metadf.columns

    # merge together on name

    mergedf = df.merge(metadf, how='left', on='name')

    #apply filter
    filtered = mergedf.query(expr=quarry_str)

    # split to df and metadf
    filter_df = filtered[df_init_cols]
    filter_metadf = filtered[metadf_init_cols]

    # set indexes
    filter_df = filter_df.set_index(df_init_idx)
    filter_metadf = filter_metadf.set_index(metadf_init_idx)

    return filter_df, filter_metadf




testdf, testmetadf = filter_data(anal.df,
                                 anal.metadf,
                                 '16.2 < temp < 18.3 & humidity < 70.0')



# def get_lc_correlation_matrix()


# # get data
# df = anal.df.reset_index()
# df = df[df['datetime'].dt.hour == hour]
# df = df[[obstype, 'name']]

# # get landcover data
# landcover_cols = [col for col in anal.metadf.columns if col.endswith(f'_{buffer_rad}m')]
# lc_df = anal.metadf[landcover_cols]


# # merge together
# df = df.merge(lc_df, how='left', left_on='name', right_index = True)




# #
# # Correlation between different variables
# #
# corr = df.corr(method='pearson')
# #
# # Set up the matplotlib plot configuration
# #
# f, ax = plt.subplots(figsize=(12, 10))
# #
# # Generate a mask for upper traingle
# #
# mask = np.triu(np.ones_like(corr, dtype=bool))
# #
# # Configure a custom diverging colormap
# #
# cmap = sns.diverging_palette(230, 20, as_cmap=True)
# #
# # Draw the heatmap
# #
# sns.heatmap(corr, annot=True, mask = mask, cmap=cmap)




# from scipy.stats import pearsonr
# import numpy as np
# rho = df.corr()
# pval = df.corr(method=lambda x, y: pearsonr(x, y)[1]) - np.eye(*rho.shape)
# p = pval.applymap(lambda x: ''.join(['*' for t in [.05, .01, .001] if x<=t]))
# rho.round(2).astype(str) + p



