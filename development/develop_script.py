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
dataset.get_lcz()




#%%
from metobs_toolkit.analysis import Analysis
an = Analysis(obsdf = dataset.df,
              metadf = dataset.metadf,
              settings = dataset.settings)


df = an.get_diurnal_statistics(refstation='vlinder08', errorbands=False, colorby='lcz')

#%%

# import pandas as pd
# import matplotlib.pyplot as plt
# import matplotlib

# names = an.df.index.get_level_values('name').to_series()

# # Define a Matplotlib colormap
# cmap = matplotlib.colormaps['viridis']

# copper = matplotlib.colormaps['copper'].resampled(names.shape[0])

# cmap2 = LinearSegmentedColormap.from_list("mycmap", list(zip(nodes, colors)))



#%%

# lczdict = {'vlinder01': 'Low plants (LCZ D)',
#  'vlinder02': 'Open midrise',
#  'vlinder03': 'Open midrise',
#  'vlinder04': 'Sparsely built',
#  'vlinder05': 'Water (LCZ G)',
#  'vlinder06': 'Scattered Trees (LCZ B)',
#  'vlinder07': 'Compact midrise',
#  'vlinder08': 'Compact midrise',
#  'vlinder09': 'Scattered Trees (LCZ B)',
#  'vlinder10': 'Compact midrise',
#  'vlinder11': 'Open lowrise',
#  'vlinder12': 'Open highrise',
#  'vlinder13': 'Compact midrise',
#  'vlinder14': 'Low plants (LCZ D)',
#  'vlinder15': 'Sparsely built',
#  'vlinder16': 'Water (LCZ G)',
#  'vlinder17': 'Scattered Trees (LCZ B)',
#  'vlinder18': 'Low plants (LCZ D)',
#  'vlinder19': 'Compact midrise',
#  'vlinder20': 'Compact midrise',
#  'vlinder21': 'Sparsely built',
#  'vlinder22': 'Low plants (LCZ D)',
#  'vlinder23': 'Low plants (LCZ D)',
#  'vlinder24': 'Dense Trees (LCZ A)',
#  'vlinder25': 'Water (LCZ G)',
#  'vlinder26': 'Open midrise',
#  'vlinder27': 'Compact midrise',
#  'vlinder28': 'Open lowrise'}

# present_lczs = list(set(lczdict.values()))
# present_lczs.sort()


# cmapname = 'tab20'
# cmapname2 = 'viridis'



# import matplotlib




def make_cat_colormapper(catlist, cmapname):
    """
    Create a dictionary {cat : color} for a list of categorical values.

    If the colormap has more colors than the catlist, optimal color distance is
    done. If a colormap has less colors than unique categories, the categories are grourped.

    Parameters
    ----------
    catlist : list
        List of categorical values.
    cmapname : str
        Matplotlib.colormaps name.

    Returns
    -------
    colordict : dict
        {cat: color} where the color is a RGBalpha tuple.

    """

    catlist = list(set(catlist)) #get unique categories

    cmap = matplotlib.colormaps[cmapname]

    # check number of colors in the cmap
    if cmap.N < len(catlist):
        print(f'Warning: colormap: {cmapname}, is not well suited to color {len(catlist)} categories. ')
        same_col_n_groups = np.ceil(len(catlist) / cmap.N)

        # group cateogries and color them by group
        colordict = {}
        col_idx = -1
        _cat_index = 0
        for cat in catlist:
            if _cat_index%same_col_n_groups == 0:
                col_idx += 1
            colordict[cat] = cmap(int(col_idx))
            _cat_index += 1
        return colordict

    # check if the colormap can be decreased (and thus increasing the colordistance)
    num_increase = np.floor(cmap.N / len(catlist))
    print('num inc: ', num_increase)
    i = 0
    colordict = {}
    for cat in catlist:
        print('cat: ', cat, ' i: ', i)
        colordict[cat] = cmap(int(i))
        i = i + num_increase
    return colordict



# test = make_cat_colormapper(an.metadf['lcz'].to_list(), 'viridis')



# print(test)


#%%


# #%% Make plot
# import matplotlib.pyplot as plt

# refstation = 'vlinder08'


# # Make title
# if title is None:
#     if start


#     if refstation is None:
#         title = "Diurnal cycle "




# fig, ax = plt.subplots
