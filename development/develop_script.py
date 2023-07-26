#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit
#
import os
import sys
from pathlib import Path
import pandas as pd
import time
import math


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))

# add testdata paths
from tests.push_test.test_data_paths import testdata




tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit



#%%


# datafile ='C:/Users/andre/Downloads/tal_netatmo_22.csv'

# template_file = 'C:/Users/andre/Downloads/est_netatmo_22_template.csv'

# dataset = metobs_toolkit.Dataset()

# dataset.update_settings(input_data_file=datafile,
#                         data_template_file=template_file)

# dataset.import_data_from_file()
# print(dataset)

#%%
# # use_dataset = 'debug_wide'
# use_dataset = 'single_netatmo_sara_station'
# use_dataset = 'vlindergent2022'
use_dataset = 'siebevlinder'
dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=None,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        data_template_file=testdata[use_dataset]['template'],
                        metadata_template_file=testdata[use_dataset]['template'],
                        )



dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

dataset.coarsen_time_resolution(freq = testdata[use_dataset]['coarsen'])
# dataset.apply_quality_control()
# dataset.get_lcz()
# dataset.update_gaps_and_missing_from_outliers()

#%%
dataset.update_qc_settings(obstype='radiation_temp', rep_max_valid_repetitions=3)
dataset.apply_quality_control(obstype='radiation_temp')

#%%
dataset.make_plot(obstype='radiation_temp', colorby='label')




#%%
dataset.get_qc_stats(obstype='radiation_temp')


df = dataset.df[['radiation_temp']]
outliers = dataset.outliersdf.xs('radiation_temp', level='obstype')




comb_df = dataset.combine_all_to_obsspace()
comb_df = comb_df.xs('radiation_temp', level='obstype')
#%%





#%%
dataset.update_gaps_and_missing_from_outliers()
dataset.fill_gaps_linear()

#%%

# vlinder01 2022-09-07 07:00:00+00:00  16.637895  gap_interpolation
#           2022-09-07 07:20:00+00:00  16.675789  gap_interpolation
#           2022-09-07 07:40:00+00:00  16.713684  gap_interpolation
#           2022-09-07 08:00:00+00:00  16.751579  gap_interpolation


#%%
test = dataset.get_analysis(add_gapfilled_values=True)


dummy = test.df

# # gapsfilled labels
# gapfill_settings =dataset.settings.gap['gaps_fill_info']
# gapfilllabels =[ val for val in gapfill_settings['label'].values()]


# # missingfilled labels
# missingfill_settings =dataset.settings.missing_obs['missing_obs_fill_info']
# missingfilllabels =[ val for val in missingfill_settings['label'].values()]

# fill_labels = gapfilllabels.copy()
# fill_labels.extend(missingfilllabels)


# #%%
# from metobs_toolkit.df_helpers import xs_save


# obstype = 'temp'
# mergedf = dataset.combine_all_to_obsspace()
# # df = xs_save(mergedf, obstype, 'obstype')

# keeplabels = ['ok', 'gap_interpolation']
# df = mergedf[mergedf['label'].isin(keeplabels)]
# df = df[['value']]
# df = df.unstack(level='obstype')
# df.droplevel(level=0, axis=1)










# ana = dataset.get_analysis()



#%%

