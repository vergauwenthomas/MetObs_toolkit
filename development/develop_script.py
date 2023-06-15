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

tmp_pickle=os.path.join(lib_folder, 'development', 'tmp', 'dev_pickle.pkl')

import metobs_toolkit


#%%

# # data
# era5_congo_file = '/home/thoverga/Downloads/era5_data_kongo.csv'
data_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/ATHTS01_all.csv'
# metadata_file = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Kobe/CONGO_meta.csv'
template_file ='/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/testdata_testday/Ian/template.csv'
#%%

dataset = metobs_toolkit.Dataset()

dataset.update_settings(output_folder=None,
                        input_data_file=data_file,
                        # input_metadata_file=metobs_toolkit.demo_metadatafile,
                        data_template_file=template_file,
                        # metadata_template_file=metobs_toolkit.demo_template,
                        )


dataset.import_data_from_file()

# dataset.coarsen_time_resolution()

# dataset.apply_quality_control()


# dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)



#%%




dataset.update_qc_settings(obstype='temp',
                           step_max_decrease_per_sec=0.5,
                           step_max_increase_per_sec=0.5)


dataset.apply_quality_control(
    obstype="temp",         # choose which observations you want to check
    gross_value=False,       # set True if you want to perform the gross value check
    persistance=False,       # set True if you want to perform the persistence check
    step=True,              # set True if you want to perform the spike check
    window_variation=False,  # set True if you want to perform the window variation check
)


# qc_statistics = dataset.get_qc_stats(
#     obstype="temp",     # Specify which observation variable you want to get the statistics for; here we choose temperature
#     stationname=None,   # None means all stations are plotted. You can also plot a specific station by saying 'station_A'
#     make_plot=True,     # Set True to make a plot
# )

#%%

dataset.make_plot(colorby='label')

#%%



#%%

# import pandas as pd
# import matplotlib.pyplot as plt

# combdf = dataset.combine_all_to_obsspace()
# combdf = combdf.xs('temp', level='obstype')
# mergedf = combdf[~combdf.index.duplicated()]
# init_idx = mergedf.index

# outl_groups = mergedf.groupby('label')

# outl_label ='ok'
# groupdf = outl_groups.get_group(outl_label)
# outl_color = 'blue'
# plot_settings = dataset.settings.app['plot_settings']

# #%%









# # add init_idx andf fill with nans (to avoid matplotlib interpolation)
# fill_idx = init_idx.to_frame().drop(groupdf.index)
# groupdf = pd.concat([groupdf, fill_idx])
# groupdf = groupdf.drop(columns=["name", "datetime"], errors="ignore")
# groupdf.sort_index()

# plotdf = groupdf.reset_index().pivot(
#     index="datetime", columns="name", values='value'
# )  # long to wide

# for col in plotdf.columns:

#     plotdf[[col]].dropna().plot(
#     kind="line",
#     color=outl_color,
#     # ax=ax,
#     legend=False,
#     zorder=plot_settings["time_series"]["linezorder"],
#     linewidth=plot_settings["time_series"]["linewidth"],
#   )






