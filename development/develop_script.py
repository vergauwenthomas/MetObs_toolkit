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


#%%


# print(dataset.outliersdf.shape)
dataset.make_plot(colorby='label')
dataset.update_gaps_and_missing_from_outliers(n_gapsize=10)


#%%

# #%%
dataset.make_plot(colorby='label')
dataset.fill_missing_obs_linear()
dataset.make_plot(colorby='label')
dataset.fill_gaps_linear()
#%%
dataset.make_plot(colorby='label')




todo: fix the timeseries plotting !!!!!
#%%

# dataset.get_station('n13').make_plot(colorby='label')

# #%%
# from metobs_toolkit.df_helpers import xs_save
# from metobs_toolkit.plotting_functions import sorting_function, _all_possible_labels_colormapper
# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd

# # sta = dataset.get_station('n13')
# mergedf = dataset.get_station('n13').combine_all_to_obsspace()
# # mergedf = mergedf.xs('temp', level='obstype')
# # print('a: ', mergedf['label'].value_counts())

# settings = dataset.settings







# #%%
# def make_line_plot_df(df, show_labels = []):
#     plotdf = df[df['label'].isin(show_labels)]

#     hide_df = df[~df['label'].isin(show_labels)]

#     hide_df.loc[hide_df.index, 'value'] = np.nan
#     return pd.concat([plotdf, hide_df]).sort_index()



# #%%

# from matplotlib.lines import Line2D


# plot_settings = settings.app["plot_settings"]
# show_legend = True
# ok_labels = ['ok']

# fill_labels = ['missing_obs_interpolation'] #verder updaten

# missing_labels = ['gap', 'missing timestamp']

# qc_labels = ['repetitions outlier']


# vlin_min = mergedf[mergedf['label'] == 'ok']['value'].min()
# vlin_max = mergedf[mergedf['label'] == 'ok']['value'].max()

# fig, ax = plt.subplots()
# col_mapper = _all_possible_labels_colormapper(settings) # get color mapper



# for sta in mergedf.index.get_level_values('name').unique():
#     stadf = xs_save(mergedf, sta, 'name')

#     # ---- ok obs ------- (green lines)
#     sta_ok_df = make_line_plot_df(df=stadf,
#                                  show_labels=ok_labels)

#     sta_ok_df.plot(
#         kind="line",
#         color=sta_ok_df['label'].map(col_mapper),
#         ax=ax,
#         legend=False,
#         zorder=plot_settings["time_series"]["linezorder"],
#         linewidth=plot_settings["time_series"]["linewidth"],
#         )

#     # ------ fill obs ------ (dashed lines)
#     sta_fill_df = make_line_plot_df(df=stadf,
#                                   show_labels=fill_labels)
#     sta_fill_df.plot(
#         kind="line",
#         style="--",
#         color=sta_fill_df['label'].map(col_mapper),
#         ax=ax,
#         legend=False,
#         zorder=plot_settings["time_series"]["linezorder"],
#         linewidth=plot_settings["time_series"]["linewidth"],
#         )

#  # ------ missing obs ------ (vertical lines)
# missing_df = mergedf[mergedf['label'].isin(missing_labels)]
# missing_df = missing_df.reset_index()
# ax.vlines(x=missing_df['datetime'].to_numpy(),
#           ymin=vlin_min,
#           ymax=vlin_max,
#           linestyle="--",
#           color=missing_df['label'].map(col_mapper),
#           zorder=plot_settings['time_series']["dashedzorder"],
#           linewidth=plot_settings['time_series']["linewidth"])


# # ------ outliers ------ (scatters)
# outlier_df = mergedf[mergedf['label'].isin(qc_labels)]
# outlier_df = outlier_df.reset_index()
# outlier_df.plot(
#     kind="scatter",
#     x="datetime",
#     y='value',
#     ax=ax,
#     color=outlier_df['label'].map(col_mapper),
#     legend=False,
#     zorder=plot_settings["time_series"]["scatterzorder"],
#     s=plot_settings["time_series"]["scattersize"],
# )


# # create legend
# if show_legend:
#     custom_handles = [] #add legend items to it
#     label_vec=[] # add type of label
#     for label in mergedf['label'].unique():
#         outl_color = col_mapper[label]

#         if label in ok_labels:
#             custom_handles.append(
#                 Line2D([0], [0], color=outl_color, label="ok", lw=4))
#             label_vec.append(1)

#         elif label in fill_labels:
#             custom_handles.append(
#                 Line2D([0],[0],
#                     color=outl_color,
#                     label=f"filled value ({label})",
#                     lw=1,
#                     linestyle="--",)
#                 )
#             label_vec.append(2)

#         elif label in missing_labels:
#             custom_handles.append(
#                  Line2D([0],[0],
#                      color=outl_color,
#                      label=f"{label}",
#                      lw=1,
#                      linestyle='--',
#                      linewidth=2,
#                      )
#                  )
#             label_vec.append(3)
#         else:
#             custom_handles.append(
#                 Line2D([0],[0], marker="o", color="w",
#                     markerfacecolor=outl_color,
#                     label=label,
#                     lw=1,)
#                 )
#             label_vec.append(4)


#     custom_handles = sorting_function(label_vec, custom_handles)
#     #ax.legend(handles=custom_handles)
#     box = ax.get_position()
#     ax.set_position([box.x0, box.y0 + box.height * 0.2,
#          box.width, box.height * 0.85])
#     ax.legend(handles=custom_handles, loc='upper center',
#         bbox_to_anchor=(0.5, -0.25),
#         fancybox=True, shadow=True,
#         ncol=plot_settings["time_series"]["legend_n_columns"])


#%%







