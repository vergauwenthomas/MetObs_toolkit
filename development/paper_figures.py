#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 11 13:16:58 2023

@author: thoverga
"""

import os
import logging
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
# handle logs to spyder console
toolkit_logger = metobs_toolkit.loggers
toolkit_logger.addHandler(logging.StreamHandler())


#%% Use the demo data and subset them

datadf = pd.read_csv(metobs_toolkit.demo_datafile, sep=';')
metadf = pd.read_csv(metobs_toolkit.demo_metadatafile, sep=',')

# Subset to regio gent

ghent_stations = [ 'vlinder24', 'vlinder25', 'vlinder05', 'vlinder27',
                  'vlinder02', 'vlinder01', 'vlinder28']


datadf = datadf[datadf['Vlinder'].isin(ghent_stations)]
metadf = metadf[metadf['Vlinder'].isin(ghent_stations)]

# subset period
datadf['dummy_dt'] = datadf['Datum'] + datadf['Tijd (UTC)']
datadf['dummy_dt'] = pd.to_datetime(datadf['dummy_dt'], format='%Y-%m-%d%H:%M:%S')

from datetime import datetime
startdt = datetime(2022, 9, 1)
enddt = datetime(2022, 9, 10)

datadf = datadf[(datadf['dummy_dt'] >= startdt) & (datadf['dummy_dt'] <= enddt)]

datadf = datadf.drop(columns=['dummy_dt'])


# Inducing outliers as demo
datadf = datadf.drop(index=datadf.iloc[180:200, :].index.tolist())







# save in paper folder
folder = '/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/paper_dataset'
datadf.to_csv(os.path.join(folder, 'datafile.csv'))
metadf.to_csv(os.path.join(folder, 'metadatafile.csv'))



#%% Importing raw data
use_dataset = 'paper_dataset'

dataset = metobs_toolkit.Dataset()


dataset.update_settings(output_folder=folder,
                        input_data_file=testdata[use_dataset]['datafile'],
                        input_metadata_file=testdata[use_dataset]['metadatafile'],
                        template_file=testdata[use_dataset]['template'],
                        )

dataset.import_data_from_file(**testdata[use_dataset]['kwargs'])

#%% Make raw data timeseries plot

import matplotlib.pyplot as plt
# ax = dataset.make_plot()


# #translate axes
# ax.set_title('Temperature for all stations')
# ax.set_ylabel('T2m in °C')

#Shape figure and resolution
def shape_timeseries_fig_for_paper(ax, filename):

    fig = ax.get_figure()
    fig.set_size_inches(12, 5, forward=True)
    fig.set_dpi(200)

    # fig.tight_layout()
    figpath = os.path.join(folder, filename)
    plt.savefig(fname=figpath)


# shape_timeseries_fig_for_paper(ax=ax, filename='raw_dataset_timeseries.png')


#%% Make gee interactive spatial plot
# dataset.make_gee_plot(gee_map='worldcover', save=True)


#%% Apply qc and make figures

# coarsen
dataset.coarsen_time_resolution(freq='20T')


#update QC settings
dataset.update_qc_settings(obstype='temp', gapsize_in_records=None,
                           dupl_timestamp_keep=None,
                           persis_time_win_to_check=None,
                           persis_min_num_obs=None,
                           rep_max_valid_repetitions=None,
                           gross_value_min_value=10.7,
                           gross_value_max_value=None,
                           win_var_max_increase_per_sec=None,
                           win_var_max_decrease_per_sec=None,
                           win_var_time_win_to_check=None,
                           win_var_min_num_obs=None,
                           step_max_increase_per_sec=5./3600.,
                           step_max_decrease_per_sec=None)



dataset.update_titan_qc_settings(obstype='temp', buddy_radius=10000,
                                   buddy_num_min=3, buddy_threshold=2.2,
                                   buddy_max_elev_diff=None,
                                   buddy_elev_gradient=None,
                                   buddy_min_std=1.0,
                                   buddy_num_iterations=None,
                                   buddy_debug=None)



dataset.apply_quality_control()
dataset.apply_titan_buddy_check(use_constant_altitude=True)

# change color for printing (avoid yellow!)
dataset.settings.app['plot_settings']['color_mapper']['gross_value'] = "#fc0303"



ax2 = dataset.make_plot(colorby='label')
#translate axes
ax2.set_title('Temperature for all stations')
ax2.set_ylabel('T2m in °C')

shape_timeseries_fig_for_paper(ax=ax2, filename='after_qc.png')


# get overview stats
dataset.get_qc_stats()

fig = plt.gcf()

# fig.tight_layout()
figpath = os.path.join(folder, 'qc_stats.png')
plt.savefig(fname=figpath)



#%% Filling gaps

# 1. Update gaps and missing from outliers
dataset.update_gaps_and_missing_from_outliers(obstype='temp', n_gapsize=6)

# 2. update settings
dataset.update_gap_and_missing_fill_settings(gap_interpolation_method=None,
                                             gap_interpolation_max_consec_fill=None,
                                             gap_debias_prefered_leading_period_hours=24,
                                             gap_debias_prefered_trailing_period_hours=4,
                                             gap_debias_minimum_leading_period_hours=24,
                                             gap_debias_minimum_trailing_period_hours=4,
                                             automatic_max_interpolation_duration_str=None,
                                             missing_obs_interpolation_method=None)

# 3. Get modeldata


# era5 = dataset.get_modeldata(modelname='ERA5_hourly',
#                       modeldata=None, obstype='temp',
#                       stations=None, startdt=None, enddt=None)
# era5.save_modeldata(outputfolder=folder, filename='era.pkl')

dummy_mod = metobs_toolkit.Modeldata('ERA5_hourly')
era5 = dummy_mod.import_modeldata(folder_path=folder,
                                  filename='era.pkl')

# 4. convert units of model
era5.convert_units_to_tlk('temp')

# 5. fill missing obs
dataset.fill_missing_obs_linear()

# 6. fill gaps

dataset.fill_gaps_era5(era5)

# 7. Make plot (of single station for clearity)

ax3 = dataset.get_station('vlinder28').make_plot(colorby='label')


#translate axes
ax3.set_title('Temperature for vlinder28')
ax3.set_ylabel('T2m in °C')

shape_timeseries_fig_for_paper(ax=ax3, filename='after_fill.png')




#%% Create analysis

dataset.get_landcover(buffers=[50, 150, 500], aggregate=True)




ana = dataset.get_analysis(add_gapfilled_values=True)

# 1. Get diurnal analysis
ax4 = ana.get_diurnal_statistics(colorby='name',
                           obstype='temp',
                           stations=None, startdt=None, enddt=None,
                           plot=True,
                           title='Hourly average temperature diurnal cycle',
                           y_label=None, legend=True,
                           errorbands=True, _return_all_stats=False)

fig = plt.gcf()
fig.set_dpi(200)

# fig.tight_layout()
figpath = os.path.join(folder, 'diurnal_cycle.png')
plt.savefig(fname=figpath)



# 2. Get corr heatmap
ana.get_lc_correlation_matrices(obstype=['temp'], groupby_labels=['hour'])
ana.plot_correlation_heatmap(groupby_value=5, title=None)

fig = plt.gcf()
fig.set_dpi(200)

# fig.tight_layout()
figpath = os.path.join(folder, 'corr_heatmap.png')
plt.savefig(fname=figpath)




# 3. evolution of corr
ana.settings.app['plot_settings']['correlation_scatter']['legend_ncols']=5
ax6 = ana.plot_correlation_variation()
fig = ax6.get_figure()
fig.set_dpi(200)

# fig.tight_layout()
figpath = os.path.join(folder, 'corr_scatter.png')
plt.savefig(fname=figpath)

