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


#metobs_toolkit.build_template_prompt(debug=True)

from datetime import datetime


dataset = metobs_toolkit.Dataset()
dataset.update_settings(input_data_file=metobs_toolkit.demo_datafile,
                        input_metadata_file=metobs_toolkit.demo_metadatafile,
                        data_template_file=metobs_toolkit.demo_template,
                        metadata_template_file=metobs_toolkit.demo_template)

dataset.import_data_from_file()
dataset.coarsen_time_resolution()


#%%


model_data = metobs_toolkit.Modeldata("ERA5")
# model_data.get_ERA5_data(metadf=dataset.metadf, startdt=datetime(2022, 9, 1), enddt=datetime(2022, 9, 2))

model_data.set_model_from_csv('/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/era5_data.csv')




#%%
# test = model_data.make_plot(dataset=dataset, starttime=datetime(2022, 9, 1), endtime=datetime(2022, 9, 2))

from metobs_toolkit.df_helpers import multiindexdf_datetime_subsetting, xs_save
from metobs_toolkit.plotting_functions import timeseries_plot, model_timeseries_plot


def make_plot(
    modeldata,
    obstype="temp",
    dataset = None,
    colorby = 'name',
    stationnames=None,
    starttime=None,
    endtime=None,
    title=None,
    show_outliers=True,
    show_filled=True,
    legend=True,
    _ax=None, #needed for GUI, not recommended use
    ):
    """
    This function creates a timeseries plot for the model data. The variable observation type
    is plotted for all stationnames from a starttime to an endtime.

    All styling attributes are extracted from the Settings.

    Parameters
    ----------

    stationnames : list, optional
        A list with stationnames to include in the timeseries. If None is given, all the stations are used, defaults to None.
    obstype : string, optional
         Fieldname to visualise. This can be an observation or station
         attribute. The default is 'temp'.
    starttime : datetime.datetime, optional
         Specifiy the start datetime for the plot. If None is given it will use the start datetime of the dataset, defaults to None.
    endtime : datetime.datetime, optional
         Specifiy the end datetime for the plot. If None is given it will use the end datetime of the dataset, defaults to None.
    title : string, optional
         Title of the figure, if None a default title is generated. The default is None.
    y_label : string, optional
         y-axes label of the figure, if None a default label is generated. The default is None.
    legend : bool, optional
         I True, a legend is added to the plot. The default is True.


    Returns
    -------
    axis : matplotlib.pyplot.axes
         The timeseries axes of the plot is returned.

    """
    # if stationnames is None:
    #     logger.info(f"Make {obstype}-timeseries plot of model data for all stations")
    # else:
    #     logger.info(f"Make {obstype}-timeseries plot of model data for {stationnames}")

    # Basic test
    if obstype not in modeldata.df.columns:
        print(f'ERROR: {obstype} is not foud in the modeldata df.')
        return
    if modeldata.df.empty:
        print('ERROR: The modeldata is empty.')
        return
    if (not dataset is None):
        if (obstype not in dataset.df.columns):
            print(f'ERROR: {obstype} is not foud in the Dataframe df.')
            return


    model_df = modeldata.df

    # ------ filter model ------------

    # Filter on obstype
    model_df = model_df[[obstype]]

    # Subset on stationnames
    if not stationnames is None:
        model_df = model_df[model_df.index.get_level_values('name').isin(stationnames)]

    # Subset on start and endtime
    model_df = multiindexdf_datetime_subsetting(model_df, starttime, endtime)


    #  -------- Filter dataset (if available) -----------
    # combine all dataframes
    mergedf = dataset.combine_all_to_obsspace()

    # subset to obstype
    mergedf = xs_save(mergedf, obstype, level='obstype')

    # Subset on stationnames
    if not stationnames is None:
        mergedf = mergedf[mergedf.index.get_level_values('name').isin(stationnames)]

    # Subset on start and endtime
    mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

    # remove outliers if required
    if not show_outliers:
        outlier_labels = [var['outlier_flag'] for var in dataset.settings.qc['qc_checks_info'].values()]
        mergedf = mergedf[~mergedf['label'].isin(outlier_labels)]

    # remove filled values if required
    if not show_filled:
        fill_labels = ['gap fill', 'missing observation fill'] #toolkit representation labels
        mergedf = mergedf[~mergedf['toolkit_representation'].isin(fill_labels)]







    # Generate ylabel
    try:
        model_true_field_name = modeldata.mapinfo[modeldata.modelname]['band_of_use'][obstype]['name']
    except KeyError:
        print (f'No model field name found for {obstype} in {modeldata}.')
        model_true_field_name = 'Unknown fieldname'
    y_label = f'{model_true_field_name}'

    if not dataset is None:
        dataset_obs_orig_name = dataset.data_template[obstype]['orig_name']
        units = dataset.data_template[obstype]['units']

        y_label = f'{y_label} \n {dataset_obs_orig_name} ({units})'




    #Generate title
    title = f'{modeldata.modelname} : {model_true_field_name}'
    if not dataset is None:
        title = f'{title} and {dataset_obs_orig_name} observations.'


    # make plot of the observations
    if not dataset is None:
        # make plot of the observations
        ax = dataset.make_plot(stationnames=stationnames,
                               colorby=colorby,
                               starttime=starttime,
                               endtime=endtime,
                               title=title,
                               y_label=y_label,
                               legend=False,
                               _ax=None)

        # Make plot of the model on the previous axes
        ax = model_timeseries_plot(
                                df=model_df,
                                obstype=obstype,
                                title=title,
                                ylabel=y_label,
                                show_legend=legend,
                                settings = modeldata._settings,
                                _ax = ax
        )
    else:

        # Make plot of model on empty axes
        ax = model_timeseries_plot(
            df=model_df,
            obstype=obstype,
            title=title,
            ylabel=y_label,
            show_legend=legend,
            settings = modeldata._settings,
            _ax = None
        )

    return ax




make_plot(modeldata = model_data,
          obstype='temp',
          dataset = dataset,
          )
