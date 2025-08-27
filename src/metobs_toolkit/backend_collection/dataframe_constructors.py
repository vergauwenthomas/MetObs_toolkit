#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DataFrame constructor methods for MetObs Toolkit classes.

This module contains centralized DataFrame construction methods that were
originally property methods in various classes throughout the toolkit.

Created on Aug 21, 2025

"""

import pandas as pd
from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.settings_collection import label_def
from metobs_toolkit.backend_collection.loggingmodule import log_entry


@log_entry
def analysis_df(analysis_instance) -> pd.DataFrame:
    """
    Analysis DataFrame constructor.

    Returns the full DataFrame without the time derivatives.

    Parameters
    ----------
    analysis_instance : Analysis
        The Analysis instance to construct the DataFrame for.

    Returns
    -------
    pd.DataFrame
        DataFrame indexed by ['datetime', 'name'] and containing
        observation columns.
    """
    return analysis_instance.fulldf.set_index(["datetime", "name"])[
        analysis_instance._df_cols
    ]


@log_entry
def dataset_df(dataset_instance) -> pd.DataFrame:
    """
    Dataset DataFrame constructor.

    Construct a DataFrame representation of all observations in the dataset.

    Parameters
    ----------
    dataset_instance : Dataset
        The Dataset instance to construct the DataFrame for.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame indexed by ['datetime', 'obstype', 'name'] with
        columns ['value', 'label'].
    """
    concatlist = []
    for sta in dataset_instance.stations:
        stadf = sta.df.reset_index()
        if stadf.empty:
            continue
        stadf["name"] = sta.name
        concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

    combdf = save_concat((concatlist))
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["value", "label"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["datetime", "obstype", "name"],
            ),
        )
    return combdf


@log_entry
def station_df(station_instance) -> pd.DataFrame:
    """
    Station DataFrame constructor.

    Construct a DataFrame representation of the observations.

    Parameters
    ----------
    station_instance : Station
        The Station instance to construct the DataFrame for.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame with a single column 'value'.
    """
    # return dataframe with ['datetime', 'obstype'] as index and 'value' as single column.
    concatdf = save_concat(
        ([sensor.df for sensor in station_instance.sensordata.values()])
    )

    # sort by datetime
    concatdf.sort_index(inplace=True)
    return concatdf


@log_entry
def gap_df(gap_instance) -> pd.DataFrame:
    """
    Gap DataFrame constructor.

    Return a DataFrame representation of the gap.

    Parameters
    ----------
    gap_instance : Gap
        The Gap instance to construct the DataFrame for.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns ['value', 'label', 'details'].
    """
    return pd.DataFrame(
        {
            "value": gap_instance.records,
            "label": gap_instance._labels,
            "details": gap_instance._extra_info,
        }
    )


@log_entry
def sensordata_df(sensordata_instance) -> pd.DataFrame:
    """
    SensorData DataFrame constructor.

    Return a DataFrame of the sensor records.

    Parameters
    ----------
    sensordata_instance : SensorData
        The SensorData instance to construct the DataFrame for.

    Returns
    -------
    pd.DataFrame
        A DataFrame indexed by ['datetime', 'obstype'] with columns ['value', 'label'].
    """
    import logging

    logger = logging.getLogger("<metobs_toolkit>")

    logger.debug(
        "Creating DataFrame from SensorData series for %s",
        sensordata_instance.stationname,
    )
    # get all records
    df = (
        sensordata_instance.series.to_frame()
        .rename(
            columns={
                sensordata_instance.obstype.name: "value",
                sensordata_instance.stationname: "value",
            }
        )
        .assign(label=label_def["goodrecord"]["label"])
    )

    outliersdf = sensordata_instance.outliersdf[["value", "label"]]

    gapsdf = sensordata_instance.gapsdf[["value", "label"]]

    # concat all together (do not change order)
    to_concat = [df]
    if not outliersdf.empty:
        to_concat.append(outliersdf)
    if not gapsdf.empty:
        to_concat.append(gapsdf)
    df = save_concat((to_concat))
    # remove duplicates
    df = df[~df.index.duplicated(keep="last")].sort_index()

    # add 'obstype' as index
    df = (
        df.assign(obstype=sensordata_instance.obstype.name)
        .reset_index()
        .set_index(["datetime", "obstype"])
    )

    return df


@log_entry
def modeltimeseries_df(modeltimeseries_instance) -> pd.DataFrame:
    """
    ModelTimeSeries DataFrame constructor.

    Return all records as a DataFrame.

    Parameters
    ----------
    modeltimeseries_instance : ModelTimeSeries
        The ModelTimeSeries instance to construct the DataFrame for.

    Returns
    -------
    pd.DataFrame
        A DataFrame with columns ['value', 'model'].
    """
    # get all records
    df = (
        modeltimeseries_instance.series.to_frame()
        .rename(
            columns={
                modeltimeseries_instance.modelobstype.name: "value",
                modeltimeseries_instance.stationname: "value",
            }
        )
        .assign(model=modeltimeseries_instance.modelname)
    )
    return df
