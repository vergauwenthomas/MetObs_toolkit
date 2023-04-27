#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:44:54 2022

@author: thoverga
"""

import sys
import pandas as pd
import numpy as np
import logging


from metobs_toolkit.df_helpers import (
    init_multiindex,
    init_multiindexdf,
    init_triple_multiindexdf,
)


logger = logging.getLogger(__name__)


# =============================================================================
# Helper functions
# =============================================================================


def make_outlier_df_for_check(
    station_dt_list, obsdf, obstype, flag, stationname=None, datetimelist=None
):
    """
    V5
    Helper function to create an outlier dataframe for the given station(s) and datetimes. This will be returned by
    a quality control check and later added to the dastes.outlierdf.

    Multiple commum inputstructures can be handles

    A multiindex dataframe with the relevant observationtypes i.e. the values_in_dict and a specific quality flag column (i.g. the labels) is returned.

    Parameters

    station_dt_list : MultiIndex or list of tuples: (name, datetime)
        The stations with corresponding datetimes that are labeled as outliers.
    values_in_dict : Dictionary
        A Dictionary {columnname: pd.Series()} with the values on which this check is applied.
    flagcolumnname : String
        Name of the labels column
    flag : String
        The label for all the outliers.
    stationname : String, optional
        It is possible to give the name of one station. The default is None.
    datetimelist : DatetimeIndex or List, optional
        The outlier timestamps for the stationname. The default is None.

    Returns

    check_outliers : pd.Dataframe
        A multiindex (name -- datetime) datatframe with one column (columname) and all
        values are the flag.
    """

    if isinstance(station_dt_list, pd.MultiIndex):
        multi_idx = station_dt_list

    elif isinstance(station_dt_list, list):  # list of tuples: (name, datetime)
        multi_idx = pd.MultiIndex.from_tuples(
            station_dt_list, names=["name", "datetime"]
        )
    elif not isinstance(stationname, type(None)):
        if isinstance(datetimelist, pd.DatetimeIndex):
            datetimelist = datetimelist.to_list()
        if isinstance(datetimelist, list):
            indexarrays = list(zip([stationname] * len(datetimelist), datetimelist))
            multi_idx = pd.MultiIndex.from_tuples(
                indexarrays, names=["name", "datetime"]
            )

        else:
            sys.exit(f"Type of datetimelist: {type(datetimelist)} is not implemented.")

    # subset outliers
    outliersdf = obsdf.loc[multi_idx]

    # make the triple multiindex
    outliersdf["obstype"] = obstype
    outliersdf = outliersdf.set_index("obstype", append=True)

    # add flag
    outliersdf["label"] = flag

    # subset columns
    outliersdf = outliersdf[[obstype, "label"]].rename(columns={obstype: "value"})

    # replace values in obsdf by Nan
    obsdf.loc[multi_idx, obstype] = np.nan

    return obsdf, outliersdf


# =============================================================================
# Quality assesment checks on data import
# =============================================================================


def invalid_input_check(df, checks_info):
    checkname = "invalid_input"

    # fast scan wich stations and obstypes have nan outliers
    groups = (
        df.reset_index()
        .groupby("name")
        .apply(lambda x: (np.isnan(x).any()) & (np.isnan(x).all() == False))
    )

    # extract all obstype that have outliers
    outl_obstypes = groups.apply(lambda x: x.any(), axis=0)
    outl_obstypes = outl_obstypes[outl_obstypes].index.to_list()

    # first loop over the smallest sample: outlier obstypes
    outl_dict = {}

    for obstype in outl_obstypes:
        # get stations that have ouliers for this obstype
        outl_stations = groups.loc[groups[obstype], obstype].index.to_list()

        outl_multiidx = init_multiindex()
        for sta in outl_stations:
            # apply check per station
            outl_idx = (
                df.xs(sta, level="name", drop_level=False)[obstype]
                .isnull()
                .loc[lambda x: x]
                .index
            )
            outl_multiidx = outl_multiidx.append(outl_idx)

        outl_dict[obstype] = outl_multiidx

    # create outliersdf for all outliers for all osbtypes
    outl_df = init_multiindexdf()
    for obstype, outliers in outl_dict.items():
        df, specific_outl_df = make_outlier_df_for_check(
            station_dt_list=outliers,
            obsdf=df,
            obstype=obstype,
            flag=checks_info[checkname]["outlier_flag"],
        )
        outl_df = pd.concat([outl_df, specific_outl_df])

    return df, outl_df


def duplicate_timestamp_check(df, checks_info, checks_settings):
    """
    V3
    Looking for duplcate timestaps per station. Duplicated records are removed by the method specified in the qc_settings.

    Parameters

    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)

    Returns

    df : pandas.DataFrame()
        The observations dataframe updated for duplicate timestamps. Duplicated timestamps are removed.

    """

    checkname = "duplicated_timestamp"

    duplicates = pd.Series(
        data=df.index.duplicated(keep=checks_settings[checkname]["keep"]),
        index=df.index,
    )

    if not df.loc[duplicates].empty:
        logging.warning(
            f" Following records are labeld as duplicates: {df.loc[duplicates]}, and are removed"
        )

    # Fill the outlierdf with the duplicates
    outliers = df[df.index.duplicated(keep=checks_settings[checkname]["keep"])]

    # convert values to nan in obsdf
    for obstype in df.columns:
        df.loc[outliers.index, obstype] = np.nan

    # ------- Create a outliersdf -----------#
    # the 'make outliersdf' function cannont be use because of duplicated indices

    outliers = outliers.rename(
        columns={col: "value_" + col for col in outliers.columns}
    )
    outliers = outliers.reset_index()
    outliers["_to_get_unique_idx"] = np.arange(outliers.shape[0])

    outliersdf = pd.wide_to_long(
        df=outliers,
        stubnames="value",
        sep="_",
        suffix=r"\w+",  # to use non-integer suffexes
        i=["name", "datetime", "_to_get_unique_idx"],
        j="obstype",
    )
    # remove the temorary level from the index
    outliersdf = outliersdf.droplevel("_to_get_unique_idx", axis=0)

    # add label column
    outliersdf["label"] = checks_info[checkname]["outlier_flag"]

    # drop duplicates in the obsdf, because this gives a lot of troubles
    # The method does not really mater because the values are set to nan in the observations
    df = df[~df.index.duplicated(keep="first")]

    return df, outliersdf


# =============================================================================
# Quality assesment checks on dataset
# =============================================================================


def gross_value_check(obsdf, obstype, checks_info, checks_settings):
    """
    Looking for values of an observation type that are not physical. These values are labeled and the physical limits are specified in the qc_settings.

    Parameters

    input_series : pandas.Series
        The observations series of the dataset object

    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.


    Returns

    updated_obs_series : pandas.Series
        The observations series updated for this check. Observations that didn't pass are removed.

    outlier_df : pandas.DataFrame
        The collection of records flagged as outliers by this check.

    """

    checkname = "gross_value"

    try:
        specific_settings = checks_settings[checkname][obstype]
    except:
        print(f"No {checkname} settings found for obstype={obstype}. Check is skipped!")
        logger.warning(
            f"No {checkname} settings found for obstype={obstype}. Check is skipped!"
        )
        return obsdf, init_multiindexdf()

    # drop outliers from the series (these are Nan's)
    input_series = obsdf[obstype].dropna()

    # find outlier observations as a list of tuples [(name, datetime), (name, datetime)]
    outl_obs = input_series.loc[
        (input_series <= specific_settings["min_value"])
        | (input_series >= specific_settings["max_value"])
    ].index.to_list()

    # make new obsdf and outlierdf
    obsdf, outlier_df = make_outlier_df_for_check(
        station_dt_list=outl_obs,
        obsdf=obsdf,
        obstype=obstype,
        flag=checks_info[checkname]["outlier_flag"],
    )

    return obsdf, outlier_df


def persistance_check(
    station_frequencies, obsdf, obstype, checks_info, checks_settings
):
    """
    V3
    Looking for values of an observation type that do not change during a timewindow. These are flagged as outliers.

    In order to perform this check, at least N observations chould be in that time window.


    Parameters

    input_series : pandas.Series
        The observations series of the dataset object

    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.


        Returns

        updated_obs_series : pandas.Series
            The observations series updated for this check. Observations that didn't pass are removed.

        outlier_df : pandas.DataFrame
            The collection of records flagged as outliers by this check.

    """

    checkname = "persistance"

    try:
        specific_settings = checks_settings[checkname][obstype]
    except:
        print(f"No {checkname} settings found for obstype={obstype}. Check is skipped!")
        logger.warning(
            f"No {checkname} settings found for obstype={obstype}. Check is skipped!"
        )
        return obsdf, init_multiindexdf()

    invalid_windows_check_df = (
        pd.to_timedelta(specific_settings["time_window_to_check"]) / station_frequencies
        < specific_settings["min_num_obs"]
    )
    invalid_stations = list(
        invalid_windows_check_df[invalid_windows_check_df == True].index
    )
    if bool(invalid_stations):
        print(
            f"The windows are too small for stations  {invalid_stations} to perform persistance check"
        )
        logger.info(
            f"The windows are too small for stations  {invalid_stations} to perform persistance check"
        )

    subset_not_used = obsdf[obsdf.index.get_level_values("name").isin(invalid_stations)]
    subset_used = obsdf[~obsdf.index.get_level_values("name").isin(invalid_stations)]

    if not subset_used.empty:
        # drop outliers from the series (these are Nan's)
        input_series = subset_used[obstype].dropna()

        # apply persistance
        def is_unique(
            window,
        ):  # comp order of N (while using the 'unique' function is Nlog(N))
            a = window.values
            a = a[~np.isnan(a)]
            return (a[0] == a).all()

        # TODO: Tis is very expensive if no coarsening is applied !!!! Can we speed this up?
        window_output = (
            input_series.reset_index(level=0)
            .groupby("name")
            .rolling(
                window=specific_settings["time_window_to_check"],
                closed="both",
                center=True,
                min_periods=specific_settings["min_num_obs"],
            )
            .apply(is_unique)
        )

        list_of_outliers = []
        outl_obs = window_output.loc[window_output[obstype] == True].index
        for outlier in outl_obs:
            outliers_list = get_outliers_in_daterange(
                input_series,
                outlier[1],
                outlier[0],
                specific_settings["time_window_to_check"],
                station_frequencies,
            )

            list_of_outliers.extend(outliers_list)

        list_of_outliers = list(set(list_of_outliers))

        # make new obsdf and outlierdf
        subset_used, outlier_df = make_outlier_df_for_check(
            station_dt_list=list_of_outliers,
            obsdf=subset_used,
            obstype=obstype,
            flag=checks_info[checkname]["outlier_flag"],
        )

        obsdf = pd.concat([subset_used, subset_not_used])

        return obsdf, outlier_df

    else:
        obsdf = pd.concat([subset_used, subset_not_used])

        return obsdf, init_multiindexdf()


def repetitions_check(obsdf, obstype, checks_info, checks_settings):
    """
    Looking for values of an observation type that are repeated at least with the frequency specified in the qc_settings. These values are labeled.

    Parameters

    input_series : pandas.Series
        The observations series of the dataset object

    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.


    Returns

    updated_obs_series : pandas.Series
        The observations series updated for this check. Observations that didn't pass are removed.

    outlier_df : pandas.DataFrame
        The collection of records flagged as outliers by this check.


    """

    checkname = "repetitions"

    try:
        specific_settings = checks_settings[checkname][obstype]
    except:
        print(f"No {checkname} settings found for obstype={obstype}. Check is skipped!")
        logger.warning(
            f"No {checkname} settings found for obstype={obstype}. Check is skipped!"
        )
        return obsdf, init_multiindexdf()

    # drop outliers from the series (these are Nan's)
    input_series = obsdf[obstype].dropna()

    # find outlier datetimes

    # add time interval between two consecutive records, group by consecutive records without missing records

    time_diff = input_series.index.get_level_values("datetime").to_series().diff()
    time_diff.index = input_series.index  # back to multiindex

    persistance_filter = ((input_series.shift() != input_series)).cumsum()

    grouped = input_series.groupby(["name", persistance_filter])
    # the above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = grouped.size()
    outlier_groups = group_sizes[
        group_sizes > specific_settings["max_valid_repetitions"]
    ]

    # add to outl_obs.
    outl_obs = []
    for group_idx in outlier_groups.index:
        groupseries = grouped.get_group(group_idx)
        if len(set(groupseries)) == 1:  # Check if all observations are equal in group
            outl_obs.extend(groupseries.index.to_list())

    # make new obsdf and outlierdf
    obsdf, outlier_df = make_outlier_df_for_check(
        station_dt_list=outl_obs,
        obsdf=obsdf,
        obstype=obstype,
        flag=checks_info[checkname]["outlier_flag"],
    )

    return obsdf, outlier_df


def step_check(obsdf, obstype, checks_info, checks_settings):
    """

    V3
    Looking for jumps of the values of an observation type that are larger than the limit specified in the qc_settings. These values are removed from
    the input series and combined in the outlier df.

    The purpose of this check is to flag observations with a value that is too much different compared to the previous (not flagged) recorded value.

    Parameters

    input_series : pandas.Series
        The observations series of the dataset object

    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.



     Returns

     updated_obs_series : pandas.Series
         The observations series updated for this check. Observations that didn't pass are removed.

     outlier_df : pandas.DataFrame
         The collection of records flagged as outliers by this check.

    """

    checkname = "step"

    try:
        specific_settings = checks_settings[checkname][obstype]
    except:
        print(f"No {checkname} settings found for obstype={obstype}. Check is skipped!")
        logger.warning(
            f"No {checkname} settings found for obstype={obstype}. Check is skipped!"
        )
        return obsdf, init_multiindexdf()

    # drop outliers from the series (these are Nan's)
    input_series = obsdf[obstype].dropna()

    list_of_outliers = []

    for name in input_series.index.droplevel("datetime").unique():
        subdata = input_series.xs(name, level="name", drop_level=False)

        time_diff = subdata.index.get_level_values("datetime").to_series().diff()
        time_diff.index = subdata.index  # back to multiindex
        # define filter
        step_filter = (
            (subdata - subdata.shift(1))
            > (
                specific_settings["max_increase_per_second"]
                * time_diff.dt.total_seconds()
            )
        ) | (
            (subdata - subdata.shift(1))
            < (
                specific_settings["max_decrease_per_second"]
                * time_diff.dt.total_seconds()
            )
        )  # &
        # (time_diff == station_frequencies[name]))
        outl_obs = step_filter[step_filter == True].index

        list_of_outliers.extend(outl_obs)

    # make new obsdf and outlierdf
    obsdf, outlier_df = make_outlier_df_for_check(
        station_dt_list=list_of_outliers,
        obsdf=obsdf,
        obstype=obstype,
        flag=checks_info[checkname]["outlier_flag"],
    )

    return obsdf, outlier_df


def window_variation_check(
    station_frequencies, obsdf, obstype, checks_info, checks_settings
):
    """

    V3
    Looking for jumps of the values of an observation type that are larger than the limit specified in the qc_settings. These values are removed from
    the input series and combined in the outlier df.

    There is a increament threshold (that is if there is a max value difference and the maximum value occured later than the minimum value occured.)
    And vice versa is there a decreament threshold.

    The check is only applied if there are at leas N observations in the time window.


    Parameters

    station_frequencies: pandas.Series
        The series containing the dataset time-resolution for all stations (as index).

    input_series : pandas.Series
        The observations series of the dataset object

    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.



     Returns

     updated_obs_series : pandas.Series
         The observations series updated for this check. Observations that didn't pass are removed.

     outlier_df : pandas.DataFrame
         The collection of records flagged as outliers by this check.

    """
    checkname = "window_variation"

    try:
        specific_settings = checks_settings[checkname][obstype]
    except:
        print(f"No {checkname} settings found for obstype={obstype}. Check is skipped!")
        logger.warning(
            f"No {checkname} settings found for obstype={obstype}. Check is skipped!"
        )
        return obsdf, init_multiindexdf()

    invalid_windows_check_df = (
        pd.to_timedelta(specific_settings["time_window_to_check"]) / station_frequencies
        < specific_settings["min_window_members"]
    )
    invalid_stations = list(
        invalid_windows_check_df[invalid_windows_check_df == True].index
    )
    if bool(invalid_stations):
        print(
            f"The windows are too small for stations  {invalid_stations} to perform window variation check"
        )
        logger.info(
            f"The windows are too small for stations  {invalid_stations} to perform window variation check"
        )

    subset_not_used = obsdf[obsdf.index.get_level_values("name").isin(invalid_stations)]
    subset_used = obsdf[~obsdf.index.get_level_values("name").isin(invalid_stations)]

    if not subset_used.empty:
        # drop outliers from the series (these are Nan's)
        input_series = subset_used[obstype].dropna()

        # Calculate window thresholds (by linear extarpolation)
        windowsize_seconds = pd.Timedelta(
            specific_settings["time_window_to_check"]
        ).total_seconds()
        max_window_increase = (
            specific_settings["max_increase_per_second"] * windowsize_seconds
        )
        max_window_decrease = (
            specific_settings["max_decrease_per_second"] * windowsize_seconds
        )

        # apply steptest
        def variation_test(window):
            if (max(window) - min(window) > max_window_increase) & (
                window.idxmax() > window.idxmin()
            ):
                return 1

            if (max(window) - min(window) > max_window_decrease) & (
                window.idxmax() < window.idxmin()
            ):
                return 1
            else:
                return 0

        window_output = (
            input_series.reset_index(level=0)
            .groupby("name")
            .rolling(
                window=specific_settings["time_window_to_check"],
                closed="both",
                center=True,
                min_periods=specific_settings["min_window_members"],
            )
            .apply(variation_test)
        )

        list_of_outliers = []
        outl_obs = window_output.loc[window_output[obstype] == 1].index

        for outlier in outl_obs:
            outliers_list = get_outliers_in_daterange(
                input_series,
                outlier[1],
                outlier[0],
                specific_settings["time_window_to_check"],
                station_frequencies,
            )

            list_of_outliers.extend(outliers_list)

        list_of_outliers = list(set(list_of_outliers))

        # make new obsdf and outlierdf
        subset_used, outlier_df = make_outlier_df_for_check(
            station_dt_list=list_of_outliers,
            obsdf=subset_used,
            obstype=obstype,
            flag=checks_info[checkname]["outlier_flag"],
        )

        obsdf = pd.concat([subset_used, subset_not_used])

        return obsdf, outlier_df

    else:
        obsdf = pd.concat([subset_used, subset_not_used])

        return obsdf, init_multiindexdf()


def get_outliers_in_daterange(input_data, date, name, time_window, station_freq):
    end_date = date + (pd.Timedelta(time_window) / 2).floor(station_freq[name])
    start_date = date - (pd.Timedelta(time_window) / 2).floor(station_freq[name])

    daterange = pd.date_range(start=start_date, end=end_date, freq=station_freq[name])

    multi_idx = pd.MultiIndex.from_arrays(
        arrays=[[name] * len(daterange), daterange.to_list()],
        sortorder=1,
        names=["name", "datetime"],
    )
    outlier_sub_df = pd.DataFrame(data=None, index=multi_idx, columns=None)

    intersection = outlier_sub_df.index.intersection(input_data.dropna().index).values

    return intersection
