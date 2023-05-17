#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of functions on dataframe that are often used.

Created on Thu Mar  2 16:00:59 2023

@author: thoverga
"""

import sys
import pandas as pd
import numpy as np
import geopandas as gpd
import itertools
from metobs_toolkit import observation_types


def init_multiindex():
    return pd.MultiIndex(
        levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
    )


def init_multiindexdf():
    return pd.DataFrame(index=init_multiindex())


def init_triple_multiindex():
    my_index = pd.MultiIndex(
        levels=[["name"], ["datetime"], ["obstype"]],
        codes=[[], [], []],
        names=["name", "datetime", "obstype"],
    )
    return my_index


def init_triple_multiindexdf():
    return pd.DataFrame(index=init_triple_multiindex())


def format_outliersdf_to_doubleidx(outliersdf):
    """
    Convert outliersdf to multiindex dataframe if needed.

    This is applied when the obstype level in the index is not relevant.


    Parameters
    ----------
    ouliersdf : Dataset.outliersdf
        The outliers dataframe to format to name - datetime index.

    Returns
    -------
    pandas.DataFrame
        The outliersdfdataframe where the 'obstype' level is dropped, if it was present.

    """

    if "obstype" in outliersdf.index.names:
        return outliersdf.droplevel("obstype")
    else:
        return outliersdf

def value_labeled_doubleidxdf_to_triple_idxdf(df,value_col_name='value', label_col_name = 'label'):
    """
    This function converts a double index dataframe with an 'obstype' column,
    and a 'obstype_final_label' column to a triple index dataframe where the
    obstype values are added to the index.

    Parameters
    ----------
    df : pd.DataFrame
        Dataframe with ['name', 'datetime'] as index and two columns: [obstype, obstype_final_label].
        Where obstype is an observation type.
    value_col_name : str, optional
        Name of the column for the values. The default is 'value'.
    label_col_name : str, optional
        Name of the column for the labels. The default is 'label'.

    Returns
    -------
    values : pd.DataFrame()
        Dataframe with a ['name', 'datetime', obstype] index and two columnd:
            [value_col_name, label_col_name]

    """
    if df.empty:
        return df

    present_obstypes = [col for col in df.columns if col in observation_types]

    # get all values in triple index form
    values = (df[present_obstypes].stack(dropna=False)
              .reset_index()
              .rename(columns={'level_2': 'obstype', 0: value_col_name})
              .set_index(['name', 'datetime', 'obstype']))


    # make a triple label dataframe
    labelsdf = pd.DataFrame()
    for obstype in present_obstypes:
        subdf = df.loc[:, [obstype+'_final_label']]
        subdf['obstype'] = obstype
        subdf = subdf.reset_index()
        subdf = subdf.set_index(['name', 'datetime', 'obstype'])
        subdf = subdf.rename(columns={obstype+'_final_label': label_col_name})

        labelsdf = pd.concat([labelsdf, subdf])

    values[label_col_name] = labelsdf[label_col_name]

    return values

def _find_closes_occuring_date(refdt, series_of_dt, where="before"):
    if where == "before":
        diff = refdt - (series_of_dt[series_of_dt < refdt])
    elif where == "after":
        diff = (series_of_dt[series_of_dt > refdt]) - refdt

    if diff.empty:
        # no occurences before of after

        return np.nan
    else:
        return min(diff).total_seconds()

def remove_outliers_from_obs(obsdf, outliersdf):
    # TODO this function can only be used with care!!!
    # because all timestamps will be removed that have an oulier in one specific obstype !!!!
    return obsdf.loc[~obsdf.index.isin(outliersdf.index)]


def conv_tz_multiidxdf(df, timezone):
    df.index = df.index.set_levels(df.index.levels[1].tz_convert(timezone), level=1)
    return df


def metadf_to_gdf(df, crs=4326):
    """
    Function to convert a dataframe with 'lat' en 'lon' columnst to a geopandas
    dataframe with a geometry column containing points.

    Special care for stations with missing coordinates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with a 'lat' en 'lon' column.
    crs : Integer, optional
        The epsg number of the coordinates. The default is 4326.

    Returns
    -------
    geodf : geopandas.GeaDataFrame
        The geodataframe equivalent of the df.

    """

    # only conver to points if coordinates are present
    coordsdf = df[(~df["lat"].isnull()) & (~df["lon"].isnull())]
    missing_coords_df = df[(df["lat"].isnull()) | (df["lon"].isnull())]

    geodf = gpd.GeoDataFrame(
        coordsdf, geometry=gpd.points_from_xy(coordsdf.lon, coordsdf.lat)
    )
    geodf = geodf.set_crs(epsg=crs)
    geodf = pd.concat([geodf, missing_coords_df])

    geodf = geodf.sort_index()
    return geodf


def multiindexdf_datetime_subsetting(df, starttime, endtime):
    "The multiindex equivalent of datetime_subsetting"
    dt_df = df.reset_index().set_index("datetime")
    subset_dt_df = datetime_subsetting(dt_df, starttime, endtime)

    # back to multiindex name-datetime
    subset_dt_df = subset_dt_df.reset_index()
    idx = pd.MultiIndex.from_frame(subset_dt_df[["name", "datetime"]])
    returndf = subset_dt_df.set_index(idx).drop(
        columns=["name", "datetime"], errors="ignore"
    )

    if returndf.empty:
        print(
            f"Warning: No observations left after subsetting datetime {starttime} -- {endtime} "
        )

    return returndf



# =============================================================================
# filters
# =============================================================================
def subset_stations(df, stationslist):
    df = df.loc[df.index.get_level_values(
                'name').isin(stationslist)]

    present_stations = df.index.get_level_values('name')
    not_present_stations = list(set(stationslist) - set(present_stations))
    if len(not_present_stations)!=0:
        print(f'WARNING: The stations: {not_present_stations} not found in the dataframe.')

    return df



def datetime_subsetting(df, starttime, endtime):
    """
    Wrapper function for subsetting a dataframe with a 'datetime' column or index with a start- and
    endtime.

    Parameters
    ----------
    df : pandas.DataFrame with datetimeindex
        The dataframe to apply the subsetting to.
    starttime : datetime.Datetime
        Starttime for the subsetting period (included).
    endtime : datetime.Datetime
        Endtime for the subsetting period (included).

    Returns
    -------
    pandas.DataFrame
        Subset of the df.

    """
    idx_names = list(df.index.names)
    df = df.reset_index()
    df = df.set_index('datetime')
    stand_format = "%Y-%m-%d %H:%M:%S"

    if isinstance(starttime, type(None)):
        startstring = None  # will select from the beginning of the df
    else:
        startstring = starttime.strftime(stand_format)
    if isinstance(endtime, type(None)):
        endstring = None
    else:
        endstring = endtime.strftime(stand_format)

    subset =  df[startstring:endstring]
    subset = subset.reset_index()
    subset = subset.set_index(idx_names)
    return subset


def conv_applied_qc_to_df(obstypes, ordered_checknames):
    if isinstance(obstypes, str):
        obstypes = [obstypes]
    if isinstance(ordered_checknames, str):
        ordered_checknames = [ordered_checknames]

    obslist = list(
        itertools.chain.from_iterable(
            itertools.repeat(item, len(ordered_checknames)) for item in obstypes
        )
    )

    checknamelist = list(
        itertools.chain.from_iterable(
            itertools.repeat(ordered_checknames, len(obstypes))
        )
    )

    df = pd.DataFrame({"obstype": obslist, "checkname": checknamelist})
    return df


# =============================================================================
# Records frequencies
# =============================================================================


def get_likely_frequency(
    timestamps, method="highest", simplify=True, max_simplify_error="2T"
):

    """
    Find the most likely observation frequency of a datetimeindex.

    Parameters
    ----------
    timestamps : pandas.Datetimeindex()
        Datetimeindex of the dataset.df.
    method : 'highest' or 'median', optional
        Select wich method to use. If 'highest', the highest apearing frequency is used.
        If 'median', the median of the apearing frequencies is used. The default is 'highest'.
    simplify : Boolean, optional
        If True, the likely frequency is converted to round hours, or round minutes.
        The "max_simplify_error' is used as a constrain. If the constrain is not met,
        the simplification is not performed.The default is True.
    max_simplify_error : datetimestring, optional
        The maximum deviation from the found frequency when simplifying. The default is '2T'.

    Returns
    -------
    assume_freq : datetime.timedelta
        The assumed (and simplified) frequency of the datetimeindex.

    """

    assert method in [
        "highest",
        "median",
    ], f"The method for frequency estimation ({method}) is not known. Use one of [highest, median]"


    try:
        pd.to_timedelta(max_simplify_error)
    except ValueError:
        sys.exit(
            f'{max_simplify_error} is not valid timeindication. Example: "5T" indicates 5 minutes.'
        )

    freqdist = abs(timestamps.to_series().diff().value_counts().index).sort_values(
        ascending=True
    )

    if method == "highest":
        assume_freq = freqdist[0]  # highest frequency

    elif method == "median":
        assume_freq = freqdist.median()

    if simplify:
        simplify_freq = None

        # try simplyfy to round hours
        trail_hour = assume_freq.ceil("H")
        lead_hour = assume_freq.floor("H")

        if (abs(lead_hour - assume_freq) <= abs(trail_hour - assume_freq)) & (
            lead_hour.total_seconds() != 0.0
        ):  # avoid assume freq of 0 seconds
            best_candidate = lead_hour
        else:
            best_candidate = trail_hour

        if abs(assume_freq - best_candidate) < pd.to_timedelta(max_simplify_error):
            simplify_freq = best_candidate

        # try simplyfy to round minutes
        if simplify_freq is None:
            trail_min = assume_freq.ceil("T")
            lead_min = assume_freq.floor("T")

            if (abs(lead_min - assume_freq) <= abs(trail_min - assume_freq)) & (
                lead_min.total_seconds() != 0.0
            ):  # avoid assume freq of 0 seconds
                best_candidate = lead_min
            else:
                best_candidate = trail_min

            if abs(assume_freq - best_candidate) < pd.to_timedelta(max_simplify_error):
                simplify_freq = best_candidate

        if simplify_freq is None:
            assume_freq = assume_freq
        else:
            assume_freq = simplify_freq


    if assume_freq == pd.to_timedelta(0):  # highly likely due to a duplicated record
        # select the second highest frequency
        assume_freq = abs(
            timestamps.to_series().diff().value_counts().index
        ).sort_values(ascending=True)[1]

    return assume_freq


def get_freqency_series(df, method="highest", simplify=True, max_simplify_error="2T"):

    """
    Find the most likely observation frequency for all stations individually
    based on the df. If an observation has less than two observations, assign
    the most commum frequency to it an raise a warning.

    Parameters
    ----------
    df : Metobs_toolkit.df
        Dataframe containing the observations.
    method : 'highest' or 'median', optional
        Select wich method to use. If 'highest', the highest apearing frequency is used.
        If 'median', the median of the apearing frequencies is used. The default is 'highest'.
    simplify : bool, optional
        If True, the likely frequency is converted to round hours, or round minutes.
        The "max_simplify_error' is used as a constrain. If the constrain is not met,
        the simplification is not performed.The default is True.
    max_simplify_error : Timedelta or str, optional
        The maximum deviation from the found frequency when simplifying. The default is '2T'.

    Returns
    -------
    freq_series : pandas.Series
        A pandas series with 'name' as index and likely frequencies as values.

    """

    problematic_stations = []
    freqs = {}
    for station in df.index.get_level_values(level="name").unique():
        subdf = df.xs(station, level="name")
        # remove rows with all obstype nans
        subdf = subdf.dropna(axis=0, how="all")

        # Check if all observations have at least two observations
        if subdf.shape[0] < 2:
            problematic_stations.append(station)
            print(
                f"WARNING! Stations {station} have to few observations to make a frequency estimate."
            )
            continue

        freqs[station] = get_likely_frequency(
            timestamps=subdf.index,
            method=method,
            simplify=simplify,
            max_simplify_error=max_simplify_error,
        )

    if len(problematic_stations) != 0:
        assign_med_freq = pd.to_timedelta(
            np.median([freq.total_seconds() for freq in freqs.values()]), unit="seconds"
        )

        print(
            f"Asigning the median of frequencies ({assign_med_freq}) to these stations {problematic_stations}."
        )
        for prob_station in problematic_stations:
            freqs[prob_station] = assign_med_freq


    return pd.Series(data=freqs)
