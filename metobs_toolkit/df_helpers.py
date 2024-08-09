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
import pytz
import logging

logger = logging.getLogger(__name__)

pd.options.mode.copy_on_write = True


def xs_save(df, key, level, drop_level=True):
    """Similar as pandas xs, but returns an empty df when key is not found."""
    try:
        return df.xs(key, level=level, drop_level=drop_level)
    except KeyError:
        # create empty df with same columns and index names
        columns = df.columns
        names = list(df.index.names)
        if drop_level:
            names.remove(level)

        levels = [[name] for name in names]
        codes = [[] for name in names]
        idx = pd.MultiIndex(
            levels=levels,
            codes=codes,
            names=names,
        )

    return pd.DataFrame(index=idx, columns=columns)


def concat_save(df_list, keep_all_nan_cols=False, **kwargs):
    """Concat dataframes row-wise without triggering the Futurwarning of concating empyt df's."""

    if all([isinstance(df, pd.DataFrame) for df in df_list]):
        if keep_all_nan_cols:
            return pd.concat(df_list, **kwargs)
        else:
            # This line will filter columns with all NAN values (so empty dfs + all NA entries are filtered out)
            return pd.concat([df.dropna(axis=1, how="all") for df in df_list], **kwargs)
    if all([isinstance(df, pd.Series) for df in df_list]):
        # This line will filter out empty series
        return pd.concat([ser for ser in df_list if not ser.empty], **kwargs)
    sys.exit("Cannot concat Dataframes and Series together")


def init_multiindex():
    """Construct a name-datetime pandas multiindex."""
    return pd.MultiIndex(
        levels=[["name"], ["datetime"]], codes=[[], []], names=["name", "datetime"]
    )


def init_multiindexdf():
    """Construct a name-datetime pandas multiindexdataframe."""
    return pd.DataFrame(index=init_multiindex())


def empty_outliers_df():
    """Construct an empty outliersdf."""
    my_index = pd.MultiIndex(
        levels=[["name"], ["obstype"], ["datetime"]],
        codes=[[], [], []],
        names=["name", "obstype", "datetime"],
    )

    outliersdf = pd.DataFrame(index=my_index, columns=["value", "label"])
    return outliersdf


def conv_tz_multiidxdf(df, timezone):
    """Convert datetime index to other timezone."""
    df.index = df.index.set_levels(df.index.levels[1].tz_convert(timezone), level=1)
    return df


def metadf_to_gdf(df, crs=4326):
    """Make geopandas dataframe.

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
    metadata_columns = geodf.columns
    geodf = concat_save([geodf, missing_coords_df])

    # Because empyt and Nan columns are skipped in the concat save, add them
    # again if needed
    for col in metadata_columns:
        if col not in geodf:
            geodf[col] = np.nan

    geodf = geodf.sort_index()
    return geodf


def multiindexdf_datetime_subsetting(df, starttime, endtime):
    """Multiindex equivalent of datetime_subsetting."""
    dt_df = df.reset_index().set_index("datetime")
    subset_dt_df = datetime_subsetting(dt_df, starttime, endtime)

    # back to multiindex name-datetime
    subset_dt_df = subset_dt_df.reset_index()
    idx = pd.MultiIndex.from_frame(subset_dt_df[["name", "datetime"]])
    returndf = subset_dt_df.set_index(idx).drop(
        columns=["name", "datetime"], errors="ignore"
    )

    if returndf.empty:
        logger.warning(
            f"No observations left after subsetting datetime {starttime} -- {endtime} "
        )

    return returndf


# =============================================================================
# filters
# =============================================================================
def subset_stations(df, stationslist):
    """Subset stations by name from a dataframe."""
    df = df.loc[df.index.get_level_values("name").isin(stationslist)]

    present_stations = df.index.get_level_values("name")
    not_present_stations = list(set(stationslist) - set(present_stations))
    if len(not_present_stations) != 0:
        logger.warning(
            f"The stations: {not_present_stations} not found in the dataframe."
        )

    return df


def datetime_subsetting(df, starttime, endtime):
    """Subset dataaframe by timeperiod.

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
    df = df.set_index("datetime")

    if isinstance(starttime, type(None)):
        starttime = df.index.min()  # will select from the beginning of the df
    else:
        if starttime.tzinfo is None:
            # set timezone when unaware
            starttime = starttime.replace(tzinfo=df.index.tzinfo)

    if isinstance(endtime, type(None)):
        endtime = df.index.max()
    else:
        if endtime.tzinfo is None:
            # set timezone when unaware
            endtime = endtime.replace(tzinfo=df.index.tzinfo)

    subset = df[(df.index >= starttime) & (df.index <= endtime)]
    subset = subset.reset_index()
    subset = subset.set_index(idx_names)
    subset = subset.sort_index()
    return subset


# =============================================================================
# Records frequencies
# =============================================================================


def _simplify_time(time, max_simplyfi_error, zero_protection=False):
    """Simplifies a time (or timedelta) to a rounded value, within a tolerance.

    time can be a timestamp of a timedelta.

    zero_protections make shure that the returned simplified time is not zero.
    This is in practice used when applied to frequencies.
    """
    # 1 try simplyfy to round hours
    hourly_candidate = time.round("h")
    if zero_protection:
        if hourly_candidate == pd.Timedelta(0):
            hourly_candidate = time.ceil(
                "h"
            )  # try the ceil as candidate (is never zero)

    if abs(time - hourly_candidate) < pd.to_timedelta(max_simplyfi_error):
        return hourly_candidate

    # 2  half_hourly_candidate
    halfhourly_candidate = time.round("30min")
    if zero_protection:
        if halfhourly_candidate == pd.Timedelta(0):
            halfhourly_candidate = time.ceil(
                "30min"
            )  # try the ceil as candidate (is never zero)
    if abs(time - halfhourly_candidate) < pd.to_timedelta(max_simplyfi_error):
        return halfhourly_candidate

    # 3. 10min fold
    tenmin_candidate = time.round("10min")
    if zero_protection:
        if tenmin_candidate == pd.Timedelta(0):
            tenmin_candidate = time.ceil(
                "10min"
            )  # try the ceil as candidate (is never zero)
    if abs(time - tenmin_candidate) < pd.to_timedelta(max_simplyfi_error):
        return tenmin_candidate

    # 4. 5min fold
    fivemin_candidate = time.round("5min")
    if zero_protection:
        if fivemin_candidate == pd.Timedelta(0):
            fivemin_candidate = time.ceil(
                "5min"
            )  # try the ceil as candidate (is never zero)
    if abs(time - fivemin_candidate) < pd.to_timedelta(max_simplyfi_error):
        return fivemin_candidate

    # 5. min fold
    min_candidate = time.round("min")
    if zero_protection:
        if min_candidate == pd.Timedelta(0):
            min_candidate = time.ceil(
                "1min"
            )  # try the ceil as candidate (is never zero)
    if abs(time - min_candidate) < pd.to_timedelta(max_simplyfi_error):
        return min_candidate

    # 6. No simplyfication posible
    if zero_protection:
        if time == pd.Timedelta(0):
            raise MetobsDfHelpersError()
    return time


def get_likely_frequency(
    timestamps, method="highest", simplify=True, max_simplify_error="2min"
):
    """Find the most likely observation frequency of a datetimeindex.

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
        The maximum deviation from the found frequency when simplifying. The default is '2min'.

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
            f'{max_simplify_error} is not valid timeindication. Example: "5min" indicates 5 minutes.'
        )

    freqs_blacklist = [pd.Timedelta(0), np.nan]  # avoid a zero frequency

    freqs = timestamps.to_series().diff()
    freqs = freqs[~freqs.isin(freqs_blacklist)]

    if method == "highest":
        assume_freq = freqs.min()  # highest frequency

    elif method == "median":
        assume_freq = freqs.median()

    if simplify:
        simplify_freq = None

        # try simplyfy to round hours
        trail_hour = assume_freq.ceil("h")
        lead_hour = assume_freq.floor("h")

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
            trail_min = assume_freq.ceil("min")
            lead_min = assume_freq.floor("min")

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


def get_freqency_series(df, method="highest", simplify=True, max_simplify_error="2min"):
    """Get the most likely frequencies of all stations.

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
        The maximum deviation from the found frequency when simplifying. The default is '2min'.

    Returns
    -------
    freq_series : pandas.Series
        A pandas series with 'name' as index and likely frequencies as values.

    """
    problematic_stations = []
    freqs = {}
    for station in df.index.get_level_values(level="name").unique():
        subdf = xs_save(df, station, level="name")
        # remove rows with all obstype nans
        subdf = subdf.dropna(axis=0, how="all")

        # Check if all observations have at least two observations
        if subdf.shape[0] < 2:
            problematic_stations.append(station)
            logger.warning(
                f"Stations {station} have to few observations to make a frequency estimate."
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

        logger.warning(
            f"Asigning the median of frequencies ({assign_med_freq}) to these stations {problematic_stations}."
        )
        for prob_station in problematic_stations:
            freqs[prob_station] = assign_med_freq

    return pd.Series(data=freqs)


class MetobsDfHelpersError(Exception):
    """Exception raised for errors in the datasetgaphandling."""

    def __init__(self):
        message = f"This error is likely due to a Bug. Please report this as an issue on the GitHub repository."
        super().__init__(message)
