import os
import logging
import concurrent.futures
from math import radians, cos, sin, asin, sqrt
from typing import Union, List, Dict, Tuple

import numpy as np
import pandas as pd

from metobs_toolkit.backend_collection.df_helpers import to_timedelta

logger = logging.getLogger("<metobs_toolkit>")


def _calculate_distance_matrix_with_haverine(
        metadf: pd.DataFrame) -> pd.DataFrame:
    """
    Calculate a distance matrix between points using the Haversine formula.

    This function computes the great circle distance between all pairs of
    geographical points in the given DataFrame using the Haversine formula.
    The distances are returned in a pandas DataFrame.

    Parameters
    ----------
    metadf : pandas.DataFrame
        DataFrame containing metadata for geographical points. Each row
        represents a point, and the `geometry` column must contain shapely
        Point objects with `x` (longitude) and `y` (latitude) attributes.

    Returns
    -------
    pandas.DataFrame
        DataFrame representing the distance matrix. The rows and columns
        are indexed by the same station identifiers as in `metadf`, and the
        values represent the distances (in meters) between the corresponding
        points.

    Notes
    -----
    The radius of the Earth is assumed to be 6,367,000 meters.
    The Haversine formula calculates the great circle distance, which is
    the shortest distance over the Earth's surface.
    """
    logger.debug("Entering _calculate_distance_matrix_with_haverine")
    if not isinstance(metadf, pd.DataFrame):
        raise TypeError("metadf must be a pandas.DataFrame")

    def haversine(lon1: float, lat1: float, lon2: float, lat2: float) -> float:
        """Calculate the great circle distance between two points."""
        # convert decimal degrees to radians
        lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2])

        # haversine formula
        dlon = lon2 - lon1
        dlat = lat2 - lat1
        a = sin(dlat / 2) ** 2 + cos(lat1) * cos(lat2) * sin(dlon / 2) ** 2
        c = 2 * asin(sqrt(a))
        r = 6367000  # Radius of earth in meter.
        return c * r

    distance_matrix = {}
    for sta1, row1 in metadf.iterrows():
        distance_matrix[sta1] = {}
        for sta2, row2 in metadf.iterrows():
            distance_matrix[sta1][sta2] = haversine(
                row1.geometry.x, row1.geometry.y,
                row2.geometry.x, row2.geometry.y
            )
    return pd.DataFrame(distance_matrix)


def synchronize_series(
    series_list: List[pd.Series], max_shift: pd.Timedelta
) -> Tuple[pd.DataFrame, Dict]:
    """
    Synchronize a list of pandas Series with datetime indexes.

    The target timestamps are defined by:
     * freq: the highest frequency present in the input series
     * origin: the earliest timestamp found, rounded down by the freq
     * closing: the latest timestamp found, rounded up by the freq.

    Parameters
    ----------
    series_list : list of pandas.Series
        List of pandas Series with datetime indexes.
    max_shift : pandas.Timedelta
        Maximum shift in time that can be applied to each timestamp
        in synchronization.

    Returns
    -------
    pandas.DataFrame
        DataFrame with synchronized Series.
    dict
        Dictionary mapping each synchronized timestamp to its
        original timestamp.
    """
    logger.debug("Entering synchronize_series")
    if not isinstance(series_list, list):
        raise TypeError("series_list must be a list of pandas.Series")
    if not all(isinstance(s, pd.Series) for s in series_list):
        raise TypeError("All elements in series_list must be pandas.Series")
    if not isinstance(max_shift, pd.Timedelta):
        raise TypeError("max_shift must be a pandas.Timedelta")

    # find highest frequency
    frequencies = [to_timedelta(s.index.inferred_freq) for s in series_list]
    trg_freq = min(frequencies)

    # find origin and closing timestamp (earliest/latest)
    origin = min([s.index.min() for s in series_list]).floor(trg_freq)
    closing = max([s.index.max() for s in series_list]).ceil(trg_freq)

    # Create target datetime axes
    target_dt = pd.date_range(start=origin, end=closing, freq=trg_freq)

    # Synchronize (merge with tolerance) series to the common index
    synchronized_series = []
    timestamp_mapping = {}
    for s in series_list:
        targetdf = (
            s.to_frame()
            .assign(orig_datetime=s.index)
            .reindex(
                index=pd.DatetimeIndex(target_dt),
                method="nearest",
                tolerance=max_shift,
                limit=1,
            )
        )

        # extract the mapping (new -> original)
        orig_timestampseries = targetdf["orig_datetime"]
        orig_timestampseries.name = "original_timestamp"
        timestamp_mapping[s.name] = orig_timestampseries

        synchronized_series.append(s)

    return pd.concat(synchronized_series, axis=1), timestamp_mapping


def _find_spatial_buddies(
    distance_df: pd.DataFrame, buddy_radius: Union[int, float]
) -> Dict:
    """
    Get neighbouring stations using buddy radius.

    Parameters
    ----------
    distance_df : pandas.DataFrame
        DataFrame containing distances between stations.
    buddy_radius : int or float
        Maximum distance (in meters) to consider as a buddy.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of its buddies.
    """
    logger.debug("Entering _find_spatial_buddies")
    if not isinstance(distance_df, pd.DataFrame):
        raise TypeError("distance_df must be a pandas.DataFrame")
    if not isinstance(buddy_radius, (int, float)):
        raise TypeError("buddy_radius must be an int or float")

    buddies = {}
    for refstation, distances in distance_df.iterrows():
        bud_stations = distances[distances <= buddy_radius].index.to_list()
        bud_stations.remove(refstation)
        buddies[refstation] = bud_stations

    return buddies


def _filter_to_altitude_buddies(
        spatial_buddies: Dict,
        altitudes: pd.Series,
        max_altitude_diff: Union[int, float]) -> Dict:
    """
    Filter neighbours by maximum altitude difference.

    Parameters
    ----------
    spatial_buddies : dict
        Dictionary mapping each station to a list of its spatial buddies.
    altitudes : pandas.Series
        Series containing altitudes for each station.
    max_altitude_diff : int or float
        Maximum allowed altitude difference.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of altitude-filtered buddies.
    """
    logger.debug("Entering _filter_to_altitude_buddies")
    if not isinstance(spatial_buddies, dict):
        raise TypeError("spatial_buddies must be a dict")
    if not isinstance(altitudes, pd.Series):
        raise TypeError("altitudes must be a pandas.Series")
    if not isinstance(max_altitude_diff, (int, float)):
        raise TypeError("max_altitude_diff must be an int or float")

    alt_buddies_dict = {}
    for refstation, buddylist in spatial_buddies.items():
        alt_diff = abs((altitudes.loc[buddylist]) - altitudes.loc[refstation])

        alt_buddies = alt_diff[alt_diff <= max_altitude_diff].index.to_list()
        alt_buddies_dict[refstation] = alt_buddies
    return alt_buddies_dict


def _filter_to_minimum_samplesize(
        buddydict: Dict,
        min_sample_size: int) -> Dict:
    """
    Filter stations that are too isolated using minimum sample size.

    Parameters
    ----------
    buddydict : dict
        Dictionary mapping each station to a list of its buddies.
    min_sample_size : int
        Minimum number of buddies required.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of buddies meeting the
        minimum sample size.
    """
    logger.debug("Entering _filter_to_minimum_samplesize")
    if not isinstance(buddydict, dict):
        raise TypeError("buddydict must be a dict")
    if not isinstance(min_sample_size, int):
        raise TypeError("min_sample_size must be an int")

    to_check_stations = {}
    for refstation, buddies in buddydict.items():
        if len(buddies) < min_sample_size:
            # not enough buddies
            to_check_stations[refstation] = []  # remove buddies
        else:
            to_check_stations[refstation] = buddies
    return to_check_stations


def create_groups_of_buddies(buddydict: Dict) -> List[Tuple]:
    """
    Create unique groups of buddies from a buddy dictionary.

    Parameters
    ----------
    buddydict : dict
        Dictionary mapping each station to a list of its buddies.

    Returns
    -------
    list of tuple
        List of tuples, each containing a group of station names.
    """
    logger.debug("Entering create_groups_of_buddies")
    if not isinstance(buddydict, dict):
        raise TypeError("buddydict must be a dict")

    grouped_stations = []
    groups = []
    for refstation, buddies in buddydict.items():
        if not bool(buddies):
            continue
        if refstation in grouped_stations:
            continue
        group = tuple([refstation, *buddies])

        grouped_stations.extend([*group])
        groups.append(group)

    return groups


def toolkit_buddy_check(
    dataset: "Dataset",  # noqa: F821
    obstype: str,
    buddy_radius: Union[int, float],
    min_sample_size: int,
    max_alt_diff: Union[int, float, None],
    min_std: Union[int, float],
    std_threshold: Union[int, float],
    N_iter: int,
    instantanious_tolerance: pd.Timedelta,  # TYPO
    lapserate: Union[float, None] = None,  # -0.0065 for temperature #TYPO
    use_mp: bool = True,
) -> Tuple[list, dict]:
    """
        Spatial buddy check.

        The buddy check compares an observation against its neighbors
        (i.e. buddies). The check loops over all the groups, which are stations
        within a radius of each other. For each group, the absolute value of
        the difference with the group mean, normalized by the standard
        deviation (with a defined minimum), is computed. If one (or more)
        exceeds the std_threshold, the most extreme (=baddest) observation of
        that group is labeled as an outlier.

        Multiple iterations of this checks can be done using the N_iter.

        A schematic step-by-step description of the buddy check:

          1. A distance matrix is constructed for all interdistances between
            the stations. This is done using the haversine approximation.
          2. Groups of buddies (neighbours) are created by using the
            buddy_radius. These groups are further filtered by:
            * removing stations from the groups that differ to much in altitude
              (based on the max_alt_diff)
            * removing groups of buddies that are too small (based on the
              min_sample_size)

          3. Observations per group are synchronized in time (using the
            max_shift as tolerance for allignment).
          4. If a lapsrate is specified, the observations are corrected for
            altitude differences.
          5. For each buddy group:
            * The mean, standard deviation (std), and sample size are computed.
            * If the std is lower than the minimum std, it is replaced by the
              minimum std.
            * Chi values are calculated for all records.
            * For each timestamp the record with the highest Chi is tested if
              it is larger then std_threshold.
              If so, that record is flagged as an outlier. It will be ignored
                in the next iteration.
            * This is repeated N_iter times.


        Parameters
        ----------
        dataset : Dataset
            The dataset to apply the buddy check on.
        obstype : str
            The observation type that has to be checked.
        buddy_radius : int or float
            The radius to define spatial neighbors in meters.
        min_sample_size : int
            The minimum sample size to calculate statistics on.
        max_alt_diff : int, float, or None
            The maximum altitude difference allowed for buddies. If None,
            no altitude filter is applied.
        min_std : int or float
            The minimum standard deviation for sample statistics. This should
            represent the accuracy of the observations.
        std_threshold : int or float
            The threshold (std units) for flagging observations as outliers.
        N_iter : int
            The number of iterations to perform the buddy check.
        instantanious_tolerance : pandas.Timedelta
            The maximum time difference allowed for synchronizing observations.
        lapserate : float or None, optional
            Describes how the obstype changes with altitude (in meters). If
            None, no altitude correction is applied. For temperature, a
            common value is -0.0065.
        use_mp : bool, optional
            Use multiprocessing to speed up the buddy check.

        Returns
        -------
        list
            A list of tuples containing the outlier station, timestamp,
            and detail message. Each tuple is in the form (station_name,
            timestamp, message).
        dict
            A dictionary mapping each synchronized timestamp to its original
              timestamp.

        Notes
        -----
        The altitude of the stations can be extracted from GEE by using the
        `Dataset.get_altitude()` method.

    """

    # -----  Part 1: construct buddy groups ------
    # compute distance metric
    metadf = dataset.metadf
    dist_matrix = _calculate_distance_matrix_with_haverine(metadf)

    # find potential buddies by distance
    buddies = _find_spatial_buddies(distance_df=dist_matrix,
                                    buddy_radius=buddy_radius)

    # filter buddies by altitude difference
    if max_alt_diff is not None:
        if metadf["altitude"].isna().any():
            raise ValueError("At least one station has a NaN \
value for 'altitude'")
        # Filter by altitude difference
        buddies = _filter_to_altitude_buddies(
            spatial_buddies=buddies,
            altitudes=metadf["altitude"],
            max_altitude_diff=max_alt_diff,
        )

    # Filter by sample size (based on the number of buddy stations)
    buddies = _filter_to_minimum_samplesize(
        buddydict=buddies, min_sample_size=min_sample_size
    )

    # create unique groups of buddies (list of tuples)
    buddygroups = create_groups_of_buddies(buddies)

    # ---- Part 2: Preparing the records  -----

    # construct a wide observation dataframe
    concatlist = []
    for sta in dataset.stations:
        if obstype in sta.sensordata.keys():
            records = sta.get_sensor(obstype).series
            records.name = sta.name
            concatlist.append(records)

    # synchronize the timestamps
    combdf, timestamp_map = synchronize_series(
        series_list=concatlist, max_shift=instantanious_tolerance
    )

    # lapse rate correction
    if lapserate is not None:
        # get altitude dataframe
        altdict = {sta.name: sta.site.altitude for sta in dataset.stations}
        altseries = pd.Series(altdict)
        altcorrectionseries = (altseries - altseries.min()) * lapserate
        combdf = combdf - altcorrectionseries  # Correct for altitude

    # ---- Part 3 : Apply buddy check on each group,
    #  rejecting the most extreme outlier

    outliersbin = []
    for i in range(N_iter):
        # convert values to NaN, if they are labeled as outlier in
        #  previous iteration
        if bool(outliersbin):
            for outlier_station, outlier_time, _msg in outliersbin:
                if outlier_station in combdf.columns:
                    combdf.loc[outlier_time, outlier_station] = np.nan

        if use_mp:
            # Use multiprocessing generator (parallelization)
            num_cpus = os.cpu_count()
            # since this check is an instantaneous check --> 
            # perfect for splitting the dataset in chunks in time
            chunks = np.array_split(combdf, num_cpus)

            # create inputargs for each buddygroup, and for each chunk in time
            inputargs = [
                (buddygroup, chunk, min_sample_size, min_std, std_threshold)
                for buddygroup in buddygroups
                for chunk in chunks
            ]

            with concurrent.futures.ProcessPoolExecutor() as executor:
                outliers = executor.map(find_buddy_group_outlier, inputargs)
            outliers = list(outliers)

        else:
            # create inputargs for each buddygroup, and for each chunk in time
            inputargs = [
                (buddygroup, combdf, min_sample_size, min_std, std_threshold)
                for buddygroup in buddygroups
            ]

            outliers = list(map(find_buddy_group_outlier, inputargs))

        # unpack double nested list
        outliers = [item for sublist in outliers for item in sublist]
        outliersbin.extend(outliers)
        i += 1

    return outliersbin, timestamp_map


def find_buddy_group_outlier(inputarg: Tuple) -> List[Tuple]:
    """
    Apply a buddy check on a group to identify outliers.

    Parameters
    ----------
    inputarg : tuple
        A tuple containing:
        - buddygroup : list
            List of station names that form the buddy group.
        - combdf : pandas.DataFrame
            DataFrame containing the combined data for all stations.
        - min_sample_size : int
            Minimum number of non-NaN values required in the buddy group for a
            valid comparison.
        - min_std : float
            Minimum standard deviation to use when calculating z-scores.
        - outlier_threshold : float
            Threshold for identifying outliers in terms of z-scores.

    Returns
    -------
    list of tuple
        Each tuple contains:
        - str : The station name of the most extreme outlier.
        - pandas.Timestamp : The timestamp of the outlier.
        - str : A detailed message describing the outlier.

    Notes
    -----
    This function performs the following steps:
    1. Subsets the data to the buddy group.
    2. Calculates the mean, standard deviation, and count of non-NaN values
      for each timestamp.
    3. Filters out timestamps with insufficient data.
    4. Replaces standard deviations below the minimum threshold with the
      minimum value.
    5. Converts station values to z-scores.
    6. Identifies timestamps with at least one outlier.
    7. Locates the most extreme outlier for each timestamp.
    8. Generates a detailed message for each outlier.
    """
    logger.debug("Entering find_buddy_group_outlier")
    if not isinstance(inputarg, tuple):
        raise TypeError("inputarg must be a tuple")
    buddygroup, combdf = inputarg[0], inputarg[1]
    min_sample_size, min_std, outlier_threshold = inputarg[2:]
    if not isinstance(buddygroup, (list, tuple)):
        raise TypeError("buddygroup must be a list or tuple")
    if not isinstance(combdf, pd.DataFrame):
        raise TypeError("combdf must be a pandas.DataFrame")
    if not isinstance(min_sample_size, int):
        raise TypeError("min_sample_size must be an int")
    if not isinstance(min_std, (int, float)):
        raise TypeError("min_std must be an int or float")
    if not isinstance(outlier_threshold, (int, float)):
        raise TypeError("outlier_threshold must be an int or float")

    # subset to the buddies
    buddydf = combdf[[*buddygroup]]

    # calculate std and mean row wise
    buddydf["mean"] = buddydf[[*buddygroup]].mean(axis=1)
    buddydf["std"] = buddydf[[*buddygroup]].std(axis=1)
    buddydf["non_nan_count"] = buddydf[[*buddygroup]].notna().sum(axis=1)

    # subset to samples with enough members (check for each timestamp 
    # specifically)
    buddydf = buddydf.loc[buddydf["non_nan_count"] >= min_sample_size]

    # replace std by minimum, if needed
    buddydf.loc[buddydf["std"] < min_std, "std"] = np.float32(min_std)

    # Convert values to sigmas
    for station in buddygroup:
        buddydf[station] = ((buddydf[station] - buddydf["mean"]).abs()
                            /
                            buddydf["std"])

    # Drop rows for which all values are smaller than the threshold 
    # (speed up the last step)
    buddydf["timestamp_with_outlier"] = buddydf[[*buddygroup]].apply(
        lambda row: any(row > outlier_threshold), axis=1
    )
    buddydf = buddydf.loc[buddydf["timestamp_with_outlier"]]

    # locate the most extreme outlier per timestamp
    buddydf["is_the_most_extreme_outlier"] = (
        buddydf[[*buddygroup]].idxmax(axis=1)
    )

    def msgcreator(row):
        retstr = f"Outlier at {row['is_the_most_extreme_outlier']}"
        retstr += f" with chi value \
{row[row['is_the_most_extreme_outlier']]:.2f},"
        retstr += f" is part of {buddygroup}, with mean: {row['mean']:.2f}, \
std: {row['std']:.2f}"
        return retstr

    # detail info string
    buddydf["detail_msg"] = buddydf.apply(
        lambda row: msgcreator(row), axis=1, result_type="reduce"
    )

    return list(
        zip(
            buddydf["is_the_most_extreme_outlier"],
            buddydf.index,
            buddydf["detail_msg"]
        )
    )
