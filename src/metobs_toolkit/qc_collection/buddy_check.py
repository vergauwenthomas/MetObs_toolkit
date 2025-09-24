import os
import logging
import concurrent.futures
from math import radians, cos, sin, asin, sqrt
from typing import Union, List, Dict, Tuple

import numpy as np
import pandas as pd

from metobs_toolkit.backend_collection.df_helpers import to_timedelta

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


def _calculate_distance_matrix_with_haverine(metadf: pd.DataFrame) -> pd.DataFrame:
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

    @log_entry
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
                row1.geometry.x, row1.geometry.y, row2.geometry.x, row2.geometry.y
            )
    return pd.DataFrame(distance_matrix)


@log_entry
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


def _find_LCZ_buddies(
    metadf: pd.DataFrame, max_dist: Union[int, float], distance_df: pd.DataFrame
) -> Dict:
    """
    Get neighbouring stations using both LCZ class and spatial distance.

    Parameters
    ----------
    metadf : pandas.DataFrame
        DataFrame containing metadata for stations. Must include an 'LCZ' column.
    max_dist : int or float
        Maximum distance (in meters) to consider as a LCZ buddy.
    distance_df : pandas.DataFrame
        DataFrame containing distances between stations.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of its LCZ buddies that are also within the specified distance.

    Notes
    -----
    - LCZ buddies are stations with the same LCZ class as the reference station.
    - The final buddies are the intersection of LCZ buddies and spatial buddies within `max_dist`.
    """

    LCZ_buddies = {}
    # Find buddies by LCZ
    for refstation in metadf.index:
        ref_LCZ = metadf.loc[refstation, "LCZ"]
        ref_buddies = metadf.loc[metadf["LCZ"] == ref_LCZ].index.to_list()
        LCZ_buddies[refstation] = ref_buddies

    # Find buddies by distance
    spatial_buddies = _find_spatial_buddies(distance_df, max_dist)

    # Crossection of both buddy defenitions
    final_buddies = {}
    for refstation in LCZ_buddies.keys():
        final_buddies[refstation] = list(
            set(LCZ_buddies[refstation]).intersection(set(spatial_buddies[refstation]))
        )

    return final_buddies


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
    buddies: Dict, altitudes: pd.Series, max_altitude_diff: Union[int, float]
) -> Dict:
    """
    Filter neighbours by maximum altitude difference.

    Parameters
    ----------
    buddies : dict
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

    alt_buddies_dict = {}
    for refstation, buddylist in buddies.items():
        alt_diff = abs((altitudes.loc[buddylist]) - altitudes.loc[refstation])

        alt_buddies = alt_diff[alt_diff <= max_altitude_diff].index.to_list()
        alt_buddies_dict[refstation] = alt_buddies
    return alt_buddies_dict


def _filter_to_minimum_samplesize(buddydict: Dict, min_sample_size: int) -> Dict:
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

    to_check_stations = {}
    for refstation, buddies in buddydict.items():
        if len(buddies) < min_sample_size:
            # not enough buddies
            to_check_stations[refstation] = []  # remove buddies
        else:
            to_check_stations[refstation] = buddies
    return to_check_stations


@log_entry
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


@log_entry
def toolkit_buddy_check(
    target_stations: list["Station"],  # noqa: F821
    metadf: pd.DataFrame,
    obstype: str,
    spatial_buddy_radius: Union[int, float],
    spatial_min_sample_size: int,
    max_alt_diff: Union[int, float, None],
    min_std: Union[int, float],
    spatial_z_threshold: Union[int, float],
    N_iter: int,
    instantaneous_tolerance: pd.Timedelta,
    # LCZ safety net
    max_LCZ_buddy_dist: Union[int, float, None],
    min_LCZ_safetynet_sample_size: Union[int, None],
    safetynet_z_threshold: Union[int, float, None],
    use_LCZ_safetynet: bool = False,
    # Technical
    lapserate: Union[float, None] = None,  # -0.0065 for temperature
    use_mp: bool = True,
) -> Tuple[list, dict]:
    """
    Spatial buddy check.

    The buddy check compares an observation against its neighbors
    (i.e. spatial buddies). The check loops over all the groups, which are stations
    within a radius of each other. For each group, the z-value of the reference
    observation is computed given the sample of spatial buddies. If one (or more)
    exceeds the `spatial_z_threshold`, the most extreme (=baddest) observation of
    that group is labeled as an outlier.

    Multiple iterations of this checks can be done using the N_iter.

    Optional a LCZ-safetynet can be applied.
    If so, the (potential) outliers, per iteration, are tested with another
    sample. This sample contains the LCZ-buddies, that are stations with the
    same LCZ as the reference station, and with a maximum distance of
    `max_LCZ_buddy_dist`. If a `max_alt_diff` is specified, a altitude-difference
    filtering is applied on these buddies aswell.  If a test is sucsesfull, that
    is if the z-value is smaller than the `safetynet_z_threshold`, then the
    outlier is saved. It will be removed from the outliers, and will pass to the
    next iteration or the end of this function.

    A schematic step-by-step description of the buddy check:

    #. A distance matrix is constructed for all interdistances between
       the stations. This is done using the haversine approximation.
    #. Groups of spatial buddies (neighbours) are created by using the
       `spatial_buddy_radius.` These groups are further filtered by:

       * removing stations from the groups that differ to much in altitude
         (based on the `max_alt_diff`)
       * removing groups of buddies that are too small (based on the
         `min_sample_size`)

    #. Observations per group are synchronized in time (using the
       `instantaneous_tolerance` for allignment).
    #. If a `lapsrate` is specified, the observations are corrected for
       altitude differences.
    #. The following steps are repeated for `N-iter` iterations:

       #. The values of outliers flaged by a previous iteration are converted to
          NaN's. Therefore they are not used in any following score or sample.
       #. For each buddy group:

          * The mean, standard deviation (std), and sample size are computed.
          * If the std is lower than the `minimum_std`, it is replaced by the
            minimum std.
          * Chi values are calculated for all records.
          * For each timestamp the record with the highest Chi is tested if
            it is larger then spatial_z_threshold.
            If so, that record is flagged as an outlier. It will be ignored
            in the next iteration.

       #. If `use_LCZ_safetynet` is True, the following steps are applied on
          the outliers flagged by the current iteration.

          * The LCZ-buddy sample is tested in size (samplesize must be bigger
            then `min_LCZ_safetynet_sample_size`). If the condition is not met,
            the safetynet test is not applied.
          * The safetynet test is applied:

            * The mean and std are computed of the LCZ-buddy sample. If the
              std is smaller then `min_std`, then the latter is used.
            * The z-value is computed for the target record (= flagged outlier).
            * If the z-value is smaller than `safetynet_z_threshold`, the
              tested outlier is "saved", and is removed from the set of outliers
              for the current iteration.

    Parameters
    ----------
    target_stations : list[Station]
        A list of Station objects to apply the buddy check on. These should be
        stations that contain the target observation type.
    metadf : pandas.DataFrame
        DataFrame containing station metadata including coordinates (geometry)
        and altitude information for all stations.
    obstype : str
        The observation type that has to be checked.
    spatial_buddy_radius : int or float
        The radius to define spatial neighbors in meters.
    spatial_min_sample_size : int
        The minimum sample size to calculate statistics on used by
        spatial-buddy samples.
    max_alt_diff : int, float, or None
        The maximum altitude difference allowed for buddies. If None,
        no altitude filter is applied.
    min_std : int or float
        The minimum standard deviation for sample statistics. This should
        represent the accuracy of the observations.
    spatial_z_threshold : int or float
        The threshold, tested with z-scores, for flagging observations as outliers.
    N_iter : int
        The number of iterations to perform the buddy check.
    instantaneous_tolerance : pandas.Timedelta
        The maximum time difference allowed for synchronizing observations.
    max_LCZ_buddy_dist:  int or float or None
        The radius to look for LCZ buddies (= same LCZ as a reference station).
        Typically this is set larger than the buddy_radius (= radius used for
        spatial buddies). This parameter is only used when use_LCZ_safetynet is
        True. The default is None.
    min_LCZ_safetynet_sample_size: int or None
        The minimum samplesize for the LCZ safetynet test. If a sample is to
        small, the safetynet test is not applied. The default is None.
    safetynet_z_threshold: int or float or None
        The threshold for a succesfull safety net test. If the z-value is
        less than `safetynet_z_threshold`, the test is succesfull and the
        outlier is "saved". It can proceed as a regular observation in the next
        iteration.
    use_LCZ_safetynet: bool
        If True, the LCZ safetynet test is applied on the stations flagged as
        outlier by the spatial buddies in each iteration. The default is False.
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

    * The altitude of the stations can be extracted from GEE by using the
      `Dataset.get_altitude()` method.
    * The LCZ of the stations can be extracted from GEE by using the
      `Dataset.get_LCZ()` method.


    """

    # -----  Part 1: construct buddy groups ------
    # compute distance metric
    dist_matrix = _calculate_distance_matrix_with_haverine(metadf)

    # find potential buddies by distance
    spatial_buddies = _find_spatial_buddies(
        distance_df=dist_matrix, buddy_radius=spatial_buddy_radius
    )

    if use_LCZ_safetynet:
        LCZ_buddies = _find_LCZ_buddies(
            metadf=metadf, max_dist=max_LCZ_buddy_dist, distance_df=dist_matrix
        )

    # filter buddies by altitude difference
    if max_alt_diff is not None:
        if metadf["altitude"].isna().any():
            raise ValueError(
                "At least one station has a NaN \
value for 'altitude'"
            )
        # Filter by altitude difference
        spatial_buddies = _filter_to_altitude_buddies(
            buddies=spatial_buddies,
            altitudes=metadf["altitude"],
            max_altitude_diff=max_alt_diff,
        )
        if use_LCZ_safetynet:
            LCZ_buddies = _filter_to_altitude_buddies(
                buddies=LCZ_buddies,
                altitudes=metadf["altitude"],
                max_altitude_diff=max_alt_diff,
            )

    # Filter by sample size (based on the number of buddy stations)
    spatial_buddies = _filter_to_minimum_samplesize(
        buddydict=spatial_buddies, min_sample_size=spatial_min_sample_size
    )

    # create unique groups of buddies (list of tuples)
    buddygroups = create_groups_of_buddies(spatial_buddies)

    # ---- Part 2: Preparing the records  -----

    # construct a wide observation dataframe
    concatlist = []
    for sta in target_stations:
        if obstype in sta.sensordata.keys():
            records = sta.get_sensor(obstype).series
            records.name = sta.name
            concatlist.append(records)

    # synchronize the timestamps
    combdf, timestamp_map = synchronize_series(
        series_list=concatlist, max_shift=instantaneous_tolerance
    )

    # lapse rate correction
    if lapserate is not None:
        # get altitude dataframe
        altdict = {sta.name: sta.site.altitude for sta in target_stations}
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
                (
                    buddygroup,
                    chunk,
                    spatial_min_sample_size,
                    min_std,
                    spatial_z_threshold,
                )
                for buddygroup in buddygroups
                for chunk in chunks
            ]

            with concurrent.futures.ProcessPoolExecutor() as executor:
                outliers = executor.map(find_buddy_group_outlier, inputargs)
            outliers = list(outliers)

        else:
            # create inputargs for each buddygroup, and for each chunk in time
            inputargs = [
                (
                    buddygroup,
                    combdf,
                    spatial_min_sample_size,
                    min_std,
                    spatial_z_threshold,
                )
                for buddygroup in buddygroups
            ]

            outliers = list(map(find_buddy_group_outlier, inputargs))

        # unpack double nested list
        outliers = [item for sublist in outliers for item in sublist]

        # apply LCZ safety net
        if use_LCZ_safetynet:
            outliers = apply_LCZ_safety_net(
                outliers=outliers,
                LCZ_buddies=LCZ_buddies,
                wideobsds=combdf,
                safety_std_threshold=safetynet_z_threshold,
                min_sample_size=min_LCZ_safetynet_sample_size,
                min_std=min_std,
            )
            # NOTE: The records that were saved by the safety net, will be tested
            # again in the following iteration. (A different result can occur
            # if the spatial-/savetynet-sample is changed in the next iteration.

        # Add reference to the iteration in the msg of the outliers
        outliers = [
            (station, timestamp, f"{msg} (iteration {i}/{N_iter})")
            for station, timestamp, msg in outliers
        ]

        outliersbin.extend(outliers)
        i += 1

    return outliersbin, timestamp_map


@log_entry
def apply_LCZ_safety_net(
    outliers: list,  # list of tuples
    LCZ_buddies: dict,
    wideobsds: pd.DataFrame,
    safety_std_threshold: Union[int, float],
    min_sample_size: int,
    min_std: Union[int, float],
) -> list:
    """
    Apply the LCZ (Local Climate Zone) safety net to outliers detected by the spatial buddy check.

    For each outlier, this function checks whether the value can be "saved" by comparison
    with its LCZ buddies (stations with the same LCZ and within a certain distance).
    If the outlier's value is within the specified z-threshold when compared to its LCZ buddies,
    it is not considered an outlier for this iteration.

    Parameters
    ----------
    outliers : list of tuple
        List of detected outliers, each as a tuple (station_name, timestamp, message).
    LCZ_buddies : dict
        Dictionary mapping each station to a list of its LCZ buddies.
    wideobsds : pandas.DataFrame
        Wide-format DataFrame with stations as columns and timestamps as index.
    safety_std_threshold : int or float
        Z-value threshold for saving an outlier using the LCZ safety net.
    min_sample_size : int
        Minimum number of LCZ buddies required to apply the safety net.
    min_std : int or float
        Minimum standard deviation to use for z-value calculation.

    Returns
    -------
    list of tuple
        List of outliers that were not saved by the LCZ safety net, each as a tuple
        (station_name, timestamp, message). Outliers that are "saved" are not included
        in the returned list.

    Notes
    -----
    - The LCZ safety net is only applied if there are enough LCZ buddies and non-NaN values.
    - Outliers from previous iterations are already set to NaN in `wideobsds` and are not considered.
    - The function appends a message to the outlier if the safety net is not applied or not passed.
    """

    checked_outliers = []
    for outl in outliers:
        outlstation, outltimestamp, outl_msg = outl

        outl_value = wideobsds.loc[outltimestamp, outlstation]
        outl_LCZ_buddies = LCZ_buddies[outlstation]

        # check if samplesize is suffcient
        if len(outl_LCZ_buddies) < min_sample_size:
            msg = f"Too few LCZ buddies to apply a safety net check ({len(outl_LCZ_buddies)} < {min_sample_size})"
            logger.debug(
                f"skip LCZ-safety net for {outlstation}, too few LCZ buddies({len(outl_LCZ_buddies)} < {min_sample_size})."
            )
            checked_outliers.append((outlstation, outltimestamp, outl_msg + msg))
            continue

        # get LCZ safetynet samples
        # NOTE: The sample is constructd using the wideobsds, thus outliers
        # from the current iteration of the spatial buddy check are not taken into account !!
        # Outliers from the previous iterations are taken into account, since
        # wideobsdf is alterd (nan's are placed at outlier records) in the beginning
        # of each iteration.
        savetynet_samples = wideobsds.loc[outltimestamp][LCZ_buddies[outlstation]]

        # Compute scores
        sample_mean = savetynet_samples.mean()
        sample_std = savetynet_samples.std()
        sample_non_nan_count = savetynet_samples.notna().sum()

        # instantanious sample size check
        if sample_non_nan_count < min_sample_size:
            msg = f"Too few non-nan LCZ buddies to apply a safety net check ({sample_non_nan_count} < {min_sample_size})"
            logger.debug(
                f"skip LCZ-safety net for {outlstation}, too few non-nan LCZ buddies({sample_non_nan_count} < {min_sample_size})."
            )
            checked_outliers.append((outlstation, outltimestamp, outl_msg + msg))
            continue

        # apply min std
        if sample_std < min_std:
            sample_std = min_std  # set minimum std

        # Check if saved
        z_value = abs(((outl_value) - sample_mean) / sample_std)
        if z_value <= safety_std_threshold:
            # is saved
            logger.debug(
                f"{outlstation} at {outltimestamp} is saved by the LCZ safetynet with a z={z_value}."
            )
            # Do not append the current outl to the checked

        else:
            # not saved by the savety net
            msg = f"LCZ safety net check applied but not saved (since z= {z_value} > {safety_std_threshold})"
            checked_outliers.append((outlstation, outltimestamp, outl_msg + msg))
            continue
    logger.debug(
        f"A total of {len(outliers) - len(checked_outliers)} records are saved by the LCZ safety net."
    )
    return checked_outliers


@log_entry
def find_buddy_group_outlier(inputarg: Tuple) -> List[Tuple]:
    """
    Apply a buddy check on a group to identify outliers.

    Parameters
    ----------
    inputarg : tuple
        A tuple containing:

        * buddygroup : list

        * buddygroup : list
            List of station names that form the buddy group.
        * combdf : pandas.DataFrame
        * combdf : pandas.DataFrame
            DataFrame containing the combined data for all stations.
        * min_sample_size : int
        * min_sample_size : int
            Minimum number of non-NaN values required in the buddy group for a
            valid comparison.
        * min_std : float
        * min_std : float
            Minimum standard deviation to use when calculating z-scores.
        * outlier_threshold : float
        * outlier_threshold : float
            Threshold for identifying outliers in terms of z-scores.

    Returns
    -------
    list of tuple
        Each tuple contains:

        * str : The station name of the most extreme outlier.
        * pandas.Timestamp : The timestamp of the outlier.
        * str : A detailed message describing the outlier.

        * str : The station name of the most extreme outlier.
        * pandas.Timestamp : The timestamp of the outlier.
        * str : A detailed message describing the outlier.

    Notes
    -----
    This function performs the following steps:


    1. Subsets the data to the buddy group.
    2. Calculates the mean, standard deviation, and count of non-NaN values
       for each timestamp.
       for each timestamp.
    3. Filters out timestamps with insufficient data.
    4. Replaces standard deviations below the minimum threshold with the
       minimum value.
       minimum value.
    5. Converts station values to z-scores.
    6. Identifies timestamps with at least one outlier.
    7. Locates the most extreme outlier for each timestamp.
    8. Generates a detailed message for each outlier.


    """

    buddygroup, combdf = inputarg[0], inputarg[1]
    min_sample_size, min_std, outlier_threshold = inputarg[2:]

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
        buddydf[station] = (buddydf[station] - buddydf["mean"]).abs() / buddydf["std"]

    # Drop rows for which all values are smaller than the threshold
    # (speed up the last step)
    buddydf["timestamp_with_outlier"] = buddydf[[*buddygroup]].apply(
        lambda row: any(row > outlier_threshold), axis=1
    )
    buddydf = buddydf.loc[buddydf["timestamp_with_outlier"]]

    # locate the most extreme outlier per timestamp
    buddydf["is_the_most_extreme_outlier"] = buddydf[[*buddygroup]].idxmax(axis=1)

    @log_entry
    def msgcreator(row):
        retstr = f"Outlier at {row['is_the_most_extreme_outlier']}"
        retstr += f" with chi value \
{row[row['is_the_most_extreme_outlier']]:.2f},"
        retstr += f" is part of {buddygroup}, with mean: {row['mean']:.2f}, \
std: {row['std']:.2f}. "
        return retstr

    # detail info string
    buddydf["detail_msg"] = buddydf.apply(
        lambda row: msgcreator(row), axis=1, result_type="reduce"
    )

    return list(
        zip(
            buddydf["is_the_most_extreme_outlier"], buddydf.index, buddydf["detail_msg"]
        )
    )
