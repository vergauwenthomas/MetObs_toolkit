from __future__ import annotations

import os
import logging
import concurrent.futures
from typing import Union, List, Dict, Tuple, TYPE_CHECKING

import numpy as np
import pandas as pd

from metobs_toolkit.backend_collection.datetime_collection import to_timedelta
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.qc_collection.distancematrix_func import generate_distance_matrix
from .whitelist import WhiteSet

if TYPE_CHECKING:
    from metobs_toolkit.station import Station

logger = logging.getLogger("<metobs_toolkit>")


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


def _validate_safety_net_configs(safety_net_configs: List[Dict]) -> None:
    """
    Validate that all required keys are present in safety_net_configs.

    Parameters
    ----------
    safety_net_configs : list of dict
        List of safety net configuration dictionaries.

    Raises
    ------
    ValueError
        If safety_net_configs is not a list or contains non-dict elements.
    KeyError
        If any required key is missing from a safety net configuration.
    """
    if safety_net_configs is None:
        return None

    required_keys = {"category", "buddy_radius", "z_threshold", "min_sample_size"}

    if not isinstance(safety_net_configs, list):
        raise ValueError(
            f"safety_net_configs must be a list, got {type(safety_net_configs).__name__}"
        )

    for i, config in enumerate(safety_net_configs):
        if not isinstance(config, dict):
            raise ValueError(
                f"Each safety net config must be a dict, but config at index {i} "
                f"is {type(config).__name__}"
            )

        missing_keys = required_keys - set(config.keys())
        if missing_keys:
            raise KeyError(
                f"Safety net config at index {i} is missing required key(s): "
                f"{', '.join(sorted(missing_keys))}. "
                f"Required keys are: {', '.join(sorted(required_keys))}"
            )


def _find_buddies_by_distance(
    distance_df: pd.DataFrame, buddy_radius: Union[int, float]
) -> Dict:
    """
    Get neighbouring stations using buddy radius (distance only).

    This is the core distance-based buddy finding function used internally
    by other buddy-finding functions.

    Parameters
    ----------
    distance_df : pandas.DataFrame
        DataFrame containing distances between stations.
    buddy_radius : int or float
        Maximum distance (in meters) to consider as a buddy.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of its buddies within the radius.
    """

    buddies = {}
    for refstation, distances in distance_df.iterrows():
        bud_stations = distances[distances <= buddy_radius].index.to_list()
        bud_stations.remove(refstation)
        buddies[refstation] = bud_stations

    return buddies


def _find_category_buddies(
    metadf: pd.DataFrame,
    category_column: str,
    max_dist: Union[int, float],
    distance_df: pd.DataFrame,
) -> Dict:
    """
    Get neighbouring stations using a categorical column and spatial distance.

    This function identifies buddy stations that share the same categorical
    value (e.g., LCZ, network, region) and are within a specified distance.

    Parameters
    ----------
    metadf : pandas.DataFrame
        DataFrame containing metadata for stations. Must include the specified
        category column.
    category_column : str
        The name of the categorical column to group stations by (e.g., 'LCZ',
        'network', 'region').
    max_dist : int or float
        Maximum distance (in meters) to consider as a category buddy.
    distance_df : pandas.DataFrame
        DataFrame containing distances between stations.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of its category buddies that
        are also within the specified distance.

    Notes
    -----
    - Category buddies are stations with the same category value as the reference
      station.
    - The final buddies are the intersection of category buddies and spatial
      buddies within `max_dist`.
    - Stations with NaN values in the category column are handled: they will not
      match with any other station (including other NaN stations).
    """
    category_buddies = {}
    # Find buddies by category
    for refstation in metadf.index:
        ref_category = metadf.loc[refstation, category_column]
        # Handle NaN values - they should not match with anything
        if pd.isna(ref_category):
            logger.warning(
                "Station %s has NaN value for category '%s' - no category buddies assigned",
                refstation,
                category_column,
            )
            category_buddies[refstation] = []
        else:
            ref_buddies = metadf.loc[
                metadf[category_column] == ref_category
            ].index.to_list()
            category_buddies[refstation] = ref_buddies

    # Find buddies by distance
    spatial_buddies = _find_buddies_by_distance(distance_df, max_dist)

    # Intersection of both buddy definitions
    final_buddies = {}
    for refstation in category_buddies.keys():
        final_buddies[refstation] = list(
            set(category_buddies[refstation]).intersection(
                set(spatial_buddies[refstation])
            )
        )

    return final_buddies


def _find_spatial_buddies(
    distance_df: pd.DataFrame,
    metadf: pd.DataFrame,
    buddy_radius: Union[int, float],
) -> Dict:
    """
    Get neighbouring stations using buddy radius (spatial distance only).

    This function is a wrapper around `_find_category_buddies` that finds
    buddies based purely on spatial distance, without any categorical
    constraints. It works by creating a dummy category column where all
    stations have the same value.

    Parameters
    ----------
    distance_df : pandas.DataFrame
        DataFrame containing distances between stations.
    metadf : pandas.DataFrame
        DataFrame containing metadata for stations. The index should be
        station names.
    buddy_radius : int or float
        Maximum distance (in meters) to consider as a buddy.

    Returns
    -------
    dict
        Dictionary mapping each station to a list of its buddies within
        the specified radius.

    See Also
    --------
    _find_category_buddies : Find buddies by category and distance.
    _find_buddies_by_distance : Core distance-based buddy finding function.
    """

    # Create a temporary metadf with a dummy category column where all
    # stations have the same value, so _find_category_buddies will not
    # filter by category
    temp_metadf = metadf.copy()
    temp_metadf["_all_same_category"] = 1

    return _find_category_buddies(
        metadf=temp_metadf,
        category_column="_all_same_category",
        max_dist=buddy_radius,
        distance_df=distance_df,
    )


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
    target_stations: list[Station],
    metadf: pd.DataFrame,
    obstype: str,
    spatial_buddy_radius: Union[int, float],
    spatial_min_sample_size: int,
    max_alt_diff: Union[int, float, None],
    min_std: Union[int, float],
    spatial_z_threshold: Union[int, float],
    N_iter: int,
    instantaneous_tolerance: pd.Timedelta,
    # Whitelist arguments
    whiteset: WhiteSet,
    # Safety nets
    safety_net_configs: List[Dict] = None,
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

    Optionally, one or more safety nets can be applied. A safety net tests
    potential outliers against a sample of stations that share a categorical
    attribute (e.g., LCZ, network). If the z-value computed using the safety
    net sample is below the specified threshold, the outlier is "saved" and
    removed from the outlier set for the current iteration.

    Safety nets are applied in the order they are specified, allowing for
    multi-level filtering (e.g., first test against LCZ buddies, then against
    network buddies).

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

       #. The values of outliers flagged by a previous iteration are converted to
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

       #. If `safety_net_configs` is provided, the following steps are applied
          on the outliers flagged by the current iteration, for each safety net
          in order:

          * Category buddies (stations sharing the same category value within
            the specified radius) are identified.
          * The safety net sample is tested in size (sample size must be at
            least `min_sample_size`). If the condition is not met, the safety
            net test is not applied.
          * The safety net test is applied:

            * The mean and std are computed of the category-buddy sample. If
              the std is smaller than `min_std`, the latter is used.
            * The z-value is computed for the target record (= flagged outlier).
            * If the z-value is smaller than the safety net's `z_threshold`,
              the tested outlier is "saved" and removed from the set of outliers
              for the current iteration.

       #. If `whiteset` contains records, any outliers that match the white-listed
          timestamps (and optionally station names) are removed from the outlier set
          for the current iteration. White-listed records participate in all buddy
          check calculations but are not flagged as outliers in the final results.

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
    whiteset : WhiteSet
        A WhiteSet instance containing records that should be excluded from
        outlier detection. Records in the WhiteSet undergo the buddy check
        iterations as regular records but are removed from the outlier set
        at the end of each iteration.
    safety_net_configs : list of dict, optional
        List of safety net configurations to apply in order. Each dict must
        contain:

        * 'category': str, the metadata column name to group by (e.g., 'LCZ',
          'network')
        * 'buddy_radius': int or float, maximum distance for category buddies
          (in meters)
        * 'z_threshold': int or float, z-value threshold for saving outliers
        * 'min_sample_size': int, minimum number of buddies required for the
          safety net test

        The default is None.
    lapserate : float or None, optional
        Describes how the obstype changes with altitude (in meters). If
        None, no altitude correction is applied. For temperature, a
        common value is -0.0065.
    use_mp : bool, optional
        Use multiprocessing to speed up the buddy check. The default is True.

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

    # Validate safety net configs if provided
    _validate_safety_net_configs(safety_net_configs)

    # -----  Part 1: construct buddy groups ------
    # compute distance metric
    logger.debug("Calculating distance matrix with Haversine formula")
    dist_matrix = generate_distance_matrix(metadf)

    # find potential buddies by distance
    logger.debug(
        "Finding spatial buddies within radius of %s meters", spatial_buddy_radius
    )
    spatial_buddies = _find_spatial_buddies(
        distance_df=dist_matrix, metadf=metadf, buddy_radius=spatial_buddy_radius
    )

    # filter buddies by altitude difference
    if max_alt_diff is not None:
        logger.debug(
            "Filtering buddies by maximum altitude difference of %s meters",
            max_alt_diff,
        )
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

    # Filter by sample size (based on the number of buddy stations)
    logger.debug(
        "Filtering buddies by minimum sample size of %s", spatial_min_sample_size
    )
    spatial_buddies = _filter_to_minimum_samplesize(
        buddydict=spatial_buddies, min_sample_size=spatial_min_sample_size
    )

    # create unique groups of buddies (list of tuples)
    logger.debug("Creating groups of buddies")
    buddygroups = create_groups_of_buddies(spatial_buddies)
    logger.debug("Number of buddy groups created: %s", len(buddygroups))

    # ---- Part 2: Preparing the records  -----

    # construct a wide observation dataframe
    logger.debug("Constructing wide observation DataFrame for obstype: %s", obstype)
    concatlist = []
    for sta in target_stations:
        if obstype in sta.sensordata.keys():
            records = sta.get_sensor(obstype).series
            records.name = sta.name
            concatlist.append(records)

    # synchronize the timestamps
    logger.debug("Synchronizing timestamps")
    combdf, timestamp_map = synchronize_series(
        series_list=concatlist, max_shift=instantaneous_tolerance
    )

    # lapse rate correction
    if lapserate is not None:
        logger.debug("Applying lapse rate correction with rate: %s", lapserate)
        # get altitude dataframe
        altdict = {sta.name: sta.site.altitude for sta in target_stations}
        altseries = pd.Series(altdict)
        altcorrectionseries = (altseries - altseries.min()) * lapserate
        combdf = combdf - altcorrectionseries  # Correct for altitude

    # ---- Part 3 : Apply buddy check on each group,
    #  rejecting the most extreme outlier

    outliersbin = []
    for i in range(N_iter):
        logger.debug("Starting iteration %s of %s", i + 1, N_iter)
        # convert values to NaN, if they are labeled as outlier in
        #  previous iteration
        if bool(outliersbin):
            logger.debug("Converting previous-iteration outliers to NaN")
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

            logger.debug("Finding outliers in each buddy group")
            outliers = list(map(find_buddy_group_outlier, inputargs))

        # unpack double nested list
        outliers = [item for sublist in outliers for item in sublist]

        # Apply safety nets (if configured)
        if safety_net_configs:
            logger.debug(
                "Applying %s safety net(s) to %s outliers",
                len(safety_net_configs),
                len(outliers),
            )
            outliers = apply_safetynets(
                outliers=outliers,
                safety_net_configs=safety_net_configs,
                wideobsds=combdf,
                metadf=metadf,
                distance_df=dist_matrix,
                max_alt_diff=max_alt_diff,
                min_std=min_std,
            )
            # NOTE: Records saved by any safety net will be tested again in
            # the following iteration. A different result can occur if the
            # spatial/safety net sample changes in the next iteration.

        # Save white-listed records
        outliers = save_whitelist_records(
            outliers=outliers,
            whiteset=whiteset,
            obstype=obstype,
        )
        # NOTE: The white-listed records are removed from the outliers at the end
        # of each iteration, similar to the safety nets. They participate in
        # the buddy check calculations but are not flagged as outliers.

        # Save white-listed records
        outliers = save_whitelist_records(
            outliers=outliers,
            whiteset=whiteset,
            obstype=obstype,
        )
        # NOTE: The white-listed records are removed from the outliers at the end
        # of each iteration, similar to the LCZ safety net. They participate in
        # the buddy check calculations but are not flagged as outliers.

        # Add reference to the iteration in the msg of the outliers
        outliers = [
            (station, timestamp, f"{msg} (iteration {i+1}/{N_iter})")
            for station, timestamp, msg in outliers
        ]

        outliersbin.extend(outliers)
        i += 1

    return outliersbin, timestamp_map


@log_entry
def apply_safety_net(
    outliers: list,
    category_buddies: dict,
    wideobsds: pd.DataFrame,
    safety_z_threshold: Union[int, float],
    min_sample_size: int,
    min_std: Union[int, float],
    category_name: str,
) -> list:
    """
    Apply a category-based safety net to outliers detected by the spatial buddy check.

    This function works with any categorical grouping (e.g., LCZ, network, region).

    For each outlier, this function checks whether the value can be "saved" by
    comparison with its category buddies (stations with the same category value
    and within a certain distance). If the outlier's value is within the specified
    z-threshold when compared to its category buddies, it is not considered an
    outlier for this iteration.

    Parameters
    ----------
    outliers : list of tuple
        List of detected outliers, each as a tuple (station_name, timestamp, message).
    category_buddies : dict
        Dictionary mapping each station to a list of its category buddies.
    wideobsds : pandas.DataFrame
        Wide-format DataFrame with stations as columns and timestamps as index.
    safety_z_threshold : int or float
        Z-value threshold for saving an outlier using the safety net.
    min_sample_size : int
        Minimum number of category buddies required to apply the safety net.
    min_std : int or float
        Minimum standard deviation to use for z-value calculation.
    category_name : str
        Name of the category being used (for logging and messages).

    Returns
    -------
    list of tuple
        List of outliers that were not saved by the safety net, each as a tuple
        (station_name, timestamp, message). Outliers that are "saved" are not
        included in the returned list.

    Notes
    -----
    - The safety net is only applied if there are enough category buddies and
      non-NaN values.
    - Outliers from previous iterations are already set to NaN in `wideobsds`
      and are not considered.
    - The function appends a message to the outlier if the safety net is not
      applied or not passed.
    """
    checked_outliers = []
    for outl in outliers:
        outlstation, outltimestamp, outl_msg = outl

        outl_value = wideobsds.loc[outltimestamp, outlstation]
        outl_category_buddies = category_buddies.get(outlstation, [])

        # Check if sample size is sufficient
        if len(outl_category_buddies) < min_sample_size:
            msg = f"Too few {category_name} buddies to apply safety net ({len(outl_category_buddies)} < {min_sample_size})."
            logger.debug(
                "Skip %s safety net for %s: too few buddies (%s < %s).",
                category_name,
                outlstation,
                len(outl_category_buddies),
                min_sample_size,
            )
            checked_outliers.append((outlstation, outltimestamp, outl_msg + msg))
            continue

        # Get safety net samples
        # NOTE: The sample is constructed using wideobsds, thus outliers
        # from the current iteration are not taken into account!
        # Outliers from previous iterations are taken into account since
        # wideobsdf is altered (NaNs placed at outlier records) at the beginning
        # of each iteration.
        safetynet_samples = wideobsds.loc[outltimestamp][outl_category_buddies]

        # Compute scores
        sample_mean = safetynet_samples.mean()
        sample_std = safetynet_samples.std()
        sample_non_nan_count = safetynet_samples.notna().sum()

        # Instantaneous sample size check
        if sample_non_nan_count < min_sample_size:
            msg = f"Too few non-NaN {category_name} buddies ({sample_non_nan_count} < {min_sample_size})."
            logger.debug(
                "Skip %s safety net for %s: too few non-NaN buddies (%s < %s).",
                category_name,
                outlstation,
                sample_non_nan_count,
                min_sample_size,
            )
            checked_outliers.append((outlstation, outltimestamp, outl_msg + msg))
            continue

        # Apply min std
        if sample_std < min_std:
            sample_std = min_std

        # Check if saved
        z_value = abs((outl_value - sample_mean) / sample_std)
        if z_value <= safety_z_threshold:
            # Is saved
            logger.debug(
                "%s at %s is saved by %s safety net with z=%.2f.",
                outlstation,
                outltimestamp,
                category_name,
                z_value,
            )
            # Do not append the current outl to checked (it's saved)
        else:
            # Not saved by the safety net
            msg = f"{category_name} safety net applied but not saved (z={z_value:.2f} > {safety_z_threshold})."
            checked_outliers.append((outlstation, outltimestamp, outl_msg + msg))
            continue

    n_saved = len(outliers) - len(checked_outliers)
    logger.debug(
        "A total of %s records are saved by the %s safety net.", n_saved, category_name
    )
    return checked_outliers


@log_entry
def apply_safetynets(
    outliers: list,
    safety_net_configs: List[Dict],
    wideobsds: pd.DataFrame,
    metadf: pd.DataFrame,
    distance_df: pd.DataFrame,
    max_alt_diff: Union[int, float, None],
    min_std: Union[int, float],
) -> list:
    """
    Apply multiple safety nets in sequence to outliers.

    Each safety net is defined by a category column, buddy radius, z-threshold,
    and minimum sample size. Outliers are tested against each safety net in order,
    and if saved by any of them, they are removed from the outlier list.

    Parameters
    ----------
    outliers : list of tuple
        List of detected outliers, each as a tuple (station_name, timestamp, message).
    safety_net_configs : list of dict
        List of safety net configurations. Each dict must contain:
        - 'category': str, the metadata column name to group by
        - 'buddy_radius': int or float, maximum distance for category buddies
        - 'z_threshold': int or float, z-value threshold for saving outliers
        - 'min_sample_size': int, minimum number of buddies required
    wideobsds : pandas.DataFrame
        Wide-format DataFrame with stations as columns and timestamps as index.
    metadf : pandas.DataFrame
        DataFrame containing station metadata.
    distance_df : pandas.DataFrame
        DataFrame containing distances between stations.
    max_alt_diff : int, float, or None
        Maximum altitude difference allowed for buddies. If None, no altitude
        filter is applied.
    min_std : int or float
        Minimum standard deviation for sample statistics.

    Returns
    -------
    list of tuple
        List of outliers that were not saved by any safety net.
    """
    if not safety_net_configs:
        return outliers

    current_outliers = outliers

    for config in safety_net_configs:
        category = config["category"]
        buddy_radius = config["buddy_radius"]
        z_threshold = config["z_threshold"]
        min_sample_size = config["min_sample_size"]

        logger.debug(
            "Applying %s safety net (radius=%s, z_threshold=%s, min_sample=%s)",
            category,
            buddy_radius,
            z_threshold,
            min_sample_size,
        )

        # Find category buddies
        category_buddies = _find_category_buddies(
            metadf=metadf,
            category_column=category,
            max_dist=buddy_radius,
            distance_df=distance_df,
        )

        # Filter by altitude difference if specified
        if max_alt_diff is not None:
            category_buddies = _filter_to_altitude_buddies(
                buddies=category_buddies,
                altitudes=metadf["altitude"],
                max_altitude_diff=max_alt_diff,
            )

        # Apply the safety net
        current_outliers = apply_safety_net(
            outliers=current_outliers,
            category_buddies=category_buddies,
            wideobsds=wideobsds,
            safety_z_threshold=z_threshold,
            min_sample_size=min_sample_size,
            min_std=min_std,
            category_name=category,
        )

    return current_outliers


@log_entry
def save_whitelist_records(
    outliers: list,
    whiteset: WhiteSet,
    obstype: str,
) -> list:
    """Remove whitelisted records from the outlier list.

    This function filters out any outliers that are present in the WhiteSet.
    Whitelisted records are known valid observations that should not be flagged
    as outliers, even if they are detected by the buddy check.

    Parameters
    ----------
    outliers : list of tuple
        List of detected outliers, each as a tuple (station_name, timestamp, message).
    whiteset : WhiteSet
        A WhiteSet instance containing records that should be excluded from outlier
        detection. The WhiteSet is converted to station-specific and obstype-specific
        SensorWhiteSet instances for each station in the outliers list.
    obstype : str
        The observation type being checked. Used to filter the whiteset for the
        target obstype.

    Returns
    -------
    list of tuple
        List of outliers excluding those that are whitelisted. Each tuple contains
        (station_name, timestamp, message).

    Notes
    -----
    * Whitelisted records undergo the buddy check iterations as if they are regular
      records.
    * Only at the end of each iteration are they filtered out from the outliers list.
    * This allows whitelisted records to still influence the statistics of their
      buddy groups.
    * The function processes each station separately by creating a SensorWhiteSet
      for each station-obstype combination.
    """

    outldf = pd.DataFrame(outliers, columns=["name", "datetime", "message"])

    for outlsta in outldf["name"].unique():
        # Create a sensorwhiteset for each station
        sensorwhiteset = whiteset.create_sensorwhitelist(
            stationname=outlsta, obstype=obstype
        )
        # get the white-listed datetimes for the station
        outliers_dts = sensorwhiteset.catch_white_records(
            outliers_idx=pd.DatetimeIndex(
                data=outldf[outldf["name"] == outlsta]["datetime"], name="datetime"
            )
        )

        # subset to the saved outliers
        outldf = outldf.drop(
            outldf[
                (outldf["name"] == outlsta) & (~outldf["datetime"].isin(outliers_dts))
            ].index
        )

    # convert back to a list of tuples (name, datetime, message)
    outliers = list(outldf.itertuples(index=False, name=None))
    return outliers


@log_entry
def find_buddy_group_outlier(inputarg: Tuple) -> List[Tuple]:
    """
    Apply a buddy check on a group to identify outliers.

    Parameters
    ----------
    inputarg : tuple
        A tuple containing:

        * buddygroup : list
            List of station names that form the buddy group.
        * combdf : pandas.DataFrame
            DataFrame containing the combined data for all stations.
        * min_sample_size : int
            Minimum number of non-NaN values required in the buddy group for a
            valid comparison.
        * min_std : float
            Minimum standard deviation to use when calculating z-scores.
        * outlier_threshold : float
            Threshold for identifying outliers in terms of z-scores.

    Returns
    -------
    list of tuple
        Each tuple contains:

        * str : The station name of the most extreme outlier.
        * pandas.Timestamp : The timestamp of the outlier.
        * str : A detailed message describing the outlier.

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
        """
        Create a detailed message describing an outlier.

        Parameters
        ----------
        row : pandas.Series
            A row from the buddy DataFrame containing outlier information,
            including 'is_the_most_extreme_outlier', 'mean', and 'std' columns.

        Returns
        -------
        str
            Formatted message describing the outlier with its z-score and
            buddy group statistics.
        """
        retstr = f"Outlier at {row['is_the_most_extreme_outlier']}"
        retstr += f" with chi value \
{row[row['is_the_most_extreme_outlier']]:.2f},"
        retstr += (
            f" is part of {sorted(buddygroup)}, with mean: {row['mean']:.2f}, "
            f"std: {row['std']:.2f}. "
        )
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
