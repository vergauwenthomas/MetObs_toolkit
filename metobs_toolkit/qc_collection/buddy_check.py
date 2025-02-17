import os
import logging
import concurrent.futures
from math import radians, cos, sin, asin, sqrt

import numpy as np
import pandas as pd


logger = logging.getLogger(__file__)


def _calculate_distance_matrix_with_haverine(metadf):

    def haversine(lon1, lat1, lon2, lat2):
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


def synchronize_series(series_list, max_shift):
    """
    Synchronize a list of pandas Series with datetime indexes.
    The target timestamps are definde by:

     * freq: the highest frequency present in the input series
     * origin: the earlyest timestamp found, rounded down by the freq
     * closing: the latest timestamp found, rounded up by the freq.

    Parameters:
    series_list (list): List of pandas Series with datetime indexes.
    max_shift (pd.Timedelta): Maximum shift in time that can be applied to each timestamp in synchronization.

    Returns:
    pd.DataFrame: DataFrame with synchronized Series.
    """
    # find highest frequency
    frequencies = [pd.Timedelta(s.index.inferred_freq) for s in series_list]
    trg_freq = min(frequencies)

    # find origin and closing timestamp (earliest/latest)
    origin = min([s.index.min() for s in series_list]).floor(trg_freq)
    closing = max([s.index.max() for s in series_list]).ceil(trg_freq)

    # Create target datetime axes
    target_dt = pd.date_range(start=origin, end=closing, freq=trg_freq)

    # Synchronize (merge with tollerance) series to the common index
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


def _find_spatial_buddies(distance_df, buddy_radius):
    """Get neighbouring stations using buddy radius."""
    buddies = {}
    for refstation, distances in distance_df.iterrows():
        bud_stations = distances[distances <= buddy_radius].index.to_list()
        bud_stations.remove(refstation)
        buddies[refstation] = bud_stations

    return buddies


# filter altitude buddies
def _filter_to_altitude_buddies(
    spatial_buddies, altitudes: pd.Series, max_altitude_diff
):
    """Filter neighbours by maximum altitude difference."""
    alt_buddies_dict = {}
    for refstation, buddylist in spatial_buddies.items():
        alt_diff = abs((altitudes.loc[buddylist]) - altitudes.loc[refstation])

        alt_buddies = alt_diff[alt_diff <= max_altitude_diff].index.to_list()
        alt_buddies_dict[refstation] = alt_buddies
    return alt_buddies_dict


def _filter_to_minimum_samplesize(buddydict, min_sample_size):
    """Filter stations that are to isolated using minimum sample size."""
    to_check_stations = {}
    for refstation, buddies in buddydict.items():
        if len(buddies) < min_sample_size:
            # not enough buddies
            to_check_stations[refstation] = []  # remove buddies
        else:
            to_check_stations[refstation] = buddies
    return to_check_stations


def create_groups_of_buddies(buddydict):

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
    dataset,
    obstype="temp",
    buddy_radius=10000.0,
    min_sample_size=4,
    max_alt_diff=None,
    min_std=1.0,
    std_threshold=3.1,
    N_iter=2,
    instantanious_tolerance=pd.Timedelta("4min"),
    lapserate=None,  # -0.0065 for temperarture
    use_mp=True,
):
    """Spatial buddy check.

    The buddy check compares an observation against its neighbors (i.e. buddies). The check loops
    over all the groups, which are stations within a radius of each other. For each group, the absolute
    value of the difference with the groupmean, normalized by the standared deviation (with a defined minimum),
    is computed. The baddest observation of that group is labeled as an outlier, if it exceeds the std_threshold.

    Multiple iterations of this checks can be done using the N_iter.

    A schematic step-by-step description of the buddy check:

      1. A distance matrix is constructed for all interdistances between the stations. This is done using the haversine approximation.
      2. A set of all (spatial) buddies per station is created by filtering out all stations that are too far.
      3. The buddies are further filtered based on altitude differences with respect to the reference station.
      4. Buddy groups are defined as sets with stations near each other.
      5. The timestamps are synchronized towards the highest frequency, earlyest and latest timestamps.
      6. For each buddy group:
        * Observations of buddies are extracted from all observations.
        * These observations are corrected for altitude differences if a lapsrate is specified.
        * For each buddy group, the mean, standard deviation (std), and sample size are computed.
        * If the std is lower than the minimum std, it is replaced by the minimum std.
        * Chi values are calculated for all records.
        * For each timestamp the record with the highest Chi is tested if it is larger then std_threshold.
        If so, that record (stationname + timestamp) is flagged as an outlier.


    Parameters
    ----------
    dataset: metobs_toolkit.Dataset
        The dataset to apply the buddy check on.
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'
    buddy_radius : numeric
        The radius to define neighbors in meters.
    min_sample_size : int
        The minimum sample size to calculate statistics on.
    max_alt_diff : numeric or None
        The maximum altitude difference allowed for buddies. I None,
        altitude is not taken into account.
    min_std : numeric
        The minimum standard deviation for sample statistics. This should
        represent the accuracy of the observations.
    std_threshold : numeric
        The threshold (std units) for flagging observations as outliers.
    haversine_approx : bool, optional
        Use the haversine approximation (earth is a sphere) to calculate
        distances between stations. The default is True.
    metric_epsg : str, optional
        EPSG code for the metric CRS to calculate distances in. Only used when
        haversine approximation is set to False. Thus becoming a better
        distance approximation but not globally applicable The default is '31370'
        (which is suitable for Belgium).
    lapserate : numeric, optional
        Describe how the obstype changes with altitude (in meters). The default is -0.0065.

    Returns
    -------
    obsdf: Pandas.DataFrame
        The dataframe containing the unflagged-observations
    outlier_df : Pandas.DataFrame
        The dataframe containing the flagged observations

    """

    # -----  Part 1: construct buddy groups ------
    # copute distance metric
    metadf = dataset.metadf
    dist_matrix = _calculate_distance_matrix_with_haverine(metadf)

    # find potential buddies by distance
    buddies = _find_spatial_buddies(distance_df=dist_matrix, buddy_radius=buddy_radius)

    # filter buddies by altitude difference
    if max_alt_diff is None:
        pass
    else:
        if metadf["altitude"].isna().any():
            raise ValueError("At least one station has a NaN value for 'altitude'")
        # Filter by altitude difference
        buddies = _filter_to_altitude_buddies(
            spatial_buddies=buddies,
            altitudes=metadf["altitude"],
            max_altitude_diff=max_alt_diff,
        )

    # Filter by samplesize (based on the number of buddy stations)
    buddies = _filter_to_minimum_samplesize(
        buddydict=buddies, min_sample_size=min_sample_size
    )

    # create unique groups of buddies (list of tuples)
    buddygroups = create_groups_of_buddies(buddies)

    # ---- Part 2: Preparing the records  -----

    # construct a wide observation dataframe
    concatlist = []
    for sta in dataset.stations:
        if obstype in sta.obsdata.keys():
            records = sta.obsdata[obstype].series
            records.name = sta.name
            concatlist.append(records)

    # synchronize the timestamps
    combdf, timestamp_map = synchronize_series(
        series_list=concatlist, max_shift=instantanious_tolerance
    )

    # lapsrate correction
    if lapserate is None:
        pass
    else:
        # get altitude dataframe
        altdict = {sta.name: sta.site.altitude for sta in dataset.stations}
        altseries = pd.Series(altdict)
        altcorrectionseries = (altseries - altseries.min()) * lapserate
        combdf = combdf - altcorrectionseries  # Correct for altitude

    # ---- Part 3 : Apply buddy check on each group, rejecting the most extreme outlier

    outliersbin = []
    for i in range(N_iter):
        # convert values to Nan, if they are labeled as outlier in previous iteration
        if bool(outliersbin):
            for outlier_station, outlier_time in outliersbin:
                if outlier_station in combdf.columns:
                    combdf.loc[outlier_time, outlier_station] = np.nan

        if use_mp:
            # Use multiprocessing generatore (parralelization)
            num_cpus = os.cpu_count()
            # since this check is an instantanious check --> perfect for splitting the dataset in chunks in time
            chunks = np.array_split(combdf, num_cpus)

            # create inputargs for each buddygroup, an for each chunk in time
            inputargs = [
                (buddygroup, chunk, min_sample_size, min_std, std_threshold)
                for buddygroup in buddygroups
                for chunk in chunks
            ]

            with concurrent.futures.ProcessPoolExecutor() as executor:
                outliers = executor.map(find_buddy_group_outlier, inputargs)
            outliers = list(outliers)

        else:
            # create inputargs for each buddygroup, an for each chunk in time
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


def find_buddy_group_outlier(inputarg):
    """apply buddy check on a group"""

    # unpack arguments
    buddygroup, combdf = inputarg[0], inputarg[1]
    min_sample_size, min_std, outlier_threshold = inputarg[2:]

    # subset to the buddies
    buddydf = combdf[[*buddygroup]]

    # calucalate std and mean row wise
    buddydf["mean"] = buddydf[[*buddygroup]].mean(axis=1)
    buddydf["std"] = buddydf[[*buddygroup]].std(axis=1)
    buddydf["non_nan_count"] = buddydf[[*buddygroup]].notna().sum(axis=1)

    # subset to samples with enough members (check for each timesstap specifically)
    buddydf = buddydf.loc[buddydf["non_nan_count"] >= min_sample_size]

    # replace std by minimum, if needed
    buddydf.loc[buddydf["std"] < min_std, "std"] = np.float32(min_std)

    # Convert values to sigmas
    for station in buddygroup:
        buddydf[station] = (buddydf[station] - buddydf["mean"]).abs() / buddydf["std"]

    # Drop rows for which all values are smaller then the treshold (speedup the last step)
    buddydf["timestamp_with_outlier"] = buddydf[[*buddygroup]].apply(
        lambda row: any(row > outlier_threshold), axis=1
    )
    buddydf = buddydf.loc[buddydf["timestamp_with_outlier"]]

    # locate the most extreme outlier per timestamp
    buddydf["is_the_most_extreme_outlier"] = buddydf[[*buddygroup]].idxmax(axis=1)

    # add (outlierstationname: timestamp) to outliers
    return list(zip(buddydf["is_the_most_extreme_outlier"], buddydf.index))


# def toolkit_buddy_check(
#     obsdf,
#     metadf,
#     obstype,
#     buddy_radius,
#     min_sample_size,
#     max_alt_diff,
#     min_std,
#     std_threshold,
#     haversine_approx=True,
#     metric_epsg="31370",
#     lapserate=-0.0065,
# ):
#     """Spatial buddy check.

#     The buddy check compares an observation against its neighbors (i.e. buddies). The check looks for
#     buddies in a neighborshood specified by a certain radius. The buddy check flags observations if the
#     (absolute value of the) difference between the observations and the average of the neighbors
#     normalized by the standard deviation in the circle is greater than a predefined threshold.

#     A schematic step-by-step description of the buddy check:

#       1. A distance matrix is constructed for all interdistances between the stations. This is done using the haversine approximation, or by first converting the Coordinate Reference System (CRS) to a metric one, specified by an EPSG code.
#       2. A set of all (spatial) buddies per station is created by filtering out all stations that are too far.
#       3. The buddies are further filtered based on altitude differences with respect to the reference station.
#       4. For each station:
#         * Observations of buddies are extracted from all observations.
#         * These observations are corrected for altitude differences by assuming a constant lapse rate.
#         * For each reference record, the mean, standard deviation (std), and sample size of the corrected buddies’ observations are computed.
#         * If the std is lower than the minimum std, it is replaced by the minimum std.
#         * Chi values are calculated for all reference records.
#         * If the Chi value is larger than the std_threshold, the record is accepted, otherwise, it is marked as an outlier.

#     Parameters
#     ----------
#     obsdf: Pandas.DataFrame
#         The dataframe containing the observations
#     metadf: Pandas.DataFrame
#         The dataframe containing the metadata (e.g. latitude, longitude...)
#     obstype: String, optional
#         The observation type that has to be checked. The default is 'temp'
#     buddy_radius : numeric
#         The radius to define neighbors in meters.
#     min_sample_size : int
#         The minimum sample size to calculate statistics on.
#     max_alt_diff : numeric
#         The maximum altitude difference allowed for buddies.
#     min_std : numeric
#         The minimum standard deviation for sample statistics. This should
#         represent the accuracy of the observations.
#     std_threshold : numeric
#         The threshold (std units) for flagging observations as outliers.
#     haversine_approx : bool, optional
#         Use the haversine approximation (earth is a sphere) to calculate
#         distances between stations. The default is True.
#     metric_epsg : str, optional
#         EPSG code for the metric CRS to calculate distances in. Only used when
#         haversine approximation is set to False. Thus becoming a better
#         distance approximation but not globally applicable The default is '31370'
#         (which is suitable for Belgium).
#     lapserate : numeric, optional
#         Describe how the obstype changes with altitude (in meters). The default is -0.0065.

#     Returns
#     -------
#     obsdf: Pandas.DataFrame
#         The dataframe containing the unflagged-observations
#     outlier_df : Pandas.DataFrame
#         The dataframe containing the flagged observations

#     """
#     outliers_idx = init_multiindex()

#     # Get spatial buddies for each station
#     if haversine_approx:
#         distance_df = _calculate_distance_matrix_with_haverine(metadf=metadf)
#     else:
#         distance_df = _calculate_distance_matrix(metadf=metadf, metric_epsg=metric_epsg)
#     buddies = _find_spatial_buddies(distance_df=distance_df, buddy_radius=buddy_radius)

#     # Filter by altitude difference
#     buddies = _filter_to_altitude_buddies(
#         spatial_buddies=buddies, metadf=metadf, max_altitude_diff=max_alt_diff
#     )

#     # Filter by samplesize
#     buddydict = _filter_to_samplesize(
#         buddydict=buddies, min_sample_size=min_sample_size
#     )

#     to_check_obsdf = xs_save(obsdf, obstype, "obstype")
#     # Apply buddy check station per station
#     for refstation, buddies in buddydict.items():
#         if len(buddies) == 0:
#             logger.debug(f"{refstation} has not enough suitable buddies.")
#             continue

#         # Get observations
#         buddies_obs = to_check_obsdf[
#             to_check_obsdf.index.get_level_values("name").isin(buddies)
#         ]["value"]
#         # Unstack
#         buddies_obs = buddies_obs.unstack(level="name")

#         # Make lapsrate correction:
#         ref_alt = metadf.loc[refstation, "altitude"]
#         buddy_correction = (
#             (metadf.loc[buddies, "altitude"] - ref_alt) * (-1.0 * lapserate)
#         ).to_dict()
#         for bud in buddies_obs.columns:
#             buddies_obs[bud] = buddies_obs[bud] - buddy_correction[bud]

#         # calucalate std and mean row wise
#         buddies_obs["mean"] = buddies_obs[buddies].mean(axis=1)
#         buddies_obs["std"] = buddies_obs[buddies].std(axis=1)
#         buddies_obs["samplesize"] = buddies_obs[buddies].count(axis=1)

#         # from titan they use std adjust which is float std_adjusted = sqrt(variance + variance / n_buddies);
#         # This is not used
#         # buddies_obs['var'] = buddies_obs[buddies].var(axis=1)
#         # buddies_obs['std_adj'] =np.sqrt(buddies_obs['var'] + buddies_obs['var']/buddies_obs['samplesize'])

#         # replace where needed with min std
#         buddies_obs["std"] = buddies_obs["std"].where(
#             cond=buddies_obs["std"] >= min_std, other=min_std
#         )

#         # Get refstation observations and merge
#         ref_obs = xs_save(to_check_obsdf, refstation, "name", drop_level=False)[
#             "value"
#         ].unstack(level="name")

#         buddies_obs = buddies_obs.merge(
#             ref_obs,
#             how="left",  # both not needed because if right, than there is no buddy sample per definition.
#             left_index=True,
#             right_index=True,
#         )
#         # Calculate sigma
#         buddies_obs["chi"] = (
#             abs(buddies_obs["mean"] - buddies_obs[refstation])
#         ) / buddies_obs["std"]

#         outliers = buddies_obs[
#             (buddies_obs["chi"] > std_threshold)
#             & (buddies_obs["samplesize"] >= min_sample_size)
#         ]

#         logger.debug(f" Buddy outlier details for {refstation}: \n {buddies}")
#         # NOTE: the outliers (above) can be interesting to pass back to the dataset??

#         # to multiindex
#         outliers["name"] = refstation
#         outliers = outliers.reset_index().set_index(["name", "datetime"]).index
#         outliers_idx = outliers_idx.append(outliers)

#     # Update the outliers and replace the obsdf
#     obsdf, outlier_df = make_outlier_df_for_check(
#         station_dt_list=outliers_idx,
#         obsdf=obsdf,
#         obstype=obstype,
#         flag=label_def["buddy_check"]["label"],
#     )

#     return obsdf, outlier_df
