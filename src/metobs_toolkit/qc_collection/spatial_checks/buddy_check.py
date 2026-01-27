from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING

import os
import numpy as np
from concurrent.futures import ProcessPoolExecutor
import itertools
import pandas as pd



from metobs_toolkit.backend_collection.datetime_collection import to_timedelta
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.qc_collection.distancematrix_func import generate_distance_matrix

from metobs_toolkit.qcresult import QCresult, flagged_cond
from .buddywrapstation import BuddyCheckStation, to_qc_labels_map, reconstruct_fractured_targets
from metobs_toolkit.settings_collection import Settings
from ..whitelist import WhiteSet
# Import methods
from . import methods as buddymethods



if TYPE_CHECKING:
    from metobs_toolkit.station import Station

logger = logging.getLogger("<metobs_toolkit>")


def _run_buddy_test(kwargs):
    #executer for mutliprocessing
    return buddymethods.buddy_test_a_station(**kwargs)

@log_entry
def toolkit_buddy_check(
    target_stations: list[Station],
    metadf: pd.DataFrame,
    obstype: str,
    spatial_buddy_radius: Union[int, float],
    spatial_min_sample_size: int,
    max_alt_diff: Union[int, float, None],
    min_sample_spread: Union[int, float],
    spatial_z_threshold: Union[int, float],
    N_iter: int,
    instantaneous_tolerance: pd.Timedelta,
    # Whitelist arguments
    whiteset: WhiteSet,
    # Safety nets
    safety_net_configs: List[Dict] = None,
    #Statistical
    use_z_robust_method: bool = True,
    # Technical
    lapserate: Union[float, None] = None,  # -0.0065 for temperature
    use_mp: bool = True,
) -> List[QCresult]:
    """
    #TODO update the docstring accordingly
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

            * If use_z_robust_method is True:

              * The mean and std are computed of the category-buddy sample. If
                the std is smaller than `min_sample_spread`, the latter is used.
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
    min_sample_spread : int or float
        The minimum sample spread for sample statistics. When use_z_robust_method is True,
        this is the equal to the minimum MAD to use (avoids division by near-zero). When
        use_z_robust_method is False, this is the standard deviation. This parameter helps
        to represent the accuracy of the observations.
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
        A dictionary mapping station names to BuddyCheckStationDetails objects
        containing detailed tracking information for each timestamp.


    Notes
    -----

    * The altitude of the stations can be extracted from GEE by using the
      `Dataset.get_altitude()` method.
    * The LCZ of the stations can be extracted from GEE by using the
      `Dataset.get_LCZ()` method.

    """
    targets = [BuddyCheckStation(station=sta) for sta in target_stations]


    # Validate safety net configs if provided
    buddymethods.validate_safety_net_configs(safety_net_configs)

    # -----  Part 1: construct buddy groups ------
    # compute distance metric
    logger.debug("Calculating distance matrix with Haversine formula")
    dist_matrix = generate_distance_matrix(metadf)

    # find potential buddies by distance
    logger.debug(
        "Finding spatial buddies within radius of %s meters", spatial_buddy_radius
    )
    
    buddymethods.assign_spatial_buddies(
        distance_df=dist_matrix,
        metadf = metadf,
        max_alt_diff=max_alt_diff,
        buddy_radius=spatial_buddy_radius,
        wrappedstations=targets,
    )


    # ---- Part 2: Preparing the records  -----

    # construct a wide observation dataframe
    
    widedf, timestamp_map = buddymethods.create_wide_obs_df(
        wrappedstations=targets,
        obstype=obstype,
        instantaneous_tolerance=instantaneous_tolerance,
    )
    
    # lapse rate correction
    widedf = buddymethods.correct_lapse_rate(widedf=widedf,
                                wrappedstations=targets,
                                lapserate=lapserate)
                       
        

    # ---- Part 3 : Apply buddy check per stationcenter,

    # valid_targets = [budsta for budsta in targets if budsta.has_enough_buddies(
    #     groupname='spatial', min_buddies = spatial_min_sample_size)]
    
    outliersbin = []
    for i in range(N_iter):
        logger.debug("Starting iteration %s of %s", i + 1, N_iter)
        # convert values to NaN, if they are labeled as outlier in
        #  previous iteration
        if bool(outliersbin):
            logger.debug("Converting previous-iteration outliers to NaN")
            for outlier_station, outlier_time in outliersbin:
                if outlier_station in widedf.columns:
                    widedf.loc[outlier_time, outlier_station] = np.nan

        if use_mp:
        
            num_cores = Settings.get('use_N_cores_for_MP')

            logger.debug(f"Running spatial buddy check with multiprocessing on {num_cores} cores")

            # split dataframe along time/index dimension
            chunks = np.array_split(widedf, num_cores)

            # build input arguments for each station and each chunk
            inputargs = [
                {
                    'centerwrapstation': sta,
                    'buddygroupname': 'spatial',
                    'widedf': chunk,
                    'min_sample_size': spatial_min_sample_size,
                    'min_sample_spread': min_sample_spread,
                    'outlier_threshold': spatial_z_threshold,
                    'iteration': i,
                    'check_type': 'spatial_check',
                    'use_z_robust_method': use_z_robust_method,
                }
                for sta in targets
                for chunk in chunks
            ]

            logger.debug(
                f"Submitting {len(inputargs)} multiprocessing tasks "
                f"({len(targets)} stations × {len(chunks)} chunks)"
            )

            # run in parallel
            with ProcessPoolExecutor(max_workers=num_cores) as executor:
                buddy_output = list(
                    executor.map(
                        _run_buddy_test,
                        inputargs
                    )
                )
            
        else:
            # create inputargs for each buddygroup, and for each chunk in time
            inputargs = [
                {
                    'centerwrapstation': sta,
                    'buddygroupname': 'spatial',
                    'widedf': widedf,
                    'min_sample_size': spatial_min_sample_size,
                    'min_sample_spread': min_sample_spread,
                    'outlier_threshold': spatial_z_threshold,
                    'iteration': i,
                    'check_type': 'spatial_check',
                    'use_z_robust_method': use_z_robust_method,
                }
                for sta in targets
            ]
            

            logger.debug("Finding outliers in each buddy group")
            buddy_output = list(map(lambda kwargs: buddymethods.buddy_test_a_station(**kwargs), inputargs))

        #buddy output is [(MultiIndex, BuddyCheckStation), ...], that needs to be unpacked
        outlier_indices, updated_stations = zip(*buddy_output)    
        #overload the Buddycheckstation 
        targets = reconstruct_fractured_targets(list(updated_stations), iteration = i)
            
        # Concatenate all outlier MultiIndices
        # Each element is a MultiIndex with (name, datetime)
        spatial_outliers = buddymethods.concat_multiindices(list(outlier_indices))
        
        # Start with spatial outliers for further processing
        current_outliers_idx = spatial_outliers
        
        # Apply safety nets (if configured)
        if safety_net_configs:
            logger.debug(
                "Applying %s safety net(s) to %s outliers",
                len(safety_net_configs),
                len(current_outliers_idx),
            )
            for safety_net_config in safety_net_configs:
                
                current_outliers_idx = buddymethods.apply_safety_net(    
                    outliers=current_outliers_idx,
                    buddycheckstations = targets,
                    buddygroupname=safety_net_config['category'],
                    metadf = metadf,
                    distance_df = dist_matrix,
                    max_distance=safety_net_config['buddy_radius'],
                    max_alt_diff=max_alt_diff, #make this configurable?
                    wideobsds=widedf,
                    safety_z_threshold=safety_net_config['z_threshold'],
                    min_sample_size=safety_net_config['min_sample_size'],
                    min_sample_spread=min_sample_spread, #make this configurable?
                    use_z_robust_method=use_z_robust_method,
                    iteration=i,
                )
            
            # NOTE: Records saved by any safety net will be tested again in
            # the following iteration. A different result can occur if the
            # spatial/safety net sample changes in the next iteration.

        # Apply whitelist filtering
        current_outliers_idx = buddymethods.save_whitelist_records(
            outliers=current_outliers_idx,
            wrappedstations=targets,
            whiteset=whiteset,
            obstype=obstype,
            iteration=i,
        )
        
        # NOTE: The white-listed records are removed from the outliers at the end
        # of each iteration, similar to the safety nets. They participate in
        # the buddy check calculations but are not flagged as outliers.

        # Convert MultiIndex to list of tuples for outliersbin
        # Format: (station_name, timestamp, message)
        for name, dt in current_outliers_idx:
            outliersbin.append((name, dt))
    
    #Prepare for output
    return_results = {} 
    qcsettings = locals()
    del qcsettings['target_stations']
    del qcsettings['metadf']
    for wrapsta in targets:
        #1. Map timestamps back to original timestamps
        wrapsta.map_timestamps(timestamp_map=timestamp_map[wrapsta.name])
    
        #2. Create final QC labels (specific for buddy check)
        final_labels = wrapsta.get_final_labels()
        
        #3 Convert these flags to default qc flags
        qcflags = final_labels.map(to_qc_labels_map)
        
        #4 Create QCresult object
        outliers = wrapsta.station.get_sensor(obstype).series.loc[qcflags[qcflags == flagged_cond].index]
        
        qcres = QCresult(
                    checkname='buddy_check',
                    checksettings=qcsettings,
                    flags=qcflags,
                    outliers=outliers,
                    detail='',
                    )
        qcres.add_details_by_series(detail_series = wrapsta.get_final_details())
        return_results[wrapsta.name] = qcres            
    
    return return_results
    
