from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")

from .findbuddies import filter_buddygroup_by_altitude
from .samplechecks import buddy_test_a_station
if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor
    from ...whitelist import WhiteSet


def save_whitelist_records(
    outliers: pd.MultiIndex,
    wrappedsensors: List[BuddyWrapSensor],
    whiteset: WhiteSet,
    obstype: str,
    iteration: int,
) -> pd.MultiIndex:
    """Apply whitelist filtering to outliers and update station details.
    
    This function filters outlier records against a whitelist. Records that
    match the whitelist are marked as 'saved' and removed from the outlier set.
    Records that don't match the whitelist remain as outliers.
    
    Parameters
    ----------
    outliers : pd.MultiIndex
        MultiIndex with levels ('name', 'datetime') containing the current 
        outlier records to filter.
    wrappedsensors : list of BuddyWrapSensor
        List of wrapped sensor objects to update with whitelist details.
    whiteset : WhiteSet
        A WhiteSet instance containing records that should be excluded from
        outlier detection.
    obstype : str
        The observation type being checked.
    iteration : int
        The current iteration number.
        
    Returns
    -------
    pd.MultiIndex
        MultiIndex with levels ('name', 'datetime') containing only the 
        outliers that were NOT saved by the whitelist.
    """
    if outliers.empty:
        logger.debug("No outliers to filter with whitelist")
        return outliers
    
    if whiteset._flag_is_empty():
        logger.debug("Whitelist is empty, no records saved")
        return outliers
    
    # Track which records are not saved
    # saved_records = pd.MultiIndex.from_tuples([], names=['name', 'datetime'])
    remaining_outliers = pd.MultiIndex.from_tuples([], names=['name', 'datetime'])
    
    # Process each station's outliers
    for wrapsta in wrappedsensors:
        if wrapsta.name not in outliers.get_level_values("name").unique():
            continue  # No outliers for this station
    
        else:
    
            # Get the outlier datetimes for this station
            sta_outlier_dts = pd.DatetimeIndex(
                outliers[outliers.get_level_values("name") == wrapsta.name].get_level_values("datetime"),
                name="datetime"
            )
            
            # Create a SensorWhiteSet for this station
            sensorwhiteset = whiteset.create_sensorwhitelist(
                stationname=wrapsta.name, obstype=obstype
            )
            
            # Filter to get remaining outliers (not whitelisted)
            remaining_dts = sensorwhiteset.catch_white_records(outliers_idx=sta_outlier_dts)
            
            # Saved records are those that were filtered out
            saved_dts = sta_outlier_dts.difference(remaining_dts)
            
            # Build MultiIndex for saved and remaining
            # if len(saved_dts) > 0:
            #     saved_idx = pd.MultiIndex.from_arrays(
            #         [[wrapsta.name] * len(saved_dts), saved_dts],
            #         names=['name', 'datetime']
            #     )
                # saved_records = saved_records.union(saved_idx)
                
                # Update the wrapped station with saved details
                
            # Create detail messages for saved records
            detail_series = pd.Series(
                [f"Saved by whitelist at iteration {iteration}" for _ in saved_dts],
                index=saved_dts
            )
            wrapsta.update_whitelist_details(
                whitelistseries=detail_series,
                iteration=iteration,
                is_saved=True
            )
            
            # Create detail messages for records not saved by whitelist
            detail_series = pd.Series(
                [f"Not saved by whitelist at iteration {iteration}" for _ in remaining_dts],
                index=remaining_dts
            )
            wrapsta.update_whitelist_details(
                whitelistseries=detail_series,
                iteration=iteration,
                is_saved=False
            )

            if len(remaining_dts) > 0:
                remaining_idx = pd.MultiIndex.from_arrays(
                    [[wrapsta.name] * len(remaining_dts), remaining_dts],
                    names=['name', 'datetime']
                )
                remaining_outliers = remaining_outliers.union(remaining_idx)
               
                    

    
    return remaining_outliers.sort_values()

