from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")

from .findbuddies import filter_buddygroup_by_altitude
from .samplechecks import buddy_test_a_station
from ..buddywrapsensor import BC_PASSED

if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor


def validate_safety_net_configs(safety_net_configs: List[Dict]) -> None:
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

    return None




def apply_safety_net(
    outliers: pd.Index,
    buddycheckstations: List[BuddyWrapSensor],
    buddygroupname: str,
    metadf: pd.DataFrame,
    distance_df: pd.DataFrame,
    max_distance: Union[int, float],
    max_alt_diff: Union[int, float, None],
    wideobsds: pd.DataFrame,
    safety_z_threshold: Union[int, float],
    min_sample_size: int,
    min_sample_spread: Union[int, float],
    use_z_robust_method: bool,
    iteration: int,
) -> pd.MultiIndex:
   
    # Track records that were saved (passed the safety net test)
    saved_records = pd.MultiIndex.from_tuples([], names=['name', 'datetime'])
    
    #create a name map of the wrappedstations
    name_map = {wrapsta.name: wrapsta for wrapsta in buddycheckstations}
    
    
    #find the categorical buddies (only for the outlier stations)
    for outlstation in outliers.get_level_values('name').unique():
        wrapsta = name_map[outlstation]
        
        ref_category = metadf.loc[wrapsta.name, buddygroupname]
        # Handle NaN values - they should not match with anything
        if pd.isna(ref_category):
            logger.warning(
                "Station %s has NaN value for category '%s' - no category buddies assigned",
                wrapsta.name,
                buddygroupname,
            )
            # Assign empty buddy list
            wrapsta.set_buddies([], groupname=buddygroupname)
        else:
            #find potential candidates
            buddy_candidates = metadf.loc[
                metadf[buddygroupname] == ref_category
            ].index.to_list()
            
            #remove self from buddy candidates
            buddy_candidates.remove(wrapsta.name)
            
            target_distances = distance_df.loc[wrapsta.name, buddy_candidates]
            #filter by distance
            ref_buddies = target_distances[target_distances <= max_distance].index.to_list()
            
            # Assign the found buddies
            wrapsta.set_buddies(ref_buddies, groupname=buddygroupname)
            
        #filter by altitude difference if needed
        if max_alt_diff is not None:
            filter_buddygroup_by_altitude(
                wrappedstation=wrapsta,
                groupname=buddygroupname,
                altitudes=metadf['altitude'],
                max_altitude_diff=max_alt_diff
            )
    
    #find outliers in the new categorical group
    # The buddy_test_a_station function updates flags/details directly
    # and returns only the outlier MultiIndex (BC_FLAGGED records)
    # We need to track BC_PASSED records to remove them from outliers
    for outlstation in outliers.get_level_values('name').unique():
        wrapsta = name_map[outlstation]
        
        # Get the timestamps for this station from the original outliers
        station_outlier_timestamps = outliers[
            outliers.get_level_values('name') == outlstation
        ].get_level_values('datetime')
        
        # Subset wideobsds to only the outlier timestamps for this station
        widedf_subset = wideobsds.loc[station_outlier_timestamps]
        
        # Run the buddy test - this updates flags/details directly on wrapsta
        # and returns outliers (BC_FLAGGED records)
        station_flagged = buddy_test_a_station(
            centerwrapstation=wrapsta,
            buddygroupname=buddygroupname,
            widedf=widedf_subset,
            min_sample_size=min_sample_size,
            min_sample_spread=min_sample_spread,
            outlier_threshold=safety_z_threshold,
            iteration=iteration,
            check_type=f'safetynet_check:{buddygroupname}',
            use_z_robust_method=use_z_robust_method,
        )
        
        # Get passed timestamps from the flags DataFrame
        # These are records where the safetynet check passed (BC_PASSED)
        check_col = f'safetynet_check:{buddygroupname}'
        if not wrapsta.flags.empty and check_col in wrapsta.flags.columns:
            # Get flags for this iteration
            iter_mask = wrapsta.flags.index.get_level_values('iteration') == iteration
            iter_flags = wrapsta.flags.loc[iter_mask, check_col]
            
            # Find passed timestamps
            passed_mask = iter_flags == BC_PASSED
            if passed_mask.any():
                passed_timestamps = iter_flags[passed_mask].index.get_level_values('datetime')
                
                # Create MultiIndex for saved records
                station_saved = pd.MultiIndex.from_arrays(
                    [[outlstation] * len(passed_timestamps), passed_timestamps],
                    names=['name', 'datetime']
                )
                saved_records = saved_records.union(station_saved)
    
    # Return original outliers minus the saved records
    remaining_outliers = outliers.difference(saved_records)
    
    return remaining_outliers.sort_values().unique()
    
   