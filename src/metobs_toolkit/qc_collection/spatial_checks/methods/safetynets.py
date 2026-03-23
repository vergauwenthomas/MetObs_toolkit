from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")

from .findbuddies import filter_buddygroup_by_altitude, subset_buddies_to_nearest
from .samplechecks import buddy_test_a_station
from ..buddywrapsensor import BC_PASSED, BC_NO_BUDDIES

if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor


def validate_safety_net_configs(safety_net_configs: List[Dict]) -> None:
    """
    Validate that all required keys are present in safety_net_configs.

    This function checks that each safety net configuration contains all required
    keys and validates optional keys when present. Required keys are: 'category',
    'buddy_radius', 'z_threshold', and 'min_sample_size'. Optional keys include
    'max_sample_size' (must be larger than 'min_sample_size') and
    'only_if_previous_had_no_buddies' (boolean, cannot be True for first safety net).

    Parameters
    ----------
    safety_net_configs : list of dict or None
        List of safety net configuration dictionaries. If None, the function
        returns without validation.

    Raises
    ------
    ValueError
        If safety_net_configs is not a list or contains non-dict elements, or
        if validation fails for max_sample_size (must be > min_sample_size),
        or if 'only_if_previous_had_no_buddies' is not a boolean or is True
        for the first safety net.
    KeyError
        If any required key is missing from a safety net configuration.
    """
    if safety_net_configs is None:
        return None

    required_keys = {"category", "buddy_radius", "z_threshold", "min_sample_size"}
    optional_keys = {"max_sample_size", "only_if_previous_had_no_buddies"}

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

        # Validate optional max_sample_size
        if "max_sample_size" in config:
            max_ss = config["max_sample_size"]
            if max_ss is not None:
                min_ss = config["min_sample_size"]
                if max_ss <= min_ss:
                    raise ValueError(
                        f"Safety net config at index {i}: 'max_sample_size' "
                        f"({max_ss}) must be larger than 'min_sample_size' "
                        f"({min_ss})."
                    )

        # Validate optional only_if_previous_had_no_buddies
        if "only_if_previous_had_no_buddies" in config:
            val = config["only_if_previous_had_no_buddies"]
            if not isinstance(val, bool):
                raise ValueError(
                    f"Safety net config at index {i}: "
                    f"'only_if_previous_had_no_buddies' must be a bool, "
                    f"got {type(val).__name__}."
                )
            if val and i == 0:
                raise ValueError(
                    f"Safety net config at index {i}: "
                    f"'only_if_previous_had_no_buddies' cannot be True for "
                    f"the first safety net because there is no previous "
                    f"safety net to fall back from."
                )

    return None


def assign_safety_net_buddies(
    wrapsta: BuddyWrapSensor,
    metadf: pd.DataFrame,
    distance_df: pd.DataFrame,
    buddygroupname: str,
    max_distance: Union[int, float],
    min_distance: Union[int, float],
    max_alt_diff: Union[int, float, None],
    max_sample_size: Union[int, None],
) -> None:
    """
    Assign category buddies to a wrapped station for safety net buddy checks.

    This function identifies buddies that share the same categorical attribute
    (e.g., LCZ, network) with the reference station, within specified distance
    constraints and altitude difference limits.

    The assigned buddy group is stored on the wrapped station and can be accessed
    using ``wrapsta.get_buddies(groupname=buddygroupname)``.

    Parameters
    ----------
    wrapsta : BuddyWrapSensor
        The wrapped station for which buddies should be assigned.
    metadf : pd.DataFrame
        DataFrame containing station metadata. Must have the category column
        specified by ``buddygroupname`` (e.g., 'LCZ', 'network') and an
        'altitude' column if ``max_alt_diff`` is not None.
    distance_df : pd.DataFrame
        Symmetric distance matrix with station names as index and columns.
        Distances should be in meters.
    buddygroupname : str
        The name of the metadata column to group by (e.g., 'LCZ', 'network').
        This is also used as the buddy group identifier on the wrapped station.
    max_distance : int or float
        Maximum distance (in meters) for category buddies. Stations farther
        than this distance will be excluded.
    min_distance : int or float
        Minimum distance (in meters) required between the station and its
        category buddies. Stations closer than this distance will be excluded.
    max_alt_diff : int, float, or None
        Maximum altitude difference (in meters) allowed for buddies. If None,
        no altitude filtering is applied.
    max_sample_size : int or None
        Maximum number of category buddies to keep per station. If not None,
        buddies are sorted by distance and only the nearest ``max_sample_size``
        are retained. If None, no limit is applied.

    Returns
    -------
    None
        The function modifies ``wrapsta`` in place by assigning the buddy group.

    Notes
    -----
    The function applies filters in the following order:

    1. Category matching: Only stations with the same category value as the
       reference station are considered.
    2. Distance filtering: Only stations within [min_distance, max_distance]
       are kept.
    3. Altitude filtering (if max_alt_diff is not None): Only stations with
       altitude difference <= max_alt_diff are kept.
    4. Sample size limiting (if max_sample_size is not None): Only the nearest
       max_sample_size buddies are kept.

    If the reference station has a NaN value for the category, an empty buddy
    list is assigned and a warning is logged.

    """

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
        # find potential candidates
        buddy_candidates = metadf.loc[
            metadf[buddygroupname] == ref_category
        ].index.to_list()

        # remove self from buddy candidates
        buddy_candidates.remove(wrapsta.name)

        target_distances = distance_df.loc[wrapsta.name, buddy_candidates]
        # filter by distance
        ref_buddies = target_distances[
            (target_distances <= max_distance) & (target_distances >= min_distance)
        ].index.to_list()

        # Assign the found buddies
        wrapsta.set_buddies(ref_buddies, groupname=buddygroupname)

    # filter by altitude difference if needed
    if max_alt_diff is not None:
        filter_buddygroup_by_altitude(
            wrappedsensor=wrapsta,
            groupname=buddygroupname,
            altitudes=metadf["altitude"],
            max_altitude_diff=max_alt_diff,
        )

    # Subset category buddies to nearest N if max_sample_size is set
    if max_sample_size is not None:
        subset_buddies_to_nearest(
            wrappedsensors=[wrapsta],
            distance_df=distance_df,
            max_sample_size=max_sample_size,
            groupname=buddygroupname,
        )


def apply_safety_net(
    outliers: pd.Index,
    buddychecksensors: List[BuddyWrapSensor],
    buddygroupname: str,
    metadf: pd.DataFrame,
    distance_df: pd.DataFrame,
    max_distance: Union[int, float],
    min_distance: Union[int, float],
    max_alt_diff: Union[int, float, None],
    wideobsds: pd.DataFrame,
    safety_z_threshold: Union[int, float],
    min_sample_size: int,
    min_sample_spread: Union[int, float],
    use_z_robust_method: bool,
    iteration: int,
    max_sample_size: Union[int, None] = None,
    only_if_previous_had_no_buddies: bool = False,
    previous_safetynet_category: Union[str, None] = None,
) -> pd.MultiIndex:
    """Run a single safety-net buddy check on the flagged outlier records.

    For each outlier record, safety-net buddies are assigned (using a
    possibly different radius/altitude filter), then the z-score buddy
    test is applied. Records that pass the safety-net test are removed
    from the outliers, and the remaining outliers are returned. Flags
    and details are updated on the BuddyWrapSensor objects in-place.

    Parameters
    ----------
    outliers : pandas.Index
        MultiIndex ``('name', 'datetime')`` of records flagged as outliers
        by the primary spatial check.
    buddychecksensors : list of BuddyWrapSensor
        All wrapped sensors involved in the buddy check.
    buddygroupname : str
        Name prefix for the safety-net buddy group (e.g. a LCZ category).
    metadf : pandas.DataFrame
        Station metadata with at least an ``'altitude'`` column and a column
        matching ``buddygroupname`` for category grouping.
    distance_df : pandas.DataFrame
        Symmetric distance matrix (metres) with station names as index
        and columns.
    max_distance : int or float
        Maximum buddy search radius (metres) for the safety net.
    min_distance : int or float
        Minimum buddy distance (metres).
    max_alt_diff : int, float, or None
        Altitude filter (metres).  None disables altitude filtering.
    wideobsds : pandas.DataFrame
        Wide-format observations DataFrame used for z-score computation.
    safety_z_threshold : int or float
        Z-score threshold for the safety-net outlier test.
    min_sample_size : int
        Minimum number of valid buddy samples required.
    min_sample_spread : int or float
        Minimum sample spread (std or MAD) to avoid near-zero division.
    use_z_robust_method : bool
        If True, use the robust (median/MAD) z-score method.
    iteration : int
        Current iteration number.
    max_sample_size : int or None, optional
        Maximum number of buddies to use.  None means no cap.
        Default is None.
    only_if_previous_had_no_buddies : bool, optional
        If True, only apply the safety net for records that had
        insufficient buddies in the previous safety-net layer.
        Default is False.
    previous_safetynet_category : str or None, optional
        Category name of the preceding safety net, required when
        ``only_if_previous_had_no_buddies`` is True.  Default is None.

    Returns
    -------
    pandas.MultiIndex
        MultiIndex ``('name', 'datetime')`` of remaining outliers after
        applying the safety net. This is the input ``outliers`` minus any
        records that passed the safety-net test (i.e., were saved from
        being flagged).
    """
    # Track records that were saved (passed the safety net test)
    saved_records = pd.MultiIndex.from_tuples([], names=["name", "datetime"])

    # create a name map of the wrappedsensors
    name_map = {wrapsens.name: wrapsens for wrapsens in buddychecksensors}

    # If only_if_previous_had_no_buddies is True, restrict outliers to only
    # those records where the previous safety net had insufficient buddies
    # (BC_NO_BUDDIES flag). This is determined by inspecting the flags
    # already stored on each BuddyWrapSensor for the current iteration.
    if only_if_previous_had_no_buddies:
        if previous_safetynet_category is None:
            raise ValueError(
                "only_if_previous_had_no_buddies is True but "
                "previous_safetynet_category is None. This should not "
                "happen -- the first safety net cannot use this option."
            )

        prev_check_col = f"safetynet_check:{previous_safetynet_category}"
        previous_no_buddies = pd.MultiIndex.from_tuples([], names=["name", "datetime"])

        for station_name in outliers.get_level_values("name").unique():
            wrapsta = name_map[station_name]
            if not wrapsta.flags.empty and prev_check_col in wrapsta.flags.columns:
                iter_mask = (
                    wrapsta.flags.index.get_level_values("iteration") == iteration
                )
                iter_flags = wrapsta.flags.loc[iter_mask, prev_check_col]
                nb_mask = iter_flags == BC_NO_BUDDIES
                if nb_mask.any():
                    nb_timestamps = iter_flags[nb_mask].index.get_level_values(
                        "datetime"
                    )
                    station_nb = pd.MultiIndex.from_arrays(
                        [
                            [station_name] * len(nb_timestamps),
                            nb_timestamps,
                        ],
                        names=["name", "datetime"],
                    )
                    previous_no_buddies = previous_no_buddies.union(station_nb)

        if previous_no_buddies.empty:
            logger.info(
                "only_if_previous_had_no_buddies is True but no records "
                "from the previous safety net ('%s') had insufficient "
                "buddies. Skipping safety net '%s' entirely.",
                previous_safetynet_category,
                buddygroupname,
            )
            return outliers

        outliers = outliers.intersection(previous_no_buddies)
        logger.info(
            "Filtering to %s outlier records that had insufficient "
            "buddies in the previous safety net ('%s').",
            len(outliers),
            previous_safetynet_category,
        )

    if outliers.empty:
        return outliers

    # find the categorical buddies (only for the outlier stations)
    for outlstation in outliers.get_level_values("name").unique():
        wrapsta = name_map[outlstation]
        assign_safety_net_buddies(
            wrapsta=wrapsta,
            metadf=metadf,
            distance_df=distance_df,
            buddygroupname=buddygroupname,
            max_distance=max_distance,
            min_distance=min_distance,
            max_alt_diff=max_alt_diff,
            max_sample_size=max_sample_size,
        )

    # find outliers in the new categorical group
    # The buddy_test_a_station function updates flags/details directly
    # and returns only the outlier MultiIndex (BC_FLAGGED records)
    # We need to track BC_PASSED records to remove them from outliers
    for outlstation in outliers.get_level_values("name").unique():
        wrapsta = name_map[outlstation]

        # Get the timestamps for this station from the original outliers
        station_outlier_timestamps = outliers[
            outliers.get_level_values("name") == outlstation
        ].get_level_values("datetime")

        # Subset wideobsds to only the outlier timestamps for this station
        widedf_subset = wideobsds.loc[station_outlier_timestamps]

        # Run the buddy test - this updates flags/details directly on wrapsta
        # and returns outliers (BC_FLAGGED records)
        station_flagged = buddy_test_a_station(
            centerwrapsensor=wrapsta,
            buddygroupname=buddygroupname,
            widedf=widedf_subset,
            min_sample_size=min_sample_size,
            min_sample_spread=min_sample_spread,
            outlier_threshold=safety_z_threshold,
            iteration=iteration,
            check_type=f"safetynet_check:{buddygroupname}",
            use_z_robust_method=use_z_robust_method,
        )

        # Get passed timestamps from the flags DataFrame
        # These are records where the safetynet check passed (BC_PASSED)
        check_col = f"safetynet_check:{buddygroupname}"
        if not wrapsta.flags.empty and check_col in wrapsta.flags.columns:
            # Get flags for this iteration
            iter_mask = wrapsta.flags.index.get_level_values("iteration") == iteration
            iter_flags = wrapsta.flags.loc[iter_mask, check_col]

            # Find passed timestamps
            passed_mask = iter_flags == BC_PASSED
            if passed_mask.any():
                passed_timestamps = iter_flags[passed_mask].index.get_level_values(
                    "datetime"
                )

                # Create MultiIndex for saved records
                station_saved = pd.MultiIndex.from_arrays(
                    [[outlstation] * len(passed_timestamps), passed_timestamps],
                    names=["name", "datetime"],
                )
                saved_records = saved_records.union(station_saved)

    # Return original outliers minus the saved records
    remaining_outliers = outliers.difference(saved_records)

    return remaining_outliers.sort_values().unique()
