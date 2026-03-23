from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor

# ------------------------------------------
#    Callables by buddy check
# ------------------------------------------


def assign_spatial_buddies(
    distance_df: pd.DataFrame,
    metadf: pd.DataFrame,
    buddy_radius: Union[int, float],
    min_buddy_distance: Union[int, float],
    wrappedsensors: List[BuddyWrapSensor],
    max_alt_diff: Union[int, float, None] = None,
) -> None:
    """Assign spatial buddy groups to all wrapped stations.

    Finds buddies for each station within ``buddy_radius`` using
    :func:`_find_buddies_by_distance`, optionally further filters by
    altitude difference, and stores the result on each
    :class:`~buddywrapsensor.BuddyWrapSensor` via :meth:`set_buddies`.

    Parameters
    ----------
    distance_df : pandas.DataFrame
        Symmetric distance matrix (metres) with station names as index and
        columns.
    metadf : pandas.DataFrame
        Station metadata; must contain an ``'altitude'`` column when
        ``max_alt_diff`` is not None.
    buddy_radius : int or float
        Maximum distance (metres) from a station for another station to be
        considered a buddy.
    min_buddy_distance : int or float
        Minimum distance (metres); stations closer than this are excluded.
    wrappedsensors : list of BuddyWrapSensor
        Wrapped sensor objects whose buddy groups are updated in-place.
    max_alt_diff : int, float, or None, optional
        If given, buddies with an altitude difference greater than this
        value (metres) are removed.  Default is None (no altitude filter).
    """
    spatial_buddies = _find_buddies_by_distance(
        distance_df=distance_df,
        buddy_radius=buddy_radius,
        min_buddy_distance=min_buddy_distance,
    )

    # update the wrapsensors
    for wrapsensor in wrappedsensors:
        wrapsensor.set_buddies(spatial_buddies[wrapsensor.name], groupname="spatial")

    if max_alt_diff is not None:
        for wrapsensor in wrappedsensors:
            filter_buddygroup_by_altitude(
                wrappedsensor=wrapsensor,
                groupname="spatial",
                altitudes=metadf["altitude"],
                max_altitude_diff=max_alt_diff,
            )


def filter_buddygroup_by_altitude(
    wrappedsensor: BuddyWrapSensor,
    groupname: str,
    altitudes: pd.Series,
    max_altitude_diff: Union[int, float],
):
    """Remove buddies whose altitude differs too much from the station's altitude.

    The filtered buddy list is stored under ``'{groupname}_filtered'`` via
    :meth:`~buddywrapsensor.BuddyWrapSensor.filter_buddies`.

    Parameters
    ----------
    wrappedsensor : BuddyWrapSensor
        The wrapped sensor whose buddy group is to be filtered.
    groupname : str
        The name of the buddy group to filter (e.g. ``'spatial'``).
    altitudes : pandas.Series
        Series mapping station names to altitude values (metres).  Must
        contain no NaN values.
    max_altitude_diff : int or float
        Maximum allowed absolute altitude difference (metres).

    Raises
    ------
    ValueError
        If ``altitudes`` contains any NaN values.
    """
    if altitudes.isnull().any():
        raise ValueError(
            "Altitude series contains NaN values. All stations must have valid altitude data for altitude filtering."
        )

    # update the filter flag
    station_altitude = altitudes.loc[wrappedsensor.name]
    alt_buddies = []
    for buddy_name in wrappedsensor.get_buddies(groupname=groupname):
        buddy_altitude = altitudes.loc[buddy_name]
        if abs(station_altitude - buddy_altitude) <= max_altitude_diff:
            alt_buddies.append(buddy_name)
    wrappedsensor.filter_buddies(groupname=groupname, filteredbuddies=alt_buddies)


def subset_buddies_to_nearest(
    wrappedsensors: List,
    distance_df: pd.DataFrame,
    max_sample_size: int,
    groupname: str,
) -> None:
    """Subset buddy groups to the nearest N stations.

    For each wrapped station, the buddies in the specified group are sorted
    by distance and only the closest ``max_sample_size`` buddies are kept.

    Parameters
    ----------
    wrappedsensors : list of BuddyWrapSensor
        The wrapped sensors whose buddy groups should be subsetted.
    distance_df : pd.DataFrame
        Symmetric distance matrix with station names as index and columns.
    max_sample_size : int
        Maximum number of buddies to keep per station.
    groupname : str
        The name of the buddy group to subset (e.g., 'spatial').
    """
    for wrapsta in wrappedsensors:
        buddies = wrapsta.get_buddies(groupname=groupname)
        if len(buddies) <= max_sample_size:
            continue
        # Sort buddies by distance and keep only the nearest N
        buddy_distances = distance_df.loc[wrapsta.name, buddies]
        nearest = buddy_distances.nsmallest(max_sample_size).index.to_list()
        wrapsta.filter_buddies(groupname=groupname, filteredbuddies=nearest)


# ------------------------------------------
#    Help functions to find buddies
# ------------------------------------------
def _find_buddies_by_distance(
    distance_df: pd.DataFrame,
    buddy_radius: Union[int, float],
    min_buddy_distance: Union[int, float] = 0,
) -> Dict:
    """Build a mapping from each station to its buddy stations within a distance range.

    Parameters
    ----------
    distance_df : pandas.DataFrame
        Symmetric distance matrix (metres) with station names as index and
        columns.
    buddy_radius : int or float
        Maximum distance (metres) for buddy inclusion.
    min_buddy_distance : int or float, optional
        Minimum distance (metres) for buddy inclusion.  Default is 0.

    Returns
    -------
    dict
        Dictionary ``{station_name: [buddy_names]}`` for every station in
        ``distance_df``.
    """
    buddies = {}
    for refstation, distances in distance_df.iterrows():
        bud_stations = distances[
            (distances <= buddy_radius) & (distances >= min_buddy_distance)
        ].index.to_list()
        if refstation in bud_stations:
            bud_stations.remove(refstation)
        buddies[refstation] = bud_stations

    return buddies
