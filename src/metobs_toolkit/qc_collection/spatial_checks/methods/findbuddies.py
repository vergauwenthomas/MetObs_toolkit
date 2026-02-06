from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from ...buddystation import BuddyCheckStation

# ------------------------------------------
#    Callables by buddy check
# ------------------------------------------

def assign_spatial_buddies(
    distance_df: pd.DataFrame,
    metadf: pd.DataFrame,
    buddy_radius: Union[int, float],
    wrappedstations: List[BuddyCheckStation],
    max_alt_diff: Union[int, float, None]=None,
) -> None:


    spatial_buddies = _find_buddies_by_distance(distance_df=distance_df,
                                                buddy_radius=buddy_radius)
    
    #update the wrapstations
    for wrapsta in wrappedstations:
        wrapsta.set_buddies(spatial_buddies[wrapsta.name], groupname='spatial')
    
    if max_alt_diff is not None:
        for wrapsta in wrappedstations:
            filter_buddygroup_by_altitude(
                wrappedstation=wrapsta,
                groupname='spatial',
                altitudes=metadf['altitude'],
                max_altitude_diff=max_alt_diff
            )

       

def filter_buddygroup_by_altitude(
    wrappedstation: BuddyCheckStation,
    groupname: str,
    altitudes: pd.Series,
    max_altitude_diff: Union[int, float]
) :
    
    if altitudes.isnull().any():
        raise ValueError("Altitude series contains NaN values. All stations must have valid altitude data for altitude filtering.")
        
    #update the filter flag
    station_altitude = altitudes.loc[wrappedstation.name]
    alt_buddies = []
    for buddy_name in wrappedstation.get_buddies(groupname=groupname):
        buddy_altitude = altitudes.loc[buddy_name]
        if abs(station_altitude - buddy_altitude) <= max_altitude_diff:
            alt_buddies.append(buddy_name)
    wrappedstation.filter_buddies(groupname=groupname, filteredbuddies=alt_buddies)
    
def subset_buddies_to_nearest(
    wrappedstations: List,
    distance_df: pd.DataFrame,
    max_sample_size: int,
    groupname: str,
) -> None:
    """Subset buddy groups to the nearest N stations.

    For each wrapped station, the buddies in the specified group are sorted
    by distance and only the closest ``max_sample_size`` buddies are kept.

    Parameters
    ----------
    wrappedstations : list of BuddyWrapSensor
        The wrapped stations whose buddy groups should be subsetted.
    distance_df : pd.DataFrame
        Symmetric distance matrix with station names as index and columns.
    max_sample_size : int
        Maximum number of buddies to keep per station.
    groupname : str
        The name of the buddy group to subset (e.g., 'spatial').
    """
    for wrapsta in wrappedstations:
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
    distance_df: pd.DataFrame, buddy_radius: Union[int, float]
) -> Dict:

    buddies = {}
    for refstation, distances in distance_df.iterrows():
        bud_stations = distances[distances <= buddy_radius].index.to_list()
        bud_stations.remove(refstation)
        buddies[refstation] = bud_stations

    return buddies

