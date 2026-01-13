from __future__ import annotations

import logging
from typing import Union, List, Dict, TYPE_CHECKING, Tuple

import pandas as pd
from metobs_toolkit.backend_collection.datetime_collection import to_timedelta

logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from ...buddystation import BuddyCheckStation


def create_wide_obs_df(wrappedstations: List[BuddyCheckStation],
                       obstype: str,
                       instantaneous_tolerance: pd.Timedelta
                       ) -> Tuple[pd.DataFrame, Dict]:
    logger.debug("Constructing wide observation DataFrame for obstype: %s", obstype)
    concatlist = []
    for wrapsta in wrappedstations:
        if obstype in wrapsta.station.sensordata.keys():
            records = wrapsta.station.get_sensor(obstype).series
            records.name = wrapsta.name
            concatlist.append(records)
            

    # synchronize the timestamps
    logger.debug("Synchronizing timestamps")
    combdf, timestamp_map = _synchronize_series(
        series_list=concatlist, max_shift=instantaneous_tolerance
    )
    
    return (combdf, timestamp_map)

def _synchronize_series(
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



def concat_multiindices(
    indices: List[pd.MultiIndex]
) -> pd.MultiIndex:
    """Concatenate a list of MultiIndex objects into a single MultiIndex.
    
    Parameters
    ----------
    indices : list of pd.MultiIndex
        List of MultiIndex objects to concatenate.
        
    Returns
    -------
    pd.MultiIndex
        Concatenated MultiIndex.
    """
    if not indices:
        return pd.MultiIndex.from_tuples([], names=['name', 'datetime'])
    
    concatenated = pd.MultiIndex.from_tuples(
        [tup for idx in indices for tup in idx],
        names=indices[0].names
    )
    
    
    # non_empty_indices = [idx for idx in outlier_indices if len(idx) > 0]
    #     if non_empty_indices:
    #         spatial_outliers = non_empty_indices[0]
    #         for idx in non_empty_indices[1:]:
    #             spatial_outliers = spatial_outliers.union(idx)
    #         spatial_outliers = spatial_outliers.drop_duplicates()
    #     else:
    #         spatial_outliers = pd.MultiIndex.from_tuples([], names=['name', 'datetime'])
    
    return concatenated 