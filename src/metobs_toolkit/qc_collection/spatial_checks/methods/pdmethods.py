from __future__ import annotations

import logging
from typing import List, Dict, TYPE_CHECKING, Tuple

import pandas as pd
from metobs_toolkit.backend_collection.datetime_collection import to_timedelta

logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor


def create_wide_obs_df(wrappedsensors: List[BuddyWrapSensor],
                       instantaneous_tolerance: pd.Timedelta
                       ) -> Tuple[pd.DataFrame, Dict]:
    """Build a wide-format observations DataFrame from wrapped sensors.

    Sensor time series are synchronised to a common regular datetime index
    using :func:`_synchronize_series` before being combined column-wise.

    Parameters
    ----------
    wrappedsensors : list of BuddyWrapSensor
        Wrapped sensors to include.  The station name is used as the column
        label.
    instantaneous_tolerance : pandas.Timedelta
        Maximum time shift allowed when merging a sensor's timestamps onto
        the common target index.

    Returns
    -------
    pandas.DataFrame
        Wide DataFrame with one column per station and a synchronised
        DatetimeIndex.
    dict
        Timestamp mapping returned by :func:`_synchronize_series`; maps
        each synchronised timestamp to the original timestamp for each
        station.
    """
    concatlist = []
    for wrapsens in wrappedsensors:
        records = wrapsens.sensor.series
        records.name = wrapsens.name
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
    
    
    return concatenated 