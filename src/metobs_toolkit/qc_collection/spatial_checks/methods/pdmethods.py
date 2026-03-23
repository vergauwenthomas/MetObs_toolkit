from __future__ import annotations

import logging
from typing import List, Dict, TYPE_CHECKING, Tuple

import pandas as pd
from metobs_toolkit.backend_collection.datetime_collection import to_timedelta

logger = logging.getLogger("<metobs_toolkit>")


if TYPE_CHECKING:
    from ..buddywrapsensor import BuddyWrapSensor


def create_wide_obs_df(
    wrappedsensors: List[BuddyWrapSensor], instantaneous_tolerance: pd.Timedelta
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
    target_dt = target_dt.rename("datetime")

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

        # Ensure each original value is used at most once.
        # When multiple target timestamps map to the same original, keep
        # only the closest match and discard the rest.
        orig_col = targetdf["orig_datetime"].dropna()
        for orig_ts, group in orig_col.groupby(orig_col):
            if len(group) > 1:
                diffs = (group.index - orig_ts).to_series(index=group.index).abs()
                keep = diffs.idxmin()
                drop = group.index.drop(keep)
                targetdf.loc[drop, s.name] = float("nan")
                targetdf.loc[drop, "orig_datetime"] = pd.NaT

        # extract the mapping (new -> original)
        orig_timestampseries = targetdf["orig_datetime"]
        orig_timestampseries.name = "original_timestamp"
        timestamp_mapping[s.name] = orig_timestampseries

        synchronized_series.append(targetdf[s.name])

    return pd.concat(synchronized_series, axis=1), timestamp_mapping


def concat_multiindices(indices: List[pd.MultiIndex]) -> pd.MultiIndex:
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
        return pd.MultiIndex.from_tuples([], names=["name", "datetime"])

    concatenated = pd.MultiIndex.from_tuples(
        [tup for idx in indices for tup in idx], names=indices[0].names
    )

    return concatenated


if __name__ == "__main__":
    print("Running _synchronize_series sanity checks...")

    # ------------------------------------------------------------------ #
    # Test 1: Same frequency, aligned timestamps
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [10.0, 11.0, 12.0],
        index=pd.date_range("2024-01-01 00:00", periods=3, freq="5min"),
        name="A",
    )
    s2 = pd.Series(
        [20.0, 21.0, 22.0],
        index=pd.date_range("2024-01-01 00:00", periods=3, freq="5min"),
        name="B",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("2min"))

    assert list(df.columns) == ["A", "B"]
    assert df.index.freq is not None
    assert len(df) == 3
    assert df["A"].tolist() == [10.0, 11.0, 12.0]
    assert df["B"].tolist() == [20.0, 21.0, 22.0]
    assert set(tmap.keys()) == {"A", "B"}
    assert tmap["A"].index.equals(df.index)

    # ------------------------------------------------------------------ #
    # Test 2: Same freq, overlapping but shifted timestamps
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [10.0, 11.0, 12.0, 13.0],
        index=pd.date_range("2024-01-01 00:00", periods=4, freq="5min"),
        name="A",
    )
    s2 = pd.Series(
        [20.0, 21.0, 22.0, 23.0],
        index=pd.date_range("2024-01-01 00:10", periods=4, freq="5min"),
        name="B",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("2min"))

    assert df.index.freq is not None
    expected_idx = pd.date_range("2024-01-01 00:00", "2024-01-01 00:25", freq="5min")
    assert df.index.equals(expected_idx)
    assert df.loc["2024-01-01 00:00", "A"] == 10.0
    assert pd.isna(df.loc["2024-01-01 00:20", "A"])
    assert pd.isna(df.loc["2024-01-01 00:00", "B"])
    assert df.loc["2024-01-01 00:10", "B"] == 20.0

    # ------------------------------------------------------------------ #
    # Test 3: Different frequencies (5min vs 10min)
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [10.0, 11.0, 12.0, 13.0, 14.0],
        index=pd.date_range("2024-01-01 00:00", periods=5, freq="5min"),
        name="high",
    )
    s2 = pd.Series(
        [20.0, 21.0, 22.0],
        index=pd.date_range("2024-01-01 00:00", periods=3, freq="10min"),
        name="low",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("1min"))

    assert df.index.freq == pd.tseries.frequencies.to_offset("5min")
    assert df.loc["2024-01-01 00:00", "low"] == 20.0
    assert df.loc["2024-01-01 00:10", "low"] == 21.0
    assert df.loc["2024-01-01 00:20", "low"] == 22.0
    assert pd.isna(df.loc["2024-01-01 00:05", "low"])
    assert pd.isna(df.loc["2024-01-01 00:15", "low"])

    # ------------------------------------------------------------------ #
    # Test 4: Slight time shifts within tolerance
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [10.0, 11.0, 12.0],
        index=pd.to_datetime([
            "2024-01-01 00:00:00",
            "2024-01-01 00:05:00",
            "2024-01-01 00:10:00",
        ]),
        name="exact",
    )
    s2 = pd.Series(
        [20.0, 21.0, 22.0],
        index=pd.to_datetime([
            "2024-01-01 00:00:30",
            "2024-01-01 00:05:30",
            "2024-01-01 00:10:30",
        ]),
        name="shifted",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("1min"))

    assert df["shifted"].notna().sum() == 3
    assert df.loc["2024-01-01 00:00", "shifted"] == 20.0
    orig = tmap["shifted"].loc[pd.Timestamp("2024-01-01 00:00")]
    assert orig == pd.Timestamp("2024-01-01 00:00:30")

    # ------------------------------------------------------------------ #
    # Test 5: Shifts exceed tolerance → NaN
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [10.0, 11.0, 12.0],
        index=pd.date_range("2024-01-01 00:00", periods=3, freq="5min"),
        name="base",
    )
    s2 = pd.Series(
        [20.0, 21.0, 22.0],
        index=pd.to_datetime([
            "2024-01-01 00:02:00",
            "2024-01-01 00:07:00",
            "2024-01-01 00:12:00",
        ]),
        name="far",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("30s"))

    assert df["base"].notna().sum() == 3
    assert df["far"].notna().sum() == 0
    assert tmap["far"].isna().all()

    # ------------------------------------------------------------------ #
    # Test 6: Each original value used at most once (deduplication)
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
        index=pd.date_range("2024-01-01 00:00", periods=6, freq="5min"),
        name="ref",
    )
    s2 = pd.Series(
        [40.0, 41.0, 42.0],
        index=pd.to_datetime([
            "2024-01-01 00:02:00",
            "2024-01-01 00:12:00",
            "2024-01-01 00:22:00",
        ]),
        name="shifted10",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("3min"))

    assert df.loc["2024-01-01 00:00", "shifted10"] == 40.0
    assert pd.isna(df.loc["2024-01-01 00:05", "shifted10"])
    assert df.loc["2024-01-01 00:10", "shifted10"] == 41.0
    assert pd.isna(df.loc["2024-01-01 00:15", "shifted10"])
    for val in [40.0, 41.0, 42.0]:
        assert (df["shifted10"] == val).sum() == 1

    # ------------------------------------------------------------------ #
    # Test 7: Single series
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [10.0, 11.0, 12.0],
        index=pd.date_range("2024-01-01", periods=3, freq="1h"),
        name="solo",
    )
    df, tmap = _synchronize_series([s1], max_shift=pd.Timedelta("5min"))

    assert df.shape == (3, 1)
    assert df.columns.tolist() == ["solo"]
    assert df.index.freq is not None

    # ------------------------------------------------------------------ #
    # Test 8: Returned index is a perfect DatetimeIndex with freq
    # ------------------------------------------------------------------ #
    s1 = pd.Series(
        [1.0, 2.0, 3.0],
        index=pd.date_range("2024-01-01 00:00", periods=3, freq="10min"),
        name="X",
    )
    s2 = pd.Series(
        [4.0, 5.0, 6.0],
        index=pd.date_range("2024-01-01 00:00", periods=3, freq="10min"),
        name="Y",
    )
    df, tmap = _synchronize_series([s1, s2], max_shift=pd.Timedelta("1min"))

    assert isinstance(df.index, pd.DatetimeIndex)
    assert df.index.freq is not None
    assert df.index[0] == pd.Timestamp("2024-01-01 00:00")
    assert df.index[-1] == pd.Timestamp("2024-01-01 00:20")

    print("✓ All sanity checks passed")

