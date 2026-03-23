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
    import numpy as np

    passed = 0
    failed = 0

    def check(condition, msg):
        global passed, failed
        if condition:
            passed += 1
            print(f"  PASS: {msg}")
        else:
            failed += 1
            print(f"  FAIL: {msg}")

    # ------------------------------------------------------------------ #
    # Test 1: Same frequency, aligned timestamps
    # ------------------------------------------------------------------ #
    print("\nTest 1: Same frequency, aligned timestamps")
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

    check(list(df.columns) == ["A", "B"], "columns are A, B")
    check(df.index.freq is not None, "result has a perfect frequency")
    check(len(df) == 3, "3 rows")
    check(df["A"].tolist() == [10.0, 11.0, 12.0], "A values preserved")
    check(df["B"].tolist() == [20.0, 21.0, 22.0], "B values preserved")
    check(set(tmap.keys()) == {"A", "B"}, "timestamp map has both keys")
    check(
        tmap["A"].index.equals(df.index),
        "timestamp map index matches df index",
    )

    # ------------------------------------------------------------------ #
    # Test 2: Same freq, overlapping but shifted timestamps
    # ------------------------------------------------------------------ #
    print("\nTest 2: Same freq, partially overlapping windows")
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

    check(df.index.freq is not None, "result has a perfect frequency")
    expected_idx = pd.date_range("2024-01-01 00:00", "2024-01-01 00:25", freq="5min")
    check(df.index.equals(expected_idx), "index spans full range")
    # A has data for 00:00–00:15, B has data for 00:10–00:25
    check(df.loc["2024-01-01 00:00", "A"] == 10.0, "A at start")
    check(pd.isna(df.loc["2024-01-01 00:20", "A"]), "A has NaN beyond its range")
    check(pd.isna(df.loc["2024-01-01 00:00", "B"]), "B has NaN before its range")
    check(df.loc["2024-01-01 00:10", "B"] == 20.0, "B at its start")

    # ------------------------------------------------------------------ #
    # Test 3: Different frequencies (5min vs 10min)
    # ------------------------------------------------------------------ #
    print("\nTest 3: Different frequencies – target freq = highest (5min)")
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

    check(df.index.freq == pd.tseries.frequencies.to_offset("5min"),
          "target frequency is 5min")
    # Low-freq values at 00:00, 00:10, 00:20 → only 00:00, 00:10, 00:20
    # map exactly; 00:05, 00:15 are >1min away → NaN
    check(df.loc["2024-01-01 00:00", "low"] == 20.0, "low-freq at 00:00")
    check(df.loc["2024-01-01 00:10", "low"] == 21.0, "low-freq at 00:10")
    check(df.loc["2024-01-01 00:20", "low"] == 22.0, "low-freq at 00:20")
    check(pd.isna(df.loc["2024-01-01 00:05", "low"]),
          "low-freq at 00:05 is NaN (no close match)")
    check(pd.isna(df.loc["2024-01-01 00:15", "low"]),
          "low-freq at 00:15 is NaN (no close match)")

    # ------------------------------------------------------------------ #
    # Test 4: Slight time shifts within tolerance
    # ------------------------------------------------------------------ #
    print("\nTest 4: Slight time shifts within tolerance")
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

    check(df["shifted"].notna().sum() == 3,
          "all shifted values matched within tolerance")
    check(df.loc["2024-01-01 00:00", "shifted"] == 20.0, "shifted at 00:00")
    # Verify timestamp mapping records the original timestamp
    orig = tmap["shifted"].loc[pd.Timestamp("2024-01-01 00:00")]
    check(orig == pd.Timestamp("2024-01-01 00:00:30"),
          "timestamp map records original shifted time")

    # ------------------------------------------------------------------ #
    # Test 5: Shifts exceed tolerance → NaN
    # ------------------------------------------------------------------ #
    print("\nTest 5: Shifts exceed tolerance → NaN")
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

    check(df["base"].notna().sum() == 3, "base values all present")
    check(df["far"].notna().sum() == 0,
          "far values all NaN (2min shift > 30s tolerance)")
    check(tmap["far"].isna().all(),
          "timestamp map all NaT for unmatched series")

    # ------------------------------------------------------------------ #
    # Test 6: Each original value used at most once
    # ------------------------------------------------------------------ #
    print("\nTest 6: Each original value used at most once (deduplication)")
    # s1 = 5min freq, s2 = 10min freq but shifted by 2min.
    # Target freq = 5min. With max_shift=3min, each s2 value is within
    # tolerance of two target timestamps – dedup must keep only the closest.
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

    # Each s2 value (at +2min) is nearest to both a 00/10/20 target (2min)
    # and a 05/15/25 target (3min). Dedup should keep 00:00, 00:10, 00:20.
    check(df.loc["2024-01-01 00:00", "shifted10"] == 40.0,
          "40.0 mapped to nearest target 00:00 (2min away)")
    check(pd.isna(df.loc["2024-01-01 00:05", "shifted10"]),
          "00:05 is NaN (40.0 already used at 00:00)")
    check(df.loc["2024-01-01 00:10", "shifted10"] == 41.0,
          "41.0 mapped to nearest target 00:10 (2min away)")
    check(pd.isna(df.loc["2024-01-01 00:15", "shifted10"]),
          "00:15 is NaN (41.0 already used at 00:10)")
    # Each value appears exactly once
    for val in [40.0, 41.0, 42.0]:
        cnt = (df["shifted10"] == val).sum()
        check(cnt == 1, f"value {val} appears exactly once (got {cnt})")

    # ------------------------------------------------------------------ #
    # Test 7: Single series
    # ------------------------------------------------------------------ #
    print("\nTest 7: Single series pass-through")
    s1 = pd.Series(
        [10.0, 11.0, 12.0],
        index=pd.date_range("2024-01-01", periods=3, freq="1h"),
        name="solo",
    )
    df, tmap = _synchronize_series([s1], max_shift=pd.Timedelta("5min"))

    check(df.shape == (3, 1), "shape is (3, 1)")
    check(df.columns.tolist() == ["solo"], "column is 'solo'")
    check(df.index.freq is not None, "frequency is set")

    # ------------------------------------------------------------------ #
    # Test 8: Returned index is a perfect DatetimeIndex with freq
    # ------------------------------------------------------------------ #
    print("\nTest 8: Index is a complete regular DatetimeIndex")
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

    check(isinstance(df.index, pd.DatetimeIndex), "index is DatetimeIndex")
    check(df.index.freq is not None, "index has a freq attribute")
    check(df.index[0] == pd.Timestamp("2024-01-01 00:00"), "starts at origin")
    check(df.index[-1] == pd.Timestamp("2024-01-01 00:20"), "ends at closing")

    # ------------------------------------------------------------------ #
    # Summary
    # ------------------------------------------------------------------ #
    print(f"\n{'=' * 60}")
    print(f"Results: {passed} passed, {failed} failed out of {passed + failed}")
    if failed:
        raise SystemExit(1)
