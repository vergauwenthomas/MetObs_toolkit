from __future__ import annotations
import logging
from typing import Literal, Union, Optional, TYPE_CHECKING

import warnings


import pandas as pd
import numpy as np

from metobs_toolkit.backend_collection.errorclasses import MetObsInternalError

if TYPE_CHECKING:
    from dateutil.tz import tzfile
    from datetime import tzinfo


def timestamps_to_datetimeindex(
    timestamps: np.ndarray | pd.Series | pd.DatetimeIndex | list,
    current_tz: Optional[Union[tzfile | tzinfo | str]] = None,
    name: Optional[str] = "datetime",
) -> pd.DatetimeIndex:
    """
    Convert timestamps to a timezone-aware DatetimeIndex.

    Parameters
    ----------
    timestamps : np.ndarray, pd.Series, pd.DatetimeIndex, or list
        The timestamps to convert.
    current_tz : timezone, tz.tzfile, tzinfo, str, or None
        The timezone of the input timestamps. Required if timestamps are
        timezone-naive. If timestamps are timezone-aware, this is used to
        validate consistency.

    Returns
    -------
    pd.DatetimeIndex
        A timezone-aware DatetimeIndex.

    Raises
    ------
    MetObsInternalError
        If timestamps are timezone-naive and current_tz is not provided.
        If timestamps are timezone-aware but timezone differs from current_tz.
    """
    # Create DatetimeIndex without explicitly putting a timezone
    dt_index = pd.DatetimeIndex(data=timestamps, name=name)

    if dt_index.empty:
        return dt_index  # no need to process empty index

    # Case 1: Timezone-naive timestamps
    if dt_index.tz is None:
        if current_tz is None:
            raise MetObsInternalError(
                "Timestamps are timezone-naive but no current_tz was provided. "
                "Please provide current_tz to make the timestamps timezone-aware."
            )
        # Localize to current_tz
        return dt_index.tz_localize(current_tz)

    # Case 2: Timezone-aware timestamps
    if current_tz is not None:
        # Check if timezone matches current_tz
        # Normalize timezone representations for comparison
        input_tz = dt_index.tz
        expected_tz = current_tz

        if str(input_tz) != str(expected_tz):
            raise MetObsInternalError(
                f"Timestamps have timezone '{input_tz}' but current_tz is "
                f"'{current_tz}'. Timezone mismatch detected."
            )

    return dt_index


def convert_timezone(
    datetimeindex: pd.DatetimeIndex,
    target_tz: Union[tzfile | tzinfo | str] = "UTC",
) -> pd.DatetimeIndex:
    """
    Convert a DatetimeIndex to a target timezone.

    Parameters
    ----------
    datetimeindex : pd.DatetimeIndex
        The DatetimeIndex to convert. If timezone-naive, it is assumed to be
        UTC and a warning is raised.
    target_tz : timezone, tz.tzfile, tzinfo, or str, optional
        The target timezone to convert to. Default is "UTC".

    Returns
    -------
    pd.DatetimeIndex
        The DatetimeIndex converted to the target timezone.
    """

    if datetimeindex.tz is None:
        raise MetObsInternalError(
            f"Connot convert a datetimeindex with no timezone to {target_tz}."
        )

    return datetimeindex.tz_convert(target_tz)


def to_timedelta(inputdelta):
    """Convert input to a pandas timedelta object."""
    if isinstance(inputdelta, pd.Timedelta):
        return inputdelta
    elif isinstance(inputdelta, str):
        # Clue: When freq is extracted for datetimeindex,
        # the string representation does not always start with
        # a number (ex: 'T' for minutes). This is not accepted by
        # pd.to_timedelta. Therefore, add a '1' in front of the string.
        if not inputdelta[0].isdigit():
            inputdelta = "1" + inputdelta
        return pd.to_timedelta(inputdelta)
    else:
        raise TypeError(f"Input {inputdelta} is not a valid timedelta.")
