from __future__ import annotations
import logging
from typing import Literal, Union, TYPE_CHECKING

import warnings


import pandas as pd
import numpy as np

if TYPE_CHECKING:
    from dateutil.tz import tzfile
    from datetime import tzinfo


def timestamps_to_datetimeindex(
    timestamps: np.ndarray | pd.Series | pd.DatetimeIndex | list,
    tz: Union[tzfile | tzinfo | str] = "UTC",
    **kwargs,
) -> pd.DatetimeIndex:
    """
    Convert timestamps to a timezone-aware DatetimeIndex.

    Parameters
    ----------
    timestamps : np.ndarray, pd.Series, pd.DatetimeIndex, or list
        The timestamps to convert.
    tz : timezone, tz.tzfile, tzinfo, or str
        The timezone to apply to the DatetimeIndex.

    Returns
    -------
    pd.DatetimeIndex
        A timezone-aware DatetimeIndex.
    """
    # Create DatetimeIndex without explicitly putting a timezone
    dt_index = pd.DatetimeIndex(data=timestamps, **kwargs)

    if dt_index.empty:
        return dt_index  # no need to convert timezone for empty index
    else:
        # Convert to target timezone
        return convert_timezone(datetimeindex=dt_index, target_tz=tz)


def convert_timezone(
    datetimeindex: pd.DatetimeIndex,
    target_tz: Union[tz.tzfile | tzinfo | str] = "UTC",
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
        warnings.warn(
            "The DatetimeIndex has no timezone information. " "Assuming UTC timezone.",
            UserWarning,
            stacklevel=2,
        )
        datetimeindex = datetimeindex.tz_localize("UTC")

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
