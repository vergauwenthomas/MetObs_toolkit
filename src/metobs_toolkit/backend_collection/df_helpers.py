#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
A collection of functions on dataframe that are often used.

Created on Thu Mar  2 16:00:59 2023

@author: thoverga
"""
import logging
import pandas as pd
import numpy as np

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger(__name__)

pd.options.mode.copy_on_write = True


@log_entry
def convert_to_numeric_series(
    arr: "array-like", datadtype: type = np.float32
) -> pd.Series:
    """Convert an array to a numeric pandas series.

    Parameters
    ----------
    arr : array-like
        Input array.
    datadtype : type, optional
        Target numeric data type, by default np.float32

    Returns
    -------
    pd.Series
        Numeric pandas series.
    """
    # First try direct conversion to check if already numeric
    series = pd.Series(arr)
    if pd.api.types.is_numeric_dtype(series):
        numseries = series
    else:
        # Convert to string, replace commas by dots, and 'nan' by np.nan
        strseries = series.astype(str).str.replace(",", ".").replace("nan", np.nan)

        # Convert to numeric, setting errors to 'coerce' (invalid parsing will be set as NaN)
        numseries = pd.to_numeric(strseries, errors="coerce")

    # convert to target numeric type
    numseries = numseries.astype(datadtype)

    return numseries


@log_entry
def save_concat(targets, **kwargs):
    if not isinstance(targets, list):
        targets = list(targets)

    isempty = [tar.empty for tar in targets]

    # if some (or none) but not all are empty
    if not all(isempty):
        return pd.concat([tar for tar in targets if not tar.empty], **kwargs)

    # if all are empty
    else:
        # retrun an empty df with columns an index the union of all targets
        all_columns = set()
        all_index = set()
        for tar in targets:
            all_columns.update(tar.columns)
            all_index.update(tar.index)
        return pd.DataFrame(columns=list(all_columns), index=list(all_index))


@log_entry
def to_timedelta(inputdelta):
    """Convert input to a pandas timedelta object."""
    if isinstance(inputdelta, pd.Timedelta):
        return inputdelta
    elif isinstance(inputdelta, str):
        # Clue: When freq is extracted for datetimeindex,
        # the string representation does nt always start with
        # a number (ex: 'T' for minutes). This is not accepted by
        # pd.to_timedelta. Therefore, add a '1' in front of the string.
        if not inputdelta[0].isdigit():
            inputdelta = "1" + inputdelta
        return pd.to_timedelta(inputdelta)
    else:
        raise TypeError(f"Input {inputdelta} is not a valid timedelta.")
