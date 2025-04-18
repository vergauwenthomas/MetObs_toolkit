import logging
from typing import Union

import numpy as np
import pandas as pd

from .common_functions import test_moving_window_condition

# Set up logging
logger = logging.getLogger("<metobs_toolkit>")


def persistence_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
) -> pd.DatetimeIndex:
    """Check if values are not constant in a moving time window.

    Perform a persistence check on a time series to identify periods where observations remain constant
    within a specified time window. If the values are constant, all records in the moving window are
    flagged as outliers.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index should be datetime-like.
    timewindow : pd.Timedelta
        The size of the rolling time window to check for persistence.
    min_records_per_window : int
        The minimum number of non-NaN records required within the time window for the check to be valid.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    Notes
    -----
    - The function uses a rolling window approach to check if all non-NaN values within the window
      are identical.
    - If the minimum number of records per window is locally not met, the function logs a warning and skips
      the persistence check.
    - This function can be computationally expensive for large datasets or small time windows.

    Warnings
    --------
    If the minimum number of records per window is not met over the full time series, a warning is logged, and the function
    returns an empty DatetimeIndex.
    """
    logger.info("Entering the function: persistence_check")

    # Validate argument types
    if not isinstance(records, pd.Series):
        raise TypeError("Argument 'records' must be of type pd.Series")
    if not isinstance(timewindow, pd.Timedelta):
        raise TypeError("Argument 'timewindow' must be of type pd.Timedelta")
    if not isinstance(min_records_per_window, int):
        raise TypeError("Argument 'min_records_per_window' must be of type int")

    # Test if the conditions for the moving window are met by the records frequency
    is_met = test_moving_window_condition(  # TYPO corrected: was 'ismet'
        records=records,
        windowsize=timewindow,
        min_records_per_window=min_records_per_window,
    )
    if not is_met:
        logger.warning(
            "The minimum number of window members for the persistence check is not met!"
        )
        return pd.DatetimeIndex(name="datetime", data=[])

    # Apply persistence
    def is_unique(window: pd.Series) -> bool:
        """
        Check if all non-NaN values in the window are identical.

        Parameters
        ----------
        window : pd.Series
            A pandas Series representing the rolling window.

        Returns
        -------
        bool
            True if all non-NaN values are identical, False otherwise.
        """
        logger.debug("Entering the helper function: is_unique")
        a = window.values
        a = a[~np.isnan(a)]
        return (a[0] == a).all()

    # TODO: This is very expensive if no coarsening is applied! Can we speed this up?
    window_is_constant = (
        records.dropna()  # Exclude outliers and gaps
        .rolling(
            window=timewindow,
            closed="both",
            center=True,
            min_periods=min_records_per_window,
        )
        .apply(is_unique)
    )
    # The returns are numeric values (0 --> False, NaN --> not checked (members/window condition not met), 1 --> outlier)
    window_is_constant = window_is_constant.map({0.0: False, np.nan: False, 1.0: True})

    logger.info("Exiting the function: persistence_check")
    return window_is_constant[window_is_constant].index
