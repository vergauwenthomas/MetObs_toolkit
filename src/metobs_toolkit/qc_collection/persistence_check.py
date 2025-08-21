import logging

import numpy as np
import pandas as pd

from .common_functions import test_moving_window_condition

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def persistence_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
) -> pd.DatetimeIndex:
    """
    Check if values are not constant in a moving time window.

    Performs a persistence check on a time series to identify periods where observations remain constant
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

    * The function uses a rolling window approach to check if all non-NaN values within the window
      are identical.
    * If the minimum number of records per window is locally not met, the function logs a warning and skips
      the persistence check.
    * This function can be computationally expensive for large datasets or small time windows.

    Warnings
    --------
    If the minimum number of records per window is not met over the full time series, a warning is logged, and the function
    returns an empty DatetimeIndex.
    """

    # Test if the conditions for the moving window are met by the records frequency
    is_met = test_moving_window_condition(
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
    @log_entry
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
        a = window.values
        a = a[~np.isnan(a)]
        return (a[0] == a).all() if len(a) > 0 else False

    # This is very expensive if no coarsening is applied! Can we speed this up?
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

    logger.debug("Exiting function persistence_check")
    return window_is_constant[window_is_constant].index
