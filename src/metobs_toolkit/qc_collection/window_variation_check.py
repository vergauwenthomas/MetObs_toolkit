import logging
from typing import Union
import pandas as pd

from .common_functions import test_moving_window_condition

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def window_variation_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
    max_increase_per_second: Union[int, float],
    max_decrease_per_second: Union[int, float],
) -> pd.DatetimeIndex:
    """
    Test if the increase or decrease in a time window exceeds a threshold.

    This function checks if the variation of observations in time does not exceed a threshold.
    It applies a moving window over the time series, defined by a duration (`timewindow`), and tests
    if the window contains at least a minimum number of records (`min_records_per_window`).

    If the observations in the window increase or decrease more than a threshold, all
    observations in the window are flagged as outliers. The threshold is defined by the
    maximum increase or decrease per second multiplied by the window size in seconds.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index
        should be datetime-like.
    timewindow : pd.Timedelta
        The duration of the moving window. This should be a pandas Timedelta object.
    min_records_per_window : int
        The minimum number of non-NaN records required within the time window for the check to be valid.
        This is dependent on the time resolution of the records.
    max_increase_per_second : int or float
        The maximum allowed increase (per second). This value is extrapolated to the window duration.
        This value must be positive.
    max_decrease_per_second : int or float
        The maximum allowed decrease (per second). This value is extrapolated to the window duration.
        This value must be negative.

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    Notes
    -----
    In general, for temperatures, the decrease threshold is set less stringent than the increase
    threshold. This is because a temperature drop is meteorologically more
    common than a sudden increase, which is often the result of a radiation error.
    A suitable value for the min_records_per_window depends on the time resolution of the records and the window size.
    This check is similar to the step check, but not identical. The step check tests a maximum allowed increase or decrease
    with respect to the previous value. The window variation check uses a moving window to test the maximum allowed variation.
    """

    # Validate argument values
    if max_decrease_per_second > 0:
        raise ValueError("max_decrease_per_second must be negative!")
    if max_increase_per_second < 0:
        raise ValueError("max_increase_per_second must be positive!")

    # Test if the conditions for the moving window are met by the records frequency
    is_met = test_moving_window_condition(
        records=records,
        windowsize=timewindow,
        min_records_per_window=min_records_per_window,
    )
    if not is_met:
        logger.warning(
            "The minimum number of window members for the window variation check is not met!"
        )
        return pd.DatetimeIndex(name="datetime", data=[])

    # Drop outliers from the series (these are NaNs)
    input_series = records.dropna()

    # Calculate window thresholds (by linear extrapolation)
    max_window_increase = (
        float(abs(max_increase_per_second)) * timewindow.total_seconds()
    )
    max_window_decrease = (
        float(abs(max_decrease_per_second)) * timewindow.total_seconds()
    )

    # Define window test (applied on the windows)
    @log_entry
    def variation_test(window: pd.Series) -> int:
        """
        Test if the variation in the window exceeds the allowed increase or decrease.

        Parameters
        ----------
        window : pd.Series
            The window of values to test.

        Returns
        -------
        int
            1 if the window contains an outlier, 0 otherwise.
        """
        if (max(window) - min(window) > max_window_increase) & (
            window.idxmax() > window.idxmin()
        ):
            return 1

        if (max(window) - min(window) > max_window_decrease) & (
            window.idxmax() < window.idxmin()
        ):
            return 1
        else:
            return 0

    # Apply rolling window
    window_outliers = input_series.rolling(
        window=timewindow,
        closed="both",
        center=True,
        min_periods=min_records_per_window,
    ).apply(variation_test)
    logger.debug("Exiting function window_variation_check")
    return window_outliers.loc[window_outliers == 1].index
