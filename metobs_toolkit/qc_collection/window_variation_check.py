import logging
import pandas as pd

from .common_functions import test_moving_window_condition

logger = logging.getLogger(__name__)


def window_variation_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
    max_increase_per_second: int | float,
    max_decrease_per_second: int | float,
) -> pd.DatetimeIndex:
    """Test if the increase/decrease in a timewindow exceeds a threshold.

    This function is used to check if the variation of observations in time,
    does not exceed a threshold. This is done by applying a moving window
    over the time series. The moving window is defined by a duration (timewindow),
    and tested if the window contains at least a minimum number of records.

    If the observations in the window increase/decrease more than a threshold, all
    observations in the window are flagged as outliers. The threshold is defined by the
    maximum increase/decrease per second multiplied by the window size in seconds.

    Parameters
    ------------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index
        should be datetime-like.
    timewindow : pd.Timedelta
        The duration of the moving window. This should be a pandas Timedelta object.
    min_records_per_window : int
        The minimum number of non-NaN records required within the time window for the check to be valid.
        This is dependent on the time resolution of the records.
    max_increase_per_second : int | float, >0
        The maximum allowed increase (per second). This value is extrapolated to the window duration.
        This value must be positive!
    max_decrease_per_second : int | float, <0
        The maximum allowed decrease (per second). This value is extrapolated to the window duration.
        This value must be negative!

    Returns
    -------
    pd.DatetimeIndex
        Timestamps of outlier records.

    Notes
    -----
    - In general, for temperatures, the decrease threshold is set less stringent than the increase
      threshold. This is because a temperature drop is meteorologically more
      common than a sudden increase which is often the result of a radiation error.
    - A suitable value for the min_records_per_window depends on the time resolution of the records and the window size.
    - This check is similar to the step check, but not identical. The step check a maximum allowed increase/decrease
      with resprect to the previous value. The window variation check uses a moving window to test the maximum allowd variation.
    """

    # Validate argument types
    if not isinstance(records, pd.Series):
        raise TypeError("Argument 'records' must be of type pd.Series")
    if not isinstance(timewindow, pd.Timedelta):
        raise TypeError("Argument 'timewindow' must be of type pd.Timedelta")
    if not isinstance(min_records_per_window, int):
        raise TypeError("Argument 'min_records_per_window' must be of type int")

    # Validate argument values
    if max_decrease_per_second > 0:
        raise ValueError("max_decrease_per_second must be negative!")
    if max_increase_per_second < 0:
        raise ValueError("max_increase_per_second must be positive!")

    # test if the conditions for the moving window are met by the records frequency
    ismet = test_moving_window_condition(
        records=records,
        windowsize=timewindow,
        min_records_per_window=min_records_per_window,
    )
    if not ismet:
        logger.warning(
            "The minimum number of window members for the window variation check is not met!"
        )
        return pd.DatetimeIndex(name="datetime", data=[])

    # drop outliers from the series (these are Nan's)
    input_series = records.dropna()

    # Calculate window thresholds (by linear extarpolation)
    max_window_increase = (
        float(abs(max_increase_per_second)) * timewindow.total_seconds()
    )
    max_window_decrease = (
        float(abs(max_decrease_per_second)) * timewindow.total_seconds()
    )

    # define window test (applied on the windows)
    def variation_test(window):
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

    # apply rolling window
    window_outliers = input_series.rolling(
        window=timewindow,
        closed="both",
        center=True,
        min_periods=min_records_per_window,
    ).apply(variation_test)
    return window_outliers.loc[window_outliers == 1].index
