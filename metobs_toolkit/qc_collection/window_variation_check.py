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
    """Test if the variation exceeds the threshold in moving time windows.

    Looking for jumps of the values of an observation type that are larger than
    the limit specified in the qc_settings. These values are removed from the
    input series and combined in the outlier df.

    There is an increament threshold (that is if there is a max value difference
    and the maximum value occurred later than the minimum value occurred.)
    And vice versa is there a decrement threshold.

    The check is only applied if there are at least N observations in the time window.

    Schematically:

    1. Find the stations that have a maximum assumed observation frequency
       that does not exceed the minimum number of records for moving window
       size. The window size is defined by a duration.
    2. Compute the maximum increase and decrease thresholds for a window.
       This is done by multiplying the maximum increase per second by the
       window size in seconds.
    3. For each station, a moving window scan is applied that validates if
       the maximum increase/decrease thresholds are exceeded. This is done
       by comparison of the minimum and maximum values inside the window. The
       validation is only applied when a sufficient amount of records are
       found in the window specified by a threshold.
    4. After the scan, *all* records found in the window that exceed one
       of these thresholds are labeled as outliers.

    Parameters
    ------------
    station_frequencies : pandas.Series
        The frequencies of all the stations. This is a column in the metadf
        attribute of the Dataset.
    obsdf : pandas.DataFrame
        The observations dataframe (Dataset.df) to check. Must have a triple
        index (name, obstype, datetime).
    obstype : str
        The observation type to check for outliers.
    checks_settings : dict
        The dictionary containing the settings for the quality control checks.


    Returns
    ----------
    obsdf : pandas.DataFrame()
        The observations dataframe updated for window-variation-outliers. These are
        represented by Nan values.
    outl_df : pandas.DataFrame
        The updated outliersdf.

    """

    assert max_decrease_per_second < 0, f"max_decrease_per_second must be negative!"

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
