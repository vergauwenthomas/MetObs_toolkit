import logging
from typing import Union
import pandas as pd
from numpy import nan

from .common_functions import test_moving_window_condition, create_qcresult_flags
from .whitelist import SensorWhiteSet
from metobs_toolkit.qcresult import QCresult
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.backend_collection.datetime_collection import (
    timestamps_to_datetimeindex,
)

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def window_variation_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
    max_increase_per_second: Union[int, float],
    max_decrease_per_second: Union[int, float],
    sensorwhiteset: SensorWhiteSet,
) -> QCresult:
    """
    Test if the increase or decrease in a time window exceeds a threshold.

    This function checks if the variation of observations in time does exceeds a threshold.
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
    sensorwhiteset : SensorWhiteSet, optional
        A SensorWhiteSet instance containing timestamps that should be excluded from outlier detection.
        Records matching the whiteset criteria will not be flagged as outliers even if they meet the
        window variation check criteria.

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

    # Drop outliers from the series (these are NaNs)
    to_check_records = records.dropna()
    
    # Test if the conditions for the moving window are met by the records frequency
    is_met = test_moving_window_condition(
        records=records, #pass records, because freq is estimated
        windowsize=timewindow,
        min_records_per_window=min_records_per_window,
    )
    if not is_met:
        logger.warning(
            "The minimum number of window members for the window variation check is not met!"
        )
        flags = create_qcresult_flags(
            all_input_records=records,
            unmet_cond_idx=to_check_records.index,
            outliers_before_white_idx=pd.DatetimeIndex([]),
            outliers_after_white_idx=pd.DatetimeIndex([]),
        )
        
        qcresult = QCresult(
            checkname="window_variation",
            checksettings=locals().pop('records', None),
            flags=flags,
            outliers=timestamps_to_datetimeindex(
                name="datetime", timestamps=[], current_tz=None
                ),
            detail="Minimum number of records per window not met.",
        )
        return qcresult

    

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
    window_flags = to_check_records.rolling(
        window=timewindow,
        closed="both",
        center=True,
        min_periods=min_records_per_window,
    ).apply(variation_test)

    # The returns are numeric values (0 --> oke, NaN --> not checked (members/window condition not met), 1 --> outlier)
    window_flags = window_flags.map(
        {0.0: 'pass', #Dummy label
         nan: 'unmet', #Dummy label
         1.0: 'flagged'}) #Dummy label     
    
    # Filter outliers
    outliers_idx = window_flags.loc[window_flags == 'flagged'].index
    # Catch the white records
    outliers_after_white_idx = sensorwhiteset.catch_white_records(outliers_idx=outliers_idx)

    #Create flags
    flags = create_qcresult_flags(
        all_input_records=records,
        unmet_cond_idx=window_flags[window_flags == 'unmet'].index,    
        outliers_before_white_idx=outliers_idx,
        outliers_after_white_idx=outliers_after_white_idx)
    
    qcresult = QCresult(
        checkname="window_variation",
        checksettings=locals().pop('records', None),
        flags=flags,
        outliers = records.loc[outliers_after_white_idx],
        detail='no details'
        )
    
    #Create and add details
    if not outliers_after_white_idx.empty:
        detailseries = pd.Series(
            data = f'Variation in {timewindow} window exceeds max increase of {max_window_increase} or max decrease of {max_window_decrease}.',
            index = outliers_after_white_idx
        )
        qcresult.add_details_by_series(detail_series = detailseries)
    
    return qcresult
