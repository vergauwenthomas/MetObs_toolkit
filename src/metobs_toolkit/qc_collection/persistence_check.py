import logging

import numpy as np
import pandas as pd

from .common_functions import test_moving_window_condition, create_qcresult_flags
from .whitelist import SensorWhiteSet
from metobs_toolkit.qcresult import (
    QCresult,
    pass_cond,
    flagged_cond,
    unmet_cond,
)
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.backend_collection.datetime_collection import (
    timestamps_to_datetimeindex,
)

logger = logging.getLogger("<metobs_toolkit>")


def _has_window_unique_values(window: pd.Series) -> bool:
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


@log_entry
def persistence_check(
    records: pd.Series,
    timewindow: pd.Timedelta,
    min_records_per_window: int,
    sensorwhiteset: SensorWhiteSet,
) -> QCresult:
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
    sensorwhiteset : SensorWhiteSet
        A SensorWhiteSet instance containing timestamps that should be excluded from outlier detection.
        Records matching the whiteset criteria will not be flagged as outliers even if they meet the
        persistence criteria.

    Returns
    -------
    QCresult
        Quality control result object containing flags, outliers, and details
        for the persistence check.

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

    to_check_records = records.dropna() # Exclude outliers and gaps
    # Test if the conditions for the moving window are met by the records frequency
    is_met = test_moving_window_condition(
        records=records, #pass records, because freq is estimated
        windowsize=timewindow,
        min_records_per_window=min_records_per_window,
    )
    if not is_met:
        logger.warning(
            "The minimum number of window members for the persistence check is not met!"
        )
        flags = create_qcresult_flags(
            all_input_records=records,
            unmet_cond_idx=to_check_records.index,
            outliers_before_white_idx=pd.DatetimeIndex([]),
            outliers_after_white_idx=pd.DatetimeIndex([]),
        )
        qcresult = QCresult(
            checkname="persistence",
            checksettings=locals().pop('records', None),
            flags=flags,
            outliers=pd.Series(index=timestamps_to_datetimeindex(
                                        name="datetime", timestamps=[], current_tz=None)),
            detail=f"Minimum number of records ({min_records_per_window}) per window ({timewindow}) not met.",
        )
        return qcresult
            

    # This is very expensive if no coarsening is applied! Can we speed this up?
    
    window_flags = (
        to_check_records  
        .rolling(
            window=timewindow,
            closed="both",
            center=True,
            min_periods=min_records_per_window,
        )
        .apply(_has_window_unique_values)
    )
    # The returns are numeric values (0 --> oke, NaN --> not checked (members/window condition not met), 1 --> outlier)
    window_flags = window_flags.map(
        {0.0: pass_cond,
         np.nan: unmet_cond,
         1.0: flagged_cond})        
    
    outliers_idx = window_flags[window_flags == flagged_cond].index

    # Catch the white records
    outliers_after_white_idx = sensorwhiteset.catch_white_records(outliers_idx=outliers_idx)
   
    #Create flags
    flags = create_qcresult_flags(
        all_input_records=records,
        unmet_cond_idx=window_flags[window_flags == unmet_cond].index,    
        outliers_before_white_idx=outliers_idx,
        outliers_after_white_idx=outliers_after_white_idx)
    
    qcresult = QCresult(
        checkname="persistence",
        checksettings=locals().pop('records', None),
        flags=flags,
        outliers = records.loc[outliers_after_white_idx],
        detail='no details'
        )
    
    #Create and add details
    if not outliers_after_white_idx.empty:
        detailseries = pd.Series(
            data = 'constant values in timewindow of ' + str(timewindow),
            index = outliers_after_white_idx
        )
        qcresult.add_details_by_series(detail_series = detailseries)
    
    return qcresult
