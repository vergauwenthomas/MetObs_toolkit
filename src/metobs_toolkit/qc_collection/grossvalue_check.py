import logging
from typing import Union
import pandas as pd

from .common_functions import create_qcresult_flags
from .whitelist import SensorWhiteSet
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.qcresult import QCresult

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def gross_value_check(
    records: pd.Series,
    lower_threshold: Union[int, float],
    upper_threshold: Union[int, float],
    sensorwhiteset: SensorWhiteSet,
) -> QCresult:
    """
    Identify outliers in a time series based on lower and upper thresholds.

    Parameters
    ----------
    records : pd.Series
        Series with a datetime-like index to check.
    lower_threshold : int or float
        Threshold below which records are flagged as outliers.
    upper_threshold : int or float
        Threshold above which records are flagged as outliers.
    sensorwhiteset : SensorWhiteSet
        A SensorWhiteSet instance containing timestamps that should be excluded from outlier detection.
        Records matching the whiteset criteria will not be flagged as outliers even if they fall
        outside the threshold range.

    Returns
    -------
    QCresult
        Quality control result object containing flags, outliers, and details
        for the gross value check.
    """

    # Drop NaN values
    to_check_records = records.dropna()
    # Identify outliers
    outliers_idx = to_check_records[
        (to_check_records < lower_threshold) | (to_check_records > upper_threshold)
    ].index

    # Exclude white records if provided
    outliers_after_white_idx = sensorwhiteset.catch_white_records(outliers_idx=outliers_idx)

    # Create QCresult flags
    flags = create_qcresult_flags(
        all_input_records=records,
        unmet_cond_idx = pd.DatetimeIndex([]),
        outliers_before_white_idx=outliers_idx,
        outliers_after_white_idx=outliers_after_white_idx,
    )

    checksettings = {
        "lower_threshold": lower_threshold,
        "upper_threshold": upper_threshold,
        "sensorwhiteset": sensorwhiteset,
    }
    
    qcresult = QCresult(
        checkname="gross_value",
        checksettings=checksettings,
        flags=flags,
        detail='no details'
        )
    
    #Create and add details
    if not outliers_after_white_idx.empty:
        detailseries = pd.Series(
            data = 'value outside gross value thresholds [' + str(lower_threshold) + ', ' + str(upper_threshold) + ']',
            index = outliers_after_white_idx
        )
        qcresult.add_details_by_series(detail_series = detailseries)
    return qcresult
