import logging
from typing import Union
import pandas as pd


from .common_functions import create_qcresult_flags
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.qcresult import QCresult
from metobs_toolkit.qc_collection.whitelist import SensorWhiteSet

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def step_check(
    records: pd.Series,
    max_increase_per_second: Union[int, float],
    max_decrease_per_second: Union[int, float],
    sensorwhiteset: SensorWhiteSet,
) -> QCresult:
    """
    Check for 'spikes' and 'dips' in a time series.

    Tests if observations produce spikes in the time series. The maximum
    allowed increase and decrease per second is set in the arguments,
    and is tested for each record (with respect to the previous record).

    If the difference between two consecutive records (i.e., the spike or dip) is larger than the
    threshold, the record is flagged as an outlier.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index should be datetime-like.
    max_increase_per_second : int or float,
        The maximum allowed increase (per second). This value is extrapolated to the time resolution of records.
        This value must be positive.
    max_decrease_per_second : int or float
        The maximum allowed decrease (per second). This value is extrapolated to the time resolution of records.
        This value must be negative.
    sensorwhiteset : SensorWhiteSet
        A SensorWhiteSet instance containing timestamps that should be excluded from outlier detection.
        Records matching the whiteset criteria will not be flagged as outliers even if they meet the
        step check criteria.

    Returns
    -------
    QCresult
        Quality control result object containing flags, outliers, and details
        for the step check.

    Notes
    -----
    In general, for temperatures, the decrease threshold is set less stringent than the increase
    threshold. This is because a temperature drop is meteorologically more
    common than a sudden increase, which is often the result of a radiation error.
    """
    checksettings = {
        "max_increase_per_second": max_increase_per_second,
        "max_decrease_per_second": max_decrease_per_second,
        "sensorwhiteset": sensorwhiteset,
    }
    
    # Validate argument values
    if max_decrease_per_second > 0:
        raise ValueError("max_decrease_per_second must be negative!")
    if max_increase_per_second < 0:
        raise ValueError("max_increase_per_second must be positive!")

    # Drop outliers from the series (these are NaNs)
    to_check_records = records.dropna()

    # Calculate timedelta between rows
    time_diff = to_check_records.index.to_series().diff()

    # Define filter
    step_filter = (
        # Step increase
        (
            (to_check_records - to_check_records.shift(1))
            > (float(max_increase_per_second) * time_diff.dt.total_seconds())
        )  # or
        |
        # Step decrease
        (
            (to_check_records - to_check_records.shift(1))
            < (max_decrease_per_second * time_diff.dt.total_seconds())
        )
    )

    outliers_idx = step_filter[step_filter].index

    # Catch the white records
    outliers_after_white_idx = sensorwhiteset.catch_white_records(outliers_idx)

    flags = create_qcresult_flags(
        all_input_records=records,
        unmet_cond_idx = pd.DatetimeIndex([]),
        outliers_before_white_idx=outliers_idx,
        outliers_after_white_idx=outliers_after_white_idx,
    )
    

    qcresult = QCresult(
        checkname="step",
        checksettings=checksettings,
        flags=flags,
        detail='no details'
        )
    
    #Create and add details
    if not outliers_after_white_idx.empty:
        detailseries = pd.Series(
            data = f'step > {max_increase_per_second:.4g} per second or step < {max_decrease_per_second:.4g} per second',
            index = outliers_after_white_idx
        )
        qcresult.add_details_by_series(detail_series = detailseries)
    return qcresult
    
    