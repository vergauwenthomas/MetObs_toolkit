import logging
import pandas as pd

from metobs_toolkit.qcresult import QCresult
from .common_functions import create_qcresult_flags
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.backend_collection.datetime_collection import (
    timestamps_to_datetimeindex,
)
from metobs_toolkit.qc_collection.whitelist import SensorWhiteSet

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def repetitions_check(
    records: pd.Series,
    max_N_repetitions: int,
    sensorwhiteset: SensorWhiteSet,
) -> QCresult:
    """
    Test if an observation changes after a number of repetitions.

    This function checks if the observation changes after a number of repetitions.
    If a value is repeated more than the specified number of times, all the repeated
    records are flagged as outliers.

    Be aware that the performance of this check depends on the `max_N_repetitions`
    and the time resolution of the observations.

    Parameters
    ----------
    records : pd.Series
        A pandas Series containing the time series data to be checked. The index should be datetime-like.
    max_N_repetitions : int
        The maximum number of repetitions allowed before the records are flagged as outliers.
        If the number of repetitions exceeds this value, all repeated records are flagged as outliers.
    sensorwhiteset : SensorWhiteSet
        A SensorWhiteSet instance containing timestamps that should be excluded from outlier detection.
        Records matching the whiteset criteria will not be flagged as outliers even if they exceed the
        max_N_repetitions threshold.

    Returns
    -------
    QCresult
        Quality control result object containing flags, outliers, and details
        for the repetitions check.

    Notes
    -----
    The repetitions check is similar to the persistence check, but not identical.
    The persistence check uses thresholds that are meteorologically based (e.g., the moving window is defined by a duration),
    in contrast to the repetitions check whose thresholds are instrumentally based (e.g., the "window" is defined by a number of records).
    """
    
    checksettings = {
        "max_N_repetitions": max_N_repetitions,
        "sensorwhiteset": sensorwhiteset,
    }
    # Drop outliers from the series (these are NaNs)
    to_check_records = records.dropna()

    # Create group definitions for repeating values that do not change
    persistence_filter = ((to_check_records.shift() != to_check_records)).cumsum()
    persdf = pd.DataFrame(
        data={"value": to_check_records, "persistgroup": persistence_filter},
        index=to_check_records.index,
    )

    # Find outlier groups
    groups = persdf.groupby(["persistgroup"])
    # The above line groups the observations which have the same value and consecutive datetimes.
    group_sizes = groups.size()
    outlier_groups = group_sizes[group_sizes > max_N_repetitions]

    
    if outlier_groups.empty:
        logger.debug("No outliers detected. Exiting repetitions_check function.")
        outliers_idx =  timestamps_to_datetimeindex(
                        timestamps=[], name="datetime", current_tz=None
                    )
        outliers = pd.Series(index=outliers_idx)
    else:
    
        # Combine all outlier groups
        outliers = pd.concat(
            [
                groups.get_group(
                    outlgroup,
                )
                for outlgroup in outlier_groups.index
            ]
        )
       

    # Catch the white records
    outliers_after_white_idx = sensorwhiteset.catch_white_records(outliers.index)
    
    # Create flags
    flags = create_qcresult_flags(
        all_input_records=records,
        unmet_cond_idx = pd.DatetimeIndex([]),
        outliers_before_white_idx=outliers.index,
        outliers_after_white_idx=outliers_after_white_idx,
    )

    qcresult = QCresult(
        checkname="repetitions",
        checksettings=checksettings,
        flags=flags,
        )
    
    #Create and add details
    if not outliers_after_white_idx.empty:
        detailseries = pd.Series(
            data = 'More than ' + str(max_N_repetitions) + ' repeated values',
            index = outliers_after_white_idx
        )
        qcresult.add_details_by_series(detail_series = detailseries)
    
    return qcresult