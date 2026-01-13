import logging

import pandas as pd


from .common_functions import create_qcresult_flags
from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.qcresult import (
    QCresult,
    pass_cond,
    flagged_cond,

)

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def duplicated_timestamp_check(records: pd.Series) -> QCresult:
    """Check for duplicated timestamps in a time series.
    
    Identifies all records that share the same timestamp. All occurrences of 
    duplicated timestamps are flagged as outliers (not just subsequent ones).
    This check is performed before invalid value checking.
    
    Parameters
    ----------
    records : pd.Series
        Series with a datetime-like index to check for duplicates.
        
    Returns
    -------
    QCresult
        Quality control result object containing flags, outliers, and details
        for the duplicated timestamp check.
        
    Notes
    -----
    * All records with duplicated timestamps are flagged, including the first occurrence.
    * Values are coerced to numeric during this check to ensure compatibility with
      downstream processing.
    * Details include a list of all values sharing each duplicated timestamp.
    """
    # find all duplicates
    duplicates = records[records.index.duplicated(keep=False)]
    
    #Drop dulicates from series, they are a mess to take along
    no_dup_records = records[~records.index.duplicated(keep='first')]
    
    #create flags (no duplicates in the index!) 
    flags = create_qcresult_flags(
        all_input_records=no_dup_records, #NO duplicates here
        unmet_cond_idx = pd.DatetimeIndex([]),
        outliers_before_white_idx=duplicates.index.unique(),
        outliers_after_white_idx=duplicates.index.unique(),
    )
    
    
    
    #Special case: this check is performed before the invalid check, so values
    #must be cast to numeric to avoid issues when combining them in the outliersdf
    duplicates = pd.to_numeric(duplicates, errors='coerce')
    
    qcresult = QCresult(
        checkname="duplicated_timestamp",
        checksettings={},
        flags=flags,
        outliers =duplicates[~duplicates.index.duplicated(keep="first")],
        detail='no details')
    
    #Create and add details
    if not duplicates.empty:
        # For each duplicated timestamp, join all values as a comma-separated string
        details = (
            duplicates.groupby(duplicates.index)
            .apply(lambda x:"duplicated timestamp with: " + ", ".join(map(str, x.values)))
        )
        qcresult.add_details_by_series(detail_series = details)
    
    
    return qcresult