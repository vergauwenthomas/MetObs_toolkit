import logging

import pandas as pd


logger = logging.getLogger("<metobs_toolkit>")



def drop_invalid_values(records: pd.Series, skip_records: pd.DatetimeIndex) -> pd.Series:
    """Remove invalid (non-numeric) values from a time series.
    
    Filters out values that could not be cast to numeric types. Invalid timestamps
    are treated as gaps and removed from the series rather than being flagged as
    outliers. This allows the gap detection mechanism to handle them appropriately.
    
    Parameters
    ----------
    records : pd.Series
        Series with a datetime-like index containing values to validate.
    skip_records : pd.DatetimeIndex
        Records to temporarily exclude from the check (typically duplicated timestamps).
        These records are preserved regardless of validity and added back after filtering.
        
    Returns
    -------
    pd.Series
        Filtered series containing only records with valid numeric values,
        plus all skipped records.
        
    Notes
    -----
    * Invalid values are interpreted as missing data (gaps) rather than outliers.
    * Skipped records are preserved to avoid interfering with prior QC checks.
    * This function does not raise an error if the check was previously applied.
    """
    skipped_data = records.loc[skip_records]
    targets = records.drop(skip_records)

    # Option 1: Create a outlier label for these invalid inputs,
    # and treat them as outliers
    # outlier_timestamps = targets[~targets.notnull()]

    # self._update_outliers(
    #     qccheckname="invalid_input",
    #     outliertimestamps=outlier_timestamps.index,
    # and treat them as outliers
    #     extra_columns={},
    #     overwrite=False,
    # )

    # Option 2: Since there is not numeric value present, these timestamps are
    # interpreted as gaps --> remove the timestamp, so that it is captured by the
    # gap finder.

    # Note: do not treat the first/last timestamps differently. That is
    # a philosiphycal choice.

    validrecords = targets[targets.notnull()]  # subset to numerical casted values
    # add the skipped records back
    validrecords = pd.concat([validrecords, skipped_data]).sort_index()
    return validrecords
        