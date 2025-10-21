import pandas as pd
import logging

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def test_moving_window_condition(
    records: pd.Series, windowsize: pd.Timedelta, min_records_per_window: int
) -> bool:
    """
    Test if the resolution of the records meets the window constraints.

    Parameters
    ----------
    records : pd.Series
        Series with a datetime-like index.
    windowsize : pd.Timedelta
        Size of the moving window.
    min_records_per_window : int
        Minimum number of records required per window.

    Returns
    -------
    bool
        True if the minimum window members condition is met, False otherwise.

    Raises
    ------
    TypeError
        If any argument is not of the expected type.
    Exception
        If the input records do not have a perfectly regular timestamp.
    """

    # Get frequency of records
    freqstr = pd.infer_freq(records.index)
    if freqstr is None:
        raise Exception("The input records do not have a perfectly regular timestamp.")
    # Convert to timedelta
    # Note: sometimes 'h' is returned, and this gives issues, so add a 1 in front
    if not freqstr[0].isdigit():
        freqstr = "1" + freqstr

    freq = pd.Timedelta(freqstr)

    # Test if minimum window members condition is met
    ismet = (windowsize / freq) >= min_records_per_window
    logger.debug("Exiting function test_moving_window_condition.")
    return ismet


def catch_white_records(outliers_idx: pd.DatetimeIndex,
                        white_records: pd.DatetimeIndex) -> pd.DatetimeIndex:
    """Remove white record timestamps from outliers index.
    
    Parameters
    ----------
    outliers_idx : pd.DatetimeIndex
        Index of outlier timestamps to filter
    white_records : pd.DatetimeIndex  
        Index of white record timestamps to exclude
        
    Returns
    -------
    pd.DatetimeIndex
        Filtered outliers index with white records removed
    """
    outliers = outliers_idx.difference(white_records)
    outliers.name = 'datetime'
    return outliers

def fmt_white_records_for_station_qc(white_records: pd.Index,
                                     trg_obstype: str,
                                     trg_station: str) -> pd.DatetimeIndex:
    """Format and filter white records for station-level QC methods.
    
    This function processes a white_records index (which may be a MultiIndex with 
    'name', 'datetime', and 'obstype' levels) and extracts only the datetime values 
    that correspond to the specified station and observation type. This is necessary 
    because station-level QC methods only need datetime information, not the full 
    MultiIndex structure.
    
    Parameters
    ----------
    white_records : pd.Index
        Input white records index. Can be a DatetimeIndex or a MultiIndex containing
        'name', 'datetime', and/or 'obstype' levels. If MultiIndex, it will be filtered
        to the target station and observation type.
    trg_obstype : str
        The target observation type to filter for (e.g., 'temp', 'humidity').
    trg_station : str
        The target station name to filter for.
        
    Returns
    -------
    pd.DatetimeIndex
        A 1-dimensional DatetimeIndex with name 'datetime' containing only the 
        timestamps relevant for the specified station and observation type.
        
    Raises
    ------
    TypeError
        If white_records is not a pd.Index.
    ValueError
        If after filtering by station and obstype, the resulting index still has
        multiple levels (indicating an unexpected structure).
        
    Notes
    -----
    The function performs the following operations:
    
    1. Filters by 'name' level if present, keeping only rows matching trg_station
    2. Filters by 'obstype' level if present, keeping only rows matching trg_obstype
    3. Drops the 'name' and 'obstype' levels if they existed
    4. Ensures the final result is a 1D DatetimeIndex named 'datetime'
    """
    if not isinstance(white_records, pd.Index):
        raise TypeError("white_records must be a pd.Index")
    
    #Filter multi-index levels if present
    if 'name' in white_records.names:
        white_records = white_records[white_records.get_level_values('name') == trg_station]
        white_records = white_records.droplevel('name')
        
    if 'obstype' in white_records.names:
        white_records = white_records[white_records.get_level_values('obstype') == trg_obstype]
        white_records = white_records.droplevel('obstype')

    #at this point, the index must be 1D representing 'datetime'
    #Only the datetime index level is relevant for stations-related QC methods
    if white_records.nlevels == 1:
        white_records.name = 'datetime'
        #TODO: make sure it is a datetime index
    else:
        raise ValueError("white_records must only contain 'datetime' (or possibly also 'name' and 'obstype' levels) after filtering.")
    
    return white_records    