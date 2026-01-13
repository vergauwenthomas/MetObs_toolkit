from __future__ import annotations
import logging
from typing import Literal, Union, TYPE_CHECKING

import numpy as np
import pandas as pd

from metobs_toolkit.backend_collection.decorators import log_entry
from metobs_toolkit.settings_collection.settings import Settings
logger = logging.getLogger("<metobs_toolkit>")


pass_cond = 'passed' #checked and successfull pass
flagged_cond = 'flagged' # checked and flagged as outlier
unmet_cond = 'condition_unmet' #not checked due to unmet specific conditions
saved_cond = 'saved' #checked and flagged but saved due to whitelist
unchecked_cond = 'unchecked' #not checked (was nan/gap before check)

class QCresult:
    """Store results of a quality control check.
    
    This class encapsulates the results of a single QC check including flags,
    detected outliers, and detailed information about the check outcome for
    all timestamps.
    
    Parameters
    ----------
    checkname : str
        Name identifying the quality control check (e.g., 'gross_value', 'persistence').
    checksettings : dict
        Dictionary of parameters and settings used for this QC check.
    flags : pd.Series
        Series with datetime index containing QC flag strings for all timestamps:
        'passed', 'flagged', 'condition_unmet', 'saved', or 'unchecked'.
    outliers : pd.Series
        Series with datetime index containing outlier values. Index should be
        a subset of the flags index.
    detail : str, optional
        Default detail string for all timestamps. Can be updated for specific
        timestamps using add_details_by_series. Default is empty string.
        
    Attributes
    ----------
    checkname : str
        Name of the QC check.
    checksettings : dict
        Settings used for the check.
    flags : pd.Series
        QC flags for all timestamps.
    outliers : pd.Series
        Detected outlier values.
    details : pd.Series
        Detailed information for each timestamp.
    """

    def __init__(
        self,
        checkname: str,
        checksettings: dict,
        flags: pd.Series,  # index: timestamps, values: 'passed', 'flagged', 'condition_unmet', 'saved'
        outliers: pd.Series,  # index: timestamps, values: outlier values
        detail: str = "",
    ):
        self.checkname = checkname
        self.checksettings = checksettings
        
        if not isinstance(flags.index, pd.DatetimeIndex):
            raise TypeError("The index of 'flags' must be a pandas.DatetimeIndex.")
        self.flags = flags
        
        
        
        if not isinstance(outliers.index, pd.DatetimeIndex):
            raise TypeError("The index of 'outliers' must be a pandas.DatetimeIndex.")
        self.outliers = outliers 
    
        #Set details (Index is Flags thus includes all timestamps!)
        self.details = pd.Series([detail] * len(flags),
                                index=flags.index)
        
        
    def __repr__(self) -> str:
        return f"QCresult(checkname={self.checkname})"
    
    
    @log_entry
    def add_details_by_series(self, detail_series: pd.Series) -> None:
        """Update the details attribute with values from a detail_series.
        
        This method updates the details attribute (a pandas Series with datetime 
        index) using index-value pairs from the provided detail_series. The 
        detail_series index must be a subset of the details attribute index.
        
        Parameters
        ----------
        detail_series : pd.Series
            A pandas Series with datetime index containing detail values to 
            update. The index should be a subset of the details attribute index.
            
        Raises
        ------
        TypeError
            If detail_series is not a pandas Series or if its index is not a 
            pandas DatetimeIndex.
        """
       
        # Update details using the index-value pairs from detail_series
        self.details.update(detail_series)
    
    def get_outlier_timestamps(self) -> pd.DatetimeIndex:
        """Return the timestamps of the outliers."""
        return self.outliers.index
        
        
    def remap_timestamps(self, mapping: dict) -> None:
        """Remap the timestamps of flags, outliers, and details using a mapping dictionary.
        
        Parameters
        ----------
        mapping : dict
            A dictionary where keys are original timestamps and values are the 
            new timestamps to map to.
        """
        self.flags.index = self.flags.index.map(lambda ts: mapping.get(ts, ts))
        self.outliers.index = self.outliers.index.map(lambda ts: mapping.get(ts, ts))
        self.details.index = self.details.index.map(lambda ts: mapping.get(ts, ts))
   
    def _flags_to_labels_map(self) -> dict:
        """Create mapping from QC flag values to display labels.
        
        Constructs a dictionary mapping internal flag values ('passed', 'flagged', etc.)
        to user-facing label strings defined in Settings. Flagged records use the
        check-specific label, while all other statuses use the 'goodrecord' label.
        
        Returns
        -------
        dict
            Mapping from flag strings to label strings.
        """
        label_mapping = {
            pass_cond: Settings.get('label_def.goodrecord.label'),
            flagged_cond: Settings.get(f'label_def.{self.checkname}.label'),
            unmet_cond: Settings.get('label_def.goodrecord.label'),
            saved_cond: Settings.get('label_def.goodrecord.label'),
            unchecked_cond: Settings.get('label_def.goodrecord.label')
        }
        return label_mapping
    def create_outliersdf(self) -> pd.DataFrame:
        """Create a DataFrame summarizing detected outliers.
        
        Constructs a DataFrame containing outlier values, their corresponding labels,
        and detailed information for each outlier timestamp. This format is compatible
        with the Dataset.outliersdf property.
        
        Returns
        -------
        pd.DataFrame
            DataFrame with datetime index and columns:
            - 'value': outlier values
            - 'label': human-readable QC check labels
            - 'details': descriptive information about each outlier
            Returns empty DataFrame with correct structure if no outliers exist.
        """
        if self.outliers.empty:
            # return empty dataframe
            return pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.DatetimeIndex([], name="datetime"),
            )
            
        labels = self.flags.loc[self.outliers.index].map(self._flags_to_labels_map())

        
        outliers_df = pd.DataFrame({
            'datetime': self.outliers.index,
            'value': self.outliers.values,
            'label': labels.values,
            'details': self.details.loc[self.outliers.index].values,
           
        })
        outliers_df.set_index('datetime', inplace=True)
        return outliers_df