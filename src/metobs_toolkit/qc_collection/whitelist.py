import pandas as pd


class SensorWhiteSet:
    def __init__(self, white_timestamps=[], all_timestamps=False):
        
        if all_timestamps:
            assert len(white_timestamps) == 0, "If all_timestamps is True, white_timestamps must be empty."
        
        self.white_timestamps = white_timestamps
        self.all_timestamps = all_timestamps
    
    def has_whites(self) -> bool:
        if self.all_timestamps:
            return True
        if len(self.white_timestamps) > 0:
            return True
        return False
    
    def all_timestamps_white(self) -> bool:
        return self.all_timestamps
    
    def get_white_timestamps(self) -> pd.DatetimeIndex:
        return pd.DatetimeIndex(data=self.white_timestamps,
                                name='datetime')
        
    def catch_white_records(self, outliers_idx: pd.DatetimeIndex) -> pd.DatetimeIndex:
        """Remove white record timestamps from outliers index.
        
        Parameters
        ----------
        outliers_idx : pd.DatetimeIndex
            Index of outlier timestamps to filter
            
        Returns
        -------
        pd.DatetimeIndex
            Filtered outliers index with white records removed
        """
        if self.has_whites():

            if self.all_timestamps_white():
                #all timestamps are white, return empty index
                outliers = pd.DatetimeIndex([], name='datetime')
            else:
                #Get the white timestamps
                white_records = self.get_white_timestamps()
                #Remove white records from outliers
                outliers = outliers_idx.difference(white_records)
        else:
            #no whites
            outliers = outliers_idx
        
        return outliers


class WhiteSet:
    def __init__(self, white_records: pd.Index=pd.Index([])):
        self.white_records = white_records
        
        #Self check
        self.test_white_records()
    
    def test_white_records(self) -> None:
        
        """Validate the structure and content of white_records index.
        
        This function checks that white_records has a valid structure for use in QC methods.
        The white_records index must contain at least one of the expected level names
        ('name', 'obstype', or 'datetime') and may not contain any unexpected levels.
        

        Raises
        ------
        ValueError
            If white_records does not contain at least one of 'name', 'obstype', or 
            'datetime' as index level names.
        ValueError
            If white_records contains index levels with names other than 'name', 
            'obstype', or 'datetime'.
            
        Notes
        -----
        
        This function is used internally to validate white_records before they are 
        processed by QC methods.
        """
        if self.is_empty():
            return
        
        if not any([idxname in self.white_records.names for idxname in ['name','obstype', 'datetime']]):
            raise ValueError("white_records must contain at least one of the following index levels: 'name', 'obstype', 'datetime'")
        if not all([idxname in ['name','obstype', 'datetime'] for idxname in self.white_records.names]):
            raise ValueError("white_records contains unexpected index levels. Only 'name', 'obstype', and 'datetime' are allowed.")
    
    def is_empty(self) -> bool:
        """Check if white_records is empty.
        
        Returns
        -------
        bool
            True if white_records is empty, False otherwise.
        """
        return self.white_records.empty

    def create_sensorwhitelist(self, trg_station: str, trg_obstype: str) -> 'SensorWhiteSet':
        """
        Create a SensorWhiteList object for a specific station and obstype
        based on the white_records MultiIndex.
        """
        if self.is_empty():
            return SensorWhiteSet(white_timestamps=[],
                                   all_timestamps=False)
        
        # Filter white_records for the target station and obstype
        trg_whitelist = self.white_records

        if 'name' in trg_whitelist.names:
            if trg_station  in trg_whitelist.get_level_values('name'):
                trg_whitelist = trg_whitelist[
                    (trg_whitelist.get_level_values('name') == trg_station)]
            else: 
                #name is specified, but no matches in whitelist
                return SensorWhiteSet(white_timestamps=[],
                                       all_timestamps=False)
            
        #filter on obstype if present
        if 'obstype' in trg_whitelist.names:
            if trg_obstype in trg_whitelist.get_level_values('obstype'):
                trg_whitelist = trg_whitelist[
                    (trg_whitelist.get_level_values('obstype') == trg_obstype)]
            else:
                #obstype is specified, but no matches in whitelist
                return SensorWhiteSet(white_timestamps=[],
                                       all_timestamps=False)
        
        if 'datetime' in trg_whitelist.names:
            #Get the white datetimes
            white_datetimes = pd.DatetimeIndex(trg_whitelist.get_level_values('datetime').unique())
            return SensorWhiteSet(white_timestamps=white_datetimes, all_timestamps=False)
        else:
            #if no datetime level is set, and name and/or obstype match, all timestamps are white
            return SensorWhiteSet(white_timestamps=[], 
                                   all_timestamps=True)