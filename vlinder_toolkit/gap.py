#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

Specific classes are created because gaps and missing obs can have out-of-sync 
datetimes wrt dataset.df. 

Created on Fri Mar  3 09:15:56 2023

@author: thoverga
"""

import pandas as pd
import numpy as np
import logging
from datetime import datetime, timedelta
import math

from vlinder_toolkit.settings import Settings

logger = logging.getLogger(__name__)




# =============================================================================
# Missing observation

# a missing observation is a missing timestamp
# =============================================================================


    

class Missingob_collection:
    def __init__(self, missing_obs_series):
        self.series = missing_obs_series
        
        
    def get_station_missingobs(self, name):
        if name in self.series.index:
            return Missingob_collection(self.series.loc[[name]])
        else:
            # return empty collection
            series = pd.Series(data=[], name='datetime', dtype=object)
            series.index.name = 'name'
            return Missingob_collection(series)
        
        
    def remove_missing_from_obs(self, obsdf):
        
        #Normally there are no missing records in the obsdf
        missing_multiidx = pd.MultiIndex.from_arrays(arrays=[self.series.index.to_list(),
                                                             self.series.to_list()],
                                                  names=[u'name', u'datetime'])
        
        obsdf = obsdf.drop(index=missing_multiidx, errors='ignore')
        
        return obsdf
    
    
    def get_missing_indx_in_obs_space(self, obsdf, resolutionseries):
        """
        Function to found which missing timestamps are expected in the observation space.
        Because of time coarsening not all missing timestamps are expected in observation space.
        
        This function handles each station seperatly because stations can have differnent resolution/timerange.
        
        

        Parameters
        ----------
        obsdf : pandas.DataFrame()
            Dataset.df.
        resolutionseries : pandas.Series() or Timedelta
            Dataset.metadf['dataset_resolution'].

        Returns
        -------
        missing_obsspace : pandas.MultiIndex
            The multiindex (name - datetime) is returned with the missing timestamps that are expexted in the observation space.

        """
                 
        missing_obsspace = pd.MultiIndex(levels=[['name'],['datetime']],
                                 codes=[[],[]],
                                 names=[u'name', u'datetime'])
        
        # per stationtion because stations can have different resolutions/timerange
        for sta in self.series.index.unique():
            
            # Get missing observations in IO space
            sta_missing = self.series.loc[sta]
            if not isinstance(sta_missing, type(pd.Series(dtype=object))):
                sta_missing = pd.Series(data=[sta_missing], index=[sta],
                                        dtype=object)
             
            
            # Get start, end and frequency of the observation in obs space
            startdt = obsdf.xs(sta, level='name').index.min()
            enddt = obsdf.xs(sta, level='name').index.max()
            obs_freq = resolutionseries.loc[sta]
            
            # Make datetimerange
            obsrange = pd.date_range(start=startdt,
                                     end=enddt,
                                     freq = obs_freq,
                                     inclusive="both")
            
            # # Look which missing timestamps appears obsspace
            sta_missing =sta_missing[sta_missing.isin(obsrange)]

            
            #Convert to multiindex
            if sta_missing.empty:
                continue
            
            sta_missing_idx =  pd.MultiIndex.from_arrays(arrays=[[sta]*len(sta_missing),
                                                                 sta_missing.to_numpy()],
                                                      names=[u'name', u'datetime'])
            
            missing_obsspace = missing_obsspace.append(sta_missing_idx)
         
            
        return missing_obsspace
         
# =============================================================================
# Gap class

# a gap is a sequence of repeting missing obs
# =============================================================================

class _Gap:
    # do not use this class outside this module
    def __init__(self, name, startdt, enddt):
        
        self.name = name
        self.startgap = startdt
        self.endgap = enddt




class Gap_collection:
    def __init__(self, gapsdf):
        self.list = [_Gap(sta, row['start_gap'], row['end_gap']) for sta, row in gapsdf.iterrows()]
        self.df = gapsdf
        
        
    
    
    def get_station_gaps(self, name):
        if name in self.df.index:
            return Gap_collection(self.df.loc[name])
        else:
            return Gap_collection(pd.DataFrame())
        
    
    
    def remove_gaps_from_obs(self, obsdf):
        
        #Create index for gaps records in the obsdf
        expanded_gabsidx = pd.MultiIndex(levels=[['name'],['datetime']],
                                 codes=[[],[]],
                                 names=[u'name', u'datetime'])
        
        for sta, row in self.df.iterrows():
        
            sta_records = obsdf.xs(sta, level='name').index #filter by name 
           
            gaps_dt = sta_records[(sta_records >= row['start_gap']) & #filter if the observations are within a gap
                                  (sta_records <= row['end_gap'])]
            
            gaps_multiidx = pd.MultiIndex.from_arrays(arrays=[[sta]*len(gaps_dt),
                                                              gaps_dt],
                                                      names=[u'name', u'datetime'])
            
            expanded_gabsidx = expanded_gabsidx.append(gaps_multiidx)
            
        #remove gaps idx from the obsdf
        obsdf = obsdf.drop(index=expanded_gabsidx)
        return obsdf
    
    
    def get_gaps_indx_in_obs_space(self, obsdf, outliersdf, resolutionseries):
        """
        test test document

        Parameters
        ----------
        obsdf : TYPE
            DESCRIPTION.
        outliersdf : TYPE
            DESCRIPTION.
        resolutionseries : TYPE
            DESCRIPTION.

        Returns
        -------
        expanded_gabsidx_obsspace : TYPE
            DESCRIPTION.

        """
            
        self._add_leading_and_trailing_obs(obsdf, outliersdf)
        
        
        expanded_gabsidx_obsspace = pd.MultiIndex(levels=[['name'],['datetime']],
                                 codes=[[],[]],
                                 names=[u'name', u'datetime'])
        
        
        for sta, row in self.df.iterrows():
            obs_freq = resolutionseries.loc[sta]
            
            gaprange = pd.date_range(start=row['leading_ob'],
                                     end=row['trailing_ob'],
                                     freq = obs_freq,
                                     inclusive="neither")
            
            gaps_multiidx = pd.MultiIndex.from_arrays(arrays=[[sta]*len(gaprange),
                                                              gaprange],
                                                      names=[u'name', u'datetime'])
            
            expanded_gabsidx_obsspace = expanded_gabsidx_obsspace.append(gaps_multiidx)
            
    
        return expanded_gabsidx_obsspace
    
    
    def _add_leading_and_trailing_obs(self, obsdf, outliersdf):
        """
        Add the leading (last obs before gap) and trailing (first obs after gap)
        as extra columns to the self.df. 
        
        The obsdf and outliersdf are both used to scan for the leading and trailing obs.

        Parameters
        ----------
        obsdf : pandas.DataFrame
            Dataset.df
        outliersdf : pandas.DataFrame
            Dataset.outliersdf

        Returns
        -------
        None.

        """
        
        leading_gap_timestamp = pd.Series(dtype=object)
        trailing_gap_timestamp = pd.Series(dtype=object)

        for sta, row in self.df.iterrows():
            # combine timestamps of observations and outliers
            sta_obs = obsdf.xs(sta, level='name').index
            sta_outl = outliersdf.xs(sta, level='name').index
            sta_comb = sta_obs.append(sta_outl)
            
            
            #find minimium timediff before
            before_diff = _find_closes_occuring_date(refdt = row['start_gap'],
                                                    series_of_dt = sta_comb,
                                                    where='before')
            
            #if no timestamps are before gap, assume gap at the start of the observations
            if math.isnan(before_diff):
                before_diff = 0.0
               
            # find minimum timediff after gap
            after_diff = _find_closes_occuring_date(refdt = row['end_gap'],
                                                    series_of_dt = sta_comb,
                                                    where='after')
            #if no timestamps are after gap, assume gap at the end of the observations
            if math.isnan(after_diff):
                after_diff = 0.0
            
            
            # get before and after timestamps
            last_ob_before = row['start_gap'] - timedelta(seconds=before_diff)
            first_ob_after = row['end_gap'] + timedelta(seconds=after_diff)

            #fill the series
            leading_gap_timestamp = pd.concat([leading_gap_timestamp, 
                                               pd.Series(index=[sta],
                                                         data=[last_ob_before])])
            
            trailing_gap_timestamp = pd.concat([trailing_gap_timestamp, 
                                               pd.Series(index=[sta],
                                                         data=[first_ob_after])])


        #add to gapsdf
        self.df['leading_ob'] = leading_gap_timestamp
        self.df['trailing_ob'] = trailing_gap_timestamp


# =============================================================================
# Find gaps and missing values
# =============================================================================
 
def get_freqency_series(df):
    freqs = {}
    for station in df.index.get_level_values(level='name').unique():
        timestamps = df.xs(station, level='name').index
        freqs[station] = get_likely_frequency(timestamps)
    return pd.Series(data=freqs)


def get_likely_frequency(timestamps):
    assume_freq = abs(timestamps.to_series().diff().value_counts().index).sort_values(ascending=True)[0]
    
    if assume_freq == pd.to_timedelta(0): #highly likely due to a duplicated record
        # select the second highest frequency
        assume_freq = abs(timestamps.to_series().diff().value_counts().index).sort_values(ascending=True)[1]
    
    return assume_freq


def _find_closes_occuring_date(refdt, series_of_dt, where='before'):
    if where=='before':
        diff = refdt - (series_of_dt[series_of_dt < refdt])
    elif where=='after':
        diff =(series_of_dt[series_of_dt > refdt]) - refdt 
        
        
    if diff.empty:
        #no occurences before of after

        return np.nan
    else:
        
        return min(diff).total_seconds()


def missing_timestamp_and_gap_check(df):
    #TODO update docstring
    """

    V3
    Looking for missing timestaps by assuming an observation frequency. The assumed frequency is the highest occuring frequency PER STATION.
    If missing observations are detected, they can be catogirized as a missing timestamp or as gap.
    
    A gap is define as a sequence of missing values with more than N repetitive missing values. N is define in the QC settings.
    


    Parameters
    ----------
    df : pandas.DataFrame
        The observations dataframe of the dataset object (Dataset.df)

    Returns
    -------

    df : pandas.DataFrame()
        The observations dataframe.
    outlier_df : pandas.DataFrame()
        The dataframe containing the missing timestamps (not gaps) with the outlier label.
    gap_df : pandas.Dataframe()
        The dataframe containing the start and end date of a specific gap.
        
    """     
    
    gaps_settings = Settings.gaps_settings
    
    gap_df = pd.DataFrame()
    gap_indices = []
    missing_timestamp_series = pd.Series(dtype=object)
    station_freqs = {}
    
    #missing timestamp per station (because some stations can have other frequencies!)

    stationnames = df.index.get_level_values(level='name').unique()
    for station in stationnames:
        
        #find missing timestamps
        timestamps = df.xs(station, level='name').index
        likely_freq = get_likely_frequency(timestamps)
     
        assert likely_freq.seconds > 0, f'The frequency is not positive!' 
        
        station_freqs[station] = likely_freq
        
        missing_datetimeseries = pd.date_range(start = timestamps.min(),
                                                end = timestamps.max(),
                                                freq=likely_freq).difference(timestamps).to_series().diff()
        
        
        if missing_datetimeseries.empty:
            continue
        
        #Check for gaps
        gap_defenition = ((missing_datetimeseries != likely_freq)).cumsum()
        consec_missing_groups = missing_datetimeseries.groupby(gap_defenition)
        group_sizes = consec_missing_groups.size()
        
        gap_groups = group_sizes[group_sizes > gaps_settings['gaps_finder']['gapsize_n']]
        
        #iterate over the gabs and fill the gapsdf
        for gap_idx in gap_groups.index:
            
            #fill the gaps df
            datetime_of_gap_records = consec_missing_groups.get_group(gap_idx).index
            gap_df = pd.concat([gap_df,
                                pd.DataFrame(data=[[datetime_of_gap_records.min(),
                                                  datetime_of_gap_records.max()]],
                                            index=[station],
                                            columns=['start_gap', 'end_gap'])])
            
            logger.debug(f'Data gap from {datetime_of_gap_records.min()} --> {datetime_of_gap_records.max()} found for {station}.')
            gap_indices.extend(list(zip([station]*datetime_of_gap_records.shape[0],
                                        datetime_of_gap_records)))
        
        # combine the missing timestams values
        missing_timestamp_groups = group_sizes[group_sizes <= gaps_settings['gaps_finder']['gapsize_n']]
        for missing_idx in missing_timestamp_groups.index:
            datetime_of_missing_records=consec_missing_groups.get_group(missing_idx).index.to_list()
            
            
            missing_timestamp_series = pd.concat([missing_timestamp_series,
                                                  pd.Series(index=[station] * len(datetime_of_missing_records),
                                                            data=datetime_of_missing_records)])
            
    
    
   
    
    

   
    df = df.sort_index()
    
   
    return df, missing_timestamp_series, gap_df