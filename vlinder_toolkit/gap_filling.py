#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:05:26 2023

@author: thoverga
"""

import pandas as pd




# =============================================================================
# Helpers
# =============================================================================


def expand_gabs_to_multiind(df, gapsdf, obstype='temp'):
    """
    Expand the gaps to timestamps, that are missing in the observations. 
    

    Parameters
    ----------
    df : dataset.df
        The observations where gaps are present as Nan values.
    gapsdf : dataset.gapsdf
        The dataframe with detaild information on the start and end of a gap.
    obstype : String, optional
        Observation type. The default is 'temp'.

    Returns
    -------
    expanded_gabsidx : pd.MultiIndex
        A Station-datetime multiindex of missing gap-records.

    """
    
    expanded_gabsidx = pd.MultiIndex(levels=[['name'],['datetime']],
                             codes=[[],[]],
                             names=[u'name', u'datetime'])
    
    
    for sta, row in gapsdf.iterrows():
    
        sta_df = df.xs(sta, level='name') #filter by name 
        posible_gaps_dt = sta_df[sta_df[obstype].isnull()].index #filter by missing observations
        gaps_dt = posible_gaps_dt[(posible_gaps_dt >= row['start_gap']) & #filter if the observations are within a gap
                                  (posible_gaps_dt <= row['end_gap'])]
        
        gaps_multiidx = pd.MultiIndex.from_arrays(arrays=[[sta]*len(gaps_dt),
                                                          gaps_dt],
                                                  names=[u'name', u'datetime'])
        
        expanded_gabsidx = expanded_gabsidx.append(gaps_multiidx)
    
    return expanded_gabsidx


# =============================================================================
# Gap fillers
# =============================================================================

def interp_gaps(obsdf, gapsdf, obstype='temp', method='time', max_consec_fill=1000):
    """
    Fill missing gap records using an interpolation of present observations. 
    The observations in the obsdf will be used to fill gaps, using the same timeresolution, 
    described in the gapsdf. 
    
    method 'time', will use datetime for interpolation.

    Parameters
    ----------
    df : dataset.df
        The observations where gaps are present as Nan values.
    gapsdf : dataset.gapsdf
        The dataframe with detaild information on the start and end of a gap.
    obstype : String, optional
        Observation type. The default is 'temp'.
    method : String, optional
        interpolation method (be awere that the datetime is the index during interpolation). The default is 'time'.
    max_consec_fill : Integer, optional
        Specify the maximum number of consecutive Nan's to fill. The default is 1000.

    Returns
    -------
    DataFrame with gap-filled values as a column and a name-datetime multiindex

    """
    #1. expand gapsdf to multiindex, at observation frequency
    expanded_gapsidx = expand_gabs_to_multiind(obsdf, gapsdf)
    
    #2. Merge gapsdf and obsdf and label the gabs
    gapsdf = pd.DataFrame(index=expanded_gapsidx).reset_index()
    gapsdf['isgap'] = True
    
    df = obsdf[obstype].reset_index()
    df = df.merge(gapsdf, how='left', on=['name', 'datetime'])
    df['isgap'] = df['isgap'].fillna(False)
    
    
    #3. subset only stations with gaps
    df = df[df['name'].isin(gapsdf['name'].unique())]
    
    #4. Drop nans, that are outside the gaps (from QC)
    
    #maybe not a good idea because now it can also be used to
    #fill in missing timestamps. 
    
    #in 6 the return is filtered by isgap, so only the gap filled values are returned.
    
    
    #5. interpolate gaps, station per stations
    interp_df = pd.DataFrame()
    for sta in gapsdf['name'].unique():
        subdf = df[df['name']==sta].set_index('datetime')
      
        subdf.interpolate(method='time',
                          limit=1000, # Maximum number of consecutive NaNs to fill. Must be greater than 0.
                          limit_area='inside', 
                          inplace=True)
        interp_df = pd.concat([interp_df, subdf])
        
    #6. Filter out the interpolated gabs
    interp_df = interp_df[interp_df['isgap']]
    interp_df = interp_df.drop(columns=['isgap'])
    interp_df.reset_index(inplace=True)
    interp_df.set_index(['name', 'datetime'], inplace=True)
    
    return interp_df

