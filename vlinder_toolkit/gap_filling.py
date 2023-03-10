#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 28 17:05:26 2023

@author: thoverga
"""
import numpy as np
import pandas as pd





# =============================================================================
# Gap fillers
# =============================================================================

def interpolate_gap(gap, obsdf, outliersdf, dataset_res, obstype,
                    method, max_consec_fill):
    
    
    #1 get trailing and leading + exploded index
    gap.update_leading_trailing_obs(obsdf, outliersdf)
    gap.update_gaps_indx_in_obs_space(obsdf, outliersdf, dataset_res)
    
    
    # initiate return value when no interpolation can be performed
    empty_interp = pd.Series(data=np.nan, index=gap.exp_gap_idx.droplevel('name'))
    empty_interp.name=obstype
    
    
    # 2 check if there is a trailing and leading gap
    if (gap.startgap == gap.leading_timestamp):
        print(f'No leading timestamp found for gap ({gap.name}): {gap.startgap} --> {gap.endgap}')
        return  empty_interp
    
    
    if (gap.endgap == gap.trailing_timestamp):
        print(f'No trailing timestamp found for gap ({gap.name}): {gap.startgap} --> {gap.endgap}')
        return  empty_interp
    
    
    
    
    # 3 check both leading and trailing are in obs, and look for alternative leading/trailing if the original is an outlier.
    sta_obs = obsdf.xs(gap.name, level='name')
    
    # leading
    if gap.leading_timestamp in sta_obs:
        # leading found in obs
        leading_dt = gap.leading_timestamp
        leading_val = sta_obs.loc[gap.leading_timestamp, obstype]
    else: 
        # look for last observation before leading timestamp
        delta_dt = (gap.leading_timestamp - sta_obs.index[sta_obs.index < gap.leading_timestamp])
        if delta_dt.empty:
            print(f'No cadidate for leading {obstype} observation found for {gap.name} with gap: {gap.startgap} --> {gap.endgap}')
            return  empty_interp
            
        leading_dt = gap.leading_timestamp - delta_dt.min()
        leading_val = sta_obs.loc[leading_dt, obstype]
    
    #trailing
    if gap.trailing_timestamp in sta_obs:
        # leading found in obs
        trailing_dt = gap.trailing_timestamp
        trailing_val = sta_obs.loc[gap.trailing_timestamp, obstype]
    else: 
        # look for last observation before leading timestamp
        delta_dt = (sta_obs.index[sta_obs.index > gap.trailing_timestamp] - gap.trailing_timestamp)
        if delta_dt.empty:
            print(f'No cadidate for trailing {obstype} observation found for {gap.name} with gap: {gap.startgap} --> {gap.endgap}')
            return  empty_interp
        # TODO: settings restrictions on maximum delta_dt ??
        trailing_dt = gap.trailing_timestamp + delta_dt.min()
        trailing_val = sta_obs.loc[trailing_dt, obstype]
    
    # Make interpolation series
    gaps_series = pd.Series(data=np.nan, index=gap.exp_gap_idx.droplevel('name'))
    gaps_series = pd.concat([gaps_series,
                             pd.Series(index=[leading_dt, trailing_dt],
                                               data=[leading_val, trailing_val])])
    gaps_series = gaps_series.sort_index()
    
    
    
    # Interpolate series
    gaps_series.interpolate(method=method,
                            limit=max_consec_fill, # Maximum number of consecutive NaNs to fill. Must be greater than 0.
                            limit_area='inside', 
                            inplace=True)
    
    # Subset only gap indixes
    gaps_series = gaps_series[gap.exp_gap_idx.droplevel('name')]
    gaps_series.name=obstype
    
    return gaps_series

