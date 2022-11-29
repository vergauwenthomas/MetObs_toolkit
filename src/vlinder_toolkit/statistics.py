# -*- coding: utf-8 -*-

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  25 13:44:54 2022

@author: thoverga
"""


import pandas as pd
import numpy as np
import logging

logger = logging.getLogger(__name__)




def get_qc_effectiveness_stats(df, obstype, observation_types, qc_labels):

    #1. Get list of relevant columns
    #extract label columns
    qc_labels_columns = [col for col in df.columns if not col in observation_types]
    #Extra savety
    qc_labels_columns = [col for col in qc_labels_columns if col.endswith('_label')]

    # for obstype in checked_obstypes:
    specific_columns = [col for col in qc_labels_columns if col.startswith(obstype+'_')] #qc applied on all obstypes not present! added later
    
    
    #2. Create mapping dict for each check to its possible labels 
    
    #make label and column mappers
    obs_labels_mappers = {}
    for col in specific_columns:
        checkname = col.replace(obstype+'_', '').replace('_label', '')
        obs_labels_mappers[col] = {'ok': qc_labels['ok'],
                                   'not checked': 'not checked',
                                   'outlier': qc_labels[checkname],
                                   'checkname': checkname}
    
    #add qc labels that are applicable on all obstypes
    if 'missing_timestamp_label' in qc_labels_columns:
        obs_labels_mappers['missing_timestamp_label'] = {
                        'ok': qc_labels['ok'],
                        'not checked': 'not checked',
                        'outlier': qc_labels['missing_timestamp'],
                        'checkname': 'missing_timestamp'}
        
        
    #3. Subset the dataframe and aggregate. Convert output to pandas.df
    df_qc = df[obs_labels_mappers.keys()]
    #TODO: maybe sort the keys first??
    
    #make counts for full dataset
    counts = df_qc.apply(lambda x: x.value_counts()).fillna(0)
    qc_countings_dict = {}
    for column, check_info  in obs_labels_mappers.items():
        try:
            ok_count = counts.loc[check_info['ok'], column]
        except KeyError:
            ok_count = 0    
        try:
            not_checked_count = counts.loc[check_info['not checked'], column]
        except KeyError:
            not_checked_count = 0    
        try:
            outlier_count = counts.loc[check_info['outlier'], column]
        except KeyError:
            outlier_count = 0    
        
        qc_countings_dict[check_info['checkname']] = {
            'ok': ok_count,
            'not checked': not_checked_count,
            'outlier': outlier_count
            }
    
    
    #Convert to df and make pieplot
    qc_counts_df = pd.DataFrame().from_dict(qc_countings_dict)
    
    # 4. Convert to persentages
    qc_percentage_df = qc_counts_df.div(qc_counts_df.sum(axis=0)) * 100.
    return qc_percentage_df


