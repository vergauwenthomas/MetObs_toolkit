#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: thoverga
"""
from numpy import nan



check_settings = {
    
    
    
    #checks on all observation types
    "duplicate_timestamp": {'keep': False}, #No numeric settings
    "missing_timestamp": {},
    
    #checks on specific observation types
    "gross_value": {'temp': {'min_value': -15.0,
                             'max_value': 39.0},
                    },
    "persistance": {'temp': {'max_valid_repetitions': 5}},
    
    "step": {'temp': {'max_value': 4}},
    
    "internal_consistency": {'temp': {'b': 18.678,
                             'c': 257.14, 'd': 234.5}} 
    }


outlier_values = {
    "missing_timestamp": nan,
    "duplicate_timestamp": nan, 
    "gross_value": nan,
    "persistance": nan,   
    "step": nan,
    "internal_consistency": nan
    }



observation_labels={
    'ok': 'ok',
    'missing_timestamp': 'missing timestamp',
    'duplicated_timestamp': 'duplicated timestamp outlier',
    'gross_value': 'gross value outlier',
    'persistance': 'persistance outlier',
    'step': 'step outlier',
    'internal_consistency': 'internal consistency outlier'

    }


#Labels are converted to numeric to compute a final quality label.
# Numeric values are arbitrary for outliers (not for ok and not checked), but
# should be invertible.
#This is done for speeding up 
numeric_label_mapper={
    'ok': 0,
    'not checked': nan,
    'missing timestamp': 1,
    'duplicated timestamp outlier': 2,
    'gross value outlier': 3,
    'persistance outlier': 4 
    
    }


