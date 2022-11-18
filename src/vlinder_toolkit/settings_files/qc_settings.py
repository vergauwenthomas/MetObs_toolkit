#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: thoverga
"""
from numpy import nan



check_settings = {
    
    #checks on all observation types
    "duplicate_timestamp": {}, #No numeric settings
    
    
    #checks on specific observation types
    "gross_value": {'temp': {'min_value': -15.0,
                             'max_value': 39.0},
                    },
    "persistance": {'temp': {'max_valid_repetitions': 5}},
    
    "step": {'temp': {'max_value': 4}},
    
    "internal_consistency": {} #No numeric settings
    }


outlier_values = {
    "duplicate_timestamp": nan, 
    "gross_value": nan,
    "persistance": nan,   
    "step": nan,
    "internal_consistency": nan
    }



observation_labels={
    'ok': 'ok',
    'duplicated_timestamp': 'duplicated timestamp outlier',
    'gross_value': 'gross value outlier',
    'persistance': 'persistance outlier',
    'step': 'step outlier',
    'internal_consistency': 'internal consistency outlier'
    }





