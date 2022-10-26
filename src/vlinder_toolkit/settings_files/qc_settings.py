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
    "persistance": {'temp': {'max_valid_repetitions': 5}}
    
    }


outlier_values = {
    "duplicate_timestamp": nan, 
    "gross_value": nan,
    "persistance": nan    
    }



observation_labels={
    'ok': 'ok',
    'duplicated_timestamp': 'duplicated timestamp outlier',
    'gross_value': 'gross value outlier',
    'persistance': 'percistance outlier'
    }





