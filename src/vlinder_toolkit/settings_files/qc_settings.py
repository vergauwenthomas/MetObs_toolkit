#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: thoverga
"""
from numpy import nan




check_settings = {
    "duplicate_timestamp": {'temp': {}}, #No numeric settings
    "gross_value": {'temp': {'min_value': -10.0,
                             'max_value': 20.0},
                    },
    "persistance": {'temp': {'max_valid_repetitions': 5}}
    
    }


outlier_values = {
    "duplicate_timestamp": 'drop', 
    "gross_value": nan,
    "persistance": nan    
    }



observation_labels={
    'ok': 'ok',
    'gross_value': 'gross value outlier',
    'presistane': 'percistance outlier'
    }





