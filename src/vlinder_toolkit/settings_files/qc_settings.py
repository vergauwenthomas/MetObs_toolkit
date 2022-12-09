#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: thoverga
"""
from numpy import nan


#Numeric settings on how the checks performs
check_settings = {
    
    "gaps_finder": {'gapsize_n': 40}, #gaps defined as n times the highest frequency on IO. 
    
    #checks on all observation types
    "duplicated_timestamp": {'keep': False}, #No numeric settings (False: drop all duplicates)

    "missing_timestamp": {},
    
    "persistance": {'temp': {'time_window_of_assumed_change': 5400,#in seconds
                             'minimum_numer': 5}}, #Minimum numer of records in window to perform check
    
    "repetitions": {'temp': {'max_valid_repetitions': 5}},
    
    #checks on specific observation types
    "gross_value": {'temp': {'min_value': -15.0,
                             'max_value': 39.0},
                    },
   
    
    "step": {'temp': {'max_increase_per_second': 8.0/3600.0, #== max 8° change in one hour
                      'max_decrease_per_second': 10.0/3600.0,
                      'min_window_members': 3,
                      'max_window_members': 5}
             },
    
    "internal_consistency": {'temp': {'b': 18.678,
                             'c': 257.14, 'd': 234.5}} 
    }



#Information on the sequence of checks and if they are applied on all observations seperatly.

#Labels are converted to numeric to compute a final quality label.
# Numeric values are arbitrary for outliers (not for ok and not checked), but
# should be invertible.
#This is done for speeding up 
checks_info={
    # 1. --> on data import
    'duplicated_timestamp':{'label_columnname': 'duplicated_timestamp_label',
                            'outlier_flag': 'duplicated timestamp outlier',
                            'numeric_flag': 1,
                            'apply_on': 'record'
                            },
    # 2(A). --> on data import
    'gaps_finder':{'label_columnname': 'gap_timestamp_label',
                            'outlier_flag': 'missing timestamp (gap)',
                            'numeric_flag': 2,
                            'apply_on': 'record'
                            },
    # 2(B). --> on data import
    'missing_timestamp':{'label_columnname': 'missing_timestamp_label',
                            'outlier_flag': 'missing timestamp',
                            'numeric_flag': 3,
                            'apply_on': 'record'
                            },
    # 3. --> on observed values
    'gross_value':{'label_columnname': 'gross_value_label', #Obstype_ is prefix
                            'outlier_flag': 'gross value outlier',
                            'numeric_flag': 4,
                            'apply_on': 'obstype'
                            },
    # 4(A). --> on observed values
    'persistance':{'label_columnname': 'persistance_label', #Obstype_ is prefix
                            'outlier_flag': 'persistance outlier',
                            'numeric_flag': 5,
                            'apply_on': 'obstype'
                            },
    # 4(B). --> on observed values
    'repetitions':{'label_columnname': 'repetitions_label', #Obstype_ is prefix
                            'outlier_flag': 'repetitions outlier',
                            'numeric_flag': 6,
                            'apply_on': 'obstype'
                            },
    # 5. --> on observed values
    'step':{'label_columnname': 'step_label', #Obstype_ is prefix
                            'outlier_flag': 'in step outlier group',
                            'numeric_flag': 7,
                            'apply_on': 'obstype'
                            },
    
    }





# observation_labels={
#     'missing_timestamp': 'missing timestamp',
#     'gap_timestamp': 'missing timestamp (gap)',
#     'duplicated_timestamp': 'duplicated timestamp outlier',
#     'gross_value': 'gross value outlier',
#     'persistance': 'persistance outlier',
#     'step': 'step outlier',
#     'internal_consistency': 'internal consistency outlier'

    # }


#Labels are converted to numeric to compute a final quality label.
# Numeric values are arbitrary for outliers (not for ok and not checked), but
# should be invertible.
#This is done for speeding up 

# numeric_label_mapper={
#     'ok': 0,
#     'not checked': nan,
#     'missing timestamp': 1,
#     'missing timestamp (gap)': 2,
#     'duplicated timestamp outlier': 3,
#     'gross value outlier': 4,
#     'persistance outlier': 5,
#     'step outlier':6,
#     }


