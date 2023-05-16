#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:47:05 2023

@author: thoverga
"""

from numpy import nan

gaps_settings = {
    "gaps_finder": {
        "gapsize_n": 40
    },  # gaps defined as n times the highest frequency on IO.
}

gaps_info = {
    "gap": {
        "label_columnname": "is_gap",
        "outlier_flag": "gap",
        "negative_flag": "no gap",
        "numeric_flag": 12,
        "apply_on": "record",
    },
    "missing_timestamp": {
        "label_columnname": "is_missing_timestamp",
        "outlier_flag": "missing timestamp",
        "negative flag": "not missing",
        "numeric_flag": 13,
        "apply_on": "record",
    },
}


# =============================================================================
#  Gap filling settings
# =============================================================================


gaps_fill_settings = {
    "linear": {"method": "time", "max_consec_fill": 100},
    "model_debias": {
        "debias_period": {
            "prefered_leading_sample_duration_hours": 48,
            "prefered_trailing_sample_duration_hours": 48,
            "minimum_leading_sample_duration_hours": 24,
            "minimum_trailing_sample_duration_hours": 24,
        }
    },
    "automatic":{'max_interpolation_duration_str': '5H'}
}


gaps_fill_info = {
    "label_columnname": "final_label",
    "label": {"linear": "gap_interpolation", "model_debias": "gap_debiased_era5"},
    "numeric_flag": 21,
}

# =============================================================================
#  Missing obs filling settings
# =============================================================================
missing_obs_fill_settings={
    'linear': {'method': 'time'}

}

missing_obs_fill_info = {
    "label_columnname": "final_label",
    "label": {"linear": "missing_obs_interpolation"},
    "numeric_flag": 23,

}