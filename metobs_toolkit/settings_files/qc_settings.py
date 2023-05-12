#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: thoverga
"""

# Numeric settings on how the checks performs
check_settings = {
    # checks on all observation types
    "duplicated_timestamp": {
        "keep": False
    },  # No numeric settings (False: drop all duplicates)
    # "missing_timestamp": {},
    "persistance": {
        "temp": {
            "time_window_to_check": "1h",  # Use this format as example: "1h20min50s"
            "min_num_obs": 5,
        }
    },  # Minimum number of records in window to perform check
    "repetitions": {"temp": {"max_valid_repetitions": 5}},
    # checks on specific observation types
    "gross_value": {
        "temp": {"min_value": -15.0, "max_value": 39.0},
    },
    "window_variation": {
        "temp": {
            "max_increase_per_second": 8.0
            / 3600.0,  # == max 8° (increase) change in 3600 seconds (==one hour)
            "max_decrease_per_second": 10.0
            / 3600.0,  # == max 10° (decrease) change in 3600 seconds (==one hour)
            "time_window_to_check": "1h",  # Use this format as example: "1h20min50s"
            "min_window_members": 3,  # Minimum numer of records in window to perform check
        }
    },
    "step": {
        "temp": {
            "max_increase_per_second": 8.0 / 3600.0,
            "max_decrease_per_second": -10.0 / 3600.0,
        }
    }
}


# Information on the sequence of checks and if they are applied on all observations seperatly.

# Labels are converted to numeric to compute a final quality label.
# Numeric values are arbitrary for outliers (not for ok and not checked), but
# should be invertible.
# This is done for speeding up
checks_info = {
    # 1. --> on data import
    "duplicated_timestamp": {
        "outlier_flag": "duplicated timestamp outlier",
        "numeric_flag": 1,
        "apply_on": "record",
    },
    # 2. --> on data import
    "invalid_input": {
        "outlier_flag": "invalid input",
        "numeric_flag": 2,
        "apply_on": "obstype",
    },
    # 3. --> on observed values
    "gross_value": {
        "outlier_flag": "gross value outlier",
        "numeric_flag": 4,
        "apply_on": "obstype",
    },
    # 4(A). --> on observed values
    "persistance": {
        "outlier_flag": "persistance outlier",
        "numeric_flag": 5,
        "apply_on": "obstype",
    },
    # 4(B). --> on observed values
    "repetitions": {
        "outlier_flag": "repetitions outlier",
        "numeric_flag": 6,
        "apply_on": "obstype",
    },
    # 5. --> on observed values
    "step": {
        "outlier_flag": "in step outlier group",
        "numeric_flag": 7,
        "apply_on": "obstype",
    },
    "window_variation": {
        "outlier_flag": "in window variation outlier group",
        "numeric_flag": 8,
        "apply_on": "obstype",
    },
}
