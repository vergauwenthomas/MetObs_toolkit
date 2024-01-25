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
    },
    "buddy_check": {
        "temp": {
            "radius": 15000,  # Search radius in meter
            "num_min": 2,  # 	int		The minimum number of buddies a station can have
            "threshold": 1.5,  # 	σ	the variance threshold for flagging a station
            "max_elev_diff": 200,  # 	 m the maximum difference in elevation for a buddy (if negative will not check for heigh difference)
            "elev_gradient": -0.0065,  # 	linear elevation gradient with height
            "min_std": 1.0,  # If the standard deviation of values in a neighborhood are less than min_std, min_std will be used instead
        }
    },
}

titan_check_settings = {
    "titan_buddy_check": {
        "temp": {
            # 'radius': 5000, #	vec	m	Search radius
            # 'num_min': 5, #	int		The minimum number of buddies a station can have
            "radius": 50000,  # 	vec	m	Search radius
            "num_min": 2,  # 	int		The minimum number of buddies a station can have
            # 'threshold': 2.5, #	float	σ	the variance threshold for flagging a station
            "threshold": 1.5,  # 	float	σ	the variance threshold for flagging a station
            "max_elev_diff": 200,  # 	float	m	the maximum difference in elevation for a buddy (if negative will not check for heigh difference)
            "elev_gradient": -0.0065,  # 	float	ou/m	linear elevation gradient with height
            "min_std": 1.0,  # 	float		If the standard deviation of values in a neighborhood are less than min_std, min_std will be used instead
            "num_iterations": 1,  # int		The number of iterations to perform
        },
    },
    "titan_sct_resistant_check": {
        "temp": {
            "num_min_outer": 3,  # 	int	    Minimal points in outer circle
            "num_max_outer": 10,  # 	int	    Maximal points in outer circle
            "inner_radius": 20000,  #  int      Radius of inner circle
            "outer_radius": 50000,  #  int      Radius of outer circle
            "num_iterations": 10,  #   int   Number of iterations
            "num_min_prof": 5,  #   int Minimum number of observations to compute vertical profile
            "min_elev_diff": 100,  #   int    Minimal elevation difference
            "min_horizontal_scale": 250,  #  int  Minimal horizontal scale
            "max_horizontal_scale": 100000,  #  int  Maximal horizontal scale
            "kth_closest_obs_horizontal_scale": 2,  #  int  Number of closest observations to consider in the adaptive estimation of the horizontal decorrelation length
            "vertical_scale": 200,  #  int  The vertical scale
            "mina_deviation": 10,  # vec Minimum admissible value
            "maxa_deviation": 10,  # vec Maximum admissible value
            "minv_deviation": 1,  # vec Minimum valid value
            "maxv_deviation": 1,  # 	vec Maximum valid value
            "eps2": 0.5,  # Ratio of observation error variance to background variance
            "tpos": 5,  # vec Positive deviation allowed
            "tneg": 8,  # vec Negative deviation allowed
            "basic": True,  #  bool  Basic mode
            "debug": False,  #  bool  Debug mode
        }
    },
}

# how to map the numeric output of titan to a 'ok' or outlier label
titan_specific_labeler = {
    "titan_buddy_check": {"ok": [0], "outl": [1]},
    "titan_sct_resistant_check": {
        "ok": [0, -999, 11, 12],  # if obs not checked, or cannot be checked assume ok
        "outl": [1],
    },
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
    "buddy_check": {
        "outlier_flag": "buddy check outlier",
        "numeric_flag": 11,
        "apply_on": "obstype",
    },
    "titan_buddy_check": {
        "outlier_flag": "titan buddy check outlier",
        "numeric_flag": 9,
        "apply_on": "obstype",
    },
    "titan_sct_resistant_check": {
        "outlier_flag": "sct resistant check outlier",
        "numeric_flag": 10,
        "apply_on": "obstype",
    },
}
