#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  2 18:47:05 2023

@author: thoverga
"""

from numpy import nan


gaps_info = {
    "gap": {
        "label_columnname": "is_gap",
        "outlier_flag": "gap",
        "negative_flag": "no gap",
        # "numeric_flag": 12,
        "apply_on": "record",
    }
}


# =============================================================================
#  Gap filling settings
# =============================================================================


gaps_fill_info = {
    "label_columnname": "final_label",
    "label": {"linear": "gap_interpolation", "model_debias": "gap_debiased_era5"},
    # "numeric_flag": 21,
}
