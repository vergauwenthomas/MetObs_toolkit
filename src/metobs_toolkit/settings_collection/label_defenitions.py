#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Label color and name defenitions for the metobs_toolkit.

@author: thoverga
"""


# =============================================================================
# Defenition of labels
# =============================================================================
# All knonw labels should be defined here:
label_def = {
    # --- Good records --------
    "goodrecord": {"label": "ok", "color": "#07f72b", "numeric_val": 0},
    "uncheckedrecord": {"label": "not checked", "color": "#f7cf07", "numeric_val": 1},
    # ------ QC labels ---------------
    "duplicated_timestamp": {
        "label": "duplicated timestamp outlier",
        "color": "#a32a1f",
        "numeric_val": 2,
    },
    # "invalid_input": {"label": "invalid input", "color": "#900357"},
    "gross_value": {
        "label": "gross value outlier",
        "color": "#f1ff2b",
        "numeric_val": 3,
    },
    "persistence": {
        "label": "persistence outlier",
        "color": "#f0051c",
        "numeric_val": 4,
    },
    "repetitions": {
        "label": "repetitions outlier",
        "color": "#056ff0",
        "numeric_val": 5,
    },
    "step": {"label": "in step outlier group", "color": "#05d4f0", "numeric_val": 6},
    "window_variation": {
        "label": "in window variation outlier group",
        "color": "#05f0c9",
        "numeric_val": 7,
    },
    "buddy_check": {
        "label": "buddy check outlier",
        "color": "#8300c4",
        "numeric_val": 8,
    },
    "buddy_check_with_LCZ_safety_net": {
        "label": "buddy check (with LCZ-safety net) outlier",
        "color": "#8300c4",
        "numeric_val": 9,
    },
    # aggregated
    "outlier": {
        "label": "QC outlier",
        "color": "#f20000",
        "agg_def": True,
        "numeric_val": 10,
    },  # specify it is an aggregated label
    # ----- Gap ----------
    "regular_gap": {"label": "gap", "color": "#f00592", "numeric_val": 11},
    # ----- Interpolation labels -----
    "interpolated_gap": {
        "label": "interpolation",
        "color": "#d406c6",
        "numeric_val": 12,
    },
    "failed_interpolation_gap": {
        "label": "failed interpolation",
        "color": "#d406c6",
        "numeric_val": 13,
    },
    # ----- raw model gapfill -----
    "raw_modeldata_fill": {
        "label": "raw modeldata fill",
        "color": "#6e1868",
        "numeric_val": 14,
    },
    "failed_raw_modeldata_fill": {
        "label": "failed raw modeldata fill",
        "color": "#6e1868",
        "numeric_val": 15,
    },
    # ----- debias model gapfill -----
    "debias_modeldata_fill": {
        "label": "debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 16,
    },
    "failed_debias_modeldata_fill": {
        "label": "failed debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 17,
    },
    "diurnal_debias_modeldata_fill": {
        "label": "diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 18,
    },
    "failed_diurnal_debias_modeldata_fill": {
        "label": "failed diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 19,
    },
    "weighted_diurnal_debias_modeldata_fill": {
        "label": "Weighted diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 20,
    },
    "failed_weighted_diurnal_debias_modeldata_fill": {
        "label": "failed Weighted diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 21,
    },
}

label_to_color_map = {group["label"]: group["color"] for group in label_def.values()}
label_to_numeric_map = {
    group["label"]: group["numeric_val"] for group in label_def.values()
}
# make lists per label-theme, so if a new thematic label is added,
# functionallity through the full toolkit is done by adding it to the
# defenition and the labelgroup

gapfill_label_group = [
    "interpolated_gap",
    "raw_modeldata_fill",
    "debias_modeldata_fill",
    "diurnal_debias_modeldata_fill",
    "weighted_diurnal_debias_modeldata_fill",
]
failed_gapfill_label_group = [
    "failed_interpolation_gap",
    "failed_raw_modeldata_fill",
    "failed_debias_modeldata_fill",
    "failed_diurnal_debias_modeldata_fill",
    "failed_weighted_diurnal_debias_modeldata_fill",
]

qc_label_group = [
    "duplicated_timestamp",
    # "invalid_input",
    "gross_value",
    "persistence",
    "repetitions",
    "step",
    "window_variation",
    "buddy_check",
    "buddy_check_with_LCZ_safety_net",
]
