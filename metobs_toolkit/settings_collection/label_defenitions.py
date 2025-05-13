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
    "goodrecord": {"label": "ok", "color": "#07f72b"},
    "uncheckedrecord": {"label": "not checked", "color": "#f7cf07"},
    # ------ QC labels ---------------
    "duplicated_timestamp": {
        "label": "duplicated timestamp outlier",
        "color": "#a32a1f",
    },
    # "invalid_input": {"label": "invalid input", "color": "#900357"},
    "gross_value": {"label": "gross value outlier", "color": "#f1ff2b"},
    "persistence": {"label": "persistence outlier", "color": "#f0051c"},
    "repetitions": {"label": "repetitions outlier", "color": "#056ff0"},
    "step": {"label": "in step outlier group", "color": "#05d4f0"},
    "window_variation": {
        "label": "in window variation outlier group",
        "color": "#05f0c9",
    },
    "buddy_check": {"label": "buddy check outlier", "color": "#8300c4"},
    "titan_buddy_check": {
            "label": "titan buddy check outlier",
            "color": "#8300c4"},
    "titan_sct_resistant_check": {
        "label": "sct resistant check outlier",
        "color": "#c17fe1",
    },
    # aggregated
    "outlier": {
        "label": "QC outlier",
        "color": "#f20000",
        "agg_def": True,
    },  # specify it is an aggregated label
    # ----- Gap ----------
    "regular_gap": {"label": "gap", "color": "#f00592"},
    # ----- Interpolation labels -----
    "interpolated_gap": {"label": "interpolation", "color": "#d406c6"},
    "failed_interpolation_gap": {
                "label": "failed interpolation",
                "color": "#d406c6"},
    # ----- raw model gapfill -----
    "raw_modeldata_fill": {"label": "raw modeldata fill", "color": "#6e1868"},
    "failed_raw_modeldata_fill": {
        "label": "failed raw modeldata fill",
        "color": "#6e1868",
    },
    # ----- debias model gapfill -----
    "debias_modeldata_fill": {
                "label": "debiased modeldata fill",
                "color": "#6e1868"},
    "failed_debias_modeldata_fill": {
        "label": "failed debiased modeldata fill",
        "color": "#6e1868",
    },
    "diurnal_debias_modeldata_fill": {
        "label": "diurnal debiased modeldata fill",
        "color": "#6e1868",
    },
    "failed_diurnal_debias_modeldata_fill": {
        "label": "failed diurnal debiased modeldata fill",
        "color": "#6e1868",
    },
    "weighted_diurnal_debias_modeldata_fill": {
        "label": "Weighted diurnal debiased modeldata fill",
        "color": "#6e1868",
    },
    "failed_weighted_diurnal_debias_modeldata_fill": {
        "label": "failed Weighted diurnal debiased modeldata fill",
        "color": "#6e1868",
    },
}

label_to_color_map = {
    group["label"]: group["color"] for group in label_def.values()
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
    "titan_buddy_check",
    "titan_sct_resistant_check",
]
