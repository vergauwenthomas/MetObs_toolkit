#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Label color and name defenitions for the metobs_toolkit.

@author: thoverga
"""
line = "line"
scatter = "scatter"
vline = "vline"


# =============================================================================
# Defenition of labels
# =============================================================================
# All knonw labels should be defined here:
label_defs = {
    # --- Good records --------
    "goodrecord": {
        "label": "ok",
        "color": "#07f72b",
        "numeric_val": 0,
        "plot_as": line,
        "ls": "-",
    },
    "uncheckedrecord": {
        "label": "not checked",
        "color": "#f7cf07",
        "numeric_val": 1,
        "plot_as": line,
        "ls": "-",
    },
    # ------ QC labels ---------------
    "duplicated_timestamp": {
        "label": "duplicated timestamp outlier",
        "color": "#a32a1f",
        "numeric_val": 2,
        "plot_as": scatter,
        "marker": ".",
    },
    # "invalid_input": {"label": "invalid input", "color": "#900357"},
    "gross_value": {
        "label": "gross value outlier",
        "color": "#f1ff2b",
        "numeric_val": 3,
        "plot_as": scatter,
        "marker": ".",
    },
    "persistence": {
        "label": "persistence outlier",
        "color": "#f0051c",
        "numeric_val": 4,
        "plot_as": scatter,
        "marker": ".",
    },
    "repetitions": {
        "label": "repetitions outlier",
        "color": "#056ff0",
        "numeric_val": 5,
        "plot_as": scatter,
        "marker": ".",
    },
    "step": {"label": "in step outlier group", "color": "#05d4f0", "numeric_val": 6},
    "window_variation": {
        "label": "in window variation outlier group",
        "color": "#05f0c9",
        "numeric_val": 7,
        "plot_as": scatter,
        "marker": ".",
    },
    "buddy_check": {
        "label": "buddy check outlier",
        "color": "#8300c4",
        "numeric_val": 8,
        "plot_as": scatter,
        "marker": ".",
    },
    "buddy_check_with_safetynets": {
        "label": "buddy check (with safetynets) outlier",
        "color": "#8300c4",
        "numeric_val": 9,
        "plot_as": scatter,
        "marker": ".",
    },
    # aggregated
    "outlier": {
        "label": "QC outlier",
        "color": "#f20000",
        "agg_def": True,
        "numeric_val": 10,
        "plot_as": scatter,
        "marker": ".",
    },  # specify it is an aggregated label
    # ----- Gap ----------
    "regular_gap": {
        "label": "gap",
        "color": "#f00592",
        "numeric_val": 11,
        "plot_as": vline,
        "ls": "-",
    },
    # ----- Interpolation labels -----
    "interpolated_gap": {
        "label": "interpolation",
        "color": "#d406c6",
        "numeric_val": 12,
        "plot_as": line,
        "ls": "--",
    },
    "failed_interpolation_gap": {
        "label": "failed interpolation",
        "color": "#d406c6",
        "numeric_val": 13,
        "plot_as": vline,
        "ls": "-",
    },
    # ----- raw model gapfill -----
    "raw_modeldata_fill": {
        "label": "raw modeldata fill",
        "color": "#6e1868",
        "numeric_val": 14,
        "plot_as": line,
        "ls": "--",
    },
    "failed_raw_modeldata_fill": {
        "label": "failed raw modeldata fill",
        "color": "#6e1868",
        "numeric_val": 15,
        "plot_as": vline,
        "ls": "-",
    },
    # ----- debias model gapfill -----
    "debias_modeldata_fill": {
        "label": "debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 16,
        "plot_as": line,
        "ls": "--",
    },
    "failed_debias_modeldata_fill": {
        "label": "failed debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 17,
        "plot_as": vline,
        "ls": "-",
    },
    "diurnal_debias_modeldata_fill": {
        "label": "diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 18,
        "plot_as": line,
        "ls": "--",
    },
    "failed_diurnal_debias_modeldata_fill": {
        "label": "failed diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 19,
        "plot_as": vline,
        "ls": "-",
    },
    "weighted_diurnal_debias_modeldata_fill": {
        "label": "Weighted diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 20,
        "plot_as": line,
        "ls": "--",
    },
    "failed_weighted_diurnal_debias_modeldata_fill": {
        "label": "failed Weighted diurnal debiased modeldata fill",
        "color": "#6e1868",
        "numeric_val": 21,
        "plot_as": vline,
        "ls": "-",
    },
}

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
    "buddy_check_with_safetynets",
]
