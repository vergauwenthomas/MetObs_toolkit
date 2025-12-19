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
        "numeric_val": 0,
        "plot_as": line,
        "plotkwargs": {
            "zorder": 1.2,
            "color": "#07f72b",
            "ls": "-",
        },
    },
    "uncheckedrecord": {
        "label": "not checked",
        "numeric_val": 1,
        "plot_as": line,
        "plotkwargs": {
            "zorder": 1.2,
            "color": "#f7cf07",
            "ls": "-",
        },
    },
    # ------ QC labels ---------------
    "duplicated_timestamp": {
        "label": "duplicated timestamp outlier",
        "numeric_val": 2,
        "plot_as": scatter,
        "plotkwargs": {"color": "#a32a1f", "marker": "o", "s": 5, "zorder": 1.5},
    },
    # "invalid_input": {"label": "invalid input", "color": "#900357"},
    "gross_value": {
        "label": "gross value outlier",
        "numeric_val": 3,
        "plot_as": scatter,
        "plotkwargs": {"color": "#f1ff2b", "marker": "o", "s": 5, "zorder": 1.5},
    },
    "persistence": {
        "label": "persistence outlier",
        "numeric_val": 4,
        "plot_as": scatter,
        "plotkwargs": {"color": "#f0051c", "marker": "o", "s": 5, "zorder": 1.5},
    },
    "repetitions": {
        "label": "repetitions outlier",
        "numeric_val": 5,
        "plot_as": scatter,
        "plotkwargs": {"color": "#056ff0", "marker": "o", "s": 5, "zorder": 1.5},
    },
    "step": {
        "label": "in step outlier group",
        "numeric_val": 6,
        "plot_as": scatter,
        "plotkwargs": {"color": "#05d4f0", "marker": "o", "s": 5, "zorder": 1.5},
    },
    "window_variation": {
        "label": "in window variation outlier group",
        "numeric_val": 7,
        "plot_as": scatter,
        "plotkwargs": {"color": "#05f0c9", "marker": "o", "s": 5, "zorder": 1.5},
    },
    "buddy_check": {
        "label": "buddy check outlier",
        "numeric_val": 8,
        "plot_as": scatter,
        "plotkwargs": {"color": "#8300c4", "marker": "o", "s": 5, "zorder": 1.5},
    },
    "buddy_check_with_safetynets": {
        "label": "buddy check (with safetynets) outlier",
        "color": "#8300c4",
        "numeric_val": 9,
        "plot_as": scatter,
        "plotkwargs": {"color": "#8300c4", "marker": "o", "s": 5, "zorder": 1.5},
    },
    # aggregated
    "outlier": {
        "label": "QC outlier",
        "agg_def": True,
        "numeric_val": 10,
        "plot_as": scatter,
        "plotkwargs": {"color": "#f20000", "marker": "o", "s": 5, "zorder": 1.5},
    },  # specify it is an aggregated label
    # ----- Gap ----------
    "regular_gap": {
        "label": "gap",
        "numeric_val": 11,
        "plot_as": vline,
        "plotkwargs": {"color": "#f00592", "ls": "-", "linewidth": 2, "zorder": 1.0},
    },
    # ----- Interpolation labels -----
    "interpolated_gap": {
        "label": "interpolation",
        "numeric_val": 12,
        "plot_as": line,
        "plotkwargs": {"color": "#d406c6", "ls": "--", "linewidth": 2, "zorder": 1.3},
    },
    "failed_interpolation_gap": {
        "label": "failed interpolation",
        "numeric_val": 13,
        "plot_as": vline,
        "plotkwargs": {"color": "#d406c6", "linewidth": 2, "zorder": 1.1},
    },
    # ----- raw model gapfill -----
    "raw_modeldata_fill": {
        "label": "raw modeldata fill",
        "numeric_val": 14,
        "plot_as": line,
        "plotkwargs": {"color": "#6e1868", "ls": "--", "linewidth": 2, "zorder": 1.3},
    },
    "failed_raw_modeldata_fill": {
        "label": "failed raw modeldata fill",
        "numeric_val": 15,
        "plot_as": vline,
        "plotkwargs": {"color": "#6e1868", "linewidth": 2, "zorder": 1.1},
    },
    # ----- debias model gapfill -----
    "debias_modeldata_fill": {
        "label": "debiased modeldata fill",
        "numeric_val": 16,
        "plot_as": line,
        "plotkwargs": {"color": "#6e1868", "ls": "--", "linewidth": 2, "zorder": 1.3},
    },
    "failed_debias_modeldata_fill": {
        "label": "failed debiased modeldata fill",
        "numeric_val": 17,
        "plot_as": vline,
        "plotkwargs": {"color": "#6e1868", "linewidth": 2, "zorder": 1.1},
    },
    "diurnal_debias_modeldata_fill": {
        "label": "diurnal debiased modeldata fill",
        "numeric_val": 18,
        "plot_as": line,
        "plotkwargs": {"color": "#6e1868", "ls": "--", "linewidth": 2, "zorder": 1.3},
    },
    "failed_diurnal_debias_modeldata_fill": {
        "label": "failed diurnal debiased modeldata fill",
        "numeric_val": 19,
        "plot_as": vline,
        "plotkwargs": {"color": "#6e1868", "linewidth": 2, "zorder": 1.1},
    },
    "weighted_diurnal_debias_modeldata_fill": {
        "label": "Weighted diurnal debiased modeldata fill",
        "numeric_val": 20,
        "plot_as": line,
        "plotkwargs": {"color": "#6e1868", "ls": "--", "linewidth": 2, "zorder": 1.3},
    },
    "failed_weighted_diurnal_debias_modeldata_fill": {
        "label": "failed Weighted diurnal debiased modeldata fill",
        "numeric_val": 21,
        "plot_as": vline,
        "plotkwargs": {"color": "#6e1868", "linewidth": 2, "zorder": 1.1},
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
