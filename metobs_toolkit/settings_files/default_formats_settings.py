#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""


plot_settings = {}
# =============================================================================
# General plot settings
# =============================================================================

# =============================================================================
# Default obs and metadata settings
# =============================================================================

# Static fields are fields (attributes and observations) that do not change in time
static_fields = [
    "network",
    "name",
    "lat",
    "lon",  # TODO make these dynamic, now used as static
    "call_name",
    "location",
    "lcz",
]

# Categorical fields are fields with values that are assumed to be categorical.
# Note: (there are static and dynamic fields that are categorical)
categorical_fields = ["wind_direction", "lcz"]


location_info = ["network", "lat", "lon", "lcz", "call_name", "location"]


default_name = "unknown_name"  # used when no station names are available

# =============================================================================
# Timeseries plots
# =============================================================================
plot_settings["time_series"] = {
    # shape
    "figsize": (15, 5),
    "colormap": "tab20",  # when colorby='name' is used
    "linewidth": 2,  #
    "linestyle_ok": "-",  # solid line
    "linestyle_fill": "--",  # dashedline
    # "linezorder": 1,  # for ok obs
    "scattersize": 4,  # for outliers
    "scatterzorder": 3,  # for outliers
    "dashedzorder": 2,  # for gapfills
    "legend_n_columns": 5,  # for the number of columns in the plot
}
# =============================================================================
# Spatial plot settings
# =============================================================================

plot_settings["spatial_geo"] = {
    # projection
    # 'proj' : 'Orthographic', #Orthographic or AlbersEqualArea
    "extent": [
        2.260609,
        49.25,
        6.118359,
        52.350618,
    ],  # used if observatioons are within
    # colors
    # 'cmap' : 'Set1',
    "cmap": "inferno_r",
    "n_for_categorical": 5,  # number of quantiles for cat data (not for LCZ)
    # shape
    "figsize": (10, 15),
    # datetime
    "fmt": "%d/%m/%Y %H:%M:%S UTC",
}

# =============================================================================
# Stats plot settings
# =============================================================================

plot_settings["pie_charts"] = {
    # shape
    "figsize": (10, 10),
    "anchor_legend_big": (-0.25, 0.75),
    "anchor_legend_small": (-3.5, 2.2),
    "radius_big": 2.0,
    "radius_small": 5.0,
    "effectiveness_colormap": {
        "ok": "green",
        "not checked": "orange",
        "outlier": "red",
    },
}

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
    "invalid_input": {"label": "invalid input", "color": "#900357"},
    "gross_value": {"label": "gross value outlier", "color": "#f1ff2b"},
    "persistence": {"label": "persistence outlier", "color": "#f0051c"},
    "repetitions": {"label": "repetitions outlier", "color": "#056ff0"},
    "step": {"label": "in step outlier group", "color": "#05d4f0"},
    "window_variation": {
        "label": "in window variation outlier group",
        "color": "#05f0c9",
    },
    "buddy_check": {"label": "buddy check outlier", "color": "#8300c4"},
    "titan_buddy_check": {"label": "titan buddy check outlier", "color": "#8300c4"},
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
    "failed_interpolation_gap": {"label": "failed interpolation", "color": "#d406c6"},
    # ----- raw model gapfill -----
    "raw_modeldata_fill": {"label": "raw modeldata fill", "color": "#6e1868"},
    "failed_raw_modeldata_fill": {
        "label": "failed raw modeldata fill",
        "color": "#6e1868",
    },
    # ----- debias model gapfill -----
    "debias_modeldata_fill": {"label": "debiased modeldata fill", "color": "#6e1868"},
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

# make lists per label-theme, so if a new thematic label is added, functionallity
# through the full toolkit is done by adding it to the defenition and the labelgroup

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
    "invalid_input",
    "gross_value",
    "persistence",
    "repetitions",
    "step",
    "window_variation",
    "buddy_check",
    "titan_buddy_check",
    "titan_sct_resistant_check",
]


# =============================================================================
# Diurnal plot settings
# =============================================================================

plot_settings["diurnal"] = {
    "figsize": (10, 10),
    "alpha_error_bands": 0.3,
    "cmap_continious": "viridis",  # if many stations are present, best to use continious rather than categorical
    "n_cat_max": 20,  # when less or equal categories are detected, use the categorical col mapping
    "cmap_categorical": "tab20",
    "legend_n_columns": 5,
}

plot_settings["annual"] = {
    "figsize": (10, 10),
    "alpha_error_bands": 0.3,
    "cmap_continious": "viridis",  # if many stations are present, best to use continious rather than categorical
    "n_cat_max": 20,  # when less or equal categories are detected, use the categorical col mapping
    "cmap_categorical": "tab20",
    "legend_n_columns": 5,
}


# =============================================================================
# correlation plot settings
# =============================================================================

plot_settings["correlation_heatmap"] = {
    "figsize": (10, 10),
    "vmin": -1,
    "vmax": 1,
    "cmap": "cool",
    "x_tick_rot": 65,
    "y_tick_rot": 0,
}

plot_settings["correlation_scatter"] = {
    "figsize": (10, 10),
    "p_bins": [0, 0.001, 0.01, 0.05, 999],  # do not change the 0.001,0.01 or 0.05
    "bins_markers": ["*", "s", "^", "X"],
    "scatter_size": 40,
    "scatter_edge_col": "black",
    "scatter_edge_line_width": 0.1,
    "ymin": -1.1,
    "ymax": 1.1,
    "cmap": "tab20",
    "legend_ncols": 3,
    "legend_text_size": 7,
}


print_settings = {"fmt_datetime": "%d/%m/%Y %H:%M:%S", "max_print_per_line": "40"}


# =============================================================================
# variables display strings
# =============================================================================
vars_display = {
    "network": "network",
    "name": "station name",
    "call_name": "pseudo name",
    "location": "region",
    "lat": "latitude",
    "lon": "longtitude",
    "temp": "temperature",
    "radiation_temp": "radiation temperature",
    "humidity": "humidity",
    "precip": "precipitation intensity",
    "precip_sum": "cummulated precipitation",
    "wind_speed": "wind speed",
    "wind_gust": "wind gust speed",
    "wind_direction": "wind direction",
    "pressure": "air pressure",
    "pressure_at_sea_level": "corrected pressure at sea level",
    "lcz": "LCZ",
}
