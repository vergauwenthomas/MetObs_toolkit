default_plot_settings = {}

# =============================================================================
# Timeseries plots
# =============================================================================
default_plot_settings["time_series"] = {
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
# Stats plot settings
# =============================================================================

default_plot_settings["pie_charts"] = {
    # shape
    "figsize": (15, 10),
    "ncols": 4,
    # "anchor_legend_big": (-0.25, 0.75),
    # "anchor_legend_small": (-3.5, 2.2),
    "radius_big": 1.0,
    "radius_small": 0.7,
    "txt_size_big_pies": 7,
    "txt_size_small_pies": 5,
}

# =============================================================================
# Diurnal plot settings
# =============================================================================

default_plot_settings["cycle_plot"] = {
    "figsize": (10, 10),
    # "alpha_error_bands": 0.3,
    # "cmap_continious": "viridis",  # if many stations are present, best to use continious rather than categorical
    # "n_cat_max": 20,  # when less or equal categories are detected, use the categorical col mapping
    "cmap_categorical": "tab20",
    "legend_n_columns": 5,
}


# ------------------------------------------
#    Default colormaps
# ------------------------------------------

default_colormaps = {
    'lcz': {
        "Compact highrise": "#8c0000",
        "Compact midrise": "#d10000",
        "Compact lowrise": "#ff0000",
        "Open highrise": "#bf4d00",
        "Open midrise": "#ff6600",
        "Open lowrise": "#ff9955",
        "Lightweight lowrise": "#faee05",
        "Large lowrise": "#bcbcbc",
        "Sparsely built": "#ffccaa",
        "Heavy industry": "#555555",
        "Dense Trees (LCZ A)": "#006a00",
        "Scattered Trees (LCZ B)": "#00aa00",
        "Bush, scrub (LCZ C)": "#648525",
        "Low plants (LCZ D)": "#b9db79",
        "Bare rock or paved (LCZ E)": "#000000",
        "Bare soil or sand (LCZ F)": "#fbf7ae",
        "Water (LCZ G)": "#6a6aff"},
    
    #SCORES
    "scores": {
        "RMSE": "#1f77b4",
        "MAE": "#ff7f0e",
        "MSE": "#2ca02c",
        "R2": "#d62728",
        "N": "#42027e",
        "modelbias": "#9467bd"}
    }