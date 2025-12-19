default_plot_settings = {}

# ------------------------------------------
#    General plotting defaults
# ------------------------------------------
default_plot_settings["coloring"] = {
    "categorical_cmap": "tab20",
    # 'continious_cmap': 'viridis',
}

# =============================================================================
# Timeseries plots
# =============================================================================
default_plot_settings["time_series"] = {
    # mpl figure settings
    "figkwargs": {"figsize": (15, 5), "tight_layout": True},
    # legend settings
    "legendkwargs": {
        "loc": "upper center",
        "bbox_to_anchor": (0.5, -0.15),
        "ncol": 5,
        "frameon": False,
    },
    # line settings
    "linekwargs": {
        "linewidth": 2,  #
    },
}

# =============================================================================
# Stats plot settings
# =============================================================================

default_plot_settings["pie_charts"] = {
    "figkwargs": {"figsize": (10, 10), "tight_layout": True},
    # grid layout settings
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
    # mpl figure settings
    "figkwargs": {"figsize": (15, 10), "tight_layout": True},
    # legend settings
    "legendkwargs": {
        "loc": "upper center",
        "bbox_to_anchor": (0.5, -0.15),
        "ncol": 5,
    },
    # line settings
    "linekwargs": {
        "linewidth": 2,  #
    },
    # horizontal line settings
    "hlinekwargs": {
        "y": 0,
        "color": "black",
        "linestyle": "--",
        "zorder": 0.9,
        "linewidth": 0.8,
    },
}
