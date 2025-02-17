import logging
from typing import Literal, Tuple
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt


from metobs_toolkit.sensordata import SensorData
import metobs_toolkit.settings_files.default_formats_settings as defaults

# ------------------------------------------
#    Label groups
# ------------------------------------------
plot_as_line_labels = present_as_line_labels = [
    defaults.label_def["goodrecord"]["label"],  #'ok']
    defaults.label_def["uncheckedrecord"]["label"],
] + defaults.gapfill_label_group  #'not checked'

plot_as_vertical_line_labels = [
    defaults.label_def["regular_gap"]["label"],  #'gap'
    defaults.label_def["duplicated_timestamp"]["label"],  # duplicated timestamp outlier
    defaults.label_def["invalid_input"]["label"],
] + defaults.failed_gapfill_label_group  # invalid input

plot_as_scatter_labels = [
    defaults.label_def[cat]["label"]
    for cat in defaults.qc_label_group
    if cat not in ["dupliacted_timestamp", "invalid_input"]
]


def create_station_color_map(catlist: list, cmapname: str = "tab20") -> dict:
    unique_elements = list(set(catlist))
    num_unique_elements = len(unique_elements)

    cmap = plt.get_cmap(cmapname)
    colors = [cmap(i / num_unique_elements) for i in range(num_unique_elements)]

    hex_colors = [matplotlib.colors.rgb2hex(color) for color in colors]

    color_map = {
        element: hex_colors[i % num_unique_elements]
        for i, element in enumerate(unique_elements)
    }

    return color_map


# ------------------------------------------
#    General plotting functions
# ------------------------------------------


def create_axes(figsize: Tuple[int, int] = (15, 5), **kwargs) -> plt.Axes:

    _fig, ax = plt.subplots(figsize=figsize, **kwargs)
    return ax


# ------------------------------------------
#    Timeseries plot layers
# ------------------------------------------


def add_lines_to_axes(
    ax: plt.Axes,
    series: pd.Series,
    legend_label: str,
    linestyle: Literal["-", "--", "-.", ":", ""] = "-",
    color: str = "navy",
    linewidth: int = 2,
    zorder: int | float = 1,
) -> plt.Axes:
    series.plot(
        ax=ax,
        color=color,
        linewidth=linewidth,
        linestyle=linestyle,
        zorder=zorder,
        label=legend_label,
    )
    return ax


def add_vertical_lines_to_axes(
    ax: plt.Axes,
    idx: pd.DatetimeIndex,
    legend_label: str,
    ymin: float,
    ymax: float,
    linestyle: Literal["-", "--", "-.", ":", ""] = "-",
    color: str = "navy",
    linewidth: int = 2,
    zorder: int | float = 1,
) -> plt.Axes:
    ax.vlines(
        x=idx.to_numpy(),
        ymin=ymin,
        ymax=ymax,
        linestyle=linestyle,
        linewidth=linewidth,
        color=color,
        zorder=zorder,
        label=legend_label,
    )
    return ax


def add_scatters_to_axes(
    ax: plt.Axes,
    series: pd.Series,
    legend_label: str,
    color: str = "navy",
    scattersize: int = 2,
    zorder: int | float = 1,
) -> plt.Axes:
    ax.scatter(
        series.index,
        series.values,
        color=color,
        s=scattersize,
        zorder=zorder,
        label=legend_label,
    )
    return ax


# ------------------------------------------
#    Styling
# ------------------------------------------
# These functions seems trivial, but by defining them in one place,
# makes it easyer to fix styling issues on all levels (sensor, station, dataset)
def set_title(ax: plt.Axes, titlestr: str) -> plt.Axes:
    ax.set_title(titlestr)
    return ax


def set_ylabel(ax: plt.Axes, ylabel: str) -> plt.Axes:
    ax.set_ylabel(ylabel)
    return ax


def set_xlabel(ax: plt.Axes, xlabel: str) -> plt.Axes:
    ax.set_xlabel(xlabel)
    return ax


def set_legend(ax: plt.Axes) -> plt.Axes:
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # this drops duplicates based on labels
    ax.legend(by_label.values(), by_label.keys())
    return ax


# ------------------------------------------
#    High level plotting
# ------------------------------------------
def plot_timeseries_color_by_label(
    sensordata: SensorData, show_gaps: bool, show_outliers: bool, ax: plt.Axes
) -> plt.Axes:
    # to use the same continious x-records for all plotlayers
    target_dt = pd.date_range(
        start=sensordata.start_datetime,
        end=sensordata.end_datetime,
        freq=sensordata.freq,
    )

    # Create labels to filter
    labels_to_plot = [
        defaults.label_def["goodrecord"]["label"],  #'ok'
        defaults.label_def["uncheckedrecord"]["label"],
    ]  #'not checked'
    if show_gaps:
        # add all labels related to gaps
        labels_to_plot += (
            [defaults.label_def["regular_gap"]["label"]]  #'gap'
            + [
                defaults.label_def[cat]["label"] for cat in defaults.gapfill_label_group
            ]  #'interpolated_gap', 'raw_modeldata_fill', ...
            + [
                defaults.label_def[cat]["label"]
                for cat in defaults.failed_gapfill_label_group
            ]  #'failed_interpolated_gap', ...
        )
    if show_outliers:
        labels_to_plot += [
            defaults.label_def[cat]["label"] for cat in defaults.qc_label_group
        ]  #'duplicated_timestamp', 'gross_value', ...

    # Get data in DataFrame style
    plotdf = (
        sensordata.df.reset_index()
        .drop(columns=["obstype"])  # is this oke when called from station/dataset obj
        .set_index("datetime")
    )
    # filter to relevant records
    plotdf = plotdf[plotdf["label"].isin(labels_to_plot)]
    # ymin, ymax are required for vertical lines
    ymin, ymax = plotdf["value"].min(), plotdf["value"].max()

    # 1. Plot data in line representation
    for label in plot_as_line_labels:
        labelseries = plotdf[plotdf["label"] == label]["value"]
        # skip if label is not present
        if labelseries.empty:
            continue
        # solid lines for good records, else dashed
        if label == defaults.label_def["goodrecord"]["label"]:
            linestyle = "-"
        else:
            linestyle = "--"
        # add the line to the axes
        ax = add_lines_to_axes(
            ax=ax,
            # Reindex to continious timestamps (timestaps without this label have nan values)
            series=labelseries.reindex(target_dt, method=None),
            legend_label=label,
            linestyle=linestyle,
            color=defaults.label_to_color_map[label],
        )

    # 2. Plot data in vertical line representation (= no numerical values)
    for label in plot_as_vertical_line_labels:
        labelseries = plotdf[plotdf["label"] == label]["value"]
        if labelseries.empty:
            continue
        ax = add_vertical_lines_to_axes(
            ax=ax,
            ymin=ymin,
            ymax=ymax,
            idx=labelseries.index,
            legend_label=label,
            color=defaults.label_to_color_map[label],
        )

    # 3. Plot data in scatter representation (=outliers with numerical values)
    for label in plot_as_scatter_labels:
        labelseries = plotdf[plotdf["label"] == label]["value"]
        if labelseries.empty:
            continue
        ax = add_scatters_to_axes(
            ax=ax,
            series=labelseries,
            legend_label=label,
            color=defaults.label_to_color_map[label],
        )

    return ax


def plot_timeseries_as_one_color(
    sensordata: SensorData,
    color: str,
    ax: plt.Axes,
    show_gaps: bool,
    show_outliers: bool,
    linestyle: str = "-",
) -> plt.Axes:

    # Get the main data in DataFrame style
    plotdf = (
        sensordata.df.reset_index()
        .drop(columns=["obstype"])  # is this okay when called from station/dataset obj
        .set_index("datetime")
    )

    # Handle gaps
    if not show_gaps:
        gaps_idx = sensordata.gapsdf.index
        plotdf.loc[gaps_idx, "value"] = np.nan

    # Handle outliers
    if not show_outliers:
        outliers_idx = sensordata.outliersdf.index
        plotdf.loc[outliers_idx, "value"] = np.nan

    # Plot the data as a single color line
    ax = add_lines_to_axes(
        ax=ax,
        series=plotdf["value"],
        legend_label=sensordata.stationname,
        linestyle=linestyle,
        color=color,
    )

    return ax
