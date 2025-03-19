import logging
from typing import Literal, Tuple

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from metobs_toolkit.sensordata import SensorData
import metobs_toolkit.settings_files.default_formats_settings as defaults

# ------------------------------------------
#    Label groups
# ------------------------------------------
plot_as_line_labels = present_as_line_labels = [
    defaults.label_def["goodrecord"]["label"],  #'ok']
    defaults.label_def["uncheckedrecord"]["label"],  #'not checked'
] + [defaults.label_def[trglab]["label"] for trglab in defaults.gapfill_label_group]


all_gap_without_val_labels = [
    defaults.label_def["regular_gap"]["label"],  #'gap'
    # defaults.label_def["duplicated_timestamp"]["label"],  # duplicated timestamp outlier
    # defaults.label_def["invalid_input"]["label"],
] + [
    defaults.label_def[trglab]["label"]
    for trglab in defaults.failed_gapfill_label_group
]


all_outlier_labels = [
    defaults.label_def[cat]["label"]
    for cat in defaults.qc_label_group
    if cat not in ["dupliacted_timestamp", "invalid_input"]
]


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
    # plotting trough pandas is problematic when the existing axes,
    # has a timerange different from the line that is added ???
    # ax = series.plot(kind='line',
    #             ax=ax,
    #             color=color,
    #             linewidth=linewidth,
    #             linestyle=linestyle,
    #             zorder=zorder,
    #             label=legend_label,
    #             )

    ax.plot(
        series.index,  # x
        series.values,  # y
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
        x=idx,
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
#    High level plotting
# ------------------------------------------


def plot_timeseries_color_by_label(
    plotdf: pd.DataFrame, show_gaps: bool, show_outliers: bool, ax: plt.Axes
) -> plt.Axes:
    # drop obstype column (not relevant for now)
    plotdf = (
        plotdf.reset_index()
        .drop(columns=["obstype"])
        .set_index(["name", "datetime"])
        .sort_index()
    )

    # Create labels to filter
    labels_to_plot = [
        defaults.label_def["goodrecord"]["label"],  #'ok'
        defaults.label_def["uncheckedrecord"]["label"],  #'not checked'
    ]

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

    # get min max values for the vertical lines
    ymin, ymax = plotdf["value"].min(), plotdf["value"].max()
    if pd.isnull(ymin):
        # ymin and ymax are nan
        ymin = 0.0
        ymax = 10

    for label in labels_to_plot:  # iterate over all labels to plot
        if label not in plotdf["label"].values:
            continue

        # 1. Plot the lines --> be aware of interpolation issues!
        if label in plot_as_line_labels:
            # solid lines for good records, else dashed
            if label == defaults.label_def["goodrecord"]["label"]:
                linestyle = "-"
            else:
                linestyle = "--"
            # iterate over stations --> to avoid interpolation over multiple stations
            for _staname, stadf in plotdf.groupby(
                plotdf.index.get_level_values("name")
            ):
                # filter to label, convert all other records to nan values (to avoid interpoltion over other labeled records)
                stalabeldf = stadf[stadf["label"] == label]
                if stalabeldf.empty:
                    continue

                # IMPORTANT!: add all other records as Nan (otherwise interpolation issue)
                stalabeldf = stalabeldf.reindex(stadf.index, method=None)

                # format the series to plot
                plotseries = stalabeldf["value"]
                plotseries.index = plotseries.index.droplevel("name")

                # add the line to the axes
                ax = add_lines_to_axes(
                    ax=ax,
                    # Reindex to continious timestamps (timestaps without this label have nan values)
                    series=plotseries,
                    legend_label=label,
                    linestyle=linestyle,
                    color=defaults.label_to_color_map[label],
                )

        # Note: no need to add it in the itergroups, no interpolation can be done
        # over different stations, sinc the plot representation (scatter/vlines)
        # do not interpolate.

        # 2. Plot data in vertical line representation (= no numerical values)
        elif label in all_gap_without_val_labels:

            # Note: a regular subset must be done since data is represented as vlines (thus no false interpolation)
            labelseries = plotdf[plotdf["label"] == label]["value"]
            if labelseries.empty:
                continue

            # format the series to plot

            labelseries.index = labelseries.index.droplevel("name")
            labelseries = labelseries.sort_index()  # testing
            # plot
            ax = add_vertical_lines_to_axes(
                ax=ax,
                ymin=ymin,
                ymax=ymax,
                idx=labelseries.index,
                legend_label=label,
                color=defaults.label_to_color_map[label],
            )

        # 3. Plot data in scatter representation (=outliers with numerical values)
        elif label in all_outlier_labels:
            # Note: a regular subset must be done since data is represented as scatters (thus no false interpolation)
            labelseries = plotdf[plotdf["label"] == label]["value"]
            if labelseries.empty:
                continue
            # format the series to plot
            labelseries.index = labelseries.index.droplevel("name")
            ax = add_scatters_to_axes(
                ax=ax,
                series=labelseries,
                legend_label=label,
                color=defaults.label_to_color_map[label],
            )
        else:
            print(f"{label} is not plotted ERROR !! ")

    return ax


def plot_timeseries_color_by_station(
    plotdf: pd.DataFrame,
    colormap: dict,  # {stationname: color}
    ax: plt.Axes,
    show_gaps: bool,
    show_outliers: bool,
    linestyle: str = "-",
    legend_prefix: str = "",
) -> plt.Axes:

    # Get the main data in DataFrame style
    # drop obstype column (not relevant for now)
    plotdf = (
        plotdf.reset_index()
        .drop(columns=["obstype"])
        .set_index(["name", "datetime"])
        .sort_index()
    )

    # Handle gaps
    if not show_gaps:
        all_gap_labels = all_gap_without_val_labels + [
            defaults.label_def[trglab]["label"]
            for trglab in defaults.gapfill_label_group
        ]
        plotdf.loc[plotdf["value"].isin(all_gap_labels), "value"] = np.nan

    # Handle outliers
    if not show_outliers:
        plotdf.loc[plotdf["value"].isin(all_outlier_labels), "value"] = np.nan

    # Plot the data as a single color line
    # iterate over stations --> to avoid interpolation over multiple stations
    for staname, stadf in plotdf.groupby(plotdf.index.get_level_values("name")):

        # format the series to plot
        plotseries = stadf["value"]
        plotseries.index = plotseries.index.droplevel("name")

        ax = add_lines_to_axes(
            ax=ax,
            series=plotseries,
            legend_label=f"{legend_prefix}{staname}",
            linestyle=linestyle,
            color=colormap[staname],
        )

    return ax
