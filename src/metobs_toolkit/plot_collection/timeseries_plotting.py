import logging
from typing import Literal, Union

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


import metobs_toolkit.settings_collection as settings

# ------------------------------------------
#    Label groups
# ------------------------------------------
plot_as_line_labels = [
    settings.label_def["goodrecord"]["label"],  # 'ok'
    settings.label_def["uncheckedrecord"]["label"],  # 'not checked'
] + [settings.label_def[trglab]["label"] for trglab in settings.gapfill_label_group]

all_gap_without_val_labels = [
    settings.label_def["regular_gap"]["label"],  # 'gap'
] + [
    settings.label_def[trglab]["label"]
    for trglab in settings.failed_gapfill_label_group
]

all_outlier_labels = [
    settings.label_def[cat]["label"]
    for cat in settings.qc_label_group
    if cat
    not in ["duplicated_timestamp", "invalid_input"]  # TYPO: 'dupliacted' corrected
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
    zorder: Union[int, float] = 1.3,
) -> plt.Axes:
    """
    Add a line plot to the given axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to which the line will be added.
    series : pandas.Series
        The data series to plot. The index represents the x-axis, and the values represent the y-axis.
    legend_label : str
        The label for the legend entry.
    linestyle : {'-', '--', '-.', ':', ''}, optional
        The style of the line. Default is '-'.
    color : str, optional
        The color of the line. Default is 'navy'.
    linewidth : int, optional
        The width of the line. Default is 2.
    zorder : int or float, optional
        The z-order of the line. Default is 1.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the line added.
    """
    logging.info(f"Entering add_lines_to_axes with legend_label={legend_label}")

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
    zorder: Union[int, float] = 1.1,
) -> plt.Axes:
    """
    Add vertical lines to the given axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to which the vertical lines will be added.
    idx : pandas.DatetimeIndex
        The x-axis positions of the vertical lines.
    legend_label : str
        The label for the legend entry.
    ymin : float
        The starting y-coordinate of the vertical lines.
    ymax : float
        The ending y-coordinate of the vertical lines.
    linestyle : {'-', '--', '-.', ':', ''}, optional
        The style of the lines. Default is '-'.
    color : str, optional
        The color of the lines. Default is 'navy'.
    linewidth : int, optional
        The width of the lines. Default is 2.
    zorder : int or float, optional
        The z-order of the lines. Default is 1.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the vertical lines added.
    """
    logging.info(
        f"Entering add_vertical_lines_to_axes with legend_label={legend_label}"
    )

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
    zorder: Union[int, float] = 1.5,
) -> plt.Axes:
    """
    Add scatter points to the given axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to which the scatter points will be added.
    series : pandas.Series
        The data series to plot. The index represents the x-axis, and the values represent the y-axis.
    legend_label : str
        The label for the legend entry.
    color : str, optional
        The color of the scatter points. Default is 'navy'.
    scattersize : int, optional
        The size of the scatter points. Default is 2.
    zorder : int or float, optional
        The z-order of the scatter points. Default is 1.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the scatter points added.
    """
    logging.info(f"Entering add_scatters_to_axes with legend_label={legend_label}")

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
    """
    Plot a timeseries with data points colored by their label.

    Parameters
    ----------
    plotdf : pandas.DataFrame
        The DataFrame containing the data to plot. It must have columns 'name', 'datetime', 'value', and 'label'.
    show_gaps : bool
        Whether to include gap-related labels in the plot.
    show_outliers : bool
        Whether to include outlier-related labels in the plot.
    ax : matplotlib.axes.Axes
        The axes on which the plot will be drawn.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the timeseries plot.
    """
    logging.info("Entering plot_timeseries_color_by_label")

    # Drop 'obstype' column and reformat DataFrame
    plotdf = (
        plotdf.reset_index()
        .drop(columns=["obstype"])
        .set_index(["name", "datetime"])
        .sort_index()
    )

    # Create labels to filter
    labels_to_plot = [
        settings.label_def["goodrecord"]["label"],  # 'ok'
        settings.label_def["uncheckedrecord"]["label"],  # 'not checked'
    ]

    if show_gaps:
        # Add all labels related to gaps
        labels_to_plot += (
            [settings.label_def["regular_gap"]["label"]]  # 'gap'
            + [
                settings.label_def[cat]["label"] for cat in settings.gapfill_label_group
            ]  # 'interpolated_gap', 'raw_modeldata_fill', ...
            + [
                settings.label_def[cat]["label"]
                for cat in settings.failed_gapfill_label_group
            ]  # 'failed_interpolated_gap', ...
        )
    if show_outliers:
        labels_to_plot += [
            settings.label_def[cat]["label"] for cat in settings.qc_label_group
        ]  # 'duplicated_timestamp', 'gross_value', ...

    # Get min and max values for the vertical lines
    ymin, ymax = plotdf["value"].min(), plotdf["value"].max()
    if pd.isnull(ymin):
        # ymin and ymax are NaN
        ymin = 0.0
        ymax = 10

    for label in labels_to_plot:  # Iterate over all labels to plot
        if label not in plotdf["label"].values:
            continue

        # 1. Plot the lines --> be aware of interpolation issues!
        if label in plot_as_line_labels:
            # Solid lines for good records, else dashed
            linestyle = (
                "-" if label == settings.label_def["goodrecord"]["label"] else "--"
            )

            # Iterate over stations --> to avoid interpolation over multiple stations
            for _staname, stadf in plotdf.groupby(
                plotdf.index.get_level_values("name")
            ):
                # Filter to label, convert all other records to NaN values (to avoid interpolation over other labeled records)
                stalabeldf = stadf[stadf["label"] == label]
                if stalabeldf.empty:
                    continue

                # IMPORTANT!: Add all other records as NaN (otherwise interpolation issue)
                stalabeldf = stalabeldf.reindex(stadf.index, method=None)

                # Format the series to plot
                plotseries = stalabeldf["value"]
                plotseries.index = plotseries.index.droplevel("name")

                # Add the line to the axes
                ax = add_lines_to_axes(
                    ax=ax,
                    # Reindex to continuous timestamps (timestamps without this label have NaN values)
                    series=plotseries,
                    legend_label=label,
                    linestyle=linestyle,
                    color=settings.label_to_color_map[label],
                )

        # 2. Plot data in vertical line representation (= no numerical values)
        elif label in all_gap_without_val_labels:
            # Note: A regular subset must be done since data is represented as vlines (thus no false interpolation)
            labelseries = plotdf[plotdf["label"] == label]["value"]
            if labelseries.empty:
                continue

            # Format the series to plot
            labelseries.index = labelseries.index.droplevel("name")
            labelseries = labelseries.sort_index()  # Sorting for consistency
            # Plot
            ax = add_vertical_lines_to_axes(
                ax=ax,
                ymin=ymin,
                ymax=ymax,
                idx=labelseries.index,
                legend_label=label,
                color=settings.label_to_color_map[label],
            )

        # 3. Plot data in scatter representation (= outliers with numerical values)
        elif label in all_outlier_labels:
            # Note: A regular subset must be done since data is represented as scatters (thus no false interpolation)
            labelseries = plotdf[plotdf["label"] == label]["value"]
            if labelseries.empty:
                continue

            # Format the series to plot
            labelseries.index = labelseries.index.droplevel("name")
            ax = add_scatters_to_axes(
                ax=ax,
                series=labelseries,
                legend_label=label,
                color=settings.label_to_color_map[label],
            )
        else:
            logging.warning(f"{label} is not plotted. ERROR!")

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
    """
    Plot a timeseries with data points colored by station.

    Parameters
    ----------
    plotdf : pandas.DataFrame
        The DataFrame containing the data to plot. It must have columns 'name', 'datetime', 'value', and 'label'.
    colormap : dict
        A dictionary mapping station names to colors.
    ax : matplotlib.axes.Axes
        The axes on which the plot will be drawn.
    show_gaps : bool
        Whether to include gap-related labels in the plot.
    show_outliers : bool
        Whether to include outlier-related labels in the plot.
    linestyle : str, optional
        The style of the line. Default is '-'.
    legend_prefix : str, optional
        A prefix to add to the legend labels. Default is an empty string.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the timeseries plot.

    Raises
    ------
    TypeError
        If any of the arguments are not of the expected type.
    """
    logging.info("Entering plot_timeseries_color_by_station")

    # Drop 'obstype' column and reformat DataFrame
    plotdf = (
        plotdf.reset_index()
        .drop(columns=["obstype"])
        .set_index(["name", "datetime"])
        .sort_index()
    )

    # Handle gaps
    if not show_gaps:
        all_gap_labels = all_gap_without_val_labels + [
            settings.label_def[trglab]["label"]
            for trglab in settings.gapfill_label_group
        ]
        plotdf.loc[plotdf["label"].isin(all_gap_labels), "value"] = np.nan

    # Handle outliers
    if not show_outliers:
        plotdf.loc[plotdf["label"].isin(all_outlier_labels), "value"] = np.nan

    # Plot the data as a single color line
    # Iterate over stations to avoid interpolation over multiple stations
    for staname, stadf in plotdf.groupby(plotdf.index.get_level_values("name")):
        # Format the series to plot
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
