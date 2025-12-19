import logging
from typing import Literal, Union

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from metobs_toolkit.backend_collection.errorclasses import MetObsInternalError
from metobs_toolkit.settings_collection import Settings

# ------------------------------------------
#    Label groups
# ------------------------------------------


def all_gap_labels() -> list[str]:
    return (
        [Settings.get("label_def.regular_gap.label")]  # 'gap'
        + [
            Settings.get(f"label_def.{cat}.label")
            for cat in Settings.get("gapfill_label_group")
        ]  # 'interpolated_gap', 'raw_modeldata_fill', ...
        + [
            Settings.get(f"label_def.{cat}.label")
            for cat in Settings.get("failed_gapfill_label_group")
        ]  # 'failed_interpolated_gap', ...
    )


def all_outlier_labels() -> list[str]:
    return [
        Settings.get(f"label_def.{cat}.label") for cat in Settings.get("qc_label_group")
    ]


# ------------------------------------------
#    Timeseries plot layers
# ------------------------------------------


def add_lines_to_axes(
    ax: plt.Axes,
    series: pd.Series,
    legend_label: str,
    **kwargs,
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
    **kwargs
        Additional keyword arguments passed to `matplotlib.axes.Axes.plot`.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the line added.
    """
    ax.plot(series.index, series.values, label=legend_label, **kwargs)  # x  # y
    return ax


def add_vertical_lines_to_axes(
    ax: plt.Axes,
    idx: pd.DatetimeIndex,
    legend_label: str,
    ymin: float,
    ymax: float,
    **kwargs,
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
    **kwargs
        Additional keyword arguments passed to `matplotlib.axes.Axes.vlines`.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the vertical lines added.
    """

    ax.vlines(x=idx, ymin=ymin, ymax=ymax, label=legend_label, **kwargs)
    return ax


def add_scatters_to_axes(
    ax: plt.Axes,
    series: pd.Series,
    legend_label: str,
    **kwargs,
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
    **kwargs
        Additional keyword arguments passed to `matplotlib.axes.Axes.scatter`.

    Returns
    -------
    matplotlib.axes.Axes
        The updated axes with the scatter points added.
    """

    ax.scatter(
        series.index,
        series.values,
        label=legend_label,
        **kwargs,
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
        Settings.get("label_def.goodrecord.label"),  # 'ok'
        Settings.get("label_def.uncheckedrecord.label"),  # 'not checked'
    ]
    if show_gaps:
        labels_to_plot += all_gap_labels()
    if show_outliers:
        labels_to_plot += all_outlier_labels()

    # Get min and max values for the vertical lines
    ymin, ymax = plotdf["value"].min(), plotdf["value"].max()
    if pd.isnull(ymin):
        # ymin and ymax are NaN
        ymin = 0.0
        ymax = 10

    label_to_checkname_map = Settings._label_to_qccheckmap()
    for label in labels_to_plot:  # Iterate over all labels to plot
        if label not in plotdf["label"].values:
            continue

        if Settings._flag_plot_as_line(label):
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
                    **Settings.get(
                        f"label_def.{label_to_checkname_map[label]}.plotkwargs", {}
                    ),
                )
        # 2. Plot data in vertical line representation (= no numerical values)
        elif Settings._flag_plot_as_vline(label):
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
                **Settings.get(
                    f"label_def.{label_to_checkname_map[label]}.plotkwargs", {}
                ),
            )

        # 3. Plot data in scatter representation (= outliers with numerical values)
        elif Settings._flag_plot_as_scatter(label):
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
                **Settings.get(
                    f"label_def.{label_to_checkname_map[label]}.plotkwargs", {}
                ),
            )
        else:
            raise MetObsInternalError(
                f"Label '{label}' is not configured to be plotted as line, vline or scatter."
            )

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
        plotdf.loc[plotdf["label"].isin(all_gap_labels()), "value"] = np.nan

    # Handle outliers
    if not show_outliers:
        plotdf.loc[plotdf["label"].isin(all_outlier_labels()), "value"] = np.nan

    linekwargs = Settings.get("plotting_settings.time_series.linekwargs", {})
    # Plot the data as a single color line
    # Iterate over stations to avoid interpolation over multiple stations
    for staname, stadf in plotdf.groupby(plotdf.index.get_level_values("name")):
        # Format the series to plot
        plotseries = stadf["value"]
        plotseries.index = plotseries.index.droplevel("name")

        # update linekwargs
        linekwargs.update({"linestyle": linestyle, "color": colormap[staname]})

        ax = add_lines_to_axes(
            ax=ax,
            series=plotseries,
            legend_label=f"{legend_prefix}{staname}",
            **linekwargs,
        )

    return ax
