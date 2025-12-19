from __future__ import annotations

import logging  # Python default package
from typing import TYPE_CHECKING

import pandas as pd
from metobs_toolkit.settings_collection import Settings
from metobs_toolkit.plot_collection import (  # Local modules
    create_categorical_color_map,
)

# Set up logging
from metobs_toolkit.backend_collection.decorators import log_entry

if TYPE_CHECKING:
    from matplotlib.pyplot import Axes

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def make_diurnal_plot(
    plotdf: pd.DataFrame,
    ax: Axes,
    colordict: dict,
    refstation: str,
) -> Axes:
    """
    Create a diurnal plot for the given data.

    Parameters
    ----------
    plotdf : pandas.DataFrame
        DataFrame containing the data to plot. Each column represents a station or category.
    ax : matplotlib.axes.Axes
        The matplotlib Axes object where the plot will be drawn.
    colordict : dict
        Dictionary mapping column names to colors. If None, a default colormap will be created.
    refstation : str
        The reference station to be plotted as a dashed line. If None, no reference is plotted.

    Returns
    -------
    matplotlib.axes.Axes
        The Axes object with the plot.

    Raises
    ------
    ValueError
        If not all present labels are in the colormap.
    """
    # Create and check colordict
    if colordict is None:
        logger.debug("Creating default colormap.")
        colmap = create_categorical_color_map(
            catlist=plotdf.columns,
            cmapname=Settings.get("plotting_settings.coloring.categorical_cmap"),
        )
    else:
        colmap = colordict

    if refstation is not None:
        logger.debug(f"Plotting reference station: {refstation}.")
        # Plot reference as dashed line
        ax.axhline(
            **Settings.get("plotting_settings.cycle_plot.hlinekwargs", {"y": 0}),
            label=f"Reference:{refstation}",
        )

    # Check if colormap is valid for the data
    if not all([col in colmap.keys() for col in plotdf.columns]):
        missing = list(set(plotdf.columns) - set(colmap.keys()))
        logger.error(f"Missing labels in colormap: {missing}.")
        raise ValueError(
            f"Not all present labels are in the colmap, these are missing: {missing}"
        )

    # Plot each column as a solid line
    for column in plotdf.columns:
        logger.debug(f"Plotting column: {column}.")
        ax.plot(
            plotdf.index,
            plotdf[column],
            label=column,
            color=colmap[column],
            linestyle="-",
            **Settings.get("plotting_settings.cycle_plot.linekwargs", {}),
        )

    logger.debug("Exiting make_diurnal_plot function.")
    return ax
