import logging  # Python default package

import matplotlib
import matplotlib.pyplot as plt  # noqa: F401  # Dependency package
import pandas as pd

from metobs_toolkit.plot_collection import (  # Local modules
    create_categorical_color_map,
    default_plot_settings,
)

# Set up logging
from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")

default_cycle_settings = default_plot_settings["cycle_plot"]


@log_entry
def make_diurnal_plot(
    plotdf: pd.DataFrame,
    ax: "matplotlib.axes.Axes",
    colordict: dict,
    refstation: str,
    figkwargs: dict,
) -> "matplotlib.axes.Axes":
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
    figkwargs : dict
        Additional keyword arguments for figure customization.

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
            catlist=plotdf.columns, cmapname=default_cycle_settings["cmap_categorical"]
        )
    else:
        colmap = colordict

    if refstation is not None:
        logger.debug(f"Plotting reference station: {refstation}.")
        # Plot reference as dashed line
        ax.axhline(
            y=0,
            color="black",
            linestyle="--",
            zorder=0.9,
            linewidth=0.8,
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
        )

    logger.debug("Exiting make_diurnal_plot function.")
    return ax
