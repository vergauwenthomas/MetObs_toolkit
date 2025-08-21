import logging
import math
from typing import Tuple

import pandas as pd
import matplotlib.pyplot as plt

from metobs_toolkit.settings_collection import label_def, label_to_color_map
from metobs_toolkit.plot_collection import default_plot_settings

# Configure logging
from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")

pieplotsettings = default_plot_settings["pie_charts"]


@log_entry
def qc_overview_pies(
    df: pd.DataFrame,
    figsize: Tuple[int, int] = pieplotsettings["figsize"],
    ncol: int = pieplotsettings["ncols"],
    radius_big: float = pieplotsettings["radius_big"],
    radius_small: float = pieplotsettings["radius_small"],
    textsize_big_pies: int = pieplotsettings["txt_size_big_pies"],
    textsize_small_pies: int = pieplotsettings["txt_size_small_pies"],
) -> plt.Figure:
    """
    Generate a quality control (QC) overview using pie charts.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing QC data. Must include columns 'N_labeled', 'N_all', and 'N_checked'.
    figsize : tuple of int, optional
        Size of the figure, by default pieplotsettings["figsize"].
    ncol : int, optional
        Number of columns in the layout, by default pieplotsettings["ncols"].
    radius_big : float, optional
        Radius of the large pie charts, by default pieplotsettings["radius_big"].
    radius_small : float, optional
        Radius of the small pie charts, by default pieplotsettings["radius_small"].
    textsize_big_pies : int, optional
        Font size for the large pie charts, by default pieplotsettings["txt_size_big_pies"].
    textsize_small_pies : int, optional
        Font size for the small pie charts, by default pieplotsettings["txt_size_small_pies"].

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure containing the QC overview pie charts.

    Raises
    ------
    TypeError
        If any of the arguments are not of the expected type.
    """

    # Validate argument types
    if not isinstance(df, pd.DataFrame):
        raise TypeError("Argument 'df' must be of type pandas.DataFrame.")
    if not isinstance(figsize, tuple):
        raise TypeError("Argument 'figsize' must be of type tuple.")
    if not isinstance(ncol, int):
        raise TypeError("Argument 'ncol' must be of type int.")
    if not isinstance(radius_big, (int, float)):
        raise TypeError("Argument 'radius_big' must be of type int or float.")
    if not isinstance(radius_small, (int, float)):
        raise TypeError("Argument 'radius_small' must be of type int or float.")
    if not isinstance(textsize_big_pies, int):
        raise TypeError("Argument 'textsize_big_pies' must be of type int.")
    if not isinstance(textsize_small_pies, int):
        raise TypeError("Argument 'textsize_small_pies' must be of type int.")

    # Define layout
    fig = plt.figure(figsize=figsize)
    fig.tight_layout()

    spec = fig.add_gridspec(4, ncol)
    ax_thl = fig.add_subplot(spec[0, :2])  # top half left
    ax_thr = fig.add_subplot(spec[0, 2:])  # top half right

    # Frequency with all
    plotdf = df
    colors = [label_to_color_map[label] for label in plotdf.index]
    plotdf.plot(
        ax=ax_thl,
        kind="pie",
        y="N_labeled",
        autopct="%1.1f%%",
        legend=False,
        colors=colors,
        radius=radius_big,
        fontsize=textsize_big_pies,
    )
    ax_thl.set_title("Label frequencies")
    ax_thl.set_ylabel("")

    # Outliers comparison
    plotdf = df[
        ~df.index.isin(
            [
                label_def["goodrecord"]["label"],  # TYPO
                label_def["regular_gap"]["label"],  # TYPO
            ]
        )
    ]

    colors = [label_to_color_map[label] for label in plotdf.index]

    if plotdf.empty:
        # No outliers --> full pie with "No QC outliers" in the color of 'ok'
        plotdf = pd.DataFrame(
            data={"N_labeled": [100]}, index=pd.Index(data=["No QC outliers"])
        )
        colors = [label_def["goodrecord"]["color"]]  # TYPO

    plotdf.plot(
        ax=ax_thr,
        kind="pie",
        y="N_labeled",
        autopct="%1.1f%%",
        legend=False,
        colors=colors,
        radius=radius_big,
        fontsize=textsize_big_pies,
    )
    ax_thr.set_title("Outlier specific frequencies")
    ax_thr.set_ylabel("")

    # Performance per check
    plotdf = df[
        ~df.index.isin(
            [
                label_def["goodrecord"]["label"],
                label_def["regular_gap"]["label"],
            ]
        )
    ]

    # Label to QC check name map
    label_too_qcname_map = {val["label"]: key for key, val in label_def.items()}

    i = 0
    for idx, row in plotdf.iterrows():
        # Target a specific axes
        subax = fig.add_subplot(spec[math.floor(i / ncol) + 1, i % ncol])

        # Construct a plot Series
        plotseries = pd.Series(
            {
                label_def["uncheckedrecord"]["label"]: row["N_all"] - row["N_checked"],
                label_def["goodrecord"]["label"]: row["N_checked"] - row["N_labeled"],
                label_def["outlier"]["label"]: row["N_labeled"],
            }
        )
        # Define colors
        colors = [label_to_color_map[label] for label in plotseries.index]

        plotseries.plot(
            ax=subax,
            kind="pie",
            autopct="%1.1f%%",
            legend=False,
            colors=colors,
            radius=radius_small,
            fontsize=textsize_small_pies,
        )

        subax.set_title(f"Effectiveness of {label_too_qcname_map[idx]}")
        subax.set_ylabel("")

        i += 1

    logger.debug("Exiting qc_overview_pies function.")
    return fig
