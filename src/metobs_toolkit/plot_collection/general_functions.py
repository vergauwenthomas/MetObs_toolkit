import logging
import locale
from typing import Tuple

import matplotlib
import matplotlib.pyplot as plt

# Set up logging
from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")

# ------------------------------------------
#    Figure/axes layouts
# ------------------------------------------


@log_entry
def create_axes(figsize: Tuple[int, int] = (15, 5), **kwargs) -> plt.Axes:
    """
    Create a matplotlib Axes object with a specified figure size.

    Parameters
    ----------
    figsize : tuple of int, optional
        Size of the figure, by default (15, 5).
    **kwargs
        Additional keyword arguments passed to plt.subplots.

    Returns
    -------
    matplotlib.axes.Axes
        The created Axes object.
    """
    _fig, ax = plt.subplots(figsize=figsize, **kwargs)
    return ax


# ------------------------------------------
#    Styling
# ------------------------------------------
# These functions seem trivial, but by defining them in one place,
# it is easier to fix styling issues on all levels (sensor, station, dataset)
@log_entry
def set_title(ax: plt.Axes, titlestr: str) -> plt.Axes:
    """
    Set the title of the axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to set the title for.
    titlestr : str
        The title string.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the title set.
    """
    ax.set_title(titlestr)
    return ax


@log_entry
def set_ylabel(ax: plt.Axes, ylabel: str) -> plt.Axes:
    """
    Set the y-axis label of the axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to set the y-label for.
    ylabel : str
        The y-axis label string.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the y-label set.
    """
    ax.set_ylabel(ylabel)
    return ax


@log_entry
def set_xlabel(ax: plt.Axes, xlabel: str) -> plt.Axes:
    """
    Set the x-axis label of the axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to set the x-label for.
    xlabel : str
        The x-axis label string.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the x-label set.
    """
    ax.set_xlabel(xlabel)
    return ax


@log_entry
def format_datetime_axes(ax: plt.Axes, set_diurnal_format: bool = False) -> plt.Axes:
    """
    Format the x-axis of the axes for datetime values.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to format.
    set_diurnal_format : bool, optional
        If True, set the format to hours and minutes ("%H:%M"), otherwise use the default date formatter.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with formatted x-axis.
    """
    # Set the locale to English (United States) for date formatting
    locale.setlocale(locale.LC_TIME, "en_US.UTF-8")

    xtick_locator = matplotlib.dates.AutoDateLocator()
    if set_diurnal_format:
        xtick_formatter = matplotlib.dates.DateFormatter("%H:%M")
    else:
        xtick_formatter = matplotlib.dates.AutoDateFormatter(xtick_locator)

    ax.xaxis.set_major_locator(xtick_locator)
    ax.xaxis.set_major_formatter(xtick_formatter)
    return ax


# ------------------------------------------
#    Legend handling
# ------------------------------------------


def _drop_cur_legend(ax: plt.Axes) -> None:
    """Remove the current legend from the axes if present."""
    if ax.get_legend() is None:
        return  # nothing to drop
    else:
        ax.get_legend().remove()


def _get_unique_handles_and_labels(ax: plt.Axes):
    """
    Get unique legend handles and labels from the axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to extract handles and labels from.

    Returns
    -------
    tuple
        Tuple of (handles, labels) with unique labels.
    """
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = {}
    for handle, label in zip(handles, labels):
        if label not in unique_labels:
            unique_labels[label] = handle

    labels = list(unique_labels.keys())
    handles = list(unique_labels.values())
    return handles, labels


def _create_main_legend_items(handles, labels):
    """
    Create main legend items, grouping record and model data.

    Parameters
    ----------
    handles : list
        List of legend handles.
    labels : list
        List of legend labels.

    Returns
    -------
    list
        List of tuples (label, handle) for main legend items.
    """
    # Get all labels related to records
    recorditems = [(lab, hand) for lab, hand in zip(labels, handles) if "@" not in lab]
    recordlabels = [item[0] for item in recorditems]

    # Get all model data labels that do not have a similar record label
    modelitems_to_show = [
        (lab, hand)
        for lab, hand in zip(labels, handles)
        if lab.split("@")[-1] not in recordlabels
    ]

    return recorditems + modelitems_to_show


@log_entry
def set_legend(ax: plt.Axes, ncols: int = 8) -> plt.Axes:
    """
    Set the legend for the axes, handling duplicate labels and model/observation distinction.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        The axes to set the legend for.
    ncols : int, optional
        Number of columns for the legend, by default 8.

    Returns
    -------
    matplotlib.axes.Axes
        The axes with the legend set.
    """
    # remove current legends from the axes
    _drop_cur_legend(ax)

    # drop duplicated labels and get labels/handles
    handles, labels = _get_unique_handles_and_labels(ax)

    # Construct main group labels (the collection of labels related to colors to plot)
    main_legenditems = _create_main_legend_items(handles, labels)

    recordlabels = [
        lab[0] for lab in main_legenditems if "@" not in lab[0]
    ]  # record label
    # Test if there are similar labels of model data and records
    similarity = any(
        [
            True if (lab.split("@")[-1] in recordlabels) & ("@" in lab) else False
            for lab in labels
        ]
    )
    # if there is similarity (reference of the same station), create a secondary legend
    if similarity:
        # get the model data info -> extract from the first label with '@' in
        for lab, hand in zip(labels, handles):
            if "@" in lab:
                modelinfo = lab.split("@")[0]
                modellinestyle = hand.get_linestyle()
                break
        secondary_items = [
            ("Observations", plt.Line2D([0], [0], color="black", linestyle="-")),
            (
                f"{modelinfo}",
                plt.Line2D([0], [0], color="black", linestyle=modellinestyle),
            ),
        ]

    # NOTE: since the main legend can change the figure shape, add the secondary axis first !!!
    # It does not work if order is switched.

    # Plot secondary legend
    if similarity:
        main_labels = [item[0] for item in secondary_items]
        handles = [item[1] for item in secondary_items]
        extra_legend = ax.legend(
            handles, main_labels
        )  # do not use ax.legend here, that will overwrite the main legend
        ax.add_artist(extra_legend)

    # add the main axes
    main_labels = [item[0] for item in main_legenditems]
    main_handles = [item[1] for item in main_legenditems]

    if len(main_labels) > 8:  # Adjust the threshold as needed
        # adjust figure
        fig = ax.get_figure()
        fig.subplots_adjust(
            bottom=0.2
        )  # Adjust the bottom margin to make space for the legend
        # create legend
        _ = ax.legend(
            main_handles,
            main_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
            ncol=ncols,
        )
    else:
        # create legend
        _ = ax.legend(main_handles, main_labels)

    return ax


# ------------------------------------------
#    Coloring
# ------------------------------------------


@log_entry
def create_categorical_color_map(catlist: list, cmapname: str = "tab20") -> dict:
    """
    Create a categorical color map for a list of categories.

    Parameters
    ----------
    catlist : list
        List of categories to assign colors to.
    cmapname : str, optional
        Name of the matplotlib colormap to use, by default "tab20".

    Returns
    -------
    dict
        Dictionary mapping each category to a hex color string.
    """
    unique_elements = list(set(catlist))
    unique_elements = sorted(
        unique_elements
    )  # sort alphabetically, so color scheme is equal in workflow
    num_unique_elements = len(unique_elements)

    cmap = plt.get_cmap(cmapname)
    colors = [cmap(i / num_unique_elements) for i in range(num_unique_elements)]

    hex_colors = [matplotlib.colors.rgb2hex(color) for color in colors]

    color_map = {
        element: hex_colors[i % num_unique_elements]
        for i, element in enumerate(unique_elements)
    }

    return color_map


@log_entry
def create_linestyle_map(
    catlist: list,
    linestyles: list = ["-", "--", "-.", ":"],
    user_linestyle_defs: dict = {},
) -> dict:
    """
    Create a linestyle map for a list of categories.

    Parameters
    ----------
    catlist : list
        List of categories to assign linestyles to.
    linestyles : list, optional
        List of matplotlib linestyle strings to cycle through, by default ['-', '--', '-.', ':'].
    user_linestyle_defs : dict, optional
        Dictionary of user-defined linestyle overrides. Keys should match categories
        from catlist, values should be matplotlib linestyle strings, by default {}.

    Returns
    -------
    dict
        Dictionary mapping each unique category to its corresponding linestyle string.
    """
    unique_elements = list(set(catlist))
    unique_elements = sorted(
        unique_elements
    )  # sort alphabetically, so color scheme is equal in workflow

    linestyle_map = {
        val: linestyles[i % len(linestyles)] for i, val in enumerate(unique_elements)
    }

    # Force the user defs
    linestyle_map.update(user_linestyle_defs)

    return linestyle_map
