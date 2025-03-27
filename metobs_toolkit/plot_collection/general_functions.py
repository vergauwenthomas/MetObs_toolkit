import logging
from typing import Literal, Tuple

import matplotlib.pyplot as plt
import matplotlib

# ------------------------------------------
#    Figure/axes layouts
# ------------------------------------------


def create_axes(figsize: Tuple[int, int] = (15, 5), **kwargs) -> plt.Axes:

    _fig, ax = plt.subplots(figsize=figsize, **kwargs)
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


def format_datetime_axes(ax: plt.Axes, set_diurnal_format: bool = False) -> plt.Axes:
    """Set the x-axes to autodateformat. Optionally set diurnal format."""
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


def _drop_cur_legend(ax):
    if ax.get_legend() is None:
        return  # nothing to drop
    else:
        ax.get_legend().remove()


def _get_unique_handles_and_labels(ax):
    handles, labels = ax.get_legend_handles_labels()
    unique_labels = {}
    for handle, label in zip(handles, labels):
        if label not in unique_labels:
            unique_labels[label] = handle

    labels = list(unique_labels.keys())
    handles = list(unique_labels.values())
    return handles, labels


def _creat_main_legend_items(handles, labels):
    # Get all labels related to records
    recorditems = [(lab, hand) for lab, hand in zip(labels, handles) if "@" not in lab]
    recordlabels = [item[0] for item in recorditems]

    # Geta all modeldata labels that do not have a similar recordlabel
    modelitems_to_show = [
        (lab, hand)
        for lab, hand in zip(labels, handles)
        if lab.split("@")[-1] not in recordlabels
    ]

    return recorditems + modelitems_to_show


def set_legend(ax: plt.Axes, ncols: int = 8) -> plt.Axes:
    # remove current legends from the axes
    _drop_cur_legend(ax)

    # drop duplicated labels and get labels/handles
    handles, labels = _get_unique_handles_and_labels(ax)

    # Construct maingroup labels (the collection of labels related to colors to plot)
    main_legenditems = _creat_main_legend_items(handles, labels)

    recordlabels = [
        lab[0] for lab in main_legenditems if "@" not in lab[0]
    ]  # record labe
    # Test if there are similar labels of modeldata and records
    similarity = any(
        [
            True if (lab.split("@")[-1] in recordlabels) & ("@" in lab) else False
            for lab in labels
        ]
    )
    # if there is similarity (reference of the same station), create a secondary legend
    if similarity:
        # get the modeldata info -> extract from the first label with '@' in
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

    # NOTE: since the main legend, can change the figure shape add the secondary axis first !!! It
    # does not work if order is switched.

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
        main_legend = ax.legend(
            main_handles,
            main_labels,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.15),
            ncol=ncols,
        )
    else:
        # create legend
        main_legend = ax.legend(main_handles, main_labels)

    return ax


# ------------------------------------------
#    Coloring
# ------------------------------------------


def create_categorical_color_map(catlist: list, cmapname: str = "tab20") -> dict:
    unique_elements = list(set(catlist))
    unique_elements = sorted(
        unique_elements
    )  # sort alphabetically, so colorscheme is equal in workflow
    num_unique_elements = len(unique_elements)

    cmap = plt.get_cmap(cmapname)
    colors = [cmap(i / num_unique_elements) for i in range(num_unique_elements)]

    hex_colors = [matplotlib.colors.rgb2hex(color) for color in colors]

    color_map = {
        element: hex_colors[i % num_unique_elements]
        for i, element in enumerate(unique_elements)
    }

    return color_map
