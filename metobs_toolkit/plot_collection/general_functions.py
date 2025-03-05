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


def set_legend(ax: plt.Axes) -> plt.Axes:
    handles, labels = ax.get_legend_handles_labels()
    by_label = dict(zip(labels, handles))  # this drops duplicates based on labels
    ax.legend(by_label.values(), by_label.keys())
    return ax


# ------------------------------------------
#    Coloring
# ------------------------------------------


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
