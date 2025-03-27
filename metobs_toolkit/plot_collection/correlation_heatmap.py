# FIXME


import matplotlib.pyplot as plt


def heatmap_plot(cor_dict, title, heatmap_settings):
    """Make a heatmap plot (i.g. matrix visualisation).

    Parameters
    ----------
    cor_dict : dict
        A dictionary of the correlations to plot.
    title : str
        The title of the figure.
    heatmap_settings : dict
        The plot settings for heatmaps.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The axes of the plot.

    """
    # make heatmap of cor
    fig, ax = plt.subplots(figsize=heatmap_settings["figsize"])
    im = ax.imshow(
        cor_dict["cor matrix"],
        interpolation="nearest",
        vmin=heatmap_settings["vmin"],
        vmax=heatmap_settings["vmax"],
        cmap=heatmap_settings["cmap"],
    )

    fig.colorbar(im, orientation="vertical", fraction=0.05)

    # Loop over data dimensions and create text annotations
    for i in range(len(cor_dict["cor matrix"].columns)):
        for j in range(len(cor_dict["cor matrix"].index)):
            ax.text(
                j,
                i,
                cor_dict["combined matrix"].to_numpy()[i, j],
                ha="center",
                va="center",
                color="black",
            )

    # styling
    # Show all ticks and label them with the dataframe column name
    ax.set_xticks(
        ticks=list(range(cor_dict["cor matrix"].shape[1])),
        labels=cor_dict["cor matrix"].columns.to_list(),
        rotation=heatmap_settings["x_tick_rot"],
    )

    ax.set_yticks(
        ticks=list(range(cor_dict["cor matrix"].shape[0])),
        labels=cor_dict["cor matrix"].index.to_list(),
        rotation=heatmap_settings["y_tick_rot"],
    )

    ax.set_title(title)

    return ax
