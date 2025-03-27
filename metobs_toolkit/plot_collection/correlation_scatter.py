import pandas as pd

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

# FIXME


def correlation_scatter(
    full_cor_dict, groupby_labels, obstypes, title, cor_scatter_settings
):
    """Plot the correlation variation as a scatterplot.

    The statistical significance is indicated by the scattertype.

    Parameters
    ----------
    full_cor_dict : dict
        A dictionary containing the 'cor matrix', and 'significance matrix'
        keys and corresponding matrices.
    groupby_labels : str or list
        The groupdefenition that is used for the x-axes label.
    obstypes : str
        The observation type to plot the correlations of.
    title : str
        The title of the figure.
    cor_scatter_settings : dict
        The specific plot settings for the correlation scatter plot.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The axes of the plot.

    """
    # combine all correlation matrices to one with multiindex
    comb_cor_df = pd.DataFrame()
    comb_p_df = pd.DataFrame()
    for key, subcordict in full_cor_dict.items():

        # if mulitple groupby are given, key is tuple --> conv to string
        if isinstance(key, tuple):
            key = str(key)
        # corelations
        subdf_cor = subcordict["cor matrix"]
        # make multi index df
        subdf_cor["group"] = key
        subdf_cor.index.name = "categories"
        subdf_cor = subdf_cor[subdf_cor.index.isin(obstypes)]
        subdf_cor = subdf_cor.reset_index().set_index(["group", "categories"])
        comb_cor_df = pd.concat([comb_cor_df, subdf_cor])

        # p values
        subdf_p = subcordict["significance matrix"]
        # make multi index df
        subdf_p["group"] = key
        subdf_p.index.name = "categories"
        subdf_p = subdf_p[subdf_p.index.isin(obstypes)]
        subdf_p = subdf_p.reset_index().set_index(["group", "categories"])
        comb_p_df = pd.concat([comb_p_df, subdf_p])

    # create plotdf structure
    plot_cor_df = comb_cor_df.unstack()
    plot_cor_df.columns = [f"{col[0]} - {col[1]}" for col in plot_cor_df.columns]
    plot_p_df = comb_p_df.unstack()
    plot_p_df.columns = [f"{col[0]} - {col[1]}" for col in plot_p_df.columns]

    # Get columns without variation (these will not be plotted)
    const_cols = plot_cor_df.columns[plot_cor_df.nunique() <= 1]
    logger.warning(
        f" The following correlations are constant for all groups and will not be included in the plot: {const_cols}"
    )

    # Subset to the columns that has to be plotted
    plot_cor_df = plot_cor_df.drop(columns=const_cols)
    plot_p_df = plot_p_df.drop(columns=const_cols)

    # make a colormap for the left over correlations
    col_mapper = make_cat_colormapper(
        catlist=plot_cor_df.columns.to_list(), cmapname=cor_scatter_settings["cmap"]
    )

    # make figure
    fig, ax = plt.subplots(figsize=cor_scatter_settings["figsize"])

    # add the zero line
    ax.axhline(y=0.0, linestyle="--", linewidth=1, color="black")

    # Define p value bins
    p_bins = cor_scatter_settings["p_bins"]  # [0, .001, 0.01, 0.05, 999]
    bins_markers = cor_scatter_settings["bins_markers"]  # ['*', 's', '^', 'x']

    # # iterate over the different corelations to plot
    custom_handles = []
    for cor_name in plot_cor_df.columns:
        to_scatter = plot_cor_df[[cor_name]]

        # convert p values to markers
        to_scatter["p-value"] = plot_p_df[cor_name]
        to_scatter["markers"] = pd.cut(
            x=to_scatter["p-value"], bins=p_bins, labels=bins_markers
        )
        to_scatter = to_scatter.reset_index()

        # plot per scatter group
        scatter_groups = to_scatter.groupby("markers", observed=True)
        for marker, markergroup in scatter_groups:
            markergroup.plot(
                x="group",
                y=cor_name,
                kind="scatter",
                ax=ax,
                s=cor_scatter_settings["scatter_size"],
                edgecolors=cor_scatter_settings["scatter_edge_col"],
                linewidth=cor_scatter_settings["scatter_edge_line_width"],
                color=col_mapper[cor_name],
                marker=marker,
                ylim=(cor_scatter_settings["ymin"], cor_scatter_settings["ymax"]),
            )

        # add legend handl for the colors
        custom_handles.append(
            Line2D([0], [0], color=col_mapper[cor_name], label=cor_name, lw=4)
        )

    # add legend handl for the scatter types
    marker_def = list(zip(p_bins[1:], bins_markers))
    for p_edge, mark in marker_def:
        custom_handles.append(
            Line2D(
                [0],
                [0],
                marker=mark,
                color="black",
                markerfacecolor="w",
                label=f"p < {p_edge}",
                lw=1,
            )
        )

    # format legend
    box = ax.get_position()
    ax.set_position([box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.85])
    ax.legend(
        handles=custom_handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.1),
        fancybox=True,
        shadow=True,
        prop={"size": cor_scatter_settings["legend_text_size"]},
        ncol=cor_scatter_settings["legend_ncols"],
    )

    # styling attributes
    ax.set_ylabel("Pearson correlation")
    ax.set_xlabel(f"Groups of {groupby_labels}")
    ax.set_title(title)

    return ax
