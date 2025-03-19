import math
import pandas as pd
import matplotlib.pyplot as plt


import metobs_toolkit.settings_files.default_formats_settings as defaults

plotsettings = defaults.plot_settings["pie_charts"]


def qc_overview_pies(
    df,
    figsize=plotsettings["figsize"],
    ncol=plotsettings["ncols"],
    radius_big=plotsettings["radius_big"],
    radius_small=plotsettings["radius_small"],
    textsize_big_pies=plotsettings["txt_size_big_pies"],
    textsize_small_pies=plotsettings["txt_size_small_pies"],
):

    # Specify rcParams
    # plt.rcParams["axes.titlelocation"] = "center"
    # plt.rcParams["axes.titlesize"] = 10
    # plt.rcParams["axes.titleweight"] = 2
    # plt.rcParams["axes.titlecolor"] = "black"

    # Define layout
    fig = plt.figure(figsize=figsize)
    fig.tight_layout()

    spec = fig.add_gridspec(
        4,
        ncol,
    )
    ax_thl = fig.add_subplot(spec[0, :2])  # top half left
    ax_thr = fig.add_subplot(spec[0, 2:])  # top half right

    # freq with all
    plotdf = df
    colors = [defaults.label_to_color_map[label] for label in plotdf.index]
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

    # outliers comparison
    plotdf = df[
        ~df.index.isin(
            [
                defaults.label_def["goodrecord"]["label"],
                defaults.label_def["regular_gap"]["label"],
            ]
        )
    ]

    colors = [defaults.label_to_color_map[label] for label in plotdf.index]

    if plotdf.empty:
        # no outliers --> full pie with "No QX outliers" in the color of 'ok
        plotdf = pd.DataFrame(
            data={"N_labeled": [100]}, index=pd.Index(data=["No QC outliers"])
        )
        colors = [defaults.label_def["goodrecord"]["color"]]

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
    ax_thr.set_title("outlier specific frequencies")
    ax_thr.set_ylabel("")

    # performance per check
    plotdf = df[
        ~df.index.isin(
            [
                defaults.label_def["goodrecord"]["label"],
                defaults.label_def["regular_gap"]["label"],
            ]
        )
    ]

    # label to qccheckname map
    labl_to_qcname_map = {val["label"]: key for key, val in defaults.label_def.items()}

    i = 0
    for idx, row in plotdf.iterrows():
        # target a specific axes
        subax = fig.add_subplot(spec[math.floor(i / ncol) + 1, i % ncol])

        # construct a plot Series
        plotseries = pd.Series(
            {
                defaults.label_def["uncheckedrecord"]["label"]: row["N_all"]
                - row["N_checked"],
                defaults.label_def["goodrecord"]["label"]: row["N_checked"]
                - row["N_labeled"],
                defaults.label_def["outlier"]["label"]: row["N_labeled"],
            }
        )
        # define colors
        colors = [defaults.label_to_color_map[label] for label in plotseries.index]

        plotseries.plot(
            ax=subax,
            kind="pie",
            autopct="%1.1f%%",
            legend=False,
            colors=colors,
            radius=radius_small,
            fontsize=textsize_small_pies,
        )

        subax.set_title(f"Effectiveness of {labl_to_qcname_map[idx]}")
        subax.set_ylabel("")

        i += 1
    return fig
