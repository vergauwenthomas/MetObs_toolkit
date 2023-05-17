#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:26:52 2022

@author: thoverga
"""


import pandas as pd
import math
import numpy as np
import geopandas as gpd
from datetime import datetime

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

from metobs_toolkit.geometry_functions import find_largest_extent
from mpl_toolkits.axes_grid1 import make_axes_locatable


# =============================================================================
# Helpers
# =============================================================================

def map_obstype(obstype, template):
    return template[obstype].to_dict()




def make_cat_colormapper(catlist, cmapname):
    """
    Create a dictionary {cat : color} for a list of categorical values.

    If the colormap has more colors than the catlist, optimal color distance is
    done. If a colormap has less colors than unique categories, the categories are grourped.

    Parameters
    ----------
    catlist : list
        List of categorical values.
    cmapname : str
        Matplotlib.colormaps name.

    Returns
    -------
    colordict : dict
        {cat: color} where the color is a RGBalpha tuple.

    """

    catlist = list(set(catlist)) #get unique categories

    cmap = matplotlib.colormaps[cmapname]

    # check number of colors in the cmap
    if cmap.N < len(catlist):
        print(f'Warning: colormap: {cmapname}, is not well suited to color {len(catlist)} categories. ')
        same_col_n_groups = np.ceil(len(catlist) / cmap.N)

        # group cateogries and color them by group
        colordict = {}
        col_idx = -1
        _cat_index = 0
        for cat in catlist:
            if _cat_index%same_col_n_groups == 0:
                col_idx += 1
            colordict[cat] = cmap(int(col_idx))
            _cat_index += 1
        return colordict

    # check if the colormap can be decreased (and thus increasing the colordistance)
    num_increase = np.floor(cmap.N / len(catlist))

    i = 0
    colordict = {}
    for cat in catlist:
        colordict[cat] = cmap(int(i))
        i = i + num_increase
    return colordict



# =============================================================================
# Plotters
# =============================================================================

def geospatial_plot(
    plotdf,
    variable,
    timeinstance,
    title,
    legend,
    vmin,
    vmax,
    plotsettings,
    categorical_fields,
    static_fields,
    display_name_mapper,
    world_boundaries_map,
):
    # Load default plot settings
    default_settings = plotsettings["spatial_geo"]

    # subset to obstype
    plotdf = plotdf[[variable, "geometry"]]

    # Subset to the stations that have coordinates
    ignored_stations = plotdf[plotdf["geometry"].isnull()]
    plotdf = plotdf[~plotdf["geometry"].isnull()]
    if plotdf.empty:
        # logger.error(f'No coordinate data found, geoplot can not be made. Plotdf: {plotdf}')
        print(f"No coordinate data found, geoplot can not be made. Plotdf: {plotdf}")
        return

    if not ignored_stations.empty:
        # logger.error(f'No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!')
        print(
            f"No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!"
        )

    # make color scheme for field
    if variable in categorical_fields:
        is_categorical = True
        if variable == "lcz":
            # use all available LCZ categories
            use_quantiles = False
        else:
            use_quantiles = True
    else:
        is_categorical = False
        use_quantiles = False

    # if observations extend is contained by default exten, use default else use obs extend
    use_extent = find_largest_extent(
        geodf=gpd.GeoDataFrame(plotdf), extentlist=default_settings["extent"]
    )

    # Style attributes
    if isinstance(title, type(None)):
        if variable in static_fields:
            title = display_name_mapper[variable]
        else:
            dtstring = datetime.strftime(timeinstance, default_settings["fmt"])
            title = display_name_mapper[variable] + " at " + dtstring

    ax = _spatial_plot(
        gdf=plotdf,
        variable=variable,
        legend=legend,
        use_quantiles=use_quantiles,
        is_categorical=is_categorical,
        k_quantiles=default_settings["n_for_categorical"],
        cmap=default_settings["cmap"],
        world_boundaries_map=world_boundaries_map,
        figsize=default_settings["figsize"],
        extent=use_extent,
        title=title,
        vmin=vmin,
        vmax=vmax,
    )
    return ax


def _spatial_plot(
    gdf,
    variable,
    legend,
    use_quantiles,
    is_categorical,
    k_quantiles,
    cmap,
    world_boundaries_map,
    figsize,
    extent,
    title,
    vmin,
    vmax,
):
    # TODO: docstring + beter positionion of the lengends + fix size of pies

    gdf = gpd.GeoDataFrame(gdf)
    fig, ax = plt.subplots(1, 1, figsize=figsize)

    # Make color scheme
    if use_quantiles:
        # maybe better to use evenly spaced intervals rather than quantiles?
        scheme = "equalinterval"
    else:
        scheme = None
        if isinstance(vmin, type(None)) | isinstance(vmax, type(None)):
            vmin = gdf[variable].min()
            vmax = gdf[variable].max()

    if is_categorical:
        legend_kwds = {"loc": "best"}
        vmin = None
        vmax = None
        cax = None
    else:
        legend_kwds = None
        divider = make_axes_locatable(ax)

        cax = divider.append_axes("right", size="5%", pad=0.1)

    # world map as underlayer
    world_boundaries = gpd.read_file(world_boundaries_map)
    world_boundaries.plot(ax=ax)

    # add observations as scatters
    gdf.plot(
        column=variable,
        scheme=scheme,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        # edgecolor='white',
        # linewidth=0.5,
        # scale='NUMBER OF PERSONS KILLED',
        # limits=(8, 24),
        categorical=is_categorical,
        legend=legend,
        # legend_var='scale',
        # legend_kwargs={'loc': 'upper left', 'markeredgecolor': 'black'},
        # legend_values=[2, 1], legend_labels=['2 Fatalities', '1 Fatality'],
        ax=ax,
        cax=cax,
        legend_kwds=legend_kwds,
    )

    # set extent
    ax.set_xlim(left=extent[0], right=extent[2])
    ax.set_ylim(bottom=extent[1], top=extent[3])

    ax.set_title(title)

    return ax


def sorting_function(label_vec, custom_handles, number_of_labels_types=3):
    # TODO: clean this up? rewrite to better code?
    sorted_vec = []            
    # group 1, 2, 3
    for i in range(1, number_of_labels_types+1): # loop over the type of labels
        for l in range(len(label_vec)): #loop over the length of the label_vec
            if label_vec[l] == i:
                sorted_vec.append(l) 
                # makes a vector of same size as label_vec
                # but with the right order of permutations. 
    sorted_handles = [custom_handles[i] for i in sorted_vec] 
    # reordering the custom handles to put 1 at the front    
    
    return sorted_handles


def timeseries_plot(
    mergedf,
    # obstype,
    title,
    xlabel,
    ylabel,
    colorby,
    show_legend,
    show_outliers,
    settings,
    _ax=None #needed for GUI, not recommended use
):

    plot_settings = settings.app["plot_settings"]


    if isinstance(_ax, type(None)):
        # init figure
        fig, ax = plt.subplots(figsize=plot_settings["time_series"]["figsize"])
    else:
        ax=_ax

    # get data ready
    mergedf = mergedf[~mergedf.index.duplicated()]
    init_idx = mergedf.index

    # define different groups (different plotting styles)

    # ok group
    ok_group_label = 'ok'

    # filled value groups
    fill_labels= [ val for val in settings.gap['gaps_fill_info']['label'].values()]
    missing_fill_labels = [ val for val in settings.missing_obs['missing_obs_fill_info']['label'].values()]
    fill_labels.extend(missing_fill_labels)

    # outlier groups
    # Catching with an else


    custom_handles = [] #add legend items to it
    label_vec=[] # add type of label
    

    if colorby == "label":
        # iterate over label groups
        col_mapper = _all_possible_labels_colormapper(settings) # get color mapper

        outl_groups = mergedf.groupby('label')
        legenddict = {}
        for outl_label, groupdf in outl_groups:
            outl_color = col_mapper[outl_label]
            

            # plot data for the 'ok' group
            if outl_label == ok_group_label:  # ok data as lines
                # add init_idx andf fill with nans (to avoid matplotlib interpolation)
                fill_idx = init_idx.to_frame().drop(groupdf.index)
                groupdf = pd.concat([groupdf, fill_idx])
                groupdf = groupdf.drop(columns=["name", "datetime"], errors="ignore")
                groupdf.sort_index()

                plotdf = groupdf.reset_index().pivot(
                    index="datetime", columns="name", values='value'
                )  # long to wide

                ax = plotdf.plot(
                    kind="line",
                    color=outl_color,
                    ax=ax,
                    legend=False,
                    zorder=plot_settings["time_series"]["linezorder"],
                    linewidth=plot_settings["time_series"]["linewidth"],
                )

                # add legend handl
                custom_handles.append(
                    Line2D([0], [0], color=outl_color, label="ok", lw=4))
                label_vec.append(1)



            # plot filled data
            elif outl_label in  fill_labels:  # fill gaps as dashed lines

                fill_idx = init_idx.to_frame().drop(groupdf.index)
                groupdf = pd.concat([groupdf, fill_idx])
                groupdf = groupdf.drop(columns=["name", "datetime"], errors="ignore")
                groupdf.sort_index()
                plotdf = groupdf.reset_index().pivot(
                    index="datetime", columns="name", values='value'
                )  # long to wide

                ax = plotdf.plot(
                    kind="line",
                    style="--",
                    color=outl_color,
                    ax=ax,
                    legend=False,
                    zorder=plot_settings['time_series']["dashedzorder"],
                    linewidth=plot_settings['time_series']["linewidth"],
                )

                # add legend handle
                custom_handles.append(
                    Line2D([0],[0],
                        color=outl_color,
                        label=f"filled value ({outl_label})",
                        lw=1,
                        linestyle="--",)
                    )
                label_vec.append(2)

            else:  # outliers as scatters
                plotdf = groupdf['value']
                plotdf.index = plotdf.index.droplevel("name")
                plotdf = plotdf.reset_index()
                ax = plotdf.plot(
                    kind="scatter",
                    x="datetime",
                    y='value',
                    ax=ax,
                    color=outl_color,
                    legend=False,
                    zorder=plot_settings["time_series"]["scatterzorder"],
                    s=plot_settings["time_series"]["scattersize"],
                )

                # add legend handle
                custom_handles.append(
                    Line2D([0],[0], marker="o", color="w",
                        markerfacecolor=outl_color,
                        label=outl_label,
                        lw=1,)
                    )
                label_vec.append(3)
            legenddict[outl_label] = outl_color

        # make legend
        if show_legend:
            # TODO: sort items
            # sort legend items
            custom_handles = sorting_function(label_vec, custom_handles)
            #ax.legend(handles=custom_handles)
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.95])
            ax.legend(handles=custom_handles, loc='upper center',
                bbox_to_anchor=(0.5, -0.2),
                fancybox=True, shadow=True,
                ncol=plot_settings["time_series"]["legend_n_columns"])
            

    elif colorby == "name":
        plotdf = mergedf.reset_index().pivot(
            index="datetime", columns="name", values='value'
        )
        #ax = plotdf.plot(kind="line", legend=show_legend, ax=ax)
        ax = plotdf.plot(kind="line", legend=False, ax=ax)
        if show_legend == True:
            box = ax.get_position()
            ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.95])
            ax.legend(plotdf.columns.values.tolist(), loc='upper center',
                bbox_to_anchor=(0.5, -0.2),
                fancybox=True, shadow=True,
                ncol=plot_settings["time_series"]["legend_n_columns"])

    # Set title
    ax.set_title(title)
    # ax.legend().set_title('')

    # Set x and y labels
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    return ax

def diurnal_plot(diurnaldf, errorbandsdf, title, tzstr, plot_settings,
                 colorby, lcz_dict, data_template, obstype,
                 show_zero_horizontal=False):
    # init figure
    fig, ax = plt.subplots(figsize=plot_settings["figsize"])

    if colorby == 'lcz':
        present_lczs = list(set(lcz_dict.values()))

        # select colormap
        if len(present_lczs) <= plot_settings['n_cat_max']:
            cmap = plot_settings['cmap_categorical']
        else:
            cmap = plot_settings['cmap_continious']

        colordict = make_cat_colormapper(catlist=present_lczs,
                                         cmapname=cmap)


        # Make plot per lcz
        custom_handles = []

        for lcz_cat  in present_lczs:
            stations = [sta for sta in diurnaldf.columns if lcz_dict[sta] == lcz_cat]
            diurnaldf[stations].plot(ax=ax, title=title, color=colordict[lcz_cat], legend=False)

            # add legend item
            custom_handles.append(
                Line2D([0], [0], color=colordict[lcz_cat], label=lcz_cat, lw=4)
            )
            
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2,
             box.width, box.height * 0.95])
        ax.legend(handles=custom_handles, loc='upper center',
            bbox_to_anchor=(0.5, -0.1),
            fancybox=True, shadow=True,
            ncol=plot_settings["legend_n_columns"])




    if colorby == 'name':
        # which colormap to use:
        if diurnaldf.shape[1] <= plot_settings['n_cat_max']:
            cmap = plot_settings['cmap_categorical']
        else:
            cmap = plot_settings['cmap_continious']

        diurnaldf.plot(ax=ax, title=title, legend=False, cmap=cmap)
        box = ax.get_position()
        ax.set_position([box.x0, box.y0 + box.height * 0.2,
                 box.width, box.height * 0.95])
        ax.legend(diurnaldf.columns.values.tolist(), loc='upper center',
                bbox_to_anchor=(0.5, -0.1),
                fancybox=True, shadow=True,
                ncol=plot_settings["legend_n_columns"])



    if not errorbandsdf is None:
        # Extract colorscheme from the plot
        col_sheme = {line.get_label(): line.get_color() for line in ax.get_lines()}

        for sta in errorbandsdf.columns:
            ax.fill_between(errorbandsdf.index,
                             diurnaldf[sta] - errorbandsdf[sta],
                             diurnaldf[sta] + errorbandsdf[sta],
                             alpha=plot_settings['alpha_error_bands'],
                             color=col_sheme[sta],
                             )

    if show_zero_horizontal:
        ax.axhline(y=0., color='black', linestyle='--')
    # Styling attributes



    ax.set_ylabel(map_obstype(obstype, data_template)['orig_name'])
    ax.xaxis.set_major_formatter('{x:.0f} h')
    ax.set_xlabel(f'Hours (timezone: {tzstr})')


    plt.show()

def _make_pie_from_freqs(
    freq_dict, colormapper, ax, plot_settings, radius, labelsize=10
):
    # To dataframe
    stats = pd.Series(freq_dict, name="freq").to_frame()

    # make color mapper
    stats["color"] = stats.index.map(colormapper)

    if (stats["freq"] == 0.0).all():
        print("No occurences in sample.")
        # add a 100% no occurences to it, so it can be plotted
        no_oc_df = pd.DataFrame(
            index=["No occurences"],
            data={"freq": [100.0], "color": [plot_settings["color_mapper"]["ok"]]},
        )
        stats = pd.concat([stats, no_oc_df])

    # Make pie
    patches, text = ax.pie(
        stats["freq"],
        colors=stats["color"],
        radius=radius,
        labels=[
            f"{l}, {s:0.1f}%"
            for l, s in zip(stats.index.to_list(), stats["freq"].to_list())
        ],
        textprops={"fontsize": labelsize},
    )

    return ax


def _outl_value_to_colormapper(plot_settings, qc_check_info):
    """Make color mapper for the outlier LABELVALUES to colors."""
    color_defenitions = plot_settings["color_mapper"]
    outl_name_mapper = {val["outlier_flag"]: key for key, val in qc_check_info.items()}
    outl_col_mapper = {
        outl_type: color_defenitions[outl_name_mapper[outl_type]]
        for outl_type in outl_name_mapper.keys()
    }
    return outl_col_mapper


def _all_possible_labels_colormapper(settings):
    """Make color mapper for all LABELVALUES to colors."""


    plot_settings = settings.app["plot_settings"]
    gap_settings = settings.gap
    qc_info_settings = settings.qc["qc_checks_info"]
    missing_obs_settings = settings.missing_obs['missing_obs_fill_info']


    color_defenitions = plot_settings["color_mapper"]

    mapper = dict()

    # get QC outlier labels

    outl_col_mapper = _outl_value_to_colormapper(
        plot_settings=plot_settings, qc_check_info=qc_info_settings
    )
    mapper.update(outl_col_mapper)
    # get 'ok' and 'not checked'
    mapper["ok"] = color_defenitions["ok"]
    mapper["not checked"] = color_defenitions["not checked"]

    # update gap and missing timestamp labels
    mapper[gap_settings["gaps_info"]["gap"]["outlier_flag"]] = color_defenitions["gap"]
    mapper[
        gap_settings["gaps_info"]["missing_timestamp"]["outlier_flag"]
    ] = color_defenitions["missing_timestamp"]

    # add fill for gaps
    for method, label in gap_settings["gaps_fill_info"]["label"].items():
        mapper[label] = color_defenitions[method]

    # add fill for missing
    for method, label in missing_obs_settings['label'].items():
        mapper[label] = color_defenitions[method]





    return mapper


def qc_stats_pie(
    final_stats, outlier_stats, specific_stats, plot_settings, qc_check_info
):
    # restore rcParams
    plt.rcParams = plt.rcParamsDefault

    # Specify rcParams

    # axestitl
    plt.rcParams["axes.titlelocation"] = "center"
    plt.rcParams["axes.titlesize"] = 10
    plt.rcParams["axes.titleweight"] = 2
    plt.rcParams["axes.titlecolor"] = "black"

    # label size
    textsize_big_pies = 10
    textsize_small_pies = 7

    color_defenitions = plot_settings["color_mapper"]
    # Define layout

    fig = plt.figure(figsize=plot_settings["pie_charts"]["figsize"])
    fig.tight_layout()
    spec = fig.add_gridspec(4, 4, wspace=10)

    ax_thl = fig.add_subplot(spec[0, :2])  # top half left
    ax_thr = fig.add_subplot(spec[0, 2:])  # top half right

    # 1. Make the finale label pieplot
    # make color mapper
    final_col_mapper = {
        "ok": color_defenitions["ok"],
        "QC outliers": color_defenitions["outlier"],
        "missing (gaps)": color_defenitions["gap"],
        "missing (individual)": color_defenitions["missing_timestamp"],
    }

    _make_pie_from_freqs(
        freq_dict=final_stats,
        colormapper=final_col_mapper,
        ax=ax_thl,
        plot_settings=plot_settings,
        radius=plot_settings["pie_charts"]["radius_big"],
        labelsize=textsize_big_pies,
    )

    ax_thl.set_title(
        label="Final label frequencies",
        y=(plot_settings["pie_charts"]["radius_big"] / 2) * 1.4,
        fontweight="bold",
    )

    # 2. Make QC overview pie
    # make color mapper
    outl_col_mapper = _outl_value_to_colormapper(plot_settings, qc_check_info)

    _make_pie_from_freqs(
        freq_dict=outlier_stats,
        colormapper=outl_col_mapper,
        ax=ax_thr,
        plot_settings=plot_settings,
        radius=plot_settings["pie_charts"]["radius_big"],
        labelsize=textsize_big_pies,
    )

    ax_thr.set_title(
        label="Outlier performance",
        y=(plot_settings["pie_charts"]["radius_big"] / 2) * 1.4,
        fontweight="bold",
    )

    # 3. Make a specific pie for each indvidual QC + gap + missing
    plt.rcParams["axes.titley"] = plot_settings["pie_charts"]["radius_small"] / 2
    # make color mapper
    spec_col_mapper = {
        "ok": color_defenitions["ok"],
        "not checked": color_defenitions["not checked"],
        "outlier": color_defenitions["outlier"],
        "gap": color_defenitions["gap"],
        "missing timestamp": color_defenitions["missing_timestamp"],
    }

    specific_df = pd.DataFrame(specific_stats)

    ncol = 4
    nrow = 4

    # create list of axes for the small pies
    axlist = []
    i = 0
    for checkname in specific_stats:
        ax = fig.add_subplot(
            spec[
                math.floor(i / ncol) + 1 : math.floor(i / ncol) + 2,
                i % nrow : i % nrow + 1,
            ]
        )

        # specific style formatting
        ax.set_title(
            label=checkname.replace("_", " "),
            y=plot_settings["pie_charts"]["radius_small"] / 2,
            fontweight="bold",
        )
        ax.yaxis.set_visible(False)  # ignore the default pandas title

        axlist.append(ax)
        i += 1

    # Make pie plots
    specific_df.plot.pie(
        subplots=True,
        labels=specific_df.index,
        legend=False,
        autopct="%1.1f%%",
        title=None,
        radius=plot_settings["pie_charts"]["radius_small"],
        textprops={"fontsize": textsize_small_pies},
        ax=axlist,
    )

    # Specific styling setings per pie
    for ax in axlist:
        # specific style formatting
        ax.yaxis.set_visible(False)  # ignore the default pandas title

    fig.subplots_adjust(hspace=0.7)
    plt.show()

    return
