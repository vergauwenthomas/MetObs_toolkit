#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 21 11:26:52 2022

@author: thoverga
"""

import sys
import pandas as pd
import math
import numpy as np
import geopandas as gpd
from datetime import datetime
import logging

import matplotlib
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.dates as mdates
from matplotlib.collections import LineCollection

import branca
import branca.colormap as brcm

import cartopy.crs as ccrs
import cartopy.feature as cfeature

import geemap.foliumap as foliumap
import folium
from folium import plugins as folium_plugins

from metobs_toolkit.geometry_functions import find_plot_extent
from mpl_toolkits.axes_grid1 import make_axes_locatable

from metobs_toolkit.landcover_functions import get_ee_obj
from metobs_toolkit.df_helpers import xs_save

logger = logging.getLogger(__name__)


def folium_plot(
    mapinfo,
    band,
    vis_params,
    labelnames,
    layername,
    basemap="SATELLITE",
    legendname=None,
    legendpos="bottomleft",
):
    """Make an interactive folium plot of an Image."""
    # get the ee.Image
    im = get_ee_obj(mapinfo, band)

    # make plot
    MAP = foliumap.Map()
    if basemap:
        MAP.add_basemap(basemap)
    MAP.add_layer(im, vis_params, layername)
    if legendname:
        MAP.add_legend(
            title=legendname,
            labels=labelnames,
            colors=vis_params.get("palette"),
            position=legendpos,
        )

    return MAP


def add_stations_to_folium_map(Map, metadf):
    """Add stations as markers to the folium map."""
    points = metadf["geometry"].to_crs("epsg:4326")
    for station, point in points.items():
        folium.Marker(
            location=[point.y, point.x], fill_color="#43d9de", popup=station, radius=8
        ).add_to(Map)

    return Map


# =============================================================================
# Helpers
# =============================================================================
def _get_init_mapcenter(gdf):
    center = gdf.dissolve().centroid.iloc[0]
    return [center.y, center.x]


def map_obstype(obstype, template):
    """Convert default obstype to the user-specific obstype."""
    return template[obstype].to_dict()


def make_cat_colormapper(catlist, cmapname):
    """Create a dictionary {cat : color} for a list of categorical values.

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
    catlist = list(set(catlist))  # get unique categories

    cmap = matplotlib.colormaps[cmapname]

    # check number of colors in the cmap
    if cmap.N < len(catlist):
        logger.warning(
            f"colormap: {cmapname}, is not well suited to color {len(catlist)} categories."
        )
        same_col_n_groups = np.ceil(len(catlist) / cmap.N)

        # group cateogries and color them by group
        colordict = {}
        col_idx = -1
        _cat_index = 0
        for cat in catlist:
            if _cat_index % same_col_n_groups == 0:
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


def make_folium_html_plot(
    gdf,
    variable_column,
    var_display_name,
    var_unit,
    label_column,
    label_col_map,
    vmin=None,
    vmax=None,
    radius=13,
    fill_alpha=0.6,
    mpl_cmap_name="viridis",
    max_fps=4,
    dt_disp_fmt="%Y-%m-%d %H:%M",
):

    # create a map
    m = folium.Map(
        location=_get_init_mapcenter(gdf),
        tiles="cartodbpositron",
        zoom_start=10,
        attr="<a href=https://github.com/vergauwenthomas/MetObs_toolkit </a>",
    )

    # add extra tiles
    folium.TileLayer("OpenStreetMap", overlay=False, name="OSM").add_to(m)
    # RIP free Stamen tiles
    # folium.TileLayer("Stamen Terrain", overlay=False, name='Terrain', show=False).add_to(m)
    # folium.TileLayer("stamentoner", overlay=False, name='Toner', show=False).add_to(m)

    # Coloring
    if vmin is None:
        vmin = gdf[variable_column].min()
    if vmax is None:
        vmax = gdf[variable_column].max()

    # Create colormap to display on the map
    norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
    mapper = matplotlib.cm.ScalarMappable(
        norm=norm, cmap=matplotlib.colormaps[mpl_cmap_name]
    )
    colormap = brcm.LinearColormap(
        colors=mapper.cmap.colors,
        index=None,
        vmin=vmin,
        vmax=vmax,
        caption=f"{var_display_name} ({var_unit}) colorbar",
    )

    # linear colorscale for values
    def map_value_to_hex(series, vmin, vmax, cmapname="viridis"):
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax, clip=True)
        mapper = matplotlib.cm.ScalarMappable(
            norm=norm, cmap=matplotlib.colormaps[cmapname]
        )

        return series.apply(lambda x: str(matplotlib.colors.to_hex(mapper.to_rgba(x))))

    gdf["value_color"] = map_value_to_hex(
        gdf[variable_column], vmin, vmax, cmapname=mpl_cmap_name
    )

    # check if all labels are defined
    if (
        len(
            [
                lab
                for lab in gdf[label_column].unique()
                if lab not in label_col_map.keys()
            ]
        )
        > 0
    ):
        sys.exit(
            f'Unmapped labels found: {[lab for lab in gdf["label"].unique() if lab not in label_col_map.keys()]}'
        )

    gdf["label_color"] = gdf[label_column].map(label_col_map)

    # Serialize Data to Features
    def make_scater_feature(row):
        dtstring = pd.to_datetime([row["datetime"]]).strftime(dt_disp_fmt)[0]
        coords = [[row["geometry"].x, row["geometry"].y]]
        popup_str = f" <b>{row['name']}</b>  <br> {'{:.1f}'.format(row[variable_column])} {var_unit} <br> {row[label_column]}"

        features_instance = {
            "type": "Feature",
            "geometry": {
                "type": "MultiPoint",
                "coordinates": coords,
            },
            "properties": {
                "times": [dtstring],
                "popup": popup_str,
                "tooltip": f'{row["name"]}',
                "id": "geenidee",
                "icon": "circle",
                "iconstyle": {
                    "fillColor": row["value_color"],
                    "fillOpacity": fill_alpha,
                    "stroke": "false",
                    "radius": radius,
                    "color": row["label_color"],
                },
            },
        }
        return features_instance

    features = gdf.apply(make_scater_feature, axis=1).to_list()

    # Add data to the map
    folium_plugins.TimestampedGeoJson(
        {
            "type": "FeatureCollection",
            "features": features,
        },
        period="PT1H",
        duration="PT1H",
        add_last_point=False,
        auto_play=False,
        loop=False,
        max_speed=max_fps,  # fps
        loop_button=True,
        date_options="YYYY/MM/DD HH:mm:ss",
        time_slider_drag_update=True,
    ).add_to(m)

    m.add_child(colormap)
    # add control
    folium.LayerControl().add_to(m)

    return m


def geospatial_plot(
    plotdf,
    variable,
    timeinstance,
    title,
    legend,
    legend_title,
    vmin,
    vmax,
    plotsettings,
    categorical_fields,
    static_fields,
    display_name_mapper,
    data_template,
    boundbox,
):
    """Make geospatial plot of a variable (matplotlib).

    Parameters
    ----------
    plotdf : geopandas.GeoDataFrame
        A geodataframe containing a geometry column and the column representing
        the variable to plot.
    variable : str
        Name of the variable to plot.
    timeinstance : datetime.datetime
        The timeinstance to plot the variable for, if the variable is
        timedependant.
    title : str
        Title of the figure.
    legend : bool
        If True the legend will be added to the figure.
    vmin : numeric
        The variable value to use the minimum-color for..
    vmax : numeric
        The variable value to use the maximum-color for.
    plotsettings : dict
        The default plotting settings.
    categorical_fields : list
        A list of variables that are interpreted to be categorical, so to use
        a categorical coloring scheme.
    static_fields : bool
        If True the variable is assumed to be time independant.
    display_name_mapper : dict
        Must contain at least {varname: varname_str_rep}, where the
        varname_str_rep is the string representation of the variable to plot.
    data_template : dict
        The dataset template for string representations.
    boundbox : shapely.box
        The boundbox to represent the spatial extend of the plot.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The plotted axes.

    """
    # Load default plot settings
    default_settings = plotsettings["spatial_geo"]

    # subset to obstype
    plotdf = plotdf[[variable, "geometry"]]

    # Subset to the stations that have coordinates
    ignored_stations = plotdf[plotdf["geometry"].isnull()]
    plotdf = plotdf[~plotdf["geometry"].isnull()]
    if plotdf.empty:
        logger.warning(
            f"No coordinate data found, geoplot can not be made. Plotdf: {plotdf}"
        )
        return

    if not ignored_stations.empty:
        # logger.error(f'No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!')
        logger.warning(
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
    use_extent = find_plot_extent(
        geodf=gpd.GeoDataFrame(plotdf),
        user_bounds=boundbox,
        default_extentlist=default_settings["extent"],
    )

    ax = _spatial_plot(
        gdf=plotdf,
        variable=variable,
        legend=legend,
        use_quantiles=use_quantiles,
        is_categorical=is_categorical,
        k_quantiles=default_settings["n_for_categorical"],
        cmap=default_settings["cmap"],
        figsize=default_settings["figsize"],
        extent=use_extent,
        title=title,
        legend_title=legend_title,
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
    figsize,
    extent,
    title,
    legend_title,
    vmin,
    vmax,
):
    # TODO: docstring + beter positionion of the lengends
    gdf = gpd.GeoDataFrame(gdf)
    gdf = gdf.to_crs("epsg:4326")

    fig, ax = plt.subplots(
        1, 1, figsize=figsize, subplot_kw={"projection": ccrs.PlateCarree()}
    )

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
        # categorical legend
        legend_kwds = {"loc": "best", "title": legend_title}
        vmin = None
        vmax = None
        cax = None
    else:
        # colorbar
        legend_kwds = {"label": legend_title}
        divider = make_axes_locatable(ax)

        cax = divider.append_axes(
            "right", size="5%", pad=0.1, axes_class=matplotlib.axes._axes.Axes
        )

    # add observations as scatters
    gdf.plot(
        column=variable,
        scheme=scheme,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        # color='black',
        edgecolor="black",
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

    ax.add_feature(cfeature.LAND)
    ax.add_feature(cfeature.BORDERS)
    ax.add_feature(cfeature.COASTLINE)

    ax.set_title(title)

    return ax


def _sorting_function(label_vec, custom_handles, number_of_labels_types=4):
    """Sort the order of legend items."""
    # TODO: clean this up? rewrite to better code?
    sorted_vec = []
    # group 1, 2, 3
    for i in range(1, number_of_labels_types + 1):  # loop over the type of labels
        for j in range(len(label_vec)):  # loop over the length of the label_vec
            if label_vec[j] == i:
                sorted_vec.append(j)
                # makes a vector of same size as label_vec
                # but with the right order of permutations.
    sorted_handles = [custom_handles[i] for i in sorted_vec]
    # reordering the custom handles to put 1 at the front

    return sorted_handles


def _format_datetime_axis(axes):
    """Set the xaxes to autodateformat."""
    xtick_locator = mdates.AutoDateLocator()
    xtick_formatter = mdates.AutoDateFormatter(xtick_locator)

    axes.xaxis.set_major_locator(xtick_locator)
    axes.xaxis.set_major_formatter(xtick_formatter)
    return axes


def _create_linecollection(
    linedf,
    colormapper,
    linestylemapper,
    plotsettings,
    const_color=None,
    value_col_name="value",
    label_col_name="label",
):

    # 1. convert datetime to numerics values
    if linedf.index.name == "datetime":
        inxval = mdates.date2num(linedf.index.to_pydatetime())
    else:
        linedf = linedf.reset_index()
        linedf = linedf.set_index("datetime")
        inxval = mdates.date2num(linedf.index.to_pydatetime())

    # 2. convert df to segments
    points = np.array([inxval, linedf[value_col_name]]).T.reshape(-1, 1, 2)
    segments = np.concatenate([points[:-1], points[1:]], axis=1)

    # 3. get styling info
    if const_color is None:
        color = linedf[label_col_name].map(colormapper).to_list()
    else:
        color = [const_color] * linedf.shape[0]
    linewidth = [plotsettings["time_series"]["linewidth"]] * linedf.shape[0]
    zorder = plotsettings["time_series"]["linezorder"]
    linestyle = linedf[label_col_name].map(linestylemapper).fillna("-").to_list()

    # 4. Make line collection
    lc = LineCollection(
        segments=segments,
        colors=color,
        linewidths=linewidth,
        zorder=zorder,
        linestyle=linestyle,
    )
    return lc


def timeseries_plot(
    mergedf,
    title,
    ylabel,
    colorby,
    show_legend,
    show_outliers,
    show_filled,
    settings,
    _ax=None,  # needed for GUI, not recommended use
    colorby_name_colordict=None,
):  # when colorscheme will be reused
    """Make a timeseries plot.

    Parameters
    ----------
    mergedf : pandas.DataFrame
        The dataframe containing the observations as a 'value'-column and
        labels to plot.
    title : str
        Title of the figure.
    ylabel : str
        The label for the vertical axes.
    colorby : "label" or "name"
        If "label", the toolkit label is used for the colorscheme. If "name",
        the name of the station is used for the colorscheme.
    show_legend : bool
        If True, the legend will be added under the plot.
    show_filled : bool
        If True, the filled values will be plotted.
    settings : dict, optional
        The default plotting settings.
    _ax : matplotlib.pyplot.axes
        An axes to plot on. If None, a new axes will be made. The
        default is None.
    colorby_name_colorscheme : dict
        A colormapper for the station names. If None, a new colormapper will
        be created. The default is None.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The plotted axes.
    colormapper : dict
        The use colormap.

    """
    plot_settings = settings.app["plot_settings"]

    if isinstance(_ax, type(None)):
        # init figure
        fig, ax = plt.subplots(figsize=plot_settings["time_series"]["figsize"])
    else:
        ax = _ax

    # get data ready
    mergedf = mergedf[~mergedf.index.duplicated()]

    # get min max datetime to set xrange
    dt_min = mergedf.index.get_level_values("datetime").min()
    dt_max = mergedf.index.get_level_values("datetime").max()

    # define different groups (different plotting styles)
    # ok group
    ok_labels = ["ok"]

    # filled value groups
    fill_labels = [val for val in settings.gap["gaps_fill_info"]["label"].values()]
    missing_fill_labels = [
        val for val in settings.missing_obs["missing_obs_fill_info"]["label"].values()
    ]
    fill_labels.extend(missing_fill_labels)

    # qc outlier labels
    qc_labels = [
        val["outlier_flag"] for key, val in settings.qc["qc_checks_info"].items()
    ]

    # no value group
    no_vals_labels = [
        settings.gap["gaps_info"]["gap"]["outlier_flag"],
        settings.gap["gaps_info"]["missing_timestamp"]["outlier_flag"],
    ]
    # duplicated timestamp and invalid input outliers do not have a known value, so add them to this group
    no_vals_labels.append(
        settings.qc["qc_checks_info"]["duplicated_timestamp"]["outlier_flag"]
    )
    no_vals_labels.append(
        settings.qc["qc_checks_info"]["invalid_input"]["outlier_flag"]
    )

    # no_vals_df = mergedf[mergedf['label'].isin(no_vals_labels)]

    if colorby == "label":

        # aggregate groups and make styling mappers

        col_mapper = _all_possible_labels_colormapper(settings)  # get color mapper

        # linestyle mapper
        line_mapper = {
            lab: plot_settings["time_series"]["linestyle_ok"] for lab in ok_labels
        }
        line_mapper.update(
            {lab: plot_settings["time_series"]["linestyle_fill"] for lab in fill_labels}
        )

        # set hight of the vertical lines for no vals
        vlin_min = mergedf[mergedf["label"] == "ok"]["value"].min()
        vlin_max = mergedf[mergedf["label"] == "ok"]["value"].max()

        # line labels
        line_labels = ["ok"]
        line_labels.extend(fill_labels)

        # ------ missing obs ------ (vertical lines)
        missing_df = mergedf[mergedf["label"].isin(no_vals_labels)]
        missing_df = missing_df.reset_index()
        ax.vlines(
            x=missing_df["datetime"].to_numpy(),
            ymin=vlin_min,
            ymax=vlin_max,
            linestyle="--",
            color=missing_df["label"].map(col_mapper),
            zorder=plot_settings["time_series"]["dashedzorder"],
            linewidth=plot_settings["time_series"]["linewidth"],
        )

        # ------ outliers ------ (scatters)
        outlier_df = mergedf[mergedf["label"].isin(qc_labels)]
        outlier_df = outlier_df.reset_index()
        outlier_df.plot(
            kind="scatter",
            x="datetime",
            y="value",
            ax=ax,
            color=outlier_df["label"].map(col_mapper),
            legend=False,
            zorder=plot_settings["time_series"]["scatterzorder"],
            s=plot_settings["time_series"]["scattersize"],
        )

        # -------- Ok and filled observation -------- (lines)
        for sta in mergedf.index.get_level_values("name").unique():
            stadf = xs_save(mergedf, sta, "name")  # subset to one station
            linedf = stadf[
                stadf["label"].isin(line_labels)
            ]  # subset all obs that are repr by lines

            # now add the other records, and convert the value to nan to avoid
            # interpolation in the plot
            stadf.loc[~stadf.index.isin(linedf.index), "value"] = np.nan
            # (WARNING): The above line converts all values in the mergedf, to
            # Nan's if the label is not in 'line_labels' !!! Thus plot all other
            # categories in advance and the line plot at the end. The zorder,
            # takes care of what is displayed on top.

            # make line collection
            sta_line_lc = _create_linecollection(
                linedf=stadf,
                colormapper=col_mapper,
                linestylemapper=line_mapper,
                plotsettings=plot_settings,
            )
            ax.add_collection(sta_line_lc)

        # create legend
        if show_legend:

            custom_handles = []  # add legend items to it
            label_vec = []  # add type of label
            for label in mergedf["label"].unique():
                outl_color = col_mapper[label]

                if label in ok_labels:
                    custom_handles.append(
                        Line2D([0], [0], color=outl_color, label="ok", lw=4)
                    )
                    label_vec.append(1)

                elif label in fill_labels:
                    custom_handles.append(
                        Line2D(
                            [0],
                            [0],
                            color=outl_color,
                            label=f"filled value ({label})",
                            lw=1,
                            linestyle="--",
                        )
                    )
                    label_vec.append(2)

                elif label in no_vals_labels:
                    custom_handles.append(
                        Line2D(
                            [0],
                            [0],
                            color=outl_color,
                            label=f"{label}",
                            lw=1,
                            linestyle="--",
                            linewidth=2,
                        )
                    )
                    label_vec.append(3)

                else:
                    custom_handles.append(
                        Line2D(
                            [0],
                            [0],
                            marker="o",
                            color="w",
                            markerfacecolor=outl_color,
                            label=label,
                            lw=1,
                        )
                    )
                    label_vec.append(4)

            custom_handles = _sorting_function(label_vec, custom_handles)

            box = ax.get_position()
            ax.set_position(
                [box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.85]
            )
            ax.legend(
                handles=custom_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.25),
                fancybox=True,
                shadow=True,
                ncol=plot_settings["time_series"]["legend_n_columns"],
            )

    elif colorby == "name":
        # subset obs to plot
        line_labels = ["ok"]
        if show_outliers:
            line_labels.extend(qc_labels)
        if show_filled:
            line_labels.extend(fill_labels)

        # all lines are solid lines
        line_style_mapper = {lab: "-" for lab in line_labels}

        # create color mapper if none is given
        if colorby_name_colordict is None:
            col_mapper = make_cat_colormapper(
                mergedf.index.get_level_values("name").unique(),
                plot_settings["time_series"]["colormap"],
            )
        else:
            col_mapper = colorby_name_colordict

        # iterate over station and make line collection to avoid interpolation
        for sta in mergedf.index.get_level_values("name").unique():
            stadf = xs_save(mergedf, sta, "name")  # subset to one station
            linedf = stadf[
                stadf["label"].isin(line_labels)
            ]  # subset all obs that are repr by lines

            # now add the other records, and convert the value to nan to avoid
            # interpolation in the plot
            stadf.loc[~stadf.index.isin(linedf.index), "value"] = np.nan

            # make line collection
            sta_line_lc = _create_linecollection(
                linedf=stadf,
                colormapper=None,
                const_color=col_mapper[sta],
                linestylemapper=line_style_mapper,
                plotsettings=plot_settings,
            )
            ax.add_collection(sta_line_lc)

        if show_legend is True:
            # create a legend item for each station
            custom_handles = []  # add legend items to it
            names = mergedf.index.get_level_values("name").unique().to_list()
            # sort legend items alphabetically
            names.sort()
            for sta in names:
                custom_handles.append(
                    Line2D([0], [0], color=col_mapper[sta], label=sta, lw=4)
                )

            box = ax.get_position()
            ax.set_position(
                [box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.88]
            )
            primary_legend = ax.legend(
                handles=custom_handles,
                loc="upper center",
                bbox_to_anchor=(0.5, -0.2),
                fancybox=True,
                shadow=True,
                ncol=plot_settings["time_series"]["legend_n_columns"],
            )
            ax.add_artist(primary_legend)

    # Set title
    ax.set_title(title)

    # datetime formatter
    ax = _format_datetime_axis(ax)

    # Set x and y labels
    ax.set_ylabel(ylabel)

    # set x,y limits
    ax.set_xlim(mdates.date2num(dt_min), mdates.date2num(dt_max))
    ax.autoscale(axis="y")

    return ax, col_mapper


def model_timeseries_plot(
    df,
    obstype,
    title,
    ylabel,
    settings,
    show_primary_legend,
    add_second_legend=True,
    _ax=None,  # needed for GUI, not recommended use
    colorby_name_colordict=None,
):
    """Make a timeseries plot for modeldata.

    The timeseries are plotted as dashed lines.

    Parameters
    ----------
    df : pandas.DataFrame
        The dataframe containing the timeseries.
    obstype : str
        The observation type to plot. Must be a column in the df.
    title : str
        Title of the figure.
    ylabel : str
        The label for the vertical axes.
    settings : dict, optional
        The default plotting settings.
    show_primary_legend : bool
        If True, all stationnames with corresponding color are presented in a
        legend.
    add_second_legend : bool, optional
        If True, a small legend is added indicating the solid lines are
        observations and the dashed lines are modeldata. The default is True.
    _ax : matplotlib.pyplot.axes
        An axes to plot on. If None, a new axes will be made. The
        default is None.
    colorby_name_colorscheme : dict
        A colormapper for the station names. If None, a new colormapper will
        be created. The default is None.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The plotted axes.
    colormapper : dict
        The use colormap.
    """
    plot_settings = settings.app["plot_settings"]

    if isinstance(_ax, type(None)):
        # init figure
        fig, ax = plt.subplots(figsize=plot_settings["time_series"]["figsize"])
    else:
        ax = _ax

    # get data ready
    df = df[~df.index.duplicated()]

    # rename and create dummy columns so that linecollection can be used
    df = df.rename(columns={obstype: "value"})
    df["label"] = "modeldata"

    # all lines are dashed lines
    line_style_mapper = {"modeldata": "--"}

    # create color mapper if none is given
    if colorby_name_colordict is None:
        col_mapper = make_cat_colormapper(
            df.index.get_level_values("name").unique(),
            plot_settings["time_series"]["colormap"],
        )
    else:
        col_mapper = colorby_name_colordict

    # iterate over station and make line collection to avoid interpolation
    for sta in df.index.get_level_values("name").unique():
        stadf = xs_save(df, sta, "name")  # subset to one station

        # make line collection
        sta_line_lc = _create_linecollection(
            linedf=stadf,
            colormapper=None,
            const_color=col_mapper[sta],
            linestylemapper=line_style_mapper,
            plotsettings=plot_settings,
        )
        ax.add_collection(sta_line_lc)

    if show_primary_legend is True:
        # create a legend item for each station
        custom_handles = []  # add legend items to it
        names = df.index.get_level_values("name").unique().to_list()
        # sort legend items alphabetically
        names.sort()
        for sta in names:
            custom_handles.append(
                Line2D(
                    [0], [0], color=col_mapper[sta], label=f"modeldata at {sta}", lw=4
                )
            )

        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.88]
        )
        primary_legend = ax.legend(
            handles=custom_handles,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.2),
            fancybox=True,
            shadow=True,
            ncol=plot_settings["time_series"]["legend_n_columns"],
        )
        ax.add_artist(primary_legend)

    if add_second_legend:
        line_solid = Line2D(
            [], [], color="black", linestyle="--", linewidth=1.5, label=r"model"
        )
        line_dashed = Line2D(
            [], [], color="black", linestyle="-", linewidth=1.5, label=r"observations"
        )
        secondary_legend = ax.legend(handles=[line_solid, line_dashed], loc="best")
        ax.add_artist(secondary_legend)

    # Set title
    ax.set_title(title)

    # datetime formatter
    ax = _format_datetime_axis(ax)

    # Set x and y labels
    ax.set_ylabel(ylabel)

    # set x lim
    # ax.set_xlim(left=dt_min, right=dt_max)
    # ax.set_ylim(bottom=y_min, top=y_max)
    ax.autoscale()

    return ax, col_mapper


def cycle_plot(
    cycledf,
    errorbandsdf,
    title,
    plot_settings,
    aggregation,
    data_template,
    obstype,
    y_label,
    legend,
    show_zero_horizontal=False,
):
    """Plot a cycle as a lineplot.


    Parameters
    ----------
    cycledf : pandas.DataFrame
        The dataframe containing the cycle values.
    errorbandsdf : pandas.dataframe
        The dataframe containing the std values.
    title : str
        Title of the plot.
    plot_settings : dict
        The cycle-specific settings.
    aggregation : list
        A list of strings to indicate the group defenition.
    data_template : dict
        The template of the dataset.
    obstype : str
        The observation type to plot.
    y_label : str
        The label for the vertical axes.
    legend : bool
        If True, a legend is added to the figure.
    show_zero_horizontal : bool, optional
        If True, a black horizontal line at y=0 is drawn. The default is False.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The axes of the plot.

    """
    # init figure
    fig, ax = plt.subplots(figsize=plot_settings["figsize"])

    # which colormap to use:
    if cycledf.shape[1] <= plot_settings["n_cat_max"]:
        cmap = plot_settings["cmap_categorical"]
    else:
        cmap = plot_settings["cmap_continious"]

    cycledf.plot(ax=ax, title=title, legend=False, cmap=cmap)
    if legend:
        box = ax.get_position()
        ax.set_position(
            [box.x0, box.y0 + box.height * 0.2, box.width, box.height * 0.88]
        )
        ax.legend(
            cycledf.columns.values.tolist(),
            loc="upper center",
            bbox_to_anchor=(0.5, -0.2),
            fancybox=True,
            shadow=True,
            ncol=plot_settings["legend_n_columns"],
        )

    if errorbandsdf is not None:
        # Extract colorscheme from the plot
        col_sheme = {line.get_label(): line.get_color() for line in ax.get_lines()}

        for sta in errorbandsdf.columns:
            ax.fill_between(
                errorbandsdf.index,
                cycledf[sta] - errorbandsdf[sta],
                cycledf[sta] + errorbandsdf[sta],
                alpha=plot_settings["alpha_error_bands"],
                color=col_sheme[sta],
            )

    if show_zero_horizontal:
        ax.axhline(y=0.0, color="black", linestyle="--")

    return ax


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


def correlation_scatter(
    full_cor_dict, groupby_labels, obstypes, title, cor_scatter_settings
):
    """Plot the correlation variation as a scatterplot.

    The statistical significance is indicate by the scattertype.

    Parameters
    ----------
    full_cor_dict : dict
        A dictionary containing the 'cor matrix', and 'significance matrix'
        keys and corresponding matrices.
    groupby_labels : str or list
        The groupdefenition that is used for the xaxes label.
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
        scatter_groups = to_scatter.groupby("markers")
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


def _make_pie_from_freqs(
    freq_dict, colormapper, ax, plot_settings, radius, labelsize=10
):
    """Make one pie for a dict of frequencies."""
    # To dataframe
    stats = pd.Series(freq_dict, name="freq").to_frame()

    # make color mapper
    stats["color"] = stats.index.map(colormapper)

    if (stats["freq"] == 0.0).all():
        # add a 100% no occurences to it, so it can be plotted
        no_oc_df = pd.DataFrame(
            index=["No occurences"],
            data={"freq": [100.0], "color": [plot_settings["color_mapper"]["ok"]]},
        )
        stats = pd.concat([stats, no_oc_df])

    # Remove zero occurence labels (they clutter up the lables in the pies)
    stats = stats[stats["freq"] != 0]
    # Make pie
    patches, text = ax.pie(
        stats["freq"],
        colors=stats["color"],
        radius=radius,
        labels=[
            f"{j}, {s:0.1f}%"
            for j, s in zip(stats.index.to_list(), stats["freq"].to_list())
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
    missing_obs_settings = settings.missing_obs["missing_obs_fill_info"]

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
    for method, label in missing_obs_settings["label"].items():
        mapper[label] = color_defenitions[method]

    return mapper


def qc_stats_pie(
    final_stats, outlier_stats, specific_stats, plot_settings, qc_check_info, title
):
    """Make overview Pie-plots for the frequency statistics of labels.

    Parameters
    ----------
    final_stats : dict
        Dictionary containing occurence frequencies for all labels.
    outlier_stats : dict
        Dictionary with frequency statistics of outlier-labels.
    specific_stats : dict
        Dictionary containing the effectiviness of quality control checks
        individually.
    plot_settings : dict
        The specific plot settings for the pie plots.
    qc_check_info : dict
        The qc info for all checks (includes the color scheme)..
    title : str
        Title of the figure.

    Returns
    -------
    None.

    """
    # restore rcParams
    plt.rcParams = plt.rcParamsDefault

    # Specify rcParams

    # axes title
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
        colors=[spec_col_mapper[col] for col in specific_df.index],
    )

    # Specific styling setings per pie
    for ax in axlist:
        # specific style formatting
        ax.yaxis.set_visible(False)  # ignore the default pandas title

    fig.subplots_adjust(hspace=0.7)
    fig.suptitle(
        title,
        # fontsize=30,
    )
    plt.show()

    return
