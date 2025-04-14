import geopandas as gpd

import matplotlib
import matplotlib.pyplot as plt


import cartopy.crs as ccrs
import cartopy.feature as cfeature

from metobs_toolkit.geometry_functions import find_plot_extent
from mpl_toolkits.axes_grid1 import make_axes_locatable

# TODO: FIXME


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
        time-dependant.
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
        If True the variable is assumed to be time independent.
    boundbox : shapely.box
        The boundbox to represent the spatial extent of the plot.

    Returns
    -------
    ax : matplotlib.pyplot.axes
        The plotted axes.

    """
    # Load default plot settings
    default_settings = plotsettings["spatial_geo"]

    # subset to obstype
    plotdf = plotdf[["plot_value", "geometry"]]

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
            vmin = gdf["plot_value"].min()
            vmax = gdf["plot_value"].max()

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
        column="plot_value",
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
