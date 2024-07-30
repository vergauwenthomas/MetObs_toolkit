#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 15:13:17 2024

@author: thoverga
"""


import logging
import os

from metobs_toolkit.plotting_functions import (
    geospatial_plot,
    timeseries_plot,
    folium_plot,
    add_stations_to_folium_map,
    make_folium_html_plot,
)

from metobs_toolkit.landcover_functions import connect_to_gee, _validate_metadf

from metobs_toolkit.df_helpers import (
    multiindexdf_datetime_subsetting,
    # fmt_datetime_argument,
    init_multiindex,
    init_multiindexdf,
    # init_triple_multiindexdf,
    metadf_to_gdf,
    get_freqency_series,
    value_labeled_doubleidxdf_to_triple_idxdf,
    xs_save,
    concat_save,
)


logger = logging.getLogger(__name__)


class DatasetVisuals:
    """Extension on the metobs_toolkit.Dataset class with visualisation methods"""

    def make_plot(
        self,
        stationnames=None,
        obstype="temp",
        colorby="name",
        starttime=None,
        endtime=None,
        title=None,
        y_label=None,
        legend=True,
        show_outliers=True,
        show_filled=True,
        _ax=None,  # needed for GUI, not recommended use
    ):
        """
        This function creates a timeseries plot for the dataset. The variable observation type
        is plotted for all stationnames from a starttime to an endtime.

        All styling attributes are extracted from the Settings.

        Parameters
        ----------

        stationnames : list, optional
            A list with stationnames to include in the timeseries. If None is given, all the stations are used, defaults to None.
        obstype : string, optional
             Fieldname to visualise. This can be an observation or station
             attribute. The default is 'temp'.
        colorby : 'label' or 'name', optional
             Indicate how colors should be assigned to the lines. 'label' will color the lines by their quality control label. 'name' will color by each station, defaults to 'name'.
        starttime : datetime.datetime, optional
             Specifiy the start datetime for the plot. If None is given it will use the start datetime of the dataset, defaults to None.
        endtime : datetime.datetime, optional
             Specifiy the end datetime for the plot. If None is given it will use the end datetime of the dataset, defaults to None.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated. The default is None.
        legend : bool, optional
             If True, a legend is added to the plot. The default is True.
        show_outliers : bool, optional
             If true the observations labeld as outliers will be included in
             the plot. This is only true when colorby == 'name'. The default
             is True.
        show_filled : bool, optional
             If true the filled values for gaps and missing observations will
             be included in the plot. This is only true when colorby == 'name'.
             The default is True.


        Returns
        -------
        axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>>
            >>> # Make plot
            >>> dataset.make_plot(stationnames=['vlinder02', 'vlinder16'],
            ...                   obstype='temp',
            ...                   colorby='label')
            <Axes: ...

        """

        if stationnames is None:
            logger.info(f"Make {obstype}-timeseries plot for all stations")
        else:
            logger.info(f"Make {obstype}-timeseries plot for {stationnames}")

        # combine all dataframes
        mergedf = self.get_full_status_df(return_as_wide=False)

        # subset to obstype
        mergedf = xs_save(mergedf, obstype, level="obstype")

        # Subset on stationnames
        if stationnames is not None:
            mergedf = mergedf[mergedf.index.get_level_values("name").isin(stationnames)]

        # Subset on start and endtime

        starttime = self._datetime_arg_check(starttime)
        endtime = self._datetime_arg_check(endtime)

        mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

        # Get plot styling attributes
        if title is None:
            if stationnames is None:
                if self._istype == "Dataset":
                    title = (
                        self.obstypes[obstype].get_orig_name() + " for all stations. "
                    )
                elif self._istype == "Station":
                    title = self.obstypes[obstype].get_orig_name() + " of " + self.name

            else:
                title = (
                    self.obstypes[obstype].get_orig_name()
                    + " for stations: "
                    + str(stationnames)
                )
        # create y label
        if y_label is None:
            y_label = self.obstypes[obstype].get_plot_y_label()
        # Make plot
        ax, _colmap = timeseries_plot(
            mergedf=mergedf,
            title=title,
            ylabel=y_label,
            colorby=colorby,
            show_legend=legend,
            show_outliers=show_outliers,
            show_filled=show_filled,
            settings=self.settings,
            _ax=_ax,
        )

        return ax

    def make_interactive_plot(
        self,
        obstype="temp",
        save=True,
        outputfile=None,
        starttime=None,
        endtime=None,
        vmin=None,
        vmax=None,
        mpl_cmap_name="viridis",
        radius=13,
        fill_alpha=0.6,
        max_fps=4,
        outlier_col="red",
        ok_col="black",
        gap_col="orange",
        fill_col="yellow",
    ):
        """Make interactive geospatial plot with time evolution.

        This function uses the folium package to make an interactive geospatial
        plot to illustrate the time evolution.



        Parameters
        ----------
        obstype : str or metobs_toolkit.Obstype, optional
            The observation type to plot. The default is 'temp'.
        save : bool, optional
            If true, the figure will be saved as an html-file. The default is True.
        outputfile : str, optional
            The path of the output html-file. The figure will be saved here, if
            save is True. If outputfile is not given, and save is True, than
            the figure will be saved in the default outputfolder (if given).
            The default is None.
        starttime : datetime.datetime, optional
             Specifiy the start datetime for the plot. If None is given it will
             use the start datetime of the dataset, defaults to None.
        endtime : datetime.datetime, optional
             Specifiy the end datetime for the plot. If None is given it will
             use the end datetime of the dataset, defaults to None.
        vmin : numeric, optional
            The value corresponding with the minimum color. If None, the
            minimum of the presented observations is used. The default is None.
        vmax : numeric, optional
            The value corresponding with the maximum color. If None, the
            maximum of the presented observations is used. The default is None.
        mpl_cmap_name : str, optional
            The name of the matplotlib colormap to use. The default is 'viridis'.
        radius : int, optional
            The radius (in pixels) of the scatters. The default is 13.
        fill_alpha : float ([0;1]), optional
            The alpha of the fill color for the scatters. The default is 0.6.
        max_fps : int (>0), optional
            The maximum allowd frames per second for the time evolution. The
            default is 4.
        outlier_col : str, optional
            The edge color of the scatters to identify an outliers. The default is 'red'.
        ok_col : str, optional
            The edge color of the scatters to identify an ok observation. The default is 'black'.
        gap_col : str, optional
            The edge color of the scatters to identify an missing/gap
            observation. The default is 'orange'.
        fill_col : str, optional
            The edge color of the scatters to identify a fillded observation.
            The default is 'yellow'.

        Returns
        -------
        m : folium.folium.map
            The interactive folium map.

        Note
        -------
        The figure will only appear when this is runned in notebooks. If you do
        not run this in a notebook, make sure to save the html file, and open it
        with a browser.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>>
            >>> # Make default interactive geospatial plot
            >>> dataset.make_interactive_plot()

        """
        # Check if obstype is known
        if isinstance(obstype, str):
            if obstype not in self.obstypes.keys():
                logger.error(
                    f"{obstype} is not found in the knonw observation types: {list(self.obstypes.keys())}"
                )
                return None
            else:
                obstype = self.obstypes[obstype]

        if save:
            if outputfile is None:
                if self.settings.IO["output_folder"] is None:
                    logger.error(
                        "No outputfile is given, and there is no default outputfolder specified."
                    )
                    return None
                else:
                    outputfile = os.path.join(
                        self.output_folder, "interactive_figure.html"
                    )
            else:
                # Check if outputfile has .html extension
                if not outputfile.endswith(".html"):
                    outputfile = outputfile + ".html"
                    logger.warning(
                        f"The .hmtl extension is added to the outputfile: {outputfile}"
                    )

        # Check if the obstype is present in the data
        if obstype.name not in self.df.columns:
            logger.error(f"{obstype.name} is not found in your the Dataset.")
            return None

        # Check if geospatial data is available
        if self.metadf["lat"].isnull().any():
            _sta = self.metadf[self.metadf["lat"].isnull()]["lat"]
            logger.error(f"Stations without coordinates detected: {_sta}")
            return None
        if self.metadf["lon"].isnull().any():
            _sta = self.metadf[self.metadf["lon"].isnull()]["lon"]
            logger.error(f"Stations without coordinates detected: {_sta}")
            return None

        # Construct dataframe
        combdf = self.get_full_status_df()
        combdf = xs_save(combdf, obstype.name, level="obstype")
        # Merge geospatial info
        combgdf = combdf.merge(
            self.metadf, how="left", left_on="name", right_index=True
        )

        # Subset on start and endtime
        starttime = self._datetime_arg_check(starttime)
        endtime = self._datetime_arg_check(endtime)

        combgdf = multiindexdf_datetime_subsetting(combgdf, starttime, endtime)
        combgdf = combgdf.reset_index()

        # to gdf
        combgdf = metadf_to_gdf(combgdf, crs=4326)

        # Make label color mapper
        label_col_map = {}
        # Ok label
        label_col_map["ok"] = ok_col
        # outlier labels
        for val in self.settings.qc["qc_checks_info"].values():
            label_col_map[val["outlier_flag"]] = outlier_col

        # missing labels (gaps and missing values)
        for val in self.settings.gap["gaps_info"].values():
            label_col_map[val["outlier_flag"]] = gap_col

        # fill labels
        for val in self.settings.missing_obs["missing_obs_fill_info"]["label"].values():
            label_col_map[val] = fill_col
        for val in self.settings.gap["gaps_fill_info"]["label"].values():
            label_col_map[val] = fill_col

        # make time estimation
        est_seconds = combgdf.shape[0] / 2411.5  # normal laptop
        logger.info(
            f'The figure will take approximatly (laptop) {"{:.1f}".format(est_seconds)} seconds to make.'
        )

        # Making the figure
        m = make_folium_html_plot(
            gdf=combgdf,
            variable_column="value",
            var_display_name=obstype.name,
            var_unit=obstype.get_standard_unit(),
            label_column="label",
            label_col_map=label_col_map,
            vmin=vmin,
            vmax=vmax,
            radius=radius,
            fill_alpha=fill_alpha,
            mpl_cmap_name=mpl_cmap_name,
            max_fps=int(max_fps),
        )
        if save:
            logger.info(f"Saving the htlm figure at {outputfile}")
            m.save(outputfile)
        return m

    def make_geo_plot(
        self,
        variable="temp",
        title=None,
        timeinstance=None,
        legend=True,
        vmin=None,
        vmax=None,
        legend_title=None,
        boundbox=[],
    ):
        """Make geospatial plot.

        This functions creates a geospatial plot for a field
        (observations or attributes) of all stations.

        If the field is timedepending, than the timeinstance is used to plot
        the field status at that datetime.

        If the field is categorical than the leged will have categorical
        values, else a colorbar is used.

        All styling attributes are extracted from the Settings.

        Parameters
        ----------
        variable : string, optional
            Fieldname to visualise. This can be an observation type or station
            or 'lcz'. The default is 'temp'.
        title : string, optional
            Title of the figure, if None a default title is generated. The default is None.
        timeinstance : datetime.datetime, optional
            Datetime moment of the geospatial plot. If None, the first occuring
            timestamp for wich most stations have a record of, is used. The default is None.
        legend : bool, optional
            I True, a legend is added to the plot. The default is True.
        vmin : numeric, optional
            The value corresponding with the minimum color. If None, the minimum of the presented observations is used. The default is None.
        vmax : numeric, optional
            The value corresponding with the maximum color. If None, the maximum of the presented observations is used. The default is None.
        legend_title : string, optional
            Title of the legend, if None a default title is generated. The default is None.
        boundbox : [lon-west, lat-south, lon-east, lat-north], optional
            The boundbox to indicate the domain to plot. The elemenst are numeric.
            If the list is empty, a boundbox is created automatically. The default
            is [].
        Returns
        -------
        axis : matplotlib.pyplot.geoaxes
            The geoaxes of the plot is returned.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>>
            >>> # Make default geospatial plot
            >>> dataset.make_geo_plot()
            <GeoAxes:...

        """
        # Load default plot settings
        # default_settings=Settings.plot_settings['spatial_geo']

        # Get timeinstance to present data of
        timeinstance = self._datetime_arg_check(timeinstance)
        if timeinstance is None:
            # If not specified, use the earlyest timestamp, for which most stations
            # are assumed to have a record of
            timeinstance = (
                self.df.dropna()
                .index.get_level_values(level="datetime")
                .value_counts()
                .idxmax()
            )

        logger.info(f"Make {variable}-geo plot at {timeinstance}")

        # check coordinates if available
        if not _validate_metadf(self.metadf):
            raise MetobsDatasetVisualisationError(
                "The coordinates of all stations could not be found in the meta data."
            )

        if bool(boundbox):
            if len(boundbox) != 4:
                logger.warning(
                    f"The boundbox ({boundbox}) does not contain 4 elements! The default boundbox is used!"
                )
                boundbox = []
        # create a plotdf with 'names' as index, 'plot_value' and 'geometry' as columns

        # situation 1: variable is an observationtype
        if variable in self._get_present_obstypes():
            plotdf = xs_save(self.df, variable, "obstype")
            plotdf = xs_save(plotdf, timeinstance, "datetime")
            plotdf = plotdf.merge(
                self.metadf[["geometry"]], how="left", left_index=True, right_index=True
            )
            plotdf = plotdf.rename(columns={"value": "plot_value"})
            if title is None:
                title = f"{self.obstypes[variable].get_orig_name()} at {timeinstance}."
            if (legend) & (legend_title is None):
                legend_title = f"{self.obstypes[variable].get_standard_unit()}"

        # situation 2: variable is an observationtype
        elif variable in self.metadf.columns:
            if self.metadf[variable].isna().all():
                raise MetobsDatasetVisualisationError(
                    f"There is no data availabel for {variable} in the metadf: \n {self.metadf[variable]}"
                )
            plotdf = self.metadf[[variable, "geometry"]]
            plotdf = plotdf.rename(columns={variable: "plot_value"})
            if title is None:
                title = f"{variable}"
                legend_title = ""
        else:
            raise MetobsDatasetVisualisationError(
                f"{variable} is not a found in the records an not found in the meta data."
            )

        axis = geospatial_plot(
            plotdf=plotdf,
            variable=variable,
            timeinstance=timeinstance,
            title=title,
            legend=legend,
            legend_title=legend_title,
            vmin=vmin,
            vmax=vmax,
            plotsettings=self.settings.app["plot_settings"],
            categorical_fields=self.settings.app["categorical_fields"],
            static_fields=self.settings.app["static_fields"],
            display_name_mapper=self.settings.app["display_name_mapper"],
            boundbox=boundbox,
        )

        return axis

    def make_gee_plot(self, gee_map, show_stations=True, save=False, outputfile=None):
        """Make an interactive plot of a google earth dataset.

        The location of the stations can be plotted on top of it.

        Parameters
        ----------
        gee_map : str, optional
            The name of the dataset to use. This name should be present in the
            settings.gee['gee_dataset_info']. If aggregat is True, an aggregation
            scheme should included as well. The default is 'worldcover'
        show_stations : bool, optional
            If True, the stations will be plotted as markers. The default is True.
        save : bool, optional
            If True, the map will be saved as an html file in the output_folder
            as defined in the settings if the outputfile is not set. The
            default is False.
        outputfile : str, optional
            Specify the path of the html file if save is True. If None, and save
            is true, the html file will be saved in the output_folder. The
            default is None.

        Returns
        -------
        Map : geemap.foliumap.Map
            The folium Map instance.


        Warning
        ---------
        To display the interactive map a graphical backend is required, which
        is often missing on (free) cloud platforms. Therefore it is better to
        set save=True, and open the .html in your browser

        """
        # Connect to GEE
        connect_to_gee()

        # get the mapinfo
        mapinfo = self.settings.gee["gee_dataset_info"][gee_map]

        # Read in covers, numbers and labels
        covernum = list(mapinfo["colorscheme"].keys())
        colors = list(mapinfo["colorscheme"].values())
        covername = [mapinfo["categorical_mapper"][covnum] for covnum in covernum]

        # create visparams
        vis_params = {
            "min": min(covernum),
            "max": max(covernum),
            "palette": colors,  # hex colors!
        }

        if "band_of_use" in mapinfo:
            band = mapinfo["band_of_use"]
        else:
            band = None

        Map = folium_plot(
            mapinfo=mapinfo,
            band=band,
            vis_params=vis_params,
            labelnames=covername,
            layername=gee_map,
            legendname=f"{gee_map} covers",
            # showmap = show,
        )

        if show_stations:
            if not _validate_metadf(self.metadf):
                logger.warning(
                    "Not enough coordinates information is provided to plot the stations."
                )
            else:
                Map = add_stations_to_folium_map(Map=Map, metadf=self.metadf)

        # Save if needed
        if save:
            if outputfile is None:
                # Try to save in the output folder
                if self.settings.IO["output_folder"] is None:
                    logger.warning(
                        "The outputfolder is not set up, use the update_settings to specify the output_folder."
                    )

                else:
                    filename = f"gee_{gee_map}_figure.html"
                    filepath = os.path.join(self.settings.IO["output_folder"], filename)
            else:
                # outputfile is specified
                # 1. check extension
                if not outputfile.endswith(".html"):
                    outputfile = outputfile + ".html"

                filepath = outputfile

            print(f"Gee Map will be save at {filepath}")
            logger.info(f"Gee Map will be save at {filepath}")
            Map.save(filepath)

        return Map


class MetobsDatasetVisualisationError(Exception):
    """Exception raised for errors in the template."""

    pass
