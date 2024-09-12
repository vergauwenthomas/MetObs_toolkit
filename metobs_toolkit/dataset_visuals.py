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
    # folium_map,
    add_stations_to_folium_map,
    make_folium_html_plot,
)

from metobs_toolkit.settings_files.default_formats_settings import (
    gapfill_label_group,
    failed_gapfill_label_group,
    qc_label_group,
    label_def,
)

from metobs_toolkit.gee_api import _validate_metadf

# from metobs_toolkit.landcover_functions import connect_to_gee, _validate_metadf
from metobs_toolkit.modeldata import GeeStaticModelData, GeeDynamicModelData
from metobs_toolkit.df_helpers import (
    multiindexdf_datetime_subsetting,
    metadf_to_gdf,
    xs_save,
    concat_save,
)


logger = logging.getLogger(__name__)


class DatasetVisuals:
    """Extension on the metobs_toolkit.Dataset class with visualization methods."""

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
        """Make a timeseries plot.

        This function creates a timeseries plot for the dataset. The variable observation type
        is plotted for all station names from a starttime to an endtime.

        All styling attributes are extracted from the Settings.

        Parameters
        ----------

        stationnames : list, optional
            A list with station names to include in the timeseries. If None is
            given, all the stations are used. The default to None.
        obstype : string, optional
             Fieldname to visualize. This can be an observation or station
             attribute. The default is 'temp'.
        colorby : 'label' or 'name', optional
             Indicate how colors should be assigned to the lines. 'label' will
             color the lines by their quality control label. 'name' will color
             by each station. The default to 'name'.
        starttime : datetime.datetime, optional
             Specify the start datetime for the plot. If None is given it will
             use the start datetime of the dataset. The default to None.
        endtime : datetime.datetime, optional
             Specify the end datetime for the plot. If None is given it will
             use the end datetime of the dataset. The default to None.
        title : string, optional
             Title of the figure, if None a default title is generated. The
             default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated.
             The default is None.
        legend : bool, optional
             If True, a legend is added to the plot. The default is True.
        show_outliers : bool, optional
             If true the observations labeled as outliers will be included in
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

        See Also
        -----------
        Dataset.make_geo_plot: geospatial plot.
        Dataset.make_gee_plot: geospatial plot of a GEE dataset.
        Dataset.make_interactive_plot: Interactive geospatial plot.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.coarsen_time_resolution(freq='15min')
            >>> print(dataset)
            Dataset instance containing:
                 *28 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *161272 observation records (not Nan's)
                 *0 records labeled as outliers
                 *8 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:45:00+00:00 (total duration:  14 days 23:45:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            We can now make a timeseries plot of the full dataset. By specifying
            `colorby='name'`, the colors indicate the stations.

            >>> dataset.make_plot(obstype='temp', colorby='name')
            <Axes: title={'center': 'Temperatuur for all stations. '}, ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            By specifying `colorby='label'` the colors indicate the different labels.
            Labels are assigned by quality control and the status of gaps. As an
            example, we apply default quality control to create some labels.

            We then plot the timeseries of a single station.

            >>> dataset.apply_quality_control(obstype='temp')
            >>> dataset.get_station('vlinder05').make_plot(colorby='label')
            <Axes: title={'center': 'Temperatuur of vlinder05'}, xlabel='datetime', ylabel='temp (Celsius)'>


        """
        # check if there is data
        self._data_is_required_check()

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
            y_label = self.obstypes[obstype]._get_plot_y_label()
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
        filename="interactive_figure.html",
        outputfolder=None,
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
        """Make an interactive geospatial plot with time evolution.

        This function uses the folium package to make an interactive geospatial
        plot to illustrate the time evolution.



        Parameters
        ----------
        obstype : str or metobs_toolkit.Obstype, optional
            The observation type to plot. The default is 'temp'.
        save : bool, optional
            If true, the figure will be saved as an HTML file. The default is True.
        filename : str, optional
            The filename for the HTML file. This is only used when save is True. If
            a filename is given, if it does not end with ".html", the postfix
            is added. If None, the map will not be saved as an HTML file.
            The default is "interactive_figure.html".
        outputfolder : str, optional
            The path of the folder where to save the HTML (given by `filename`), if
            save is True. If outputfoler is None, and save is True, the
            default outputfolder (see `Dataset.settings.IO['output_folder']`)
            is used. The default is None.
        starttime : datetime.datetime, optional
             Specify the start datetime for the plot. If None is given it will
             use the start datetime of the dataset. The default to None.
        endtime : datetime.datetime, optional
             Specify the end datetime for the plot. If None is given it will
             use the end datetime of the dataset. The default to None.
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
            The maximum allowed frames per second for the time evolution. The
            default is 4.
        outlier_col : str, optional
            The edge color of the scatters to identify an outliers. The default is 'red'.
        ok_col : str, optional
            The edge color of the scatters to identify an ok observation. The default is 'black'.
        gap_col : str, optional
            The edge color of the scatters to identify gaps. The default is 'orange'.
        fill_col : str, optional
            The edge color of the scatters to identify a filled observation.
            The default is 'yellow'.

        Returns
        -------
        m : folium.folium.map
            The interactive folium map.

        See Also
        -----------
        Dataset.make_plot: plot timeseries.
        Dataset.make_geo_plot: geospatial plot.
        Dataset.make_gee_plot: geospatial plot of a GEE dataset.

        Note
        -------
        The figure will only appear when this is run in notebooks. If you do
        not run this in a notebook, make sure to save the HTML file, and open it
        with a browser.

        Examples
        --------
        We start by creating a Dataset, and importing data.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> print(dataset)
        Dataset instance containing:
             *28 stations
             *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
             *483828 observation records (not Nan's)
             *0 records labeled as outliers
             *8 gaps
             *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
             *time zone of the records: UTC
             *Coordinates are available for all stations.
             *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

        We apply (default) quality control.

        >>> dataset.apply_quality_control(obstype='temp')

        To create an interactive plot, we use the `Dataset.make_interactive_plot()`
        method. We must specify a target path for the HTML file, as an example
        we save it to the current working directory (`os.getcwd()`).

        >>> import os
        >>> dataset.make_interactive_plot(obstype="temp",
        ...                              save=True,
        ...                              filename='interactive_temp_plot.html',
        ...                              outputfolder=os.getcwd(),
        ...                              starttime=None,
        ...                              endtime=None)
        <folium.folium.Map object at ...

        (You can open an HTML file with a browser.)

        """
        # check if there is data
        self._data_is_required_check()

        # Check if obstype is known
        if isinstance(obstype, str):
            if obstype not in self.obstypes.keys():
                raise MetobsDatasetVisualisationError(
                    f"{obstype} is not found in the known observation types: {list(self.obstypes.keys())}"
                )
            else:
                obstype = self.obstypes[obstype]

        if save:
            # Checkout outputfolder
            if outputfolder is None:
                if self.settings.IO["output_folder"] is None:
                    # no argument, no defualt
                    raise MetobsDatasetVisualisationError(
                        "No outputfolder is specified (as argment or as default)."
                    )
                else:  # use default
                    trg_outputdir = self.settings.IO["output_folder"]
            else:  # outputfolder argumen
                if not os.path.isdir(outputfolder):
                    raise MetobsDatasetVisualisationError(
                        f"{outputfolder} is not an existing directory."
                    )
                else:
                    trg_outputdir = outputfolder

            # Check if outputfile has .html extension
            if not filename.endswith(".html"):
                filename = filename + ".html"
                logger.warning(
                    f"The .hmtl extension is added to the outputfile: {filename}"
                )
            trg_path = os.path.join(trg_outputdir, filename)

        # Check if the obstype is present in the data
        if obstype.name not in self._get_present_obstypes():
            raise MetobsDatasetVisualisationError(
                f"{obstype.name} is not found in your the Dataset."
            )

        # Check if geospatial data is available
        if self.metadf["lat"].isnull().any():
            _sta = self.metadf[self.metadf["lat"].isnull()]["lat"]
            raise MetobsDatasetVisualisationError(
                f"Stations without coordinates detected: {_sta}"
            )

        if self.metadf["lon"].isnull().any():
            _sta = self.metadf[self.metadf["lon"].isnull()]["lon"]
            raise MetobsDatasetVisualisationError(
                f"Stations without coordinates detected: {_sta}"
            )

        starttime = self._datetime_arg_check(starttime)
        endtime = self._datetime_arg_check(endtime)

        # Construct dataframe
        combdf = self.get_full_status_df()[obstype.name].reset_index()
        # Merge geospatial info
        # TODO: Are all columns of metadf usable?? if not, do not merge them
        combgdf = combdf.merge(
            self.metadf, how="left", left_on="name", right_index=True
        )

        # Subset on start and endtime
        combgdf = multiindexdf_datetime_subsetting(combgdf, starttime, endtime)
        combgdf = combgdf.reset_index()

        # to gdf
        combgdf = metadf_to_gdf(combgdf, crs=4326)

        # Make label color mapper
        label_col_map = {}
        # Ok label
        label_col_map[label_def["goodrecord"]["label"]] = ok_col

        # outlier labels
        for qclab in [label_def[check]["label"] for check in qc_label_group]:
            label_col_map[qclab] = outlier_col

        # gaps and failed gapfill
        label_col_map[label_def["regular_gap"]["label"]] = gap_col
        for gaplab in [
            label_def[check]["label"] for check in failed_gapfill_label_group
        ]:
            label_col_map[gaplab] = gap_col

        # fill labels
        for gaplab in [label_def[check]["label"] for check in gapfill_label_group]:
            label_col_map[gaplab] = fill_col

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
            logger.info(f"Saving the htlm figure at {trg_path}")
            m.save(trg_path)
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

        This function creates a geospatial plot for a field
        (observations or attributes) of all stations.

        If the field is time-depending, than the timeinstance is used to plot
        the field status at that datetime.

        If the field is categorical then the legend will have categorical
        values, or else a colorbar is used.

        All styling attributes are extracted from the Settings.

        Parameters
        ----------
        variable : string, optional
            Fieldname to visualize. This can be an observation type or station
            or 'lcz'. The default is 'temp'.
        title : string, optional
            Title of the figure, if None a default title is generated. The default is None.
        timeinstance : datetime.datetime, optional
            Datetime moment of the geospatial plot. If None, the first occurring
            timestamp for which most stations have a record, is used. The default is None.
        legend : bool, optional
            I True, a legend is added to the plot. The default is True.
        vmin : numeric, optional
            The value corresponding with the minimum color. If None, the minimum of the presented observations is used. The default is None.
        vmax : numeric, optional
            The value corresponding with the maximum color. If None, the maximum of the presented observations is used. The default is None.
        legend_title : string, optional
            Title of the legend, if None a default title is generated. The default is None.
        boundbox : [lon-west, lat-south, lon-east, lat-north], optional
            The bound box indicates the domain to plot. The elements are numeric.
            If the list is empty, a boundbox is created automatically. The default
            is [].

        Returns
        -------
        axis : matplotlib.pyplot.geoaxes
            The geoaxes of the plot is returned.

        See Also
        -----------
        Dataset.make_plot: plot timeseries.
        Dataset.make_gee_plot: geospatial plot of a GEE dataset.
        Dataset.make_interactive_plot: Interactive geospatial plot.

        Note
        --------
        If a timezone unaware datetime is given as an argument, it is interpreted
        as if it has the same timezone as the observations.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> print(dataset)
            Dataset instance containing:
                 *28 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *483828 observation records (not Nan's)
                 *0 records labeled as outliers
                 *8 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            To create a spatial plot, we use the `Dataset.make_geo_plot()`
            method.

            >>> dataset.make_geo_plot(variable="temp")
            <GeoAxes: title={'center': 'Temperatuur at 2022-09-01 00:00:00+00:00.'}>

        .. plot::
            :context: close-figs

            If you want a different extend, timeinstance of colorscale, you can do
            that like so

            >>> import datetime
            >>> dataset.make_geo_plot(variable="temp",
            ...                       title=None,
            ...                       timeinstance=datetime.datetime(2022,9,4,16),
            ...                       legend=True,
            ...                       vmin=18,
            ...                       vmax=32,
            ...                       legend_title=None,
            ...                       boundbox=[2.17, 50.67,
            ...                                5.89, 51.42])
            <GeoAxes: title={'center': 'Temperatuur at 2022-09-04 16:00:00+00:00.'}>

        .. plot::
            :context: close-figs

        """

        # check if there is data
        self._data_is_required_check()  # not strickly a data-only method

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

    def make_gee_static_spatialplot(
        self,
        Model="lcz",
        outputfolder=None,
        filename=None,
        vmin=None,
        vmax=None,
        overwrite=False,
    ):
        """Make an interactive spatial plot of the GEE dataset and the stations.

        This method will create an interactive plot of the GEE dataset. If
        metadata is present, it will be displayed as markers on the map.

        The interactive map can be saved as an HTML file, by specifying the
        target path.


        Parameters
        ----------
        Model : str or metobs_toolkit.GeeStaticModelData, optional
            The GeeStaticModelData to plot. If a string is given, it is assumed
            to be the name of a known GeeStaticModelData. The default is 'lcz'.
        outputfolder : str or None, optional
            Path to the folder to save the HTML file. If None, the map will
            not be saved as an HTML file. The default is None.
        filename : str or None, optional
            The filename for the HTML file. If a filename is given, if it does
            not end with ".html", the prefix is added. If None, the map will not
            be saved as an HTML file. The default is None.
        vmin : num or None, optional
            If the dataset is not categorical, vmin is the minimum value
            assigned to the colormap. If None, vmin is computed by extracting
            the values at the locations of the stations. If no metadata is
            available, and vmin is None then vmin is set to 0. The default is None.
        vmax : num or None, optional
            If the dataset is not categorical, vmax is the minimum value
            assigned to the colormap. If None, vmax is computed by extracting
            the values at the locations of the stations. If no metadata is
            available, and vmax is None then vmax is set to 1. The default is None.
        overwrite : bool, optional
            If True, and if the target file exists, then it will be overwritten.
            Else, an error will be raised when the target file already exists.
            The default is False.

         Returns
         -------
         MAP : geemap.foliummap.Map
             The interactive map of the GeeStaticModelData.

        See Also
        -----------
        GeeStaticModelData : Gee Modeldata dataset without time dimension.
        Dataset.make_gee_dynamic_spatialplot: Gee interactive spatial plot of GeeDynamicModelData.

        Warning
        ---------
        To display the interactive map a graphical interactive backend is
        required, which could be missing. You can recognise this when no
        map is displayed, but the Python console prints out a message similar
        to `<geemap.foliumap.Map at 0x7ff7586b8d90>`.

        In that case, you can specify a `outputfolder` and `outputfile`, save the map as a HTML file, and
        open it with a browser.

        Examples
        -----------
        As an example, we will make a plot of the LCZ map, which is a default `GeeStaticModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> lcz_model = dataset.gee_datasets['lcz']

        If you want your stations present in the map, then you must add the
        metadf to the Modeldata. This step is not required.

        >>> #we will use the demo metadata
        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)

        We will save the map as an (HTML) file. You can specify where to save it,
        for this example we will store it in the current working directory
        (`os.getcwd()`)

        >>> import os
        >>> map = dataset.make_gee_static_spatialplot(
        ...                   Model=lcz_model,
        ...                  outputfolder = os.getcwd(),
        ...                  filename = 'LCZ_map.html',
        ...                  overwrite=True)

        """

        if isinstance(Model, str):
            if Model in self.gee_datasets.keys():
                Model = self.gee_datasets[Model]
            else:
                raise MetobsDatasetVisualisationError(
                    f"{Model} is not a known Model (name), these are the knonw Gee Modeldata: {self.gee_datasets}."
                )

        if isinstance(Model, GeeStaticModelData):
            pass
        else:
            raise MetobsDatasetVisualisationError(
                f"{Model} is not a GeeStaticModelData, but a {type(Model)}."
            )

        # Set metadata
        Model.set_metadf(self.metadf)

        return Model.make_gee_plot(
            outputfolder=outputfolder,
            filename=filename,
            vmin=vmin,
            vmax=vmax,
            overwrite=overwrite,
        )

    def make_gee_dynamic_spatialplot(
        self,
        timeinstance,
        Model="ERA5-land",
        modelobstype="temp",
        outputfolder=None,
        filename=None,
        vmin=None,
        vmax=None,
        overwrite=False,
    ):
        """Make an interactive spatial plot of the GEE dataset and the stations.

        This method will create an interactive plot of the GEE dataset at an
        instance in time. If metadata is present, it will be diplayed as
        markers on the map.

        The interactive map can be saved as an HTML file, by specifying the
        target path.


        Parameters
        ----------
        timeinstance : datetime.datetime or pandas.Timestamp
            The timeinstance to plot the GEE dataset of. If a timezone naive
            timeinstance is given, it is asumed to be in UTC. This timestamp is
            rounded down with the time resolution (.time_res). The timeinstance,
            is interpreted as UTC.
        Model : str of GeeDynamicModelData, optional
            The GeeDynamicModelData to plot. If a string is given, it is assumed
            to be the name of a known GeeDynamicModelData. The default is 'ERA5-land'.
        modelobstype : str, optional
            The name of the ModelObstype to plot. The modelobstype name must be
            known (--> not the same as an Obstype!). The default is "temp".
        outputfolder : str or None, optional
            Path to the folder to save the HTML file. If None, the map will
            not be saved as an HTML file. The default is None.
        filename : str or None, optional
            The filename for the HTML file. If a filename is given, if it does
            not end with ".html", the prefix is added. If None, the map will not
            be saved as an HTML file. The default is None.
        vmin : num or None, optional
            vmin is the minimum value assigned to the colormap. If None, and
            metadata is set, vmin is computed by computing the minimum
            modelvalue in a boundbox defined by the locations of the stations.
            If no metadata is available, and vmin is None then vmin is set to
            0. The default is None.
        vmax : num or None, optional
            vmax is the minimum value assigned to the colormap. If None, and
            metadata is set, vmax is computed by computing the minimum
            modelvalue in a boundbox defined by the locations of the stations.
            If no metadata is available, and vmax is None then vmax is set to
            1. The default is None.
        overwrite : bool, optional
            If True, and if the target file exists, then it will be overwritten.
            Else, an error will be raised when the target file already exists.
            The default is False.

        Returns
        -------
        MAP : geemap.foliummap.Map
            The interactive map of the GeeStaticModelData.

        See Also
        -----------
        GeeDynamicModelData : Gee Modeldata dataset with a time dimension.
        Dataset.make_gee_static_spatialplot: Gee interactive spatial plot of GeeStaticModelData.

        Warning
        ---------
        To display the interactive map a graphical interactive backend is
        required, which could be missing. You can recognice this when no
        map is displayed, but the python console prints out a message similar
        to `<geemap.foliumap.Map at 0x7ff7586b8d90>`.

        In that case, you can specify a `outputfolder` and `outputfile`, save the map as a HTML file, and
        open in with a browser.

        Examples
        -----------
        As an example we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()` and we will plot the 2m-temperature field.

        >>> import metobs_toolkit
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        For more details (i.g. see which modelobstypes are knonw), use the
        `GeeDynamicModelData.get_info()` method.

        >>> era5.get_info()
        Empty GeeDynamicModelData instance of ERA5-land
        ------ Details ---------
        <BLANKLINE>
         * name: ERA5-land
         * location: ECMWF/ERA5_LAND/HOURLY
         * value_type: numeric
         * scale: 2500
         * is_static: False
         * is_image: False
         * is_mosaic: False
         * credentials:
         * time res: 1h
        <BLANKLINE>
         -- Known Modelobstypes --
        <BLANKLINE>
         * temp : ModelObstype instance of temp (linked to band: temperature_2m)
            (conversion: Kelvin --> Celsius)
         * pressure : ModelObstype instance of pressure (linked to band: surface_pressure)
            (conversion: pa --> pa)
         * wind : ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)
            vectorfield that will be converted to:
              * wind_speed
              * wind_direction
            (conversion: m/s --> m/s)
        <BLANKLINE>
         -- Metadata --
        <BLANKLINE>
        No metadf is set.
        <BLANKLINE>
         -- Modeldata --
        <BLANKLINE>
        No model data is set.


        We will add metadata (station locations) in the Dataset, so that the locations appear on the map. We use the demo
        dataset's metadata for this.

        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)

        We can no use the `Dataset.make_gee_dynamic_spatialplot()` method to make
        an interactive spatial plot of the era5 Modeldata.

        We will save the output as a (HTML) file and store it in the
        current working directory (`os.getcwd`) as illustration.

        We specify a timeinstance, which is rounded-down, respecting
        the time resolution of the Gee dataset.

        >>> import os
        >>> import datetime
        >>> dt = datetime.datetime(2006,11,18, 20, 15)
        >>> str(dt)
        '2006-11-18 20:15:00'

        >>> dataset.make_gee_dynamic_spatialplot(
        ...     timeinstance=dt, #will be rounded down to 18/11/2006 20:00:00
        ...     Model=era5,
        ...     modelobstype="temp", #which modelobstype to plot
        ...     outputfolder=os.getcwd(),
        ...     filename=f'era5_temp_at_{dt}.html',
        ...     overwrite=True)
        <geemap.foliumap.Map object at ...
        """
        dt = self._datetime_arg_check(timeinstance)

        if isinstance(Model, str):
            if Model in self.gee_datasets.keys():
                Model = self.gee_datasets[Model]
            else:
                raise MetobsDatasetVisualisationError(
                    f"{Model} is not a known Model (name), these are the knonw Gee Modeldata: {self.gee_datasets}."
                )

        if isinstance(Model, GeeDynamicModelData):
            pass
        else:
            raise MetobsDatasetVisualisationError(
                f"{Model} is not a GeeDynamicModelData, but a {type(Model)}."
            )

        # Set metadata
        Model.set_metadf(self.metadf)

        return Model.make_gee_plot(
            timeinstance=dt,
            modelobstype=modelobstype,
            outputfolder=outputfolder,
            filename=filename,
            vmin=vmin,
            vmax=vmax,
            overwrite=overwrite,
        )


# =============================================================================
# Errors
# =============================================================================
class MetobsDatasetVisualisationError(Exception):
    """Exception raised for errors in the template."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
