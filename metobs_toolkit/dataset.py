#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""

import os
import copy
from datetime import datetime
import pytz
import logging
import pandas as pd
import numpy as np


from metobs_toolkit.settings import Settings
from metobs_toolkit.data_import import (
    import_data_from_csv,
    import_data_from_db,
    template_to_package_space,
    import_metadata_from_csv,
)

from metobs_toolkit.printing import print_dataset_info
from metobs_toolkit.landcover_functions import (
    connect_to_gee,
    lcz_extractor,
    height_extractor,
    lc_fractions_extractor,
)

from metobs_toolkit.plotting_functions import (
    geospatial_plot,
    timeseries_plot,
    qc_stats_pie,
)

from metobs_toolkit.qc_checks import (
    gross_value_check,
    persistance_check,
    repetitions_check,
    duplicate_timestamp_check,
    step_check,
    window_variation_check,
    invalid_input_check,
)


from metobs_toolkit.qc_statistics import get_freq_statistics
from metobs_toolkit.writing_files import write_dataset_to_csv

from metobs_toolkit.missingobs import Missingob_collection

from metobs_toolkit.gap import (
    Gap,
    remove_gaps_from_obs,
    missing_timestamp_and_gap_check,
    get_gaps_indx_in_obs_space,
    get_station_gaps,
    apply_interpolate_gaps,
    make_gapfill_df,
    apply_debias_era5_gapfill,
    gaps_to_df,
)


from metobs_toolkit.df_helpers import (
    multiindexdf_datetime_subsetting,
    remove_outliers_from_obs,
    init_multiindex,
    init_multiindexdf,
    init_triple_multiindexdf,
    metadf_to_gdf,
    conv_applied_qc_to_df,
    get_freqency_series,
    value_labeled_doubleidxdf_to_triple_idxdf,
)

from metobs_toolkit.analysis import Analysis
from metobs_toolkit.modeldata import Modeldata

from metobs_toolkit import observation_types


logger = logging.getLogger(__name__)


# =============================================================================
# Dataset class
# =============================================================================


class Dataset:
    def __init__(self):
        """
        Constructs all the necessary attributes for Dataset object.

        """
        logger.info("Initialise dataset")

        # Dataset with 'good' observations
        self.df = pd.DataFrame()

        # Dataset with outlier observations
        self.outliersdf = init_triple_multiindexdf()

        self.missing_obs = None  # becomes a Missingob_collection after import
        self.gaps = None  # becomes a list of gaps

        self.gapfilldf = init_multiindexdf()
        self.missing_fill_df = init_multiindexdf()

        # Dataset with metadata (static)
        self.metadf = pd.DataFrame()
        # dataframe containing all information on the description and mapping
        self.data_template = pd.DataFrame()

        self._istype = "Dataset"
        self._freqs = pd.Series(dtype=object)

        self._applied_qc = pd.DataFrame(columns=["obstype", "checkname"])
        self._qc_checked_obstypes = []  # list with qc-checked obstypes

        self.settings = copy.deepcopy(Settings())



    def __str__(self):
        if self.df.empty:
            return f"Empty instance of a Dataset."
        add_info = ''
        n_stations = self.df.index.get_level_values('name').unique().shape[0]
        n_obs_tot = self.df.shape[0]
        n_outl = self.outliersdf.shape[0]

        if ((not self.metadf['lat'].isnull().all()) &
            (not self.metadf['lon'].isnull().all())):
            add_info += '     *Coordinates are available for all stations. \n'


        return (f"Dataset instance containing: \n \
    *{n_stations} stations \n \
    *{n_obs_tot} observation records \n \
    *{n_outl} records labeled as outliers \n \
    *{len(self.gaps)} gaps \n \
    *{self.missing_obs.series.shape[0]} missing observations \n" + add_info)

    def __repr__(self):
        return self.__str__()


    def show_settings(self):
        """
        A function that prints out all the settings, structured per thematic.

        Returns
        -------
        None.

        """

        self.settings.show()

    def get_station(self, stationname):
        """
        Extract a metobs_toolkit.Station object from the dataset by name.

        Parameters
        ----------
        stationname : string
            The name of the station.

        Returns
        -------
        metobs_toolkit.Station
            The station object.

        """
        from metobs_toolkit.station import Station

        logger.info(f"Extract {stationname} from dataset.")

        # important: make shure all station attributes are of the same time as dataset.
        # so that all methods can be inherited.

        try:
            sta_df = self.df.xs(stationname, level="name", drop_level=False)
            sta_metadf = self.metadf.loc[stationname].to_frame().transpose()
        except KeyError:
            logger.warning(f"{stationname} not found in the dataset.")
            print(f"{stationname} not found in the dataset.")
            return None

        try:
            sta_outliers = self.outliersdf.xs(
                stationname, level="name", drop_level=False
            )
        except KeyError:
            sta_outliers = init_multiindexdf()

        sta_gaps = get_station_gaps(self.gaps, stationname)
        sta_missingobs = self.missing_obs.get_station_missingobs(stationname)

        try:
            sta_gapfill = self.gapfilldf.xs(stationname, level="name", drop_level=False)
        except KeyError:
            sta_gapfill = init_multiindexdf()

        try:
            sta_missingfill = self.missing_fill_df.xs(stationname, level="name", drop_level=False)
        except KeyError:
            sta_missingfill = init_multiindexdf()

        return Station(
            name=stationname,
            df=sta_df,
            outliersdf=sta_outliers,
            gaps=sta_gaps,
            missing_obs=sta_missingobs,
            gapfilldf=sta_gapfill,
            missing_fill_df = sta_missingfill,
            metadf=sta_metadf,
            data_template=self.data_template,
            settings=self.settings,
            _qc_checked_obstypes=self._qc_checked_obstypes,
            _applied_qc=self._applied_qc,
        )

    def show(self):
        """
        A function to print out some overview information about the Dataset.

        Returns
        -------
        None.

        """

        logger.info("Show basic info of dataset.")

        try:
            gapsdf = self.gaps.to_df()
        except:
            gapsdf = init_multiindexdf()

        if self.missing_obs is None:
            missing_obs_series = pd.Series(dtype=object)
        else:
            missing_obs_series = self.missing_obs.series

        print_dataset_info(
            self.df,
            self.outliersdf,
            gapsdf,
            missing_obs_series,
            self.settings.app["print_fmt_datetime"],
        )

    def make_plot(
        self,
        stationnames=None,
        obstype="temp",
        colorby="name",
        starttime=None,
        endtime=None,
        title=None,
        legend=True,
        show_outliers=True,
        show_filled = True,
        _ax=None, #needed for GUI, not recommended use
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
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.
        show_outliers : bool, optional
             If true the observations labeld as outliers will be included in the plot, defaults to True
        show_filled : bool, optional
             If true the filled values for gaps and missing observations will
             be included in the plot, defaults to True


        Returns
        -------
        axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

        """

        logger.info(f"Make {obstype}-timeseries plot for {stationnames}")

        # combine all dataframes
        mergedf = self.combine_all_to_obsspace()

        # subset to obstype
        mergedf = mergedf.xs(obstype, level='obstype')

        # Subset on stationnames
        if not stationnames is None:
            mergedf = mergedf.reset_index()
            mergedf = mergedf.loc[mergedf['name'].isin(stationnames)]
            mergedf = mergedf.set_index(['name', 'datetime'])


        # Subset on start and endtime
        mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

        # remove outliers if required
        if not show_outliers:
            outlier_labels = [var['outlier_flag'] for var in self.settings.qc['qc_checks_info'].values()]
            mergedf = mergedf[~mergedf['label'].isin(outlier_labels)]

        # remove filled values if required
        if not show_filled:
            fill_labels = ['gap fill', 'missing observation fill'] #toolkit representation labels
            mergedf = mergedf[~mergedf['toolkit_representation'].isin(fill_labels)]

        # Get plot styling attributes
        if title is None:
            if stationnames is None:
                if self._istype == "Dataset":
                    title = (
                        self.settings.app["display_name_mapper"][obstype]
                        + " for all stations. "
                    )
                elif self._istype == "Station":
                    title = (
                        self.settings.app["display_name_mapper"][obstype]
                        + " of "
                        + self.name
                    )

            else:
                title = (
                    self.settings.app["display_name_mapper"][obstype]
                    + " for stations: "
                    + str(stationnames)
                )


        # Make plot
        ax = timeseries_plot(
            mergedf=mergedf,
            title=title,
            xlabel="Timestamp",
            ylabel=self.data_template[obstype]["orig_name"],
            colorby=colorby,
            show_legend=legend,
            show_outliers=show_outliers,
            settings = self.settings,
            _ax = _ax
        )

        return ax

    def make_geo_plot(
        self,
        obstype="temp",
        title=None,
        timeinstance=None,
        legend=True,
        vmin=None,
        vmax=None,
    ):
        """
        This functions creates a geospatial plot for a field
        (observations or attributes) of all stations.

        If the field is timedepending, than the timeinstance is used to plot
        the field status at that datetime.

        If the field is categorical than the leged will have categorical
        values, else a colorbar is used.

        All styling attributes are extracted from the Settings.

        Parameters
        ----------
        obstype : string, optional
            Fieldname to visualise. This can be an observation or station
            attribute. The default is 'temp'.
        title : string, optional
            Title of the figure, if None a default title is generated. The default is None.
        timeinstance : datetime.datetime, optional
            Datetime moment of the geospatial plot. If None, the first available datetime is used. The default is None.
        legend : bool, optional
            I True, a legend is added to the plot. The default is True.
        vmin : numeric, optional
            The value corresponding with the minimum color. If None, the minimum of the presented observations is used. The default is None.
        vmax : numeric, optional
            The value corresponding with the maximum color. If None, the maximum of the presented observations is used. The default is None.

        Returns
        -------
        axis : matplotlib.pyplot.geoaxes
            The geoaxes of the plot is returned.

        """

        # Load default plot settings
        # default_settings=Settings.plot_settings['spatial_geo']

        # get first timeinstance of the dataset if not given
        if timeinstance is None:
            timeinstance = self.df.index.get_level_values("datetime").min()

        logger.info(f"Make {obstype}-geo plot at {timeinstance}")

        # subset to timeinstance
        plotdf = self.df.xs(timeinstance, level="datetime")

        # merge metadata
        plotdf = plotdf.merge(
            self.metadf, how="left", left_index=True, right_index=True
        )

        axis = geospatial_plot(
            plotdf=plotdf,
            variable=obstype,
            timeinstance=timeinstance,
            title=title,
            legend=legend,
            vmin=vmin,
            vmax=vmax,
            plotsettings=self.settings.app["plot_settings"],
            categorical_fields=self.settings.app["categorical_fields"],
            static_fields=self.settings.app["static_fields"],
            display_name_mapper=self.settings.app["display_name_mapper"],
            world_boundaries_map=self.settings.app["world_boundary_map"],
        )

        return axis

    # =============================================================================
    #   Gap Filling
    # =============================================================================
    def get_modeldata(
        self, modelname="ERA5_hourly", stations=None, startdt=None, enddt=None
    ):
        """
        Make a metobs_toolkit.Modeldata object with modeldata at the locations
        of the stations present in the dataset.

        Parameters
        ----------
        modelname : 'ERA5_hourly', optional
            Which dataset to download timeseries from. The default is 'ERA5_hourly'.
        stations : string or list of strings, optional
            Stationnames to subset the modeldata to. If None, all stations will be used. The default is None.
        startdt : datetime.datetime, optional
            Start datetime of the model timeseries. If None, the start datetime of the dataset is used. The default is None.
        enddt : datetime.datetime, optional
            End datetime of the model timeseries. If None, the last datetime of the dataset is used. The default is None.

        Returns
        -------
        Modl : metobs_toolkit.Modeldata
            The extracted modeldata for period and a set of stations.

        NOTE
        ------
        Only 2mT extraction of ERA5 is implemented at the moment.

        """

        Modl = Modeldata(modelname)

        # Filters
        if startdt is None:
            startdt = self.df.index.get_level_values("datetime").min()
        if enddt is None:
            enddt = self.df.index.get_level_values("datetime").max()
        if not stations is None:
            if isinstance(stations, str):
                metadf = self.metadf.loc[[stations]]
            if isinstance(stations, list):
                metadf = self.metadf.iloc[self.metadf.index.isin(stations)]
        else:
            metadf = self.metadf

        # Convert to UTC
        startdt_utc = startdt.astimezone(pytz.utc)
        enddt_utc = enddt.astimezone(pytz.utc)

        # fill modell with data
        if modelname == "ERA5_hourly":
            Modl.get_ERA5_data(metadf, startdt_utc, enddt_utc)

            return Modl
        else:
            print(f"{modelname} for set_modeldata is not implemented yet")
            return None

    def update_gaps_and_missing_from_outliers(self, obstype='temp', n_gapsize=None):
        """
        Interpret the outliers as missing observations. If there is a sequence
        of these outliers for a station, larger than n_gapsize than this will
        be interpreted as a gap.

        The outliers are not removed.

        Parameters
        ----------
        obstype : str, optional
            Use the outliers on this observation type to update the gaps and
            missing timestamps.The obstype should be an element of
            metobs_toolkit.observation_types. The default is 'temp'.
        n_gapsize : int, optional
            The minimum number of consecutive missing observations to define
            as a gap. If None, n_gapsize is taken from the settings defenition
            of gaps. The default is None.

        Returns
        -------
        None.

        Note
        -------
        Gaps and missing observations resulting from an outlier on a specific
        obstype, are assumed to be gaps/missing observation for all obstypes.

        Note
        ------
        Be aware that n_gapsize is used for the current resolution of the Dataset,
        this is different from the gap check applied on the inported data, if
        the dataset is coarsend.

        """
        if n_gapsize is None:
            n_gapsize = self.settings.gap['gaps_settings']['gaps_finder']['gapsize_n']
            if not self.metadf["assumed_import_frequency"].eq(self.metadf['dataset_resolution']).all():
                print(f'The defenition of the gapsize (n_gapsize = {n_gapsize}) \
                               will have another effect on the update of the gaps and missing \
                                   timestamps because coarsening is applied and the defenition \
                                   of the gapsize is not changed.')



        # combine to one dataframe
        mergedf = self.combine_all_to_obsspace()
        mergedf = mergedf.xs(obstype, level='obstype')


        # ignore labels
        possible_outlier_labels = [vals['outlier_flag'] for vals in self.settings.qc['qc_checks_info'].values()]


        # create groups when the final label changes
        persistance_filter = ((mergedf['label'].shift() != mergedf['label'])).cumsum()
        grouped = mergedf.groupby(['name', persistance_filter])

        #locate new gaps by size of consecutive the same final label per station
        group_sizes = grouped.size()
        outlier_groups = group_sizes[
            group_sizes > n_gapsize
        ]

        # find only groups with final label as an outlier
        gaps = []
        # new_gapsdf = pd.DataFrame()
        new_gaps_idx = init_multiindex()
        for group_idx in outlier_groups.index:
            groupdf = grouped.get_group(group_idx)
            group_final_label = groupdf['label'].iloc[0]
            if not group_final_label in possible_outlier_labels:
                #no gap candidates
                continue
            else:
                gap =Gap(name=groupdf.index.get_level_values('name')[0],
                         startdt=groupdf.index.get_level_values('datetime').min(),
                         enddt=groupdf.index.get_level_values('datetime').max())

                gaps.append(gap)
                # new_gapsdf = pd.concat([new_gapsdf,
                #                         pd.DataFrame(data={'start_gap': [groupdf.index.get_level_values('datetime').min()],
                #                                            'end_gap': [groupdf.index.get_level_values('datetime').max()]},
                #                                      index=[groupdf.index.get_level_values('name')[0]])])

                new_gaps_idx = new_gaps_idx.union(groupdf.index, sort=False)



        # add all the outliers, that are not in the new gaps to the new missing obs
        new_missing_obs = mergedf[mergedf['label'].isin(possible_outlier_labels)].index
        new_missing_obs = new_missing_obs.drop(new_gaps_idx.to_numpy(), errors='ignore')


        # to series
        missing_obs_series = new_missing_obs.to_frame().reset_index(drop=True).set_index('name')['datetime']
        # Create missing obs
        new_missing_collection = Missingob_collection(missing_obs_series)



        # update self
        self.gaps.extend(gaps)
        self.missing_obs = self.missing_obs + new_missing_collection



    # =============================================================================
    #   Gap Filling
    # =============================================================================

    # def fill_gaps_automatic(self, modeldata, obstype='temp',
    #                         max_interpolate_duration_str=None,
    #                         overwrite=True):
    #     """
    #     Fill gaps using an automatic decision on which method to use to fill.
    #     For small gaps (gap_duration <= max_interpolate_duration_str) interpolation
    #     is applied, for larger gaps the model debias method is used.

    #     Parameters
    #     ----------
    #     modeldata : metobs_toolkit.Modeldata
    #         The ERA5 Modeldata instance containing observations for the periods
    #         gaps and the leading/trailing periods..
    #     obstype : str, optional
    #         Observation type to fill the gaps for. The default is 'temp'.
    #     max_interpolate_duration_str : timedelta or timedeltastring, optional
    #         A time indication (like '5H') to indicate the maximum gapsize to use
    #         interpolation for. If None, the default settings will be used. The default is None.
    #     overwrite : bool, optional
    #         If True, present gapfill values will be overwritten. The default is True.

    #     Returns
    #     -------
    #     comb_df : TYPE
    #         DESCRIPTION.

    #     """

    #     # Validate input
    #     # check if modeldata is available
    #     if modeldata is None:
    #         print(
    #             "The dataset has no modeldate. Use the set_modeldata() function to add modeldata."
    #         )
    #         return None

    #     # check if obstype is present in eramodel
    #     assert (
    #         obstype in modeldata.df.columns
    #     ), f"{obstype} is not present in the modeldate: {modeldata}"

    #     # check if all station are present in eramodeldata
    #     stations = self.gaps.to_df().index.unique().to_list()
    #     assert all(
    #         [sta in modeldata.df.index.get_level_values("name") for sta in stations]
    #     ), f"Not all stations with gaps are in the modeldata!"

    #     if not self.gapfilldf.empty:
    #         if overwrite:
    #             print("Gapfilldf will be overwritten!")
    #             self.gapfilldf = init_multiindexdf()
    #         else:
    #             print('Gapfilldf is not empty, set "overwrite=True" to overwrite it!')
    #             print("CANCEL gap fill with ERA5")
    #             return

    #     if max_interpolate_duration_str is None:
    #         max_interpolate_duration_str = self.settings.gap["gaps_fill_settings"]["automatic"]["max_interpolation_duration_str"]

    #     fill_info = self.settings.gap["gaps_fill_info"]

    #     # select the method to apply gapfill per gap
    #     interpolate_gaps = []
    #     debias_gaps = []

    #     for gap in self.gaps.list:
    #         if gap.duration <= pd.to_timedelta(max_interpolate_duration_str):
    #             interpolate_gaps.append(gap)
    #         else:
    #             debias_gaps.append(gap)


    #     # convert to Gap_collection
    #     interpolate_gap_collection = _gap_collection_from_list_of_gaps(interpolate_gaps)
    #     debias_gap_collection =_gap_collection_from_list_of_gaps(debias_gaps)


    #     #1  Fill by interpolation

    #     filldf_interp = init_multiindexdf()

    #     fill_settings_interp = self.settings.gap["gaps_fill_settings"]["linear"]


    #     filldf_interp[obstype] = interpolate_gap_collection.apply_interpolate_gaps(
    #                         obsdf=self.df,
    #                         outliersdf=self.outliersdf,
    #                         dataset_res=self.metadf["dataset_resolution"],
    #                         obstype=obstype,
    #                         method=fill_settings_interp["method"],
    #                         max_consec_fill=fill_settings_interp["max_consec_fill"],
    #                         )

    #     # add label column
    #     filldf_interp[obstype + "_" + fill_info["label_columnname"]] = fill_info["label"]["linear"]


    #     #2 Fill by debias
    #     filldf_debias = init_multiindexdf()

    #     fill_settings_debias = self.settings.gap["gaps_fill_settings"]["model_debias"]


    #     apply_debias_era5_gapfill(gapslist = self.gaps
    #                               dataset=self,
    #                               eraModelData=modeldata,
    #                               obstype=obstype,
    #                               debias_settings=fill_settings_debias)

    #     # add label column
    #     filldf_debias[obstype + "_" + fill_info["label_columnname"]] = fill_info["label"]["model_debias"]


    #     # combine both fill df's
    #     comb_df = pd.concat([filldf_interp, filldf_debias])
    #     if overwrite:
    #         self.gapfilldf = comb_df
    #     return comb_df


    def fill_gaps_automatic(self, modeldata, obstype='temp',
                            max_interpolate_duration_str=None,
                            overwrite=True):

        # Validate input
        # check if modeldata is available
        if modeldata is None:
            print(
                "The dataset has no modeldate. Use the set_modeldata() function to add modeldata."
            )
            return None

        # check if obstype is present in eramodel
        assert (
            obstype in modeldata.df.columns
        ), f"{obstype} is not present in the modeldate: {modeldata}"

        # check if all station are present in eramodeldata
        # stations = self.gaps.to_df().index.unique().to_list()
        stations = list(set([gap.name for gap in self.gaps]))
        assert all(
            [sta in modeldata.df.index.get_level_values("name") for sta in stations]
        ), f"Not all stations with gaps are in the modeldata!"

        if not self.gapfilldf.empty:
            if overwrite:
                print("Gapfilldf will be overwritten!")
                self.gapfilldf = init_multiindexdf()
            else:
                print('Gapfilldf is not empty, set "overwrite=True" to overwrite it!')
                print("CANCEL gap fill with ERA5")
                return

        if max_interpolate_duration_str is None:
            max_interpolate_duration_str = self.settings.gap["gaps_fill_settings"]["automatic"]["max_interpolation_duration_str"]

        fill_info = self.settings.gap["gaps_fill_info"]

        # select the method to apply gapfill per gap
        interpolate_gaps = []
        debias_gaps = []

        for gap in self.gaps:
            if gap.duration <= pd.to_timedelta(max_interpolate_duration_str):
                interpolate_gaps.append(gap)
            else:
                debias_gaps.append(gap)

        #1  Fill by interpolation

        fill_settings_interp = self.settings.gap["gaps_fill_settings"]["linear"]


        apply_interpolate_gaps(
                            gapslist=interpolate_gaps,
                            obsdf=self.df,
                            outliersdf=self.outliersdf,
                            dataset_res=self.metadf["dataset_resolution"],
                            gapfill_settings=self.settings.gap['gaps_fill_info'],
                            obstype=obstype,
                            method=fill_settings_interp["method"],
                            max_consec_fill=fill_settings_interp["max_consec_fill"],
                            )

        filldf_interp = make_gapfill_df(interpolate_gaps)

        #2 Fill by debias

        fill_settings_debias = self.settings.gap["gaps_fill_settings"]["model_debias"]


        apply_debias_era5_gapfill(gapslist=debias_gaps,
                                        dataset=self,
                                        eraModelData=modeldata,
                                        obstype=obstype,
                                        debias_settings=fill_settings_debias)

        # add label column
        filldf_debias = make_gapfill_df(debias_gaps)


        # combine both fill df's
        comb_df = pd.concat([filldf_interp, filldf_debias])
        if overwrite:
            self.gapfilldf = comb_df
        return comb_df


    def fill_gaps_linear(self, obstype="temp", overwrite=True):
        """
        Fill the gaps using linear interpolation.

        Parameters
        ----------
        obstype : string, optional
            Fieldname to visualise. This can be an observation or station
            attribute. The default is 'temp'.

        Returns
        -------
        None.

        """
        # TODO logging
        fill_settings = self.settings.gap["gaps_fill_settings"]["linear"]
        fill_info = self.settings.gap["gaps_fill_info"]

        # fill gaps
        apply_interpolate_gaps(
            gapslist=self.gaps,
            obsdf=self.df,
            outliersdf=self.outliersdf,
            dataset_res=self.metadf["dataset_resolution"],
            gapfill_settings=self.settings.gap['gaps_fill_info'],
            obstype=obstype,
            method=fill_settings["method"],
            max_consec_fill=fill_settings["max_consec_fill"],
        )

        # get gapfilldf
        gapfilldf = make_gapfill_df(self.gaps)

        if overwrite:
            self.gapfilldf = gapfilldf

        return gapfilldf


    def fill_missing_obs_linear(self, obstype='temp'):
        # TODO logging
        fill_settings = self.settings.missing_obs['missing_obs_fill_settings']['linear']
        fill_info = self.settings.missing_obs['missing_obs_fill_info']



        # fill missing obs
        self.missing_obs.interpolate_missing(
                                            obsdf=self.df,
                                            resolutionseries=self.metadf["dataset_resolution"],
                                            obstype=obstype,
                                            method=fill_settings["method"],
        )
        missing_fill_df = self.missing_obs.fill_df
        missing_fill_df[obstype+'_' + fill_info["label_columnname"]] = fill_info["label"]["linear"]

        # Update attribute

        self.missing_fill_df = missing_fill_df

    def get_gaps_df(self):
        """
        List all gaps into an overview dataframe.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with stationnames as index, and the start, end and duretion
            of the gaps as columns.

        """
        return gaps_to_df(self.gaps)



    def get_analysis(self):
        """
        Create a MetObs_toolkit.Analysis instance from the Dataframe

        Returns
        -------
        metobs_toolkit.Analysis
            The Analysis instance of the Dataset.

        """

        return Analysis(obsdf = self.df,
                        metadf = self.metadf,
                        settings = self.settings,
                        data_template=self.data_template)


    def fill_gaps_era5(
        self, modeldata, method="debias", obstype="temp", overwrite=True
    ):
        """
        Fill the gaps using a metobs_toolkit.Modeldata object.


        Parameters
        ----------
        modeldata : metobs_toolkit.Modeldata
            The modeldata to use for the gapfill. This model data should the required
            timeseries to fill all gaps present in the dataset.
        method : 'debias', optional
            Specify which method to use. The default is 'debias'.
        obstype : TYPE, optional
            Fieldname to visualise. This can be an observation or station
            attribute. The default is 'temp'.
        overwrite : bool, optional
            If True, the Dataset.Gapfilldf will be overwritten. The default is True.

        Returns
        -------
        Gapfilldf : pandas.DataFrame
            A dataframe containing all gap filled values and the use method.

        """

        fill_info = self.settings.gap["gaps_fill_info"]

        # check if modeldata is available
        if modeldata is None:
            print(
                "The dataset has no modeldate. Use the set_modeldata() function to add modeldata."
            )
            return None
        # check if obstype is present in eramodel
        assert (
            obstype in modeldata.df.columns
        ), f"{obstype} is not present in the modeldate: {modeldata}"
        # check if all station are present in eramodeldata
        # stations = self.gaps.to_df().index.unique().to_list()
        stations = list(set([gap.name for gap in self.gaps]))
        assert all(
            [sta in modeldata.df.index.get_level_values("name") for sta in stations]
        ), f"Not all stations with gaps are in the modeldata!"


        if method == "debias":

            fill_settings_debias = self.settings.gap["gaps_fill_settings"]["model_debias"]


            apply_debias_era5_gapfill(gapslist=self.gaps,
                                            dataset=self,
                                            eraModelData=modeldata,
                                            obstype=obstype,
                                            debias_settings=fill_settings_debias)

            # get fill df
            filldf = make_gapfill_df(self.gaps)
        else:
            print("not implemented yet")

        if overwrite:
            self.gapfilldf = filldf
        return filldf

    def write_to_csv(
        self,
        obstype=None,
        filename=None,
        include_outliers=True,
        include_fill_values=True,
        add_final_labels=True,
        use_tlk_obsnames=True,
        overwrite_outliers_by_gaps_and_missing=True,
        seperate_metadata_file = True
    ):
        """
        Write the dataset to a file where the observations, metadata and
        (if available) the quality labels per observation type are merged
        together.

        A final qualty control label for each
        quality-controlled-observation type can be added in the outputfile.

        The file will be writen to the outputfolder specified in the settings.

        Parameters
        ----------
        obstype : string, optional
            Specify an observation type to subset all observations to. If None,
            all available observation types are writen to file. The default is
            None.
        filename : string, optional
            The name of the output csv file. If none, a standard-filename
            is generated based on the period of data. The default is None.
        include_outliers : bool, optional
            If True, the outliers will be present in the csv file. The default is True.
        include_fill_values : bool, optional
            If True, the filled gap and missing observation values will be
            present in the csv file. The default is True.
        add_final_labels : bool, optional
            If True, a column is added containing the final label of an observation. The default is True.
        use_tlk_obsnames : bool, optional
            If True, the standard naming of the metobs_toolkit is used, else
            the original names for obstypes is used. The default is True.
        overwrite_outliers_by_gaps_and_missing : bool, optional
            If the gaps and missing observations are updated using outliers,
            interpret these records as gaps/missing outliers if True. Else these
            will be interpreted as outliers. The default is True.
        seperate_metadata_file : bool, optional
            If true, the metadat is writen to a seperate file, else the metadata
            is merged to the observation in one file. The default is True.
        Returns
        -------
        None.

        """

        logger.info("Writing the dataset to a csv file")

        assert (
            not self.settings.IO["output_folder"] is None
        ), "Specify Settings.output_folder in order to export a csv."

        assert os.path.isdir(
            self.settings.IO["output_folder"]
        ), f'The outputfolder: \
            {self.settings.IO["output_folder"]} is not found. '

        # combine all dataframes
        mergedf = self.combine_all_to_obsspace(
            overwrite_outliers_by_gaps_and_missing=overwrite_outliers_by_gaps_and_missing)  # with outliers
        # Unstack mergedf
        # remove duplicates
        mergedf = mergedf[~mergedf.index.duplicated(keep='first')]

        # drop outliers if required
        if not include_outliers:
            outlier_labels = [var['outlier_flag'] for var in self.settings.qc['qc_checks_info']]
            mergedf = mergedf[~mergedf['label'].isin(outlier_labels)]


        # drop fill values if required
        if not include_fill_values:
            fill_labels = ['gap fill', 'missing observation fill'] #toolkit representation labels
            mergedf = mergedf[~mergedf['toolkit_representation'].isin(fill_labels)]

        if not obstype is None:
            mergedf = mergedf.xs(obstype, level='obstype', drop_level=False)



        # Map obstypes columns
        if not use_tlk_obsnames:
            mapper = self.data_template.transpose()['orig_name'].to_dict()
            mergedf = mergedf.reset_index()
            mergedf['new_names'] = mergedf['obstype'].map(mapper)
            mergedf = mergedf.drop(columns=['obstype'])
            mergedf = mergedf.rename(columns={'new_names': 'obstype'})
            mergedf = mergedf.set_index(['name', 'datetime', 'obstype'])


        mergedf = mergedf.unstack('obstype')

        # to one level for the columns
        mergedf.columns = [' : '.join(col).strip() for col in mergedf.columns.values]



        # columns to write
        write_dataset_to_csv(
            df=mergedf,
            metadf=self.metadf,
            filename=filename,
            outputfolder=self.settings.IO["output_folder"],
            location_info=self.settings.app["location_info"],
            seperate_metadata_file=seperate_metadata_file,
        )


    # =============================================================================
    #     Quality control
    # =============================================================================

    def apply_quality_control(
        self,
        obstype="temp",
        gross_value=True,
        persistance=True,
        repetitions=True,
        step=True,
        window_variation=True,
        # internal_consistency=True,
    ):
        """
         Apply quality control methods to the dataset.

         The default settings are used, and can be changed in the
         settings_files/qc_settings.py

         The checks are performed in a sequence: gross_vallue -->
         persistance --> ..., Outliers by a previous check are ignored in the
         following checks!

         The dataset is updated inline.

         Parameters
         ----------
         obstype : String, optional
             Name of the observationtype you want to apply the checks on. The
             default is 'temp'.
         gross_value : Bool, optional
             If True the gross_value check is applied if False not. The default
             is True.
         persistance : Bool, optional
            If True the persistance check is applied if False not. The default
            is True.. The default is True.
         step : Bool, optional
            If True the step check is applied if False not. The default is True.
         internal_consistency : Bool, optional
            If True the internal consistency check is applied if False not. The
            default is True.
         qc_info: Bool, optional
            If True info about the quality control is printed if False not. The
            default is True.
         ignore_val : numeric, optional
             Values to ignore in the quality checks. The default is np.nan.

         Returns
         ---------

         None.

        """

        if repetitions:
            print("Applying the repetitions-check on all stations.")
            logger.info("Applying repetitions check on the full dataset")

            obsdf, outl_df = repetitions_check(
                obsdf=self.df,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["qc_check_settings"],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.outliersdf = pd.concat([self.outliersdf, outl_df])

            # add this check to the applied checks
            self._applied_qc = pd.concat(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames="repetitions"
                    ),
                ],
                ignore_index=True,
            )

        if gross_value:
            print("Applying the gross-value-check on all stations.")
            logger.info("Applying gross value check on the full dataset")

            obsdf, outl_df = gross_value_check(
                obsdf=self.df,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["qc_check_settings"],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.outliersdf = pd.concat([self.outliersdf, outl_df])

            # add this check to the applied checks
            self._applied_qc = pd.concat(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames="gross_value"
                    ),
                ],
                ignore_index=True,
            )

        if persistance:
            print("Applying the persistance-check on all stations.")
            logger.info("Applying persistance check on the full dataset")

            obsdf, outl_df = persistance_check(
                station_frequencies=self.metadf["dataset_resolution"],
                obsdf=self.df,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["qc_check_settings"],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.outliersdf = pd.concat([self.outliersdf, outl_df])

            # add this check to the applied checks
            self._applied_qc = pd.concat(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames="persistance"
                    ),
                ],
                ignore_index=True,
            )

        if step:
            print("Applying the step-check on all stations.")
            logger.info("Applying step-check on the full dataset")

            obsdf, outl_df = step_check(
                obsdf=self.df,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["qc_check_settings"],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.outliersdf = pd.concat([self.outliersdf, outl_df])

            # add this check to the applied checks
            self._applied_qc = pd.concat(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(obstypes=obstype, ordered_checknames="step"),
                ],
                ignore_index=True,
            )

        if window_variation:
            print("Applying the window variation-check on all stations.")
            logger.info("Applying window variation-check on the full dataset")

            obsdf, outl_df = window_variation_check(
                station_frequencies=self.metadf["dataset_resolution"],
                obsdf=self.df,
                obstype=obstype,
                checks_info=self.settings.qc["qc_checks_info"],
                checks_settings=self.settings.qc["qc_check_settings"],
            )

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.outliersdf = pd.concat([self.outliersdf, outl_df])

            # add this check to the applied checks
            self._applied_qc = pd.concat(
                [
                    self._applied_qc,
                    conv_applied_qc_to_df(
                        obstypes=obstype, ordered_checknames="window_variation"
                    ),
                ],
                ignore_index=True,
            )

        self._qc_checked_obstypes.append(obstype)
        self._qc_checked_obstypes = list(set(self._qc_checked_obstypes))
        self.outliersdf = self.outliersdf.sort_index()



    def combine_all_to_obsspace(self, repr_outl_as_nan=False,
                                overwrite_outliers_by_gaps_and_missing=True):
        
        """
         Combine all observations, outliers, missing observations and gaps into
         one Dataframe. All observation types are combined an a label is added
         in a serperate column. 

         When gaps and missing records are updated from outliers one has to choice
         to represent these records as outliers or gaps. There can not be duplicates
         in the return dataframe. 

         By default the observation values of the outliers are saved, one can 
         choice to use these values or NaN's.
         following checks!

        

         Parameters
         ----------
         repr_outl_as_nan : bool, optional
             If True, Nan's are use for the values of the outliers. The
             default is False.
         overwrite_outliers_by_gaps_and_missing : Bool, optional
             If True, records that are labeld as gap/missing and outlier are
             labeled as gaps/missing. This has only effect when the gaps/missing
             observations are updated from the outliers. The default
             is True.

         Returns
         ---------
         combdf : pandas.DataFrame()
            A dataframe containing a continious time resolution of records, where each
            record is labeld. 

        """




        # TODO: label values from settings not hardcoding

        # =============================================================================
        # Stack outliers
        # =============================================================================

        outliersdf = self.outliersdf
        outliersdf['toolkit_representation'] = 'outlier'
        # TODO: use the repr_outl_as_nan argumenten here
        # =============================================================================
        # Stack observations
        # =============================================================================
        df = self.df
        # better save than sorry
        present_obstypes = [col for col in df if col in observation_types]
        df = df[present_obstypes]

        # to tripple index
        df = df.stack(dropna=False).reset_index().rename(columns={'level_2': 'obstype', 0: 'value'}).set_index(['name', 'datetime', 'obstype'])

        df['label'] = 'ok'
        df['toolkit_representation'] = 'observation'

        # remove outliers from the observations
        df = df[~ df.index.isin(outliersdf.index)]


        # =============================================================================
        # Stack gaps
        # =============================================================================


        # add gapfill and remove the filled records from gaps
        gapsfilldf = self.gapfilldf.copy()

        # to triple index
        gapsfilldf = value_labeled_doubleidxdf_to_triple_idxdf(gapsfilldf)
        gapsfilldf['toolkit_representation'] = 'gap fill'

        gapsidx = get_gaps_indx_in_obs_space(gapslist=self.gaps,
                                             obsdf = self.df,
                                             outliersdf = self.outliersdf,
                                             resolutionseries=self.metadf["dataset_resolution"])

        gapsdf = pd.DataFrame(index=gapsidx, columns=present_obstypes)
        gapsdf = gapsdf.stack(dropna=False).reset_index().rename(columns={'level_2': 'obstype', 0: 'value'}).set_index(['name', 'datetime', 'obstype'])

        gapsdf['label'] = self.settings.gap['gaps_info']['gap']['outlier_flag']
        gapsdf['toolkit_representation'] = 'gap'


        # Remove gaps from df
        df = df[~ df.index.isin(gapsdf.index)]

        if overwrite_outliers_by_gaps_and_missing:
            outliersdf = outliersdf.drop(index=gapsdf.index, errors='ignore')


        # Remove gapfill values records from the gaps
        gapsdf = gapsdf.drop(index=gapsfilldf.index)



        # =============================================================================
        # Stack missing
        # =============================================================================

        missingfilldf = self.missing_fill_df.copy()
        missingfilldf = value_labeled_doubleidxdf_to_triple_idxdf(missingfilldf)
        missingfilldf['toolkit_representation'] = 'missing observation fill'


        # add missing observations if they occure in observation space
        missingidx = self.missing_obs.get_missing_indx_in_obs_space(
            self.df, self.metadf["dataset_resolution"]
        )

        missingdf = pd.DataFrame(index=missingidx, columns=present_obstypes)


        missingdf = missingdf.stack(dropna=False).reset_index().rename(columns={'level_2': 'obstype', 0: 'value'}).set_index(['name', 'datetime', 'obstype'])

        missingdf['label'] = self.settings.gap['gaps_info']['missing_timestamp']['outlier_flag']
        missingdf['toolkit_representation'] = 'missing observation'

        # Remove missing from df
        df = df[~ df.index.isin(missingdf.index)]

        if overwrite_outliers_by_gaps_and_missing:
            outliersdf = outliersdf.drop(index=missingdf.index, errors='ignore')

        # Remove missingfill values records from the missing
        missingdf = missingdf.drop(index=missingfilldf.index)


        # =============================================================================
        # combine all
        # =============================================================================

        combdf = pd.concat([df, outliersdf, gapsdf, gapsfilldf, missingdf, missingfilldf]).sort_index()

        # To be shure?
        combdf = combdf[~combdf.index.duplicated(keep='first')]
        return combdf




    def get_qc_stats(self, obstype="temp", stationnames=None, make_plot=True):
        """
        Compute frequency statistics on the qc labels for an observationtype.
        The output is a dataframe containing the frequency statistics presented
        as percentages.

        These frequencies can also be presented as a collection of piecharts
        per check.

        With stationnames you can subset the data to one ore multiple stations.

        Parameters

        obstype : Str, optional
            Observation type to analyse the QC labels on. The default is
            'temp'.
        stationnames : List, Str, optional
            Stationname(s) to subset the quality labels on. If None, all
            stations are used. The default is None.
        make_plot : Bool, optional
            If True, a plot with piecharts is generated. The default is True.

        Returns

        dataset_qc_stats : pandas.DataFrame
            A table containing the label frequencies per check presented
            as percentages0.

        """

        # cobmine all and get final label
        comb_df = self.combine_all_to_obsspace()

        # subset to relevant columnt
        comb_df = comb_df.xs(obstype, level='obstype')[['label']]


        # compute freq statistics
        final_freq, outl_freq, specific_freq = get_freq_statistics(
            comb_df=comb_df,
            obstype=obstype,
            checks_info=self.settings.qc["qc_checks_info"],
            gaps_info=self.settings.gap["gaps_info"],
            applied_qc_order=self._applied_qc,
        )

        if any([stat is None for stat in [final_freq, outl_freq, specific_freq]]):
            return None

        if make_plot:
            # make pie plots
            qc_stats_pie(
                final_stats=final_freq,
                outlier_stats=outl_freq,
                specific_stats=specific_freq,
                plot_settings=self.settings.app["plot_settings"],
                qc_check_info=self.settings.qc["qc_checks_info"],
            )

        return (final_freq, outl_freq, specific_freq)

    def update_outliersdf(self, add_to_outliersdf):
        """V5"""

        self.outliersdf = pd.concat([self.outliersdf, add_to_outliersdf])

    # =============================================================================
    #     importing data
    # =============================================================================

    def coarsen_time_resolution(
        self, origin=None, origin_tz=None, freq=None, method=None, limit=None
    ):
        """
        Resample the observations to coarser timeresolution. The assumed
        dataset resolution (stored in the metadf attribute) will be updated.

        Parameters
        ----------
        origin : datetime.datetime, optional
            Define the origin (first timestamp) for the obervations. The origin
            is timezone naive, and is assumed to have the same timezone as the
            obervations. If None, the earliest occuring timestamp is used as
            origin. The default is None.
        origin_tz : str, optional
            Timezone string of the input observations. Element of
            pytz.all_timezones. If None, the timezone from the settings is
            used. The default is None.
        freq : DateOffset, Timedelta or str, optional
            The offset string or object representing target conversion.
            Ex: '15T' is 15 minuts, '1H', is one hour. If None, the target time
            resolution of the dataset.settings is used. The default is None.
        method : 'nearest' or 'bfill', optional
            Method to apply for the resampling. If None, the resample method of
            the dataset.settings is used. The default is None.
        limit : int, optional
            Limit of how many values to fill with one original observations. If
            None, the target limit of the dataset.settings is used. The default
            is None.

        Returns
        -------
        None.

        """

        if freq is None:
            freq = self.settings.time_settings["target_time_res"]
        if method is None:
            method = self.settings.time_settings["resample_method"]
        if limit is None:
            limit = int(self.settings.time_settings["resample_limit"])
        if origin_tz is None:
            origin_tz = self.settings.time_settings["timezone"]

        logger.info(
            f"Coarsening the timeresolution to {freq} using \
                    the {method}-method (with limit={limit})."
        )
        # TODO: implement buffer method
        # TODO: implement startdt point
        df = self.df.reset_index()

        if origin is None:
            # find earlyest timestamp, if it is on the hour, use it else use the following hour
            tstart = df["datetime"].min()

            if tstart.minute != 0 or tstart.second != 0 or tstart.microsecond != 0:
                # Round up to nearest hour
                tstart = tstart.ceil(freq=freq)
        else:
            origin_tz_aware = pytz.timezone(origin_tz).localize(origin)
            tstart = origin_tz_aware.astimezone(
                pytz.timezone(self.settings.time_settings["timezone"])
            )

        # Coarsen timeresolution

        if method == "nearest":
            df = (
                df.set_index("datetime")
                .groupby("name")
                .resample(freq, origin=tstart)
                .nearest(limit=limit)
            )

        elif method == "bfill":
            df = (
                df.set_index("datetime")
                .groupby("name")
                .resample(freq, origin=tstart)
                .bfill(limit=limit)
            )

        else:
            print(f"The coarsening method: {method}, is not implemented yet.")
            df = df.set_index(["name", "datetime"])

        if "name" in df.columns:
            df = df.drop(columns=["name"])

        # Update resolution info in metadf
        self.metadf["dataset_resolution"] = pd.to_timedelta(freq)
        # update df
        self.df = df

        # Remove gaps and missing from the observatios
        # most gaps and missing are already removed but when increasing timeres,
        # some records should be removed as well.
        self.df = remove_gaps_from_obs(gaplist = self.gaps, obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)


    def sync_observations(self, tollerance, verbose=True):
        """
        Simplify and syncronize the observation timestamps along different stations.

        To simplify the resolution (per station), a tollerance is use to shift timestamps. The tollerance indicates the
        maximum translation in time that can be applied to an observation.

        The sycronisation tries to group stations that have an equal simplified resolution, and syncronize them. The origin
        of the sycronized timestamps will be set to round hours, round 10-minutes or round-5 minutes if possible given the tollerance.

        The observations present in the input file are used.

        After syncronization, the IO outliers, missing observations and gaps are recomputed.

        Parameters
        ----------
        tollerance, Timedelta or str
            The tollerance string or object representing the maximum translation in time.
            Ex: '5T' is 5 minuts, '1H', is one hour.
        verbose : bool, optional
            If True, a dataframe illustrating the mapping from original datetimes to simplified and syncronized is returned. The default is True.

        Note
        --------
        Keep in mind that this method will overwrite the df, outliersdf, missing timestamps and gaps.

        Note
        --------
        Because the used observations are from the input file, previously coarsend timeresolutions are ignored.


        Returns
        -------
        pandas.DataFrame (if verbose is True)
            A dataframe containing the original observations with original timestamps and the corresponding target timestamps.

        """



        df = self.input_df

        self.df = init_multiindexdf()
        self.outliersdf = init_triple_multiindexdf()
        self.gapfilldf = init_multiindexdf()
        self.missing_obs = None
        self.gaps = None

        # find simplified resolution
        simplified_resolution = get_freqency_series(
            df=df, method="median", simplify=True, max_simplify_error=tollerance
        )

        occuring_resolutions = simplified_resolution.unique()

        df = df.reset_index()

        def find_simple_origin(tstart, tollerance):
            if tstart.minute == 0 and tstart.second == 0 and tstart.microsecond == 0:
                return tstart  # already a round hour

            # try converting to a round hour
            tstart_round_hour = tstart.round("60min")
            if abs(tstart - tstart_round_hour) <= pd.to_timedelta(tollerance):
                return tstart_round_hour

            # try converting to a tenfold in minutes
            tstart_round_tenfold = tstart.round("10min")
            if abs(tstart - tstart_round_tenfold) <= pd.to_timedelta(tollerance):
                return tstart_round_tenfold

            # try converting to a fivefold in minutes
            tstart_round_fivefold = tstart.round("5min")

            if abs(tstart - tstart_round_fivefold) <= pd.to_timedelta(tollerance):
                return tstart_round_fivefold

            # no suitable conversion found
            return tstart

        merged_df = pd.DataFrame()
        _total_verbose_df = pd.DataFrame()
        for occur_res in occuring_resolutions:
            group_stations = simplified_resolution[
                simplified_resolution == occur_res
            ].index.to_list()
            print(
                f" Grouping stations with simplified resolution of {pd.to_timedelta(occur_res)}: {group_stations}"
            )
            groupdf = df[df["name"].isin(group_stations)]

            tstart = groupdf["datetime"].min()
            tend = groupdf["datetime"].max()

            # find a good origin point
            origin = find_simple_origin(tstart=tstart, tollerance=tollerance)

            # Create records index
            target_records = pd.date_range(
                start=origin, end=tend, freq=pd.Timedelta(occur_res)
            ).to_series()

            target_records.name = "target_datetime"
            # convert records to new target records, station per station

            for sta in group_stations:
                stadf = groupdf[groupdf["name"] == sta]
                # Drop all nan values! these will be added later from the outliersdf
                stadf = stadf.set_index(["name", "datetime"])
                stadf = stadf.dropna(axis=0, how="all")
                stadf = stadf.reset_index()

                mergedstadf = pd.merge_asof(
                    left=stadf.sort_values("datetime"),
                    right=target_records.to_frame(),
                    right_on="target_datetime",
                    left_on="datetime",
                    direction="nearest",
                    tolerance=pd.Timedelta(tollerance),
                )

                # possibility 1: record is mapped crrectly
                correct_mapped = mergedstadf[~mergedstadf["target_datetime"].isnull()]


                # possibility2: records that ar not mapped to target
                # not_mapped_records =mergedstadf[mergedstadf['target_datetime'].isnull()]


                # possibilyt 3 : no suitable candidates found for the target
                # these will be cached by the missing and gap check
                # no_record_candidates = target_records[~target_records.isin(mergedstadf['target_datetime'])].values


                merged_df = pd.concat([merged_df, correct_mapped])
                if verbose:
                    _total_verbose_df = pd.concat([_total_verbose_df, mergedstadf])

        # overwrite the df with the synced observations
        merged_df = (
            merged_df.rename(
                columns={"datetime": "original_datetime", "target_datetime": "datetime"}
            )
            .set_index(["name", "datetime"])
            .drop(["original_datetime"], errors="ignore", axis=1)
            .sort_index()
        )
        # self.df = merged_df

        # Recompute the dataset attributes, apply qc, gap and missing searches, etc.
        self._construct_dataset(
            df=merged_df,
            freq_estimation_method="highest",
            freq_estimation_simplify=False,
            freq_estimation_simplify_error=None,
            fixed_freq_series=simplified_resolution,
            update_full_metadf=False,
        )  # Do not overwrite full metadf, only the frequencies

        if verbose:
            _total_verbose_df = _total_verbose_df.rename(
                columns={"datetime": "original_datetime", "target_datetime": "datetime"}
            ).set_index(["name", "datetime"])
            return _total_verbose_df

    def import_data_from_file(
        self,
        long_format=True,
        obstype=None,
        freq_estimation_method=None,
        freq_estimation_simplify=None,
        freq_estimation_simplify_error=None,
    ):

        """
        Read observations from a csv file as defined in the
        Settings.input_file. The input file columns should have a template
        that is stored in Settings.template_list.

        If the metadata is stored in a seperate file, and the
        Settings.input_metadata_file is correct, than this metadata is also
        imported (if a suitable template is in the Settings.template_list.)

        The dataset is by default assumed to be in long-format (each column represent an observation type, one column indicates the stationname).
        Wide-format can be used if 'long_format' is set to False and if the observation type is specified by obstype.

        An estimation of the observational frequency is made per station. This is used
        to find missing observations and gaps.


        The Dataset attributes are set and the following checks are executed:
                * Duplicate check
                * Invalid input check
                * Find missing observations
                * Find gaps


        Parameters
        ----------
        long_format : bool, optional
            True if the inputdata has a long-format, False if it has a wide-format. The default is True.
        obstype : str, optional
            If the dataformat is wide, specify which observation type the
            observations represent. The obstype should be an element of
            metobs_toolkit.observation_types. The default is None.
        freq_estimation_method : 'highest' or 'median', optional
            Select wich method to use for the frequency estimation. If
            'highest', the highest apearing frequency is used. If 'median', the
            median of the apearing frequencies is used. If None, the method
            stored in the
            Dataset.settings.time_settings['freq_estimation_method'] is used.
            The default is None.
        freq_estimation_simplify : bool, optional
            If True, the likely frequency is converted to round hours, or round minutes.
            The "freq_estimation_simplify_error' is used as a constrain. If the constrain is not met,
            the simplification is not performed. If None, the method
            stored in the
            Dataset.settings.time_settings['freq_estimation_simplify'] is used.
            The default is None.
        freq_estimation_simplify_error : Timedelta or str, optional
            The tollerance string or object representing the maximum translation in time to form a simplified frequency estimation.
            Ex: '5T' is 5 minuts, '1H', is one hour. If None, the method
            stored in the
            Dataset.settings.time_settings['freq_estimation_simplify_error'] is
            used. The default is None.

        Returns
        -------
        None.

        """

        print("Settings input data file: ", self.settings.IO["input_data_file"])
        logger.info(f'Importing data from file: {self.settings.IO["input_data_file"]}')

        if freq_estimation_method is None:

            freq_estimation_method = self.settings.time_settings[
                "freq_estimation_method"
            ]
        if freq_estimation_simplify is None:
            freq_estimation_simplify = self.settings.time_settings[
                "freq_estimation_simplify"
            ]
        if freq_estimation_simplify_error is None:
            freq_estimation_simplify_error = self.settings.time_settings[
                "freq_estimation_simplify_error"
            ]

        # check if obstype is valid
        if not obstype is None:
            assert (
                obstype in observation_types
            ), f'{obstype} is not a default obstype. Use one of: {self.settings.app["observation_types"]}'


        # Read observations into pandas dataframe
        df, template = import_data_from_csv(
            input_file=self.settings.IO["input_data_file"],
            template_file=self.settings.templates["data_template_file"],
            long_format=long_format,
            obstype=obstype,  # only relevant in wide format
        )


        # Set timezone information
        df.index = df.index.tz_localize(
            tz=self.settings.time_settings["timezone"],
            ambiguous="infer",
            nonexistent="shift_forward",
        )

        logger.debug(
            f'Data from {self.settings.IO["input_data_file"]} \
                     imported to dataframe.'
        )

        # drop Nat datetimes if present
        df = df.loc[pd.notnull(df.index)]

        if not "name" in df.columns:
            logger.warning(
                f'No station names find in the observations! \
                           Assume the dataset is for ONE station with the \
                         default name: {self.settings.app["default_name"]}.'
            )
            df["name"] = str(self.settings.app["default_name"])

        if self.settings.IO["input_metadata_file"] is None:
            print(
                "WARNING: No metadata file is defined.\
                  Add your settings object."
            )
            logger.warning(
                "No metadata file is defined,\
                    no metadata attributes can be set!"
            )
        else:
            logger.info(
                f'Importing metadata from file:\
                        {self.settings.IO["input_metadata_file"]}'
            )
            meta_df = import_metadata_from_csv(
                input_file=self.settings.IO["input_metadata_file"],
                template_file=self.settings.templates["metadata_template_file"],
            )

            # merge additional metadata to observations
            meta_cols = [
                colname for colname in meta_df.columns if not colname.startswith("_")
            ]
            additional_meta_cols = list(set(meta_cols).difference(df.columns))

            if bool(additional_meta_cols):
                logger.debug(
                    f"Merging metadata ({additional_meta_cols})\
                             to dataset data by name."
                )
                additional_meta_cols.append("name")  # merging on name
                # merge deletes datetime index somehow? so add it back.
                df_index = df.index
                df = df.merge(
                    right=meta_df[additional_meta_cols], how="left", on="name"
                )
                df.index = df_index

        # update dataset object
        self.data_template = pd.DataFrame().from_dict(template)

        # convert dataframe to multiindex (datetime - name)
        df = df.set_index(["name", df.index])

        # dataframe with all data of input file
        self.input_df = df.sort_index()

        self._construct_dataset(
            df=df,
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify=freq_estimation_simplify,
            freq_estimation_simplify_error=freq_estimation_simplify_error,
        )

    def import_data_from_database(
        self, start_datetime=None, end_datetime=None, coarsen_timeres=False
    ):
        """
        Function to import data directly from the framboos database and
        updating the network and station objects.


        Parameters
        ----------

        start_datetime : datetime, optional
            Start datetime of the observations. The default is None and using
            yesterday's midnight.
        end_datetime : datetime, optional
            End datetime of the observations. The default is None and using
            todays midnight.
        coarsen_timeres : Bool, optional
            If True, the observations will be interpolated to a coarser
            time resolution as is defined in the Settings. The default
            is False.

        Returns
        ----------

        None.

        Note
        ----------
        A Ugent VPN connection must be present, as well as the username and password
        stored in the settings.

        """
        if start_datetime is None:
            start_datetime = datetime.date.today() - datetime.timedelta(days=1)
        if end_datetime is None:
            end_datetime = datetime.date.today()

        # Read observations into pandas dataframe
        df = import_data_from_db(
            self.settings.db, start_datetime=start_datetime, end_datetime=end_datetime
        )

        if df.empty:  # No data has, probably connection error
            return

        # Make data template
        self.data_template = pd.DataFrame().from_dict(
            template_to_package_space(self.settings.db["vlinder_db_obs_template"])
        )

        # convert dataframe to multiindex (datetime - name)
        df = df.set_index(["name", df.index])
        df = df.sort_index()

        # If an ID has changed or not present in the metadatafile,
        # the stationname and metadata is Nan
        # These observations will be removed
        unknown_obs = df[df.index.get_level_values("name").isnull()]
        if not unknown_obs.empty:
            logger.warning(
                "There is an unknown station in the dataset \
                           (probaply due to an ID that is not present in \
                           the metadata file). This will be removed from the dataset."
            )
            df = df[~df.index.get_level_values("name").isnull()]
        self._construct_dataset(df)


    def _construct_dataset(
        self,
        df,
        freq_estimation_method,
        freq_estimation_simplify,
        freq_estimation_simplify_error,
        fixed_freq_series=None,
        update_full_metadf=True,
    ):

        """
        Helper function to construct the Dataset class from a IO dataframe.

        The df, metadf, outliersdf, gaps and missing timestamps attributes are set.

        Qc on IO is applied (duplicated check and invalid check) + gaps and missing
        values are defined by assuming a frequency per station.

        Parameters
        ----------
        df : pandas.dataframe
            The dataframe containing the input observations and metadata.
        freq_estimation_method : 'highest' or 'median'
            Select wich method to use for the frequency estimation. If
            'highest', the highest apearing frequency is used. If 'median', the
            median of the apearing frequencies is used.
        freq_estimation_simplify : bool
            If True, the likely frequency is converted to round hours, or round minutes.
            The "freq_estimation_simplify_error' is used as a constrain. If the constrain is not met,
            the simplification is not performed.
        freq_estimation_simplify_error : Timedelta or str, optional
            The tollerance string or object representing the maximum translation in time to form a simplified frequency estimation.
            Ex: '5T' is 5 minuts, '1H', is one hour.
        fixed_freq_series : pandas.series or None, optional
            If you do not want the frequencies to be recalculated, one can pass the
            frequency series to update the metadf["dataset_resolution"]. If None, the frequencies will be estimated. The default is None.
        update_full_metadf : bool, optional
            If True, the full Dataset.metadf will be updated. If False, only the frequency columns in the Dataset.metadf will be updated. The default is True.


        Returns
        -------
        None.

        """


        # Convert dataframe to dataset attributes
        self._initiate_df_attribute(dataframe=df, update_metadf=update_full_metadf)

        # Apply quality control on Import resolution
        self._apply_qc_on_import()


        if fixed_freq_series is None:
            freq_series = get_freqency_series(
                df=self.df,
                method=freq_estimation_method,
                simplify=freq_estimation_simplify,
                max_simplify_error=freq_estimation_simplify_error,
            )

            freq_series_import = freq_series

        else:
            if "assumed_import_frequency" in self.metadf.columns:
                freq_series_import = self.metadf[
                    "assumed_import_frequency"
                ]  # No update
            else:
                freq_series_import = fixed_freq_series
            freq_series = fixed_freq_series


        # add import frequencies to metadf (after import qc!)
        self.metadf["assumed_import_frequency"] = freq_series_import

        self.metadf["dataset_resolution"] = freq_series

        # Remove gaps and missing from the observations AFTER timecoarsening
        self.df = remove_gaps_from_obs(gaplist = self.gaps, obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)


    def _initiate_df_attribute(self, dataframe, update_metadf=True):
        logger.info(
            f"Updating dataset by dataframe with shape:\
                    {dataframe.shape}."
        )

        # Create dataframe with fixed order of observational columns
        obs_col_order = [col for col in observation_types if col in dataframe.columns]

        self.df = dataframe[obs_col_order].sort_index()

        if update_metadf:
            # create metadataframe with fixed number and order of columns
            metadf = dataframe.reindex(columns=self.settings.app["location_info"])
            metadf.index = metadf.index.droplevel("datetime")  # drop datetimeindex
            # drop dubplicates due to datetime
            metadf = metadf[~metadf.index.duplicated(keep="first")]

            self.metadf = metadf_to_gdf(metadf)


    def _apply_qc_on_import(self):
        # find missing obs and gaps, and remove them from the df
        self.missing_obs, self.gaps = missing_timestamp_and_gap_check(
            df=self.df,
            gapsize_n=self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"],
        )

        # Create gaps and missing obs objects
        # self.gaps = gaps_list
        # self.missing_obs = Missingob_collection(missing_obs)

        # Perform QC checks on original observation frequencies
        self.df, dup_outl_df = duplicate_timestamp_check(
            df=self.df,
            checks_info=self.settings.qc["qc_checks_info"],
            checks_settings=self.settings.qc["qc_check_settings"],
        )
        if not dup_outl_df.empty:
            self.update_outliersdf(add_to_outliersdf=dup_outl_df)

        self.df, nan_outl_df = invalid_input_check(
            self.df, checks_info=self.settings.qc["qc_checks_info"]
        )
        if not nan_outl_df.empty:
            self.update_outliersdf(nan_outl_df)

        self.outliersdf = self.outliersdf.sort_index()

        # update the order and which qc is applied on which obstype
        checked_obstypes = [obs for obs in self.df.columns if obs in observation_types]

        checknames = ["duplicated_timestamp", "invalid_input"]  # KEEP order

        self._applied_qc = pd.concat(
            [
                self._applied_qc,
                conv_applied_qc_to_df(
                    obstypes=checked_obstypes, ordered_checknames=checknames
                ),
            ],
            ignore_index=True,
        )

    # =============================================================================
    # Physiography extractions
    # =============================================================================

    def get_lcz(self):
        """
        Function to extract the Local CLimate zones (LCZ) from the 
        wudapt global LCZ map on the Google engine for all stations.

        A 'LCZ' column will be added to the metadf, and series is returned.

        Returns
        -------
        lcz_series : pandas.Series()
            A series with the stationnames as index and the LCZ as values.

        """





        # connect to gee
        connect_to_gee()

        # Extract LCZ for all stations
        lcz_series = lcz_extractor(
            metadf=self.metadf,
            mapinfo=self.settings.gee["gee_dataset_info"]["global_lcz_map"],
        )

        # drop column if it was already present
        if "lcz" in self.metadf:
            self.metadf = self.metadf.drop(columns=["lcz"])

        # update metadata
        self.metadf = self.metadf.merge(
            lcz_series.to_frame(), how="left", left_index=True, right_index=True
        )
        return lcz_series

    def get_altitude(self):
        # connect to gee
        connect_to_gee()

        # Extract LCZ for all stations
        altitude_series = height_extractor(
            metadf=self.metadf, mapinfo=self.settings.gee["gee_dataset_info"]["DEM"]
        )

        # drop column if it was already present
        if "altitude" in self.metadf:
            self.metadf = self.metadf.drop(columns=["altitude"])

        # update metadata
        self.metadf = self.metadf.merge(
            altitude_series.to_frame(), how="left", left_index=True, right_index=True
        )
        return altitude_series

    def get_landcover(self, buffers=[100], aggregate=True, overwrite=True, gee_map='worldcover'):
        """
        Extract the landcover fractions in a buffer with a specific radius for
        all stations. If an aggregation scheme is define, one can choose to
        aggregate the landcoverclasses.

        The landcover fractions will be added to the Dataset.metadf if overwrite
        is True. Presented as seperate columns where each column represent the
        landcovertype and corresponding buffer.




        Parameters
        ----------
        buffers : num, optional
            The list of buffer radia in dataset units (meters for ESA worldcover) . The default is 100.
        aggregate : bool, optional
            If True, the classes will be aggregated with the corresponding
            aggregation scheme. The default is True.
        overwrite : bool, optional
            If True, the Datset.metadf will be updated with the generated
            landcoverfractions. The default is True.
        gee_map : str, optional
            The name of the dataset to use. This name should be present in the
            settings.gee['gee_dataset_info']. If aggregat is True, an aggregation
            scheme should included as well. The default is 'worldcover'

        Returns
        -------
        frac_df : pandas.DataFrame
            A Dataframe with index: name, buffer_radius and the columns are the
            fractions.

        """
        # connect to gee
        connect_to_gee()

        df_list = []
        for buffer in buffers:

            print(f'Extracting landcover from {gee_map} with buffer radius = {buffer}')
            # Extract landcover fractions for all stations
            lc_frac_df, buffer = lc_fractions_extractor(
                metadf=self.metadf,
                mapinfo=self.settings.gee["gee_dataset_info"][gee_map],
                buffer=buffer,
                agg=aggregate,
            )

            # add buffer to the index
            lc_frac_df['buffer_radius'] = buffer
            lc_frac_df = lc_frac_df.reset_index().set_index(['name', 'buffer_radius'])
            lc_frac_df = lc_frac_df.sort_index()

            # add to the list
            df_list.append(lc_frac_df)


        # concat all df for different buffers to one
        frac_df = pd.concat(df_list)
        frac_df = frac_df.sort_index()


        if overwrite:

            for buf in frac_df.index.get_level_values('buffer_radius').unique():
                buf_df = frac_df.xs(buf, level='buffer_radius')
                buf_df.columns= [col + f'_{int(buf)}m' for col in buf_df.columns]

                # overwrite the columns or add them if they did not exist
                self.metadf[buf_df.columns] = buf_df

        return frac_df
