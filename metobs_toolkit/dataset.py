#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""

import os
from datetime import datetime
from pytz import all_timezones, common_timezones
import logging
import pandas as pd
import numpy as np


from metobs_toolkit.settings import Settings
from metobs_toolkit.data_import import (import_data_from_csv,
                                         import_data_from_db,
                                         template_to_package_space,
                                         import_metadata_from_csv)

from metobs_toolkit.printing import print_dataset_info
from metobs_toolkit.landcover_functions import (connect_to_gee,
                                                 lcz_extractor,
                                                 height_extractor,
                                                 lc_fractions_extractor,
                                                 )

from metobs_toolkit.plotting_functions import (geospatial_plot,
                                                timeseries_plot,
                                                qc_stats_pie)

from metobs_toolkit.qc_checks import (gross_value_check,
                                       persistance_check,
                                       repetitions_check,
                                       duplicate_timestamp_check,
                                       step_check,
                                       window_variation_check,
                                       invalid_input_check)


from metobs_toolkit.qc_statistics import get_freq_statistics
from metobs_toolkit.writing_files import write_dataset_to_csv

from metobs_toolkit.missingobs import Missingob_collection

from metobs_toolkit.gap import (Gap_collection,
                                 missing_timestamp_and_gap_check,
                                 get_freqency_series)


from metobs_toolkit.df_helpers import (add_final_label_to_outliersdf,
                                        multiindexdf_datetime_subsetting,
                                        remove_outliers_from_obs,
                                        init_multiindexdf,
                                        metadf_to_gdf)


from metobs_toolkit.modeldata import Modeldata


logger = logging.getLogger(__name__)




# =============================================================================
# Dataset class
# =============================================================================

class Dataset:
    def __init__(self):
        """
        Constructs all the necessary attributes for Dataset object.

        """
        logger.info('Initialise dataset')

        # Dataset with 'good' observations
        self.df = pd.DataFrame()

        # Dataset with outlier observations
        self.outliersdf = init_multiindexdf()

        self.missing_obs = None  # becomes a Missingob_collection after import
        self.gaps = None  # becomes a gap_collection after import

        self.gapfilldf = init_multiindexdf()


        # Dataset with metadata (static)
        self.metadf = pd.DataFrame()
        # dataframe containing all information on the description and mapping
        self.data_template = pd.DataFrame()

        self._istype = 'Dataset'
        self._freqs = pd.Series(dtype=object)


        self.settings = Settings()


    def update_settings(self, output_folder=None, input_data_file=None,
                  input_metadata_file=None, data_template_file=None,
                  metadata_template_file=None):

        """
        Update the most common input-output (IO) settings.
        (This should be applied before importing the observations.)

        When an update value is None, the specific setting will not be updated.


        :param output_folder: a directory to store the output to, defaults to None.
        :type output_folder: String, optional
        :param input_data_file: Path to the input data file with observations, defaults to None.
        :type input_data_file: String, optional
        :param input_metadata_file: Path to the input metadata file, defaults to None
        :type input_metadata_file: String, optional
        :param data_template_file: Path to the mapper-template csv file to be used on the observations.
        If not given, the default template is used.
        :type data_template_file: String, optional
        :param metadata_template_file: Path to the mapper-template csv file to be used on the metadata.
        If not given, the default template is used.
        :type metadata_template_file: String, optional
        :return: No return
        :rtype: No return
        """

        self.settings.update_IO(output_folder=output_folder,
                                input_data_file=input_data_file,
                                input_metadata_file=input_metadata_file,
                                data_template_file=data_template_file,
                                metadata_template_file=metadata_template_file)

    def update_timezone(self, timezonestr):
        """
        Change the timezone of the input data. By default the Brussels timezone is assumed.
        A valid timezonestring is an element of the pytz.all_timezones.

        :param timezonestr: Timezone string of the input observations.
        :type timezonestr: String
        :return: None
        :rtype: None
        """
        self.settings.update_timezone(timezonestr)

    def update_default_name(self, default_name):
        """
        Update the default name (the name of the station). This name will be
        used when no names are found in the observational dataset.

        (All observations are assumed to come from one station.)

        :param default_name: Default name to use when no names are present in the data.
        :type default_name: String
        :return: None
        :rtype: None

        """
        self.settings.app['default_name'] = str(default_name)

    def show_settings(self):
        """
        A function that prints out all the settings, structured per thematic.
        :return: No return
        :rtype: No return

        """
        self.settings.show()

    def get_station(self, stationname):
        """
        Extract a station object from the dataset by name.

        :param stationname: The name of the station
        :type stationname: String
        :return: The station object.
        :rtype: metobs_toolkit.Station

        """

        logger.info(f'Extract {stationname} from dataset.')

        # important: make shure all station attributes are of the same time as dataset.
        # so that all methods can be inherited.

        try:
            sta_df = self.df.xs(stationname, level='name', drop_level=False)
            sta_metadf = self.metadf.loc[stationname].to_frame().transpose()
        except KeyError:
            logger.warning(f'{stationname} not found in the dataset.')
            print(f'{stationname} not found in the dataset.')
            return None

        try:
            sta_outliers = self.outliersdf.xs(
                stationname, level='name', drop_level=False)
        except KeyError:
            sta_outliers = init_multiindexdf()

        sta_gaps = self.gaps.get_station_gaps(stationname)
        sta_missingobs = self.missing_obs.get_station_missingobs(stationname)

        try:
            sta_gapfill=self.gapfilldf.xs(
                stationname, level='name', drop_level=False)
        except KeyError:
            sta_gapfill = init_multiindexdf()

        return Station(name=stationname,
                       df=sta_df,
                       outliersdf=sta_outliers,
                       gaps=sta_gaps,
                       missing_obs=sta_missingobs,
                       gapfilldf=sta_gapfill,
                       metadf=sta_metadf,
                       data_template=self.data_template,
                       settings = self.settings)

    def show(self):
        """
        A function to print out some overview information about the Dataset.

        :return: No return
        :rtype: No return
        """

        logger.info('Show basic info of dataset.')

        try:
            gapsdf = self.gaps.to_df()
        except:
            gapsdf=init_multiindexdf()


        print_dataset_info(self.df, self.outliersdf, gapsdf,
                           self.settings.app['print_fmt_datetime'])

    def make_plot(self, stationnames=None, obstype='temp', colorby='name',
                  starttime=None, endtime=None,
                  title=None, legend=True, show_outliers=True):
        """
        This function creates a timeseries plot for the dataset. The variable observation type
        is plotted for all stationnames from a starttime to an endtime.


        :param stationnames: A list with stationnames to include in the timeseries. If None is given, all the stations are used, defaults to None.
        :type stationnames: list or None, optional
        :param obstype: The observation type to plot, defaults to 'temp'.
        :type variable: String, optional
        :param colorby: Indicate how colors should be assigned to the lines. 'label' will color the lines by their quality control label. 'name' will color by each station, defaults to 'name'.
        :type colorby: 'label' or 'name', optional
        :param starttime: Specifiy the start datetime for the plot. If None is given it will use the start datetime of the dataset, defaults to None.
        :type starttime: datetime.datetime, optional
        :param endtime: Specifiy the end datetime for the plot. If None is given it will use the end datetime of the dataset, defaults to None.
        :type endtime: datetime.datetime, optional
        :param title: The title for the plot. If None is given a default title is used, defaults to None.
        :type title: String, optional
        :param legend: When thrue a legend is added to the plot, defaults to True
        :type legend: Boolean, optional
        :param show_outliers: If true the observations labeld as outliers will be included in the plot, defaults to True
        :type show_outliers: boolean, optional
        :return: The plot axis
        :rtype: matplotlib.pyplot.axis

        """

        logger.info(f'Make {obstype}-timeseries plot for {stationnames}')

        # combine all dataframes
        mergedf = self.combine_all_to_obsspace()

        # Subset on stationnames
        if not isinstance(stationnames, type(None)):
            mergedf = mergedf.loc[mergedf.index.get_level_values(
                'name').isin(stationnames)]

        # Subset on start and endtime
        mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

        # Get plot styling attributes
        if isinstance(title, type(None)):
            if isinstance(stationnames, type(None)):
                if self._istype == 'Dataset':
                    title = self.settings.app['display_name_mapper'][obstype] + \
                        ' for all stations. '
                elif self._istype == 'Station':
                    title = self.settings.app['display_name_mapper'][obstype] + \
                        ' of ' + self.name

            else:
                title = self.settings.app['display_name_mapper'][obstype] + \
                    ' for stations: ' + str(stationnames)

        if ((obstype+'_final_label' not in mergedf.columns) &
            ((colorby == 'label') | (show_outliers))):
            # user whant outier information but no QC is applied on this obstype
            print(f' No quality control is applied on {obstype}! \
                  No outlier information is available.')
            print('Colorby is set to "name" and show_outliers \
                  is set to False.')
            colorby = 'name'
            show_outliers = False

        # Make plot
        ax = timeseries_plot(mergedf=mergedf,
                             obstype=obstype,
                             title=title,
                             xlabel='Timestamp',
                             ylabel=self.data_template[obstype]['orig_name'],
                             colorby=colorby,
                             show_legend=legend,
                             show_outliers=show_outliers,
                             plot_settings = self.settings.app['plot_settings'],
                             gap_settings = self.settings.gap,
                             qc_info_settings=self.settings.qc['qc_checks_info'])

        return ax

    def make_geo_plot(self, obstype='temp', title=None,
                      timeinstance=None, legend=True,
                      vmin=None, vmax=None):
        """
        This functions creates a geospatial plot for a field
        (observations or attributes) of all stations.

        If the field is timedepending, than the timeinstance is used to plot
        the field status at that datetime.

        If the field is categorical than the leged will have categorical
        values, else a colorbar is used.

        All styling attributes are extracted from the Settings.


        Parameters

        variable : String, optional
            Fieldname to visualise. This can be an observation or station
            attribute. The default is 'temp'.
        title : String, optional
            Title of the figure, if None a default title is generated. The
            default is None.
        timeinstance : datetime, optional
            Datetime moment of the geospatial plot. The default is None and the
            first datetime available is used.
        legend : Bool, optional
            Add legend to the figure. The default is True.
        vmin : float, optional
            The minimum value corresponding to the minimum color.
            The default is None and the minimum of the variable is used.
        vmax : float, optional
           The maximum value corresponding to the minimum color. The default is
           None and the maximum of the variable is used.

        Returns

        ax : Geoaxes
            The geoaxes is returned.

        """

        # Load default plot settings
        # default_settings=Settings.plot_settings['spatial_geo']

        # get first timeinstance of the dataset if not given
        if isinstance(timeinstance, type(None)):
            timeinstance = self.df.index.get_level_values('datetime').min()

        logger.info(f'Make {obstype}-geo plot at {timeinstance}')

        # subset to timeinstance
        plotdf = self.df.xs(timeinstance, level='datetime')

        # merge metadata
        plotdf = plotdf.merge(self.metadf, how='left',
                              left_index=True, right_index=True)

        axis = geospatial_plot(plotdf=plotdf,
                             variable=obstype,
                             timeinstance=timeinstance,
                             title=title,
                             legend=legend,
                             vmin=vmin,
                             vmax=vmax,
                             plotsettings=self.settings.app['plot_settings'],
                             categorical_fields=self.settings.app['categorical_fields'],
                             static_fields = self.settings.app['static_fields'],
                             display_name_mapper=self.settings.app['display_name_mapper'],
                             world_boundaries_map = self.settings.app['world_boundary_map'])

        return axis
    # =============================================================================
    #   Gap Filling
    # =============================================================================
    def get_modeldata(self, modelname='ERA5_hourly', stations=None, startdt=None, enddt=None):
        Modl = Modeldata(modelname)

        #Filters
        if isinstance(startdt, type(None)):
            startdt=self.df.index.get_level_values('datetime').min()
        if isinstance(enddt, type(None)):
            enddt=self.df.index.get_level_values('datetime').max()
        if not isinstance(stations, type(None)):
            if isinstance(stations, str):
                metadf=self.metadf.loc[[stations]]
            if isinstance(stations, list):
                metadf = self.metadf.iloc[self.metadf.index.isin(stations)]
        else:
            metadf = self.metadf


        # fill modell with data
        if modelname == 'ERA5_hourly':
            Modl.get_ERA5_data(metadf, startdt, enddt)

            return Modl
        else:
            print(f"{modelname} for set_modeldata is not implemented yet")
            return None

    # =============================================================================
    #   Gap Filling
    # =============================================================================

    def fill_gaps_linear(self, obstype='temp'):
        #TODO logging
        fill_settings =self.settings.gap['gaps_fill_settings']['linear']
        fill_info = self.settings.gap['gaps_fill_info']


        #fill gaps
        self.gapfilldf[obstype] = self.gaps.apply_interpolate_gaps(
                                    obsdf = self.df,
                                    outliersdf = self.outliersdf,
                                    dataset_res=self.metadf['dataset_resolution'],
                                    obstype=obstype,
                                    method=fill_settings['method'],
                                    max_consec_fill=fill_settings['max_consec_fill'])

        #add label column
        self.gapfilldf[obstype + '_' + fill_info['label_columnname']] = fill_info['label']['linear']


    def fill_gaps_era5(self, modeldata, method='debias', obstype='temp', overwrite=True):

        fill_info = self.settings.gap['gaps_fill_info']

        # check if modeldata is available
        if isinstance(modeldata, type(None)):
            print('The dataset has no modeldate. Use the set_modeldata() function to add modeldata.')
            return None
        # check if obstype is present in eramodel
        assert obstype in modeldata.df.columns, f'{obstype} is not present in the modeldate: {modeldata}'
        # check if all station are present in eramodeldata
        stations = self.gaps.to_df().index.unique().to_list()
        assert all([sta in modeldata.df.index.get_level_values('name') for sta in stations]),\
            f'Not all stations with gaps are in the modeldata!'

        if not self.gapfilldf.empty:
            if overwrite:
                print('Gapfilldf will be overwritten!')
                self.gapfilldf = init_multiindexdf()
            else:
                print('Gapfilldf is not empty, set "overwrite=True" to overwrite it!')
                print('CANCEL gap fill with ERA5')
                return




        if method=='debias':
            test = self.gaps.apply_debias_era5_gapfill(dataset=self,
                                                       eraModelData=modeldata,
                                                       obstype=obstype,
                                                       debias_settings=self.settings.gap['gaps_fill_settings']['model_debias'])

            self.gapfilldf[obstype] = test
            #add label column
            self.gapfilldf[obstype + '_' + fill_info['label_columnname']] = fill_info['label']['model_debias']
        else:
            print('not implemented yet')



    def write_to_csv(self, filename=None, include_outliers=True,
                     include_gapfill=True,
                     add_final_labels=True, use_tlk_obsnames=True):
        """
            Write the dataset to a file where the observations, metadata and
            (if available) the quality labels per observation type are merged
            together.

            A final qualty controll label for each
            quality-controlled-observation type can be added in the outputfile.

            The file will be writen to the Settings.outputfolder.

            Parameters

            filename : string, optional
                The name of the output csv file. If none, a standard-filename
                is generated based on the period of data. The default is None.
            add_final_labels : Bool, optional
                If True, a final qualty control label per observation type
                is added as a column. The default is True.

            Returns

            None

            """

        logger.info('Writing the dataset to a csv file')

        assert not isinstance(self.settings.IO['output_folder'], type(None)), 'Specify \
            Settings.output_folder in order to export a csv.'
        assert os.path.isdir(self.settings.IO['output_folder']), f'The outputfolder: \
            {self.settings.IO["output_folder"]} is not found. '

        # combine all dataframes
        mergedf = self.combine_all_to_obsspace()  # with outliers

        # select which columns to keep
        if include_outliers:
            if not add_final_labels:
                _fin_cols = [col for col in mergedf.columns
                             if col.endswith('_final_label')]
                mergedf = mergedf.drop(columns=_fin_cols)

        else:  # exclude outliers
            if add_final_labels:
                cols_to_keep = [col for col in mergedf.columns
                                if col in self.settings.app['observation_types']]
                cols_to_keep.extend([col for col in mergedf.columns
                                     if col.endswith('_final_label')])
                mergedf = mergedf[cols_to_keep]

        if not include_gapfill:
            # locate all filled values
            filled_df =  init_multiindexdf()
            final_columns = [col for col in mergedf.columns if col.endswith('_final_label')]
            for final_column in final_columns:
                filled_df = pd.concat([filled_df,
                                       mergedf.loc[mergedf[final_column] ==
                                                   self.settings.gaps['gaps_fill_info']['label']]])

            # drop filled values from mergedf
            mergedf = mergedf.drop(filled_df.index, errors='ignore')

            #fill with numpy nan
            nan_columns = {col: np.nan for col in mergedf.columns if col in self.settings.app['observation_types']}
            filled_df = filled_df.assign(**nan_columns)
            # rename label
            filled_df = filled_df.replace({self.settings.gaps['gaps_fill_info']['label']: self.settings.gaps['gaps_info']['gap']['outlier_flag']})
            #add to mergedf
            mergedf = pd.concat([mergedf, filled_df]).sort_index()






        # Map obstypes columns
        if not use_tlk_obsnames:
            # TODO
            print('not implemented yet')

        # TODO Convert units if needed.

        # columns to write
        write_dataset_to_csv(df=mergedf,
                             metadf=self.metadf,
                             filename=filename,
                             outputfolder = self.settings.IO['output_folder'],
                             location_info = self.settings.app['location_info'],
                             observation_types=self.settings.app['observation_types']
                             )


    # =============================================================================
    #     Quality control
    # =============================================================================

    def apply_quality_control(self, obstype='temp',
                              gross_value=True,
                              persistance=True,
                              repetitions=True,
                              step=True,
                              window_variation=True,
                              # internal_consistency=True,
                              ):
        """
        Apply quality control methods to the dataset. The default settings
        are used, and can be changed in the settings_files/qc_settings.py

        The checks are performed in a sequence: gross_vallue -->
        persistance --> ..., Outliers by a previous check are ignored in the
        following checks!

        The dataset is updated inline.

        Parameters

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

        None.

        """

        if repetitions:

            print('Applying the repetitions-check on all stations.')
            logger.info('Applying repetitions check on the full dataset')

            obsdf, outl_df = repetitions_check(
                                            obsdf=self.df,
                                            obstype=obstype,
                                            checks_info=self.settings.qc['qc_checks_info'],
                                            checks_settings = self.settings.qc['qc_check_settings'])

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.update_outliersdf(outl_df)

        if gross_value:
            print('Applying the gross-value-check on all stations.')
            logger.info('Applying gross value check on the full dataset')

            obsdf, outl_df = gross_value_check(
                                            obsdf=self.df,
                                            obstype=obstype,
                                            checks_info=self.settings.qc['qc_checks_info'],
                                            checks_settings = self.settings.qc['qc_check_settings'])

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.update_outliersdf(outl_df)

        if persistance:
            print('Applying the persistance-check on all stations.')
            logger.info('Applying persistance check on the full dataset')

            obsdf, outl_df = persistance_check(
                        station_frequencies=self.metadf['dataset_resolution'],
                        obsdf=self.df,
                        obstype=obstype,
                        checks_info=self.settings.qc['qc_checks_info'],
                        checks_settings = self.settings.qc['qc_check_settings'])

            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.update_outliersdf(outl_df)

        if step:
            print('Applying the step-check on all stations.')
            logger.info('Applying step-check on the full dataset')

            obsdf, outl_df = step_check(obsdf=self.df,
                                        obstype=obstype,
                                        checks_info=self.settings.qc['qc_checks_info'],
                                        checks_settings = self.settings.qc['qc_check_settings'])


            # update the dataset and outliers
            self.df = obsdf
            if  not outl_df.empty:
                self.update_outliersdf(outl_df)

        if window_variation:
            print('Applying the window variation-check on all stations.')
            logger.info('Applying window variation-check on the full dataset')

            obsdf, outl_df = window_variation_check(
                        station_frequencies=self.metadf['dataset_resolution'],
                        obsdf=self.df,
                        obstype=obstype,
                        checks_info=self.settings.qc['qc_checks_info'],
                        checks_settings = self.settings.qc['qc_check_settings'])


            # update the dataset and outliers
            self.df = obsdf
            if not outl_df.empty:
                self.update_outliersdf(outl_df)


        self.outliersdf = self.outliersdf.sort_index()



    def combine_all_to_obsspace(self):
        """
        Combine observations, outliers, gaps and missing timesteps to one
        dataframe in the resolution of the dataset. Final quality labels are
        calculated for all checked obstypes.

        If an observation value exist for an outlier, it will be used in the
        corresponding obstype column.

        Returns

        comb_df : pandas.DataFrame()
            Multi index dataframe with observations and labels.

        """


        outliersdf = self.outliersdf

        # if outliersdf is empty, create columns with 'not checked'
        if outliersdf.empty:
            outliercolumns = [col['label_columnname']
                              for col in self.settings.qc['qc_checks_info'].values()]
            for column in outliercolumns:
                outliersdf[column] = 'not checked'
        # get final label
        outliersdf = add_final_label_to_outliersdf(
                        outliersdf=self.outliersdf,
                        data_res_series=self.metadf['dataset_resolution'],
                        observation_types=self.settings.app['observation_types'],
                        checks_info=self.settings.qc['qc_checks_info'])

        #remove duplicate indixes (needed for update)
        outliersdf = outliersdf[~outliersdf.index.duplicated(keep='first')]


        # add gaps observations and fill with default values
        gapsidx = self.gaps.get_gaps_indx_in_obs_space(
            self.df, self.outliersdf, self.metadf['dataset_resolution'])
        gapsdf = gapsidx.to_frame()

        # add missing observations if they occure in observation space
        missingidx = self.missing_obs.get_missing_indx_in_obs_space(
            self.df, self.metadf['dataset_resolution'])
        missingdf = missingidx.to_frame()

        # add gapfill and remove the filled records from gaps
        gapsfilldf = self.gapfilldf.copy()

        gapsdf = gapsdf.drop(gapsfilldf.index, errors='ignore')
        # missingdf = missingdf.drop(gapsfilldf, errors='ignore') #for future when missing are filled


        # get observations
        df = self.df
        #fill QC outlier columns with custom values
        label_cols = [col for col in outliersdf.columns if col.endswith('_label')]
        for col in label_cols:
            df[col] = 'ok'

        # Merge obs and outliers, where obs values will be updated by outliers
        df.update(other=outliersdf,
                         join='left',
                         overwrite=True,
                         errors='ignore')


        # initiate default values
        for col in df.columns:
            if col in self.settings.app['observation_types']:
                default_value_gap = np.nan  # nan for observations
                default_value_missing = np.nan

            elif col.endswith('_final_label'):
                # 'gap' for final label
                default_value_gap = self.settings.gap['gaps_info']['gap']['outlier_flag']
                # 'is_missing_timestamp' for final label
                default_value_missing = self.settings.gap['gaps_info']['missing_timestamp']['outlier_flag']

            else:
                default_value_gap = 'not checked'
                default_value_missing = 'not checked'
                gapsfilldf[col] = 'not checked'

            gapsdf[col] = default_value_gap
            missingdf[col] = default_value_missing

        # sort columns
        gapsdf = gapsdf[list(df.columns)]
        missingdf = missingdf[list(df.columns)]


        # Merge all together
        comb_df = pd.concat([df, gapsdf, missingdf, gapsfilldf]).sort_index()

        return comb_df


    def get_qc_stats(self, obstype='temp', stationnames=None, make_plot=True):
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

        # drop observation columns that are not obstype
        ignore_obstypes = self.settings.app['observation_types'].copy()
        ignore_obstypes.remove(obstype)
        comb_df = comb_df.drop(columns=ignore_obstypes)

        # drop label columns not applicable on obstype
        relevant_columns = [col for col in comb_df.columns if col.startswith(obstype)]

        # add all columns of checks applied on records (i.g. without obs prefix like duplicate timestamp)
        record_check = {key: item['label_columnname'] for key, item in self.settings.qc['qc_checks_info'].items() if item['apply_on'] == 'record'}
        relevant_columns.extend(list(record_check.values()))
        relevant_columns = [col for col in relevant_columns if col in comb_df.columns]
        # filter relevant columns
        comb_df = comb_df[relevant_columns]


        # compute freq statistics
        final_freq, outl_freq, specific_freq = get_freq_statistics(
            comb_df = comb_df,
            obstype=obstype,
            checks_info=self.settings.qc['qc_checks_info'],
            gaps_info =self.settings.gap['gaps_info'],
            )

        if any([isinstance(stat, type(None)) for stat in [final_freq,
                                                          outl_freq,
                                                          specific_freq]]):

            return None

        if make_plot:
            # make pie plots
            qc_stats_pie(final_stats=final_freq,
                         outlier_stats=outl_freq,
                         specific_stats=specific_freq,
                         plot_settings=self.settings.app['plot_settings'],
                         qc_check_info=self.settings.qc['qc_checks_info'])

        return (final_freq, outl_freq, specific_freq)

    def update_outliersdf(self, add_to_outliersdf):

        # Get the flag column labels and find the newly added columnlabelname
        previous_performed_checks_columns = [
            col for col in self.outliersdf.columns if col.endswith('_label')]
        new_performed_checks_columns = list(set([
                    col for col in add_to_outliersdf.columns
                    if col.endswith('_label')])
                    - set(previous_performed_checks_columns))

        # add to the outliersdf
        self.outliersdf = pd.concat([self.outliersdf, add_to_outliersdf])

        # Fix labels
        self.outliersdf[previous_performed_checks_columns] = self.outliersdf[
                        previous_performed_checks_columns].fillna(value='ok')

        self.outliersdf[new_performed_checks_columns] = self.outliersdf[
            new_performed_checks_columns].fillna(value='not checked')


    # =============================================================================
    #     importing data
    # =============================================================================

    def coarsen_time_resolution(self, freq='1H', method='nearest', limit=1):
        logger.info(f'Coarsening the timeresolution to {freq} using \
                    the {method}-method (with limit={limit}).')
        # TODO: implement buffer method
        # TODO: implement startdt point
        # Coarsen timeresolution
        df = self.df.reset_index()
        if method == 'nearest':
            df = df.set_index('datetime').groupby(
                'name').resample(freq).nearest(limit=limit)

        elif method == 'bfill':
            df = df.set_index('datetime').groupby(
                'name').resample(freq).bfill(limit=limit)

        else:
            print(f'The coarsening method: {method}, is not implemented yet.')
            df = df.set_index(['name', 'datetime'])

        if 'name' in df.columns:
            df = df.drop(columns=['name'])

        # Update resolution info in metadf
        self.metadf['dataset_resolution'] = pd.to_timedelta(freq)
        # update df
        self.df = df

    def import_data_from_file(self, network='vlinder', coarsen_timeres=False):
        """
        Read observations from a csv file as defined in the
        Settings.input_file. The input file columns should have a template
        that is stored in Settings.template_list.

        If the metadata is stored in a seperate file, and the
        Settings.input_metadata_file is correct, than this metadata is also
        imported (if a suitable template is in the Settings.template_list.)


        It is possible to apply a
        resampling (downsampling) of the observations as defined in the settings.

        After the import there is always a call to Dataset.update_dataset_by_df, that
        sets up the dataset with the observations and applies some sanity checks.

        Parameters

        network : String, optional
            The name of the network for these observations. The default
            is 'vlinder'.
        coarsen_timeres : Bool, optional
            If True, the observations will be interpolated to a coarser
            time resolution as is defined in the Settings. The default
            is False.

        Returns

        None.

        """
        print('Settings input data file: ', self.settings.IO['input_data_file'])
        logger.info(f'Importing data from file: {self.settings.IO["input_data_file"]}')

        # Read observations into pandas dataframe
        df, template = import_data_from_csv(
                            input_file=self.settings.IO['input_data_file'],
                            template_file=self.settings.templates['data_template_file'])

        # Set timezone information
        df.index = df.index.tz_localize(tz=self.settings.time_settings['timezone'],
                                        ambiguous='infer',
                                        nonexistent='shift_forward')

        logger.debug(f'Data from {self.settings.IO["input_data_file"]} \
                     imported to dataframe.')

        # drop Nat datetimes if present
        df = df.loc[pd.notnull(df.index)]

        if not 'name' in df.columns:
            logger.warning(f'No station names find in the observations! \
                           Assume the dataset is for ONE station with the \
                         default name: {self.settings.app["default_name"]}.')
            df['name'] =str(self.settings.app["default_name"])

        if isinstance(self.settings.IO['input_metadata_file'], type(None)):
            print('WARNING: No metadata file is defined.\
                  Add your settings object.')
            logger.warning('No metadata file is defined,\
                    no metadata attributes can be set!')
        else:
            logger.info(f'Importing metadata from file:\
                        {self.settings.IO["input_metadata_file"]}')
            meta_df = import_metadata_from_csv(
                        input_file=self.settings.IO["input_metadata_file"],
                        template_file=self.settings.templates['metadata_template_file'])

            # merge additional metadata to observations
            meta_cols = [colname for colname in meta_df.columns
                         if not colname.startswith('_')]
            additional_meta_cols = list(set(meta_cols).difference(df.columns))

            if bool(additional_meta_cols):
                logger.debug(f'Merging metadata ({additional_meta_cols})\
                             to dataset data by name.')
                additional_meta_cols.append('name')  # merging on name
                # merge deletes datetime index somehow? so add it back.
                df_index = df.index
                df = df.merge(right=meta_df[additional_meta_cols],
                              how='left',
                              on='name')
                df.index = df_index

        # update dataset object
        self.data_template = pd.DataFrame().from_dict(template)

        # convert dataframe to multiindex (datetime - name)
        df = df.set_index(['name', df.index])

        # dataframe with all data of input file
        self.input_df = df

        self.update_dataset_by_df(dataframe=df,
                                  coarsen_timeres=coarsen_timeres)

    def import_data_from_database(self,
                                  start_datetime=None,
                                  end_datetime=None,
                                  coarsen_timeres=False):
        """
        Function to import data directly from the framboos database and
        updating the network and station objects.


        Parameters

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

        None.

        """
        if isinstance(start_datetime, type(None)):
            start_datetime = datetime.date.today() - datetime.timedelta(days=1)
        if isinstance(end_datetime, type(None)):
            end_datetime = datetime.date.today()

        # Read observations into pandas dataframe
        df = import_data_from_db(self.settings.db,
                                start_datetime=start_datetime,
                                end_datetime=end_datetime)

        if df.empty: #No data has, probably connection error
            return

        # Make data template
        self.data_template = pd.DataFrame().from_dict(
            template_to_package_space(self.settings.db['vlinder_db_obs_template']))

        # convert dataframe to multiindex (datetime - name)
        df = df.set_index(['name', df.index])
        df = df.sort_index()

        # If an ID has changed or not present in the metadatafile,
        # the stationname and metadata is Nan
        # These observations will be removed
        unknown_obs = df[df.index.get_level_values('name').isnull()]
        if not unknown_obs.empty:
            logger.warning('There is an unknown station in the dataset \
                           (probaply due to an ID that is not present in \
                           the metadata file). This will be removed.')
            df = df[~df.index.get_level_values('name').isnull()]

        self.update_dataset_by_df(
            dataframe=df, coarsen_timeres=coarsen_timeres)

    def update_dataset_by_df(self, dataframe, coarsen_timeres=False):

        """
        Update the dataset object by a dataframe.


        When filling the observations, there is an automatic check for
        missing timestamps and duplicating timestamps. If a missing timestamp
        is detected, the timestamp is created with Nan values for all
        observation types.

        If metadata is present and a LCZ-tiff file availible, than the LCZ's
        of the stations are computed.

        Parameters

        dataframe : pandas.DataFrame
        A dataframe that has an datetimeindex and following columns:
            'name, temp, radiation_temp, humidity, ...'


        Returns

        None.

        """


        logger.info(f'Updating dataset by dataframe with shape:\
                    {dataframe.shape}.')

        # Create dataframe with fixed number and order of observational columns
        df = dataframe.reindex(columns=self.settings.app['observation_types'])
        self.df = df

        # create metadataframe with fixed number and order of columns
        metadf = dataframe.reindex(columns=self.settings.app['location_info'])
        metadf.index = metadf.index.droplevel('datetime')  # drop datetimeindex
        # drop dubplicates due to datetime
        metadf = metadf[~metadf.index.duplicated(keep='first')]

        self.metadf = metadf_to_gdf(metadf)

        # add import frequencies to metadf
        self.metadf['assumed_import_frequency'] = get_freqency_series(self.df)
        self.df = df.sort_index()
        self.original_df = df.sort_index()

        # find missing obs and gaps, and remove them from the df
        self.df, missing_obs, gaps_df = missing_timestamp_and_gap_check(
            df=self.df,
            gapsize_n=self.settings.gap['gaps_settings']['gaps_finder']['gapsize_n'])

        # Create gaps and missing obs objects
        self.gaps = Gap_collection(gaps_df)
        self.missing_obs = Missingob_collection(missing_obs)

        # Perform QC checks on original observation frequencies

        self.df, dup_outl_df = duplicate_timestamp_check(df=self.df,
                                                         checks_info=self.settings.qc['qc_checks_info'],
                                                         checks_settings = self.settings.qc['qc_check_settings'])
        if not dup_outl_df.empty:
            self.update_outliersdf(dup_outl_df)

        self.df, nan_outl_df = invalid_input_check(self.df,
                                                   checks_info=self.settings.qc['qc_checks_info'])
        if not nan_outl_df.empty:
            self.update_outliersdf(nan_outl_df)


        if coarsen_timeres:
            self.coarsen_time_resolution(freq=self.settings.time_settings['target_time_res'],
                                          method=self.settings.time_settings['resample_method'],
                                          limit=self.settings.time_settings['resample_limit'])


        else:
            self.metadf['dataset_resolution'] = self.metadf['assumed_import_frequency']

        # Remove gaps and missing from the observations AFTER timecoarsening
        self.df = self.gaps.remove_gaps_from_obs(obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)

    # =============================================================================
    # Physiography extractions
    # =============================================================================

    def get_lcz(self):
        # connect to gee
        connect_to_gee()

        # Extract LCZ for all stations
        lcz_series = lcz_extractor(metadf = self.metadf,
                                   mapinfo=self.settings.gee['gee_dataset_info']['global_lcz_map'])

        #drop column if it was already present
        if 'lcz' in self.metadf:
            self.metadf = self.metadf.drop(columns=['lcz'])

        # update metadata
        self.metadf = self.metadf.merge(lcz_series.to_frame(),
                                        how='left',
                                        left_index=True, right_index=True)
        return lcz_series


    def get_altitude(self):
        # connect to gee
        connect_to_gee()

        # Extract LCZ for all stations
        altitude_series = height_extractor(metadf = self.metadf,
                                   mapinfo=self.settings.gee['gee_dataset_info']['DEM'])

        #drop column if it was already present
        if 'altitude' in self.metadf:
            self.metadf = self.metadf.drop(columns=['altitude'])

        # update metadata
        self.metadf = self.metadf.merge(altitude_series.to_frame(),
                                        how='left',
                                        left_index=True, right_index=True)
        return altitude_series


    def get_landcover(self, buffer=100, aggregate=True):


        # connect to gee
        connect_to_gee()

        # Extract landcover fractions for all stations
        lc_frac_df = lc_fractions_extractor(metadf = self.metadf,
                                   mapinfo=self.settings.gee['gee_dataset_info']['worldcover'],
                                   buffer=buffer,
                                   agg=aggregate)

        # remove columns if they exists
        self.metadf = self.metadf.drop(columns=list(lc_frac_df.columns),
                                       errors='ignore')

        # update metadf
        self.metadf = self.metadf.merge(lc_frac_df,
                                        how='left',
                                        left_index=True, right_index=True)

        return lc_frac_df



# =============================================================================
# Class stations (inherit all methods from dataset)
# =============================================================================

class Station(Dataset):
    def __init__(self, name, df, outliersdf, gaps, missing_obs, gapfilldf,
                 metadf, data_template, settings):
        self.name = name
        self.df = df
        self.outliersdf = outliersdf
        self.gaps = gaps
        self.missing_obs = missing_obs
        self.gapfilldf = gapfilldf
        self.metadf = metadf
        self.data_template = data_template
        self.settings=settings


        self._istype = 'Station'




def loggin_nan_warnings(df):
    """
    Function to feed the logger if Nan values are found in the df

    """
    for column in df.columns:
        bool_series = df[column].isnull()
        if bool_series.values.any():
            if bool_series.values.all():
                logger.warning(f'No values for {column}, they are initiated\
                               as Nan.')
            else:
                outliers = bool_series[bool_series]
                logger.warning(f'No values for stations: \
                               {outliers.index.to_list()}, for {column},\
                                   they are initiated as Nan.')

