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
                                        init_triple_multiindexdf,
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
        self.outliersdf = init_triple_multiindexdf()

        self.missing_obs = None  # becomes a Missingob_collection after import
        self.gaps = None  # becomes a gap_collection after import

        self.gapfilldf = init_multiindexdf()

        # Dataset with metadata (static)
        self.metadf = pd.DataFrame()
        # dataframe containing all information on the description and mapping
        self.data_template = pd.DataFrame()

        self._istype = 'Dataset'
        self._freqs = pd.Series(dtype=object)

        self._qc_checked_obstypes = [] #list with qc-checked obstypes

        self.settings = Settings()

    def update_settings(self, output_folder=None, input_data_file=None,
                  input_metadata_file=None, data_template_file=None,
                  metadata_template_file=None):
        """
        Update the most common input-output (IO) settings.
        (This should be applied before importing the observations.)

        When an update value is None, the specific setting will not be updated.

        Parameters
        ----------
        output_folder : string, optional
            A directory to store the output to. The default is None.
        input_data_file : string, optional
            Path to the input data file with observations. The default is None.
        input_metadata_file : string, optional
            Path to the input metadata file. The default is None.
        data_template_file : string, optional
            Path to the mapper-template csv file to be used on the observations.. The default is None.
        metadata_template_file : string, optional
            Path to the mapper-template csv file to be used on the metadata.. The default is None.

        Returns
        -------
        None.

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

        Parameters
        ----------
        timezonestr : string
            Timezone string of the input observations. Element of pytz.all_timezones.

        Returns
        -------
        None.

        """

        self.settings.update_timezone(timezonestr)

    def update_default_name(self, default_name):
        """
        Update the default name (the name of the station). This name will be
        used when no names are found in the observational dataset.

        (All observations are assumed to come from one station.)

        Parameters
        ----------
        default_name : string
            Default name to use when no names are present in the data.

        Returns
        -------
        None.

        """

        self.settings.app['default_name'] = str(default_name)

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
            sta_gapfill = self.gapfilldf.xs(
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
                       settings=self.settings)

    def show(self):
        """
        A function to print out some overview information about the Dataset.

        Returns
        -------
        None.

        """

        logger.info('Show basic info of dataset.')

        try:
            gapsdf = self.gaps.to_df()
        except:
            gapsdf = init_multiindexdf()

        print_dataset_info(self.df, self.outliersdf, gapsdf,
                           self.settings.app['print_fmt_datetime'])

    def make_plot(self, stationnames=None, obstype='temp', colorby='name',
                  starttime=None, endtime=None,
                  title=None, legend=True, show_outliers=True):
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


         Returns
         -------
         axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

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

        """

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
        None.

        """



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

        A final qualty control label for each
        quality-controlled-observation type can be added in the outputfile.

        The file will be writen to the outputfolder specified in the settings.

        Parameters
        ----------
        filename : string, optional
            The name of the output csv file. If none, a standard-filename
            is generated based on the period of data. The default is None.
        include_outliers : bool, optional
            If True, the outliers will be present in the csv file. The default is True.
        include_gapfill : bool, optional
            If True, the filled gap values will be present in the csv file. The default is True.
        add_final_labels : bool, optional
            If True, a column is added containing the final label of an observation. The default is True.
        use_tlk_obsnames : bool, optional
            If True, the standard naming of the metobs_toolkit is used, else
            the original names for obstypes is used. The default is True.

        Returns
        -------
        None.

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

        self._qc_checked_obstypes.append(obstype)
        self.outliersdf = self.outliersdf.sort_index()



    def combine_all_to_obsspace(self, repr_outl_as_nan=False):
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

        # =============================================================================
        # Unstack outliers to regular multiindex
        # =============================================================================
        outliersdf = self.outliersdf


        #remove duplicate indixes (needed for update)
        outliersdf = outliersdf[~outliersdf.index.duplicated(keep='first')]


        if not outliersdf.empty:
            outliersdf_values = outliersdf['value'].unstack() # for later use
            # convert to wide df with labels
            outliersdf = outliersdf['label'].unstack()

            # convert to final label names for columns
            outliersdf = outliersdf.rename(columns={col: col+'_final_label' for col in outliersdf.columns})


        else:
            outliersdf = init_multiindexdf()
            outliersdf_values = init_multiindexdf() # for later use
            outliercolumns = [col+'_final_label' for col in self.df if col in self.settings.app['observation_types']]
            for column in outliercolumns:
                outliersdf[column] = 'not checked'



        # =============================================================================
        # Combine observations and outliers
        # =============================================================================
        # get observations
        df = self.df


        # 1. Merge the label columns
        df_and_outl = df.merge(outliersdf, how='outer', left_index=True, right_index=True)

        # 2. fill the missing labels

        # split between obstype that are checked by qc and obstypes that are not checked

        checked_cols = [col+'_final_label' for col in self._qc_checked_obstypes]
        not_checked_cols = [col for col in df_and_outl.columns if ((col.endswith('_final_label')) and (not col in checked_cols))]

        # if obstype checked, and value is nan --> label ok
        df_and_outl[checked_cols] = df_and_outl[checked_cols].fillna('ok')
        # if obstype is not checked and label is missing --> label 'not checked'
        df_and_outl[not_checked_cols] = df_and_outl[not_checked_cols].fillna('not checked')


        # 3. Update the values if needed
        if not repr_outl_as_nan:

        # Merge obs and outliers, where obs values will be updated by outliers
            df_and_outl.update(other=outliersdf_values,
                         join='left',
                         overwrite=True,
                         errors='ignore')


        # =============================================================================
        # Make gaps, gapsfill and missing dataframes
        # =============================================================================

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




        # initiate default values
        for col in df_and_outl.columns:
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
        gapsdf = gapsdf[list(df_and_outl.columns)]
        missingdf = missingdf[list(df_and_outl.columns)]


        # Merge all together
        comb_df = pd.concat([df_and_outl, gapsdf, missingdf, gapsfilldf]).sort_index()

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
        """ V5 """

        self.outliersdf = pd.concat([self.outliersdf, add_to_outliersdf])



    # =============================================================================
    #     importing data
    # =============================================================================

    def coarsen_time_resolution(self, freq='1H', method='nearest', limit=1):
        """
        Resample the observations to coarser timeresolution. The assumed
        dataset resolution (stored in the metadf attribute) will be updated.


        Parameters
        ----------
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
        if isinstance(freq, type(None)):
            freq = self.settings.time_settings['target_time_res']
        if isinstance(method, type(None)):
            method = self.settings.time_settings['resample_method']
        if isinstance(freq, type(None)):
            limit = int(self.settings.time_settings['resample_limit'])


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

        # Remove gaps and missing from the observatios
        # most gaps and missing are already removed but when increasing timeres,
        # some records should be removed as well.
        self.df = self.gaps.remove_gaps_from_obs(obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)



    def import_data_from_file(self, coarsen_timeres=False):
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
        ----------
        coarsen_timeres : Bool, optional
            If True, the observations will be interpolated to a coarser
            time resolution as is defined in the Settings. The default
            is False.

        Returns
        ----------

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



        # Convert dataframe to dataset attributes
        self._initiate_df_attribute(dataframe=df)

        # Apply quality control on Import resolution
        self._apply_qc_on_import()

        # Remove gaps and missing from the observations AFTER timecoarsening
        self.df = self.gaps.remove_gaps_from_obs(obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)



    def import_data_from_database(self,
                                  start_datetime=None,
                                  end_datetime=None,
                                  coarsen_timeres=False):
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
                           the metadata file). This will be removed from the dataset.')
            df = df[~df.index.get_level_values('name').isnull()]


        # Convert dataframe to dataset attributes
        self._initiate_df_attribute(dataframe=df)

        # Apply quality control on Import resolution
        self._apply_qc_on_import()

        # Remove gaps and missing from the observations AFTER timecoarsening
        self.df = self.gaps.remove_gaps_from_obs(obsdf=self.df)
        self.df = self.missing_obs.remove_missing_from_obs(obsdf=self.df)


    def _initiate_df_attribute(self, dataframe):

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
        self.metadf['dataset_resolution'] = self.metadf['assumed_import_frequency']
        self.df = df.sort_index()
        self.original_df = df.sort_index()

    def _apply_qc_on_import(self):
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
            self.update_outliersdf(add_to_outliersdf=dup_outl_df)

        self.df, nan_outl_df = invalid_input_check(self.df,
                                                    checks_info=self.settings.qc['qc_checks_info'])
        if not nan_outl_df.empty:
            self.update_outliersdf(nan_outl_df)




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

