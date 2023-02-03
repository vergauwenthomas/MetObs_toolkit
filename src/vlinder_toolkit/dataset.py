#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""
# from collections.abc import Iterable

import pandas as pd
import geopandas as gpd
import numpy as np
import os
from datetime import datetime
import logging

from .settings import Settings
from .data_import import import_data_from_csv, import_data_from_database, template_to_package_space, import_metadata_from_csv
# from .data_import import coarsen_time_resolution
from .landcover_functions import geotiff_point_extraction
from .geometry_functions import find_largest_extent
from .plotting_functions import spatial_plot, timeseries_plot, timeseries_comp_plot, qc_stats_pie

from .qc_checks import gross_value_check, persistance_check, repetitions_check, duplicate_timestamp_check
from .qc_checks import step_check, missing_timestamp_and_gap_check, get_freqency_series
from .qc_checks import init_outlier_multiindexdf, gaps_to_outlier_format, window_variation_check

from .statistics import get_qc_effectiveness_stats

from .station import Station

logger = logging.getLogger(__name__)

# =============================================================================
# field classifiers
# =============================================================================

#Static fields are fields (attributes and observations) that do not change in time
static_fields = ['network', 'name', 
                'lat', 'lon', #TODO make these dynamic, now used as static 
                'call_name', 'location',
                'lcz']

#Categorical fields are fields with values that are assumed to be categorical.
#Note: (there are static and dynamic fields that are categorical)
categorical_fields = ['wind_direction', 'lcz']


observation_types = ['temp', 'radiation_temp', 'humidity', 'precip',
                     'precip_sum', 'wind_speed', 'wind_gust', 'wind_direction',
                     'pressure', 'pressure_at_sea_level']

location_info = ['network', 'lat', 'lon', 'lcz', 'call_name', 'location' ]



# =============================================================================
# Dataset class
# =============================================================================

class Dataset:
    def __init__(self):
        logger.info('Initialise dataset')

        self.df = pd.DataFrame() #Dataset with 'good' observations
           
        self.outliersdf = init_outlier_multiindexdf() #Dataset with outlier observations
        
        
        self.gapsdf = pd.DataFrame(columns=['start_gap', 'end_gap']) #Dataframe containig the gaps per station
        #TODO: a static metadf excludes moving observations :s
        self.metadf = pd.DataFrame() #Dataset with metadata (static)
        self.data_template = pd.DataFrame() #dataframe containing all information on the description and mapping
        
        self._freqs = pd.Series()
        
    
        
    def get_station(self, stationname):
        
        """
        Extract a station object from the dataset.

        Parameters
        ----------
        stationname : String
            Name of the station, example 'vlinder16'

        Returns
        -------
        station_obj : vlinder_toolkit.station.Station
            

        """
        logger.info(f'Extract {stationname} from dataset.')
        
        try:
            sta_df = self.df.xs(stationname, level='name')
            sta_meta_series = self.metadf.loc[stationname]
        except KeyError:
            logger.warning(f'{stationname} not found in the dataset.')
            print(f'{stationname} not found in the dataset.')
            return None
        
        try:
            sta_outliers = self.outliersdf.xs(stationname, level='name')
        except KeyError:
            sta_outliers = init_outlier_multiindexdf()
        
        try:
            sta_gaps = self.gapsdf.loc[stationname]
        except KeyError:
            sta_gaps = pd.DataFrame()
        

            
        return Station(name=stationname, df=sta_df,
                       outliersdf=sta_outliers,
                       gapsdf = sta_gaps,
                       meta_series=sta_meta_series,
                       data_template=self.data_template)
       
        
        
    
    def show(self):
        """
        Print basic information about the dataset.

        Returns
        -------
        None.

        """
        logger.info('Show basic info of dataset.')
        if self.df.empty:
            print("This dataset is empty!")
            # logger.error('The dataset is empty!')
        else: 
            print('\n','--------  General ---------', '\n')
            print(f' .... ')
            
            
            
            print('\n','--------  Observations ---------', '\n')
            starttimestr = datetime.strftime(min(self.df.index.get_level_values(level='datetime')),
                                             Settings.print_fmt_datetime)
            endtimestr = datetime.strftime(max(self.df.index.get_level_values(level='datetime')), Settings.print_fmt_datetime)
            
            stations_available = list(self.df.index.get_level_values(level='name').unique())
            print(f'Observations found for period: {starttimestr} --> {endtimestr}')
            # logger.debug(f'Observations found for period: {starttimestr} --> {endtimestr}')
            print(f'Following stations are in dataset: {stations_available}')
            # logger.debug(f'Following stations are in dataset: {stations_available}')
            
            print('\n', '--------  Outliers ---------', '\n')
            print(f'There are {self.outliersdf.shape[0]} flagged observations found in total. They occure in these stations: {list(self.outliersdf.index.get_level_values("name").unique())}')
            
            print('\n', '--------  Gaps ---------', '\n')
            print(f'There are {self.gapsdf.shape[0]} gaps found in total. They occure in these stations: {list(self.gapsdf.index.unique())}')
            
            
            
        
    def make_plot(self, stationnames=None, variable='temp',
                                   starttime=None, endtime=None,
                                   title=None, legend=True):
        """
        This function create a timeseries plot for the dataset. The variable observation type
        is plotted for all stationnames from a starttime to an endtime.

        Parameters
        ----------
        stationnames : List, Iterable
            Iterable of stationnames to plot. If None, all available stations ar plotted. The default is None.
        variable : String, optional
            The name of the observation type to plot. The default is 'temp'.
        starttime : datetime, optional
            The starttime of the timeseries to plot. The default is None and all observations 
            are used.
        endtime : datetime, optional
            The endtime of the timeseries to plot. The default is None and all observations 
            are used..
        title : String, optional
            Title of the figure, if None a default title is generated. The default is None.
        legend : Bool, optional
           Add legend to the figure. The default is True.
        Returns
        -------
        ax : matplotlib.axes
            The plot axes is returned.

        """
        logger.info(f'Make {variable}-timeseries plot for {stationnames}')
        
        default_settings=Settings.plot_settings['time_series']
        
        
        #Subset on obseravtion type
        plotdf = self.df[variable]
        
        #Unstack dataframe on name
        plotdf = plotdf.unstack('name')
        
        #Subset on stationnames
        if not isinstance(stationnames, type(None)):
            plotdf = plotdf[stationnames]
       
        #Subset on start and endtime
        plotdf = datetime_subsetting(plotdf, starttime, endtime)
        
        
        #plotdf is a dataframe with this structure:
            #datatime --> stationnameA, stationnameB, stationnameC
          # 2022-09-15 ...    valueA        valueB       valueC
        
       
        #Get plot styling attributes
        if isinstance(title, type(None)):
            if isinstance(stationnames, type(None)):
                title=Settings.display_name_mapper[variable] + ' for all stations. '
            else:
                title=Settings.display_name_mapper[variable] + ' for stations: ' + str(stationnames)
        
        
        #make plot
        ax = timeseries_comp_plot(plotdf=plotdf,
                                  title=title,
                                  xlabel='',
                                  ylabel=Settings.display_name_mapper[variable],
                                  figsize=default_settings['figsize'])
        

        return ax
        
        
        
        
        
    def make_geo_plot(self, variable='temp', title=None, timeinstance=None, legend=True,
                      vmin=None, vmax=None):
        """
        This functions creates a geospatial plot for a field (observations or attributes) of all stations.
        
        If the field is timedepending, than the timeinstance is used to plot the field status at that datetime.
        If the field is categorical than the leged will have categorical values, else a colorbar is used. 
        
        All styling attributes are extracted from the Settings.
        

        Parameters
        ----------
        variable : String, optional
            Fieldname to visualise. This can be an observation or station attribute. The default is 'temp'.
        title : String, optional
            Title of the figure, if None a default title is generated. The default is None.
        timeinstance : datetime, optional
            Datetime moment of the geospatial plot. The default is None and the first datetime available
            is used.
        legend : Bool, optional
            Add legend to the figure. The default is True.
        vmin : float, optional
            The minimum value corresponding to the minimum color. The default is None and 
            the minimum of the variable is used.
        vmax : float, optional
           The maximum value corresponding to the minimum color. The default is None and 
           the maximum of the variable is used.

        Returns
        -------
        ax : Geoaxes
            The geoaxes is returned.

        """
        
        
        
        #Load default plot settings
        default_settings=Settings.plot_settings['spatial_geo']
        
        #get first timeinstance of the dataset if not given
        if isinstance(timeinstance, type(None)):
            timeinstance=self.df.index.get_level_values('datetime').min()
            
        logger.info(f'Make {variable}-geo plot at {timeinstance}')
        
        #subset to timeinstance
        plotdf = self.df.xs(timeinstance, level='datetime')
        
        #merge metadata
        plotdf = plotdf.merge(self.metadf, how='left',
                              left_index=True, right_index=True )
        
        #subset to obstype
        plotdf = plotdf[[variable, 'geometry']]
        
        #Subset to the stations that have coordinates
        ignored_stations = plotdf[plotdf['geometry'].isnull()]
        plotdf = plotdf[~plotdf['geometry'].isnull()]
        if plotdf.empty:
            logger.error(f'No coordinate data found, geoplot can not be made. Plotdf: {plotdf}')
            print(f'No coordinate data found, geoplot can not be made. Plotdf: {plotdf}')
            return
        
        if not ignored_stations.empty:
            logger.error(f'No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!')
            print(f'No coordinate found for following stations: {ignored_stations.index.to_list()}, these will be ignored in the geo-plot!')
        


        #make color scheme for field
        if variable in categorical_fields:
            is_categorical=True
            if variable == 'lcz':
                #use all available LCZ categories
                use_quantiles=False
            else:
                use_quantiles=True
        else:
            is_categorical=False
            use_quantiles=False
     
        
        #if observations extend is contained by default exten, use default else use obs extend
        use_extent=find_largest_extent(geodf=gpd.GeoDataFrame(plotdf),
                                       extentlist=default_settings['extent'])
        
        
        #Style attributes
        if isinstance(title, type(None)):
            if variable in static_fields:
                title = Settings.display_name_mapper[variable]
            else:
                dtstring = datetime.strftime(timeinstance, default_settings['fmt'])
                title = Settings.display_name_mapper[variable] + ' at ' + dtstring
        
        ax = spatial_plot(gdf=plotdf,
                          variable=variable,
                          legend=legend,
                          use_quantiles=use_quantiles,
                          is_categorical=is_categorical,
                          k_quantiles=default_settings['n_for_categorical'],
                          cmap = default_settings['cmap'],
                          world_boundaries_map=Settings.world_boundary_map,
                          figsize=default_settings['figsize'],
                          extent=use_extent,
                          title=title,
                          vmin=vmin,
                          vmax=vmax
                          )
        

        return ax
    
    
    def write_to_csv(self, filename=None, include_outliers=True, add_final_labels=True):
        """
            Write the dataset to a file where the observations, metadata and (if available)
            the quality labels per observation type are merged together. 
            
            A final qualty controll label for each quality-controlled-observation type
            can be added in the outputfile.
            
            The file will be writen to the Settings.outputfolder.
    
            Parameters
            ----------
            filename : string, optional
                The name of the output csv file. If none, a standard-filename is generated
                based on the period of data. The default is None.
            add_final_labels : Bool, optional
                If True, a final qualty control label per observation type
                is added as a column. The default is True.
    
            Returns
            -------
            None
    
            """
        
        logger.info('Writing the dataset to a csv file')
        assert not isinstance(Settings.output_folder, type(None)), 'Specify Settings.output_folder in order to export a csv.'
        
        #add final quality control labels per observation type
        if add_final_labels:
            outliersdf = add_final_label_to_outliersdf(outliersdf=self.outliersdf,
                                              gapsdf= self.gapsdf,
                                              data_res_series = self.metadf['dataset_resolution'])
        else:
            outliersdf = self.outliersdf
         
       
        #Get observations and metadata columns in the right order
        logger.debug('Merging data and metadata')
        
        
        #make column ordering
        df_columns = observation_types.copy() #observations
        df_columns.extend(location_info) #metadata
        qc_columns = [col for col in outliersdf if col.endswith('_label')] #add qc labels
        df_columns.extend(qc_columns)
        df_columns.insert(0, 'datetime') # timestamp as first column
        df_columns.insert(1, 'name') #station name as second column
        
        
    
        
        #unstack observations and merge with metadf
        df = self.df 
        df[qc_columns] = 'ok'
        df = pd.concat([df, outliersdf])
        df = df.reset_index()
        
        metadf = self.metadf.reset_index()
        df = df.merge(metadf, how='left', on='name')
        
        #sort and subset columns
        df = df[df_columns]
        
                
        #find observation type that are not present
        ignore_obstypes = [col for col in observation_types if df[col].isnull().all()]
        
        df = df.drop(columns=ignore_obstypes)
        
        logger.debug(f'Skip quality labels for obstypes: {ignore_obstypes}.')
        
        df = df.sort_values(['name', 'datetime'])
       
        #make filename
        if isinstance(filename, type(None)):
            startstr = self.df.index.min().strftime('%Y%m%d') 
            endstr = self.df.index.max().strftime('%Y%m%d') 
            filename= 'dataset_' + startstr + '_' + endstr
        else:
            if filename.endswith('.csv'):
                filename = filename[:-4] #to avoid two times .csv.csv
            
        filepath = os.path.join(Settings.output_folder, filename + '.csv')
        
        #write to csv in output folder
        logger.info(f'write dataset to file: {filepath}')
        df.to_csv(path_or_buf=filepath,
                       sep=';',
                       na_rep='NaN',
                       index=True)        
        
    

    
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
        V3
        Apply quality control methods to the dataset. The default settings are used, and can be changed
        in the settings_files/qc_settings.py
        
        The checks are performed in a sequence: gross_vallue --> persistance --> ...,
        Outliers by a previous check are ignored in the following checks!
        
        The dataset is updated inline.

        Parameters
        ----------
        obstype : String, optional
            Name of the observationtype you want to apply the checks on. The default is 'temp'.
        gross_value : Bool, optional
            If True the gross_value check is applied if False not. The default is True.
        persistance : Bool, optional
           If True the persistance check is applied if False not. The default is True.. The default is True.
        step : Bool, optional
           If True the step check is applied if False not. The default is True.
       internal_consistency : Bool, optional	        
           If True the internal consistency check is applied if False not. The default is True.	          
       qc_info: Bool, optional
           If True info about the quality control is printed if False not. The default is True.
        ignore_val : numeric, optional
            Values to ignore in the quality checks. The default is np.nan.

        Returns
        -------
        None.

        """
        

        if repetitions:
           print('Applying the repetitions-check on all stations.')
           logger.info('Applying repetitions check on the full dataset')
           
           checked_series, outl_df = repetitions_check(input_series=self.df[obstype],
                                                     obstype=obstype)
           
           #update the dataset and outliers
           self.df[obstype] = checked_series
           self.update_outliersdf(outl_df)
         
            
        if gross_value:
            print('Applying the gross-value-check on all stations.')
            logger.info('Applying gross value check on the full dataset')
    
            checked_series, outl_df = gross_value_check(input_series = self.df[obstype],
                                                        obstype=obstype)
            
            #update the dataset and outliers
            self.df[obstype] = checked_series
            self.update_outliersdf(outl_df)
            
            
        if persistance:
            print('Applying the persistance-check on all stations.')
            logger.info('Applying persistance check on the full dataset')
          
            checked_series, outl_df = persistance_check(station_frequencies=self.metadf['dataset_resolution'], input_series=self.df[obstype],
                                                        obstype=obstype)

            #update the dataset and outliers
            self.df[obstype] = checked_series
            self.update_outliersdf(outl_df)
            
        if step:
            print('Applying the step-check on all stations.')
            logger.info('Applying step-check on the full dataset')
           
            checked_series, outl_df = step_check(input_series=self.df[obstype],
                                                 obstype=obstype)
                                                      
            
            #update the dataset and outliers
            self.df[obstype] = checked_series
            self.update_outliersdf(outl_df)
            
        if window_variation:
            print('Applying the window variation-check on all stations.')
            logger.info('Applying window variation-check on the full dataset')
           
            checked_series, outl_df = window_variation_check(station_frequencies=self.metadf['dataset_resolution'], input_series=self.df[obstype],
                                                 obstype=obstype)
                                                      
            
            #update the dataset and outliers
            self.df[obstype] = checked_series
            self.update_outliersdf(outl_df)
        
        self.outliersdf = self.outliersdf.sort_index()
        


        # if internal_consistency:
        #     print('Applying the internal-concsistency-check on all stations.')
        #     logger.info('Applying step-check on the full dataset')
           
        #     checked_obs, qc_flags = internal_consistency_check(input_series=self.df[obstype],
        #                                                        humidity_series=self.df['humidity'],
        #                                                        obstype=obstype,
        #                                                        ignore_val=ignore_val)

        #     #update the dataset
        #     self.df[obstype] = checked_obs
        #     label_column_name = qc_flags.name
        #     self.df[label_column_name] = qc_flags
        
    
    def get_qc_stats(self, obstype='temp', stationnames=None, make_plot=True, coarsen_timeres=False):
        """
        Compute frequency statistics on the qc labels for an observationtype.
        The output is a dataframe containing the frequency statistics presented
        as percentages. 
        
        These frequencies can also be presented as a collection of piecharts per check.
        
        With stationnames you can subset the data to one ore multiple stations.

        Parameters
        ----------
        obstype : Str, optional
            Observation type to analyse the QC labels on. The default is 'temp'.
        stationnames : List, Str, optional
            Stationname(s) to subset the quality labels on. If None, all stations are used. The default is None.
        make_plot : Bool, optional
            If True, a plot with piecharts is generated. The default is True.

        Returns
        -------
        dataset_qc_stats : pandas.DataFrame
            A table containing the label frequencies per check presented as percentages0.

        """
        
        outliersdf = self.outliersdf
        #Add gaps to the outliers for comuting scores
        outliersdf = pd.concat([outliersdf,
                                gaps_to_outlier_format(gapsdf=self.gapsdf,
                                                       dataset_res_series=self.metadf['dataset_resolution'])])
        
        if coarsen_timeres:
            gaps_1 = outliersdf[outliersdf["gap_timestamp_label"] == "missing timestamp (gap)"]
            outliersdf.drop(gaps_1.index, inplace=True)
            outliersdf = outliersdf.loc[outliersdf.index.intersection(self.df.index)]
            outliersdf = pd.concat([outliersdf, gaps_1], sort=True)
            
        outliersdf = outliersdf.fillna(value='ok')
        
        if isinstance(stationnames, type(None)):
            df = self.df
            qc_labels_df = self.get_final_qc_labels()
  
            title=f'Frequency for {obstype}-qc-checks on all stations.'
        else: 
            title=f'Frequency for {obstype}-qc-checks on {stationnames}.'
            if isinstance(stationnames, str):
                stationnames = [stationnames]
            
            qc_labels_df = self.get_final_qc_labels().loc[self.get_final_qc_labels().index.get_level_values(level='name').isin(stationnames)]
            df = self.df.loc[self.df.index.get_level_values(level='name').isin(stationnames)]
            outliersdf = outliersdf.loc[outliersdf.index.get_level_values(level='name').isin(stationnames)]
        
        
        #Do not include the 'final_qc_label' in the statis
        if obstype + '_final_label' in outliersdf.columns:
            logger.debug(f'The {obstype}_final_label, is ingored for the QC-stats. ')
            outliersdf = outliersdf.loc[:, outliersdf.columns != obstype+'_final_label']
            
       
        
        #stats on datset level
        qc_labels = {key: val['outlier_flag'] for key, val in Settings.qc_checks_info.items()}
     
        if coarsen_timeres:
            for row in self.gapsdf.iterrows():
                df = df.iloc[df.index.get_level_values('name') == row[0]]
                df.reset_index(level='name', inplace=True)
                indices_to_remove = df.loc[(df.index <= row[1]['end_gap']) & (df.index >= row[1]['start_gap'])].index
                df.drop(indices_to_remove, inplace=True)
                df.set_index(['name', df.index], inplace=True)
        
        dataset_qc_stats = get_qc_effectiveness_stats(outliersdf = outliersdf,
                                                      df =df,
                                                      obstype=obstype,
                                                      observation_types = observation_types,
                                                      qc_labels=qc_labels)
        
       
        if coarsen_timeres:
            gaps_2 = qc_labels_df[qc_labels_df["gap_timestamp_label"] == "missing timestamp (gap)"]
            qc_labels_df.drop(gaps_2.index, inplace=True)
            qc_labels_df = qc_labels_df.loc[qc_labels_df.index.intersection(self.df.index)]
            qc_labels_df = pd.concat([qc_labels_df, gaps_2], sort=True)
        
        valid_records_df = df.drop(qc_labels_df.index.intersection(df.index).dropna())           
        
        
        if make_plot:
            qc_stats_pie(valid_records_df, qc_labels_df, qc_stats=dataset_qc_stats,
                         figsize=Settings.plot_settings['qc_stats']['figsize'],
                         title=title)
               
        
        return dataset_qc_stats
    
    
    def update_outliersdf(self, add_to_outliersdf):
        
        #Get the flag column labels and find the newly added columnlabelname
        previous_performed_checks_columns =  [col for col in self.outliersdf.columns if col.endswith('_label')]
        new_performed_checks_columns = list(set([col for col in add_to_outliersdf.columns if col.endswith('_label')]) - set(previous_performed_checks_columns))

        #add to the outliersdf
        self.outliersdf = pd.concat([self.outliersdf, add_to_outliersdf])
        
        #Fix labels
        self.outliersdf[previous_performed_checks_columns] = self.outliersdf[previous_performed_checks_columns].fillna(value='ok')
        self.outliersdf[new_performed_checks_columns] = self.outliersdf[new_performed_checks_columns].fillna(value='not checked')
    
    def get_final_qc_labels(self):
        """
       
        """
        #add final quality labels
    
        outldf = add_final_label_to_outliersdf(outliersdf=self.outliersdf,
                                               gapsdf=self.gapsdf,
                                               data_res_series=self.metadf['dataset_resolution'])
        
        return outldf
        
        
    # =============================================================================
    #     importing data        
    # =============================================================================
      
    def coarsen_time_resolution(self, freq='1H', method='nearest', limit=1):
        logger.info(f'Coarsening the timeresolution to {freq} using the {method}-method (with limit={limit}).')
        #TODO: implement buffer method
        
        #Coarsen timeresolution
        df = self.df.reset_index()
        if method == 'nearest':
            df = df.set_index('datetime').groupby('name').resample(freq).nearest(limit=limit)
            
        elif method=='bfill':
            df = df.set_index('datetime').groupby('name').resample(freq).bfill(limit=limit)
           
        else: 
            print(f'The coarsening method: {method}, is not implemented yet.')
            df = df.set_index(['name', 'datetime'])
            
        if 'name' in df.columns:
            df = df.drop(columns=['name'])
            
        #Update resolution info in metadf
        self.metadf['dataset_resolution'] = pd.to_timedelta(freq)
        #update df
        self.df = df
        
    
    def import_data_from_file(self, network='vlinder', coarsen_timeres=False):
        """
        Read observations from a csv file as defined in the Settings.input_file. 
        The input file columns should have a template that is stored in Settings.template_list.
        
        If the metadata is stored in a seperate file, and the Settings.input_metadata_file is correct,
        than this metadata is also imported (if a suitable template is in the Settings.template_list.)
        
        
        It is possible to apply a 
        resampling (downsampling) of the observations as defined in the settings.
        
        After the import there is always a call to Dataset.update_dataset_by_df, that 
        sets up the dataset with the observations and applies some sanity checks.
        
        Parameters
        ----------
        network : String, optional
            The name of the network for these observationsThe default is 'vlinder'.
        coarsen_timeres : Bool, optional
            If True, the observations will be interpolated to a coarser time resolution
            as is defined in the Settings. The default is False.

        Returns
        -------
        None.

        """
        print('Settings input data file: ', Settings.input_data_file)
        logger.info(f'Importing data from file: {Settings.input_data_file}')
        
        # Read observations into pandas dataframe
        df, template = import_data_from_csv(input_file = Settings.input_data_file,
                                  file_csv_template=Settings.input_csv_template,
                                  template_list = Settings.template_list)
        
        logger.debug(f'Data from {Settings.input_data_file} imported to dataframe.')

        #drop Nat datetimes if present
        df = df.loc[pd.notnull(df.index)]
        
        
        if isinstance(Settings.input_metadata_file, type(None)):
            print('WARNING: No metadata file is defined. Add your settings object.')
            logger.warning('No metadata file is defined, no metadata attributes can be set!')
        else:
            logger.info(f'Importing metadata from file: {Settings.input_metadata_file}')
            meta_df = import_metadata_from_csv(input_file=Settings.input_metadata_file,
                                               file_csv_template=Settings.input_metadata_template,
                                               template_list = Settings.template_list)
            
            #merge additional metadata to observations
            meta_cols = [colname for colname in meta_df.columns if not colname.startswith('_')]
            additional_meta_cols = list(set(meta_cols).difference(df.columns))
            if bool(additional_meta_cols):
                logger.debug(f'Merging metadata ({additional_meta_cols}) to dataset data by name.')
                additional_meta_cols.append('name') #merging on name
                df_index = df.index #merge deletes datetime index somehow? so add it back on the merged df
                df = df.merge(right=meta_df[additional_meta_cols],
                              how='left', 
                              on='name')
                df.index = df_index
        
        
        #update dataset object
        self.data_template = pd.DataFrame().from_dict(template)
        

        #convert dataframe to multiindex (datetime - name)
        df = df.set_index(['name', df.index])
        
        #dataframe with all data of input file
        self.input_df = df
        
        self.update_dataset_by_df(dataframe = df, 
                                  coarsen_timeres=coarsen_timeres)

        
    
    def import_data_from_database(self,
                              start_datetime=None,
                              end_datetime=None,
                              coarsen_timeres=False):
        """
        Function to import data directly from the framboos database and updating the 
        network and station objects. 
        

        Parameters
        ----------
        start_datetime : datetime, optional
            Start datetime of the observations. The default is None and using 
            yesterday's midnight.
        end_datetime : datetime, optional
            End datetime of the observations. The default is None and using todays
            midnight.
        coarsen_timeres : Bool, optional
            If True, the observations will be interpolated to a coarser time resolution
            as is defined in the Settings. The default is False.

        Returns
        -------
        None.

        """
        if isinstance(start_datetime, type(None)):
            start_datetime=datetime.date.today() - datetime.timedelta(days=1)
        if isinstance(end_datetime, type(None)):
            end_datetime=datetime.date.today()
       
            
        # Read observations into pandas dataframe
        df = import_data_from_database(Settings,
                                       start_datetime=start_datetime,
                                       end_datetime=end_datetime)
        
        
        #Make data template
        self.data_template =  pd.DataFrame().from_dict(template_to_package_space(Settings.vlinder_db_obs_template))
        
        # if coarsen_timeres:
        #     df = coarsen_time_resolution(df=df,
        #                                   freq=Settings.target_time_res,
        #                                   method=Settings.resample_method)
            
        #convert dataframe to multiindex (datetime - name)
        df = df.set_index(['name', df.index])
        df = df.sort_index()
        
        #If an ID has changed or not present in the metadatafile, the stationname and metadata is Nan
        #These observations will be removed
        unknown_obs = df[df.index.get_level_values('name').isnull()]
        if not unknown_obs.empty:
            logger.warning('There is an unknown station in the dataset (probaply due to an ID that is not present in the metadata file). This will be removed.')
            df = df[~df.index.get_level_values('name').isnull()]
        
        self.update_dataset_by_df(dataframe=df, coarsen_timeres=coarsen_timeres)
        
        
   


    
        
    
    def update_dataset_by_df(self, dataframe, coarsen_timeres=False):
        """
        Update the dataset object by a dataframe.
        
        
        When filling the observations, there is an automatic check for missing timestamps and duplicating timestamps. 
        If a missing timestamp is detected, the timestamp is created with Nan values for all observation types.
        
        If metadata is present and a LCZ-tiff file availible, than the LCZ's of the stations are computed.

        Parameters
        ----------
        dataframe : pandas.DataFrame
        A dataframe that has an datetimeindex and following columns: 'name, temp, radiation_temp, humidity, ...'
            

        Returns
        -------
        None.

        """
       
        logger.info(f'Updating dataset by dataframe with shape: {dataframe.shape}.')
        
        #Create dataframe with fixed number and order of observational columns
        df = dataframe.reindex(columns = observation_types)
        self.df = df
        
        #create metadataframe with fixed number and order of columns
        metadf = dataframe.reindex(columns = location_info)
        metadf.index = metadf.index.droplevel('datetime') #drop datetimeindex
        metadf = metadf[~metadf.index.duplicated(keep='first')]#drop dubplicates due to datetime
        
        self.metadf = metadf_to_gdf(metadf)
        

        #add import frequencies to metadf
        self.metadf['assumed_import_frequency'] = get_freqency_series(self.df)
        
        #TODO: How to implement the choise to apply QC on import freq or on coarsened frequency
        self.df = df.sort_index()
        
        self.df, missing_outl_df, self.gapsdf, station_freqs= missing_timestamp_and_gap_check(df=self.df)
        self.df, dup_outl_df = duplicate_timestamp_check(df=self.df)
        #print(self.df.iloc[:100,].to_string())
        
        #update outliersdf
        self.outliersdf = pd.concat([self.outliersdf, dup_outl_df, missing_outl_df])
       
        if coarsen_timeres:
            self.coarsen_time_resolution(freq=Settings.target_time_res,
                                          method=Settings.resample_method,
                                          limit=Settings.resample_limit)
            
        else:
            self.metadf['dataset_resolution'] = self.metadf['assumed_import_frequency']
            
        #get LCZ values (if coords are availible)
        self.metadf =  get_lcz(self.metadf)
        

def metadf_to_gdf(df, crs=4326):
    """
    Function to convert a dataframe with 'lat' en 'lon' columnst to a geopandas 
    dataframe with a geometry column containing points.
    
    Special care for stations with missing coordinates.

    Parameters
    ----------
    df : pandas.DataFrame
        Dataframe with a 'lat' en 'lon' column.
    crs : Integer, optional
        The epsg number of the coordinates. The default is 4326.

    Returns
    -------
    geodf : geopandas.GeaDataFrame
        The geodataframe equivalent of the df.

    """
    
    # only conver to points if coordinates are present
    coordsdf = df[(~df['lat'].isnull()) & (~df['lon'].isnull())]
    missing_coords_df =  df[(df['lat'].isnull()) | (df['lon'].isnull())]
    
    geodf = gpd.GeoDataFrame(coordsdf,
                              geometry=gpd.points_from_xy(coordsdf.lon,
                                                          coordsdf.lat)) 
    geodf = geodf.set_crs(epsg = crs)
    geodf = pd.concat([geodf, missing_coords_df])
    
    geodf = geodf.sort_index()
    return geodf


def get_lcz(metadf):
    """
    Function to extract the LCZ's for all locations in the metadf.
    A column 'lcz' is added tot the metadf.
    
    All information on the LCZ is extracted from the Setting object (class mapper, location)
    ----------
    metadf : geopandas.GeoDataFrame
        Geodataframe with the coordinates present as geometry.

    Returns
    -------
    metadf : geopandas.GeoDataFrame
        The metadf with the added 'lcz'-column.

    """
    logger.debug('Extract LCZs')

    
    
    if metadf['geometry'].x.isnull().values.all():
        logger.info('Extract LCZs is not possible because no longtitude information is found.')
        metadf['lcz'] = 'Location unknown'
        return metadf
    if metadf['geometry'].y.isnull().values.all():
        logger.info('Extract LCZs is not possible because no latitude information is found.')
        metadf['lcz'] = 'Location unknown'
        return metadf
    

    
    geo_templates = Settings.geo_datasets_templates
    lcz_file = Settings.geo_lcz_file
    
    if isinstance(lcz_file, type(None)):
        print('No lcz tif location in the settings. Update settings: ')
        print('settings_obj.update_settings(geotiff_lcz_file="...."')
        logger.error('Extracting LCZ but no geotiff file specified!')
        metadf['lcz'] = 'Unkown lcz file'
        return metadf
    
    lcz_templates = [geo_templ for geo_templ in geo_templates if geo_templ['usage']=='LCZ']
    
    assert len(lcz_templates)==1, 'More (or no) lcz template found!'
    
    lcz_template = lcz_templates[0]
    
    human_mapper = {num: lcz_template['covers'][num]['cover_name'] 
                    for num in lcz_template['covers'].keys()}
    

   

    lcz = geotiff_point_extraction(geodf=metadf,
                                   geotiff_location=lcz_file,
                                   geotiff_crs=lcz_template['epsg'],
                                   class_to_human_mapper=human_mapper)
    metadf['lcz'] = lcz
        
    
    return metadf
          
def loggin_nan_warnings(df):
    """
    Function to feed the logger if Nan values are found in the df
    
    """
    for column in df.columns:
        bool_series = df[column].isnull()
        if bool_series.values.any():
            if bool_series.values.all():
                logger.warning(f'No values for {column}, they are initiated as Nan.')
            else:
                
                outliers = bool_series[bool_series==True]
                logger.warning(f'No values for stations: {outliers.index.to_list()}, for {column}, they are initiated as Nan.')
        


def datetime_subsetting(df, starttime, endtime):
    """
    Wrapper function for subsetting a dataframe with datetimeindex with a start- and 
    endtime. 

    Parameters
    ----------
    df : pandas.DataFrame with datetimeindex
        The dataframe to apply the subsetting to.
    starttime : datetime.Datetime
        Starttime for the subsetting period (included).
    endtime : datetime.Datetime
        Endtime for the subsetting period (included).

    Returns
    -------
    pandas.DataFrame
        Subset of the df.

    """
    
    stand_format = '%Y-%m-%d %H:%M:%S'
    
    if isinstance(starttime, type(None)):
        startstring = None #will select from the beginning of the df
    else:
        startstring = starttime.strftime(stand_format)
    if isinstance(endtime, type(None)):
        endstring = None
    else: 
        endstring = endtime.strftime(stand_format)

    return df[startstring: endstring]




def add_final_label_to_outliersdf(outliersdf, gapsdf, data_res_series):
    """
    V3
        This function creates a final label based on de individual qc labels. The final label will be that of the individual qc-label
        which rejected the obseration.
       
        This functions converts labels to numeric values, algebra to get final label, and inversly 
        convert to labels. This is faster than looping over the rows.

        Parameters
        ----------
        outliersdf : pandas.DataFrame 
            The dataset outliers dataframe containing the observations and QC labels. 
        gapsdf : pandas.Dataframe
            The dataset gaps dataframe. The gaps will be exploded and added to the outliersdf.
        data_res_series : Pandas.Series
            The series that contain the dataset resolution (values) per station (index). This 
            is stored in the dataset.metadf as column 'dataset_resolution'. These are used to explode the gaps.

        Returns
        -------
        outliersdf : pd.DataFrame
            The outliersdf with extra columns indicated by example 'temp_final_label' and 'humid_final_label'.

        """ 
       

    #combine gaps with outliers
    gaps_exploded = gaps_to_outlier_format(gapsdf, data_res_series)
    outliersdf = pd.concat([outliersdf, gaps_exploded])
    #fill Nan's by 'ok'
    if Settings.qc_checks_info['gaps_finder']['label_columnname'] in outliersdf.columns:
        outliersdf[Settings.qc_checks_info['gaps_finder']['label_columnname']].fillna(value='ok',
                                                                                      inplace=True)
    
    # order columns
    labels_columns = [column for column in outliersdf.columns if not column in observation_types]
    checked_obstypes = [obstype for obstype in observation_types if any([qc_column.startswith(obstype+'_') for qc_column in labels_columns])]
    columns_on_record_lvl = [info['label_columnname'] for checkname, info in Settings.qc_checks_info.items() if info['apply_on'] == 'record']
   
    
    # Construct numeric mapper
    labels_to_numeric_mapper = {info['outlier_flag']:info['numeric_flag'] for info in Settings.qc_checks_info.values()}
    # add 'ok' and 'not checked' labels
    labels_to_numeric_mapper['ok'] = 0
    labels_to_numeric_mapper['not checked'] = np.nan
    #invert numeric mapper
    inv_label_to_num = {v: k for k, v in labels_to_numeric_mapper.items()}
    
    
    #generete final label per obstype
    for obstype in checked_obstypes:
        # logger.debug(f'Generating final QC labels for {obstype}.')
        #Get qc column namse specific for this obstype
        specific_columns = [col for col in labels_columns if col.startswith(obstype+'_')]
        #add qc labels that are applicable on all obstypes
        specific_columns.extend(columns_on_record_lvl)
        
        #Drop columns that are not present
        specific_columns = [colmn for colmn in specific_columns if colmn in outliersdf.columns]
        
        
        
        #get labels dataframe
        qc_df = outliersdf[specific_columns]
        num_qc_df = pd.DataFrame()
        num_qc_df = qc_df.applymap(labels_to_numeric_mapper.get )
       
    
        outliersdf[obstype+'_final_label'] = num_qc_df.sum(axis=1, skipna=True).map(inv_label_to_num)
       
   
    return outliersdf
# def final_qc_label_maker(df, label_to_numeric_mapper):
#     """
#     This function creates a final label based on de individual qc labels. If all labels
#     are ok, the final label is ok. Else the final label will be that of the individual qc-label
#     which rejected the obseration.
    
#     This functions converts labels to numeric values, algebra to get final label, and inversly 
#     convert to labels. This is faster than looping over the rows.

#     Parameters
#     ----------
#     qc_df : pandas.DataFrame 
#         the specific qc_label_df with the datetimeindex, the first column the observations,
#         and labels for each QC check per column.
#     label_to_numeric_mapper : dict
#         The dictionary that maps qc-labels to numeric values (for speedup).

#     Returns
#     -------
#     final_labels : pd.Series
#         A series with the final labels and the same index as the qc_df.

#     """ 
        
#     #invert numeric mapper
#     inv_label_to_num = {v: k for k, v in label_to_numeric_mapper.items()}
    
    
#     #extract label columns
#     qc_labels_columns = [col for col in df.columns if not col in observation_types]
#     #Extra savety
#     qc_labels_columns = [col for col in qc_labels_columns if col.endswith('_label')]
    
    
#     # find the observation types on which QC is applied
#     checked_obstypes = [obstype for obstype in observation_types if any([qc_column.startswith(obstype+'_') for qc_column in qc_labels_columns])]
    
    
#     #generete final label per obstype
#     for obstype in checked_obstypes:
#         logger.debug(f'Generating final QC labels for {obstype}.')
#         #Get qc column namse specific for this obstype
#         specific_columns = [col for col in qc_labels_columns if col.startswith(obstype+'_')]
#         #add qc labels that are applicable on all obstypes
#         if 'missing_timestamp_label' in qc_labels_columns:
#             specific_columns.append('missing_timestamp_label')
        
        
#         #get labels dataframe
#         qc_df = df[specific_columns]
#         num_qc_df = pd.DataFrame()
#         num_qc_df = qc_df.applymap(label_to_numeric_mapper.get )
       
    
#         df[obstype+'_final_label'] = num_qc_df.sum(axis=1, skipna=True).map(inv_label_to_num)
    
#     #return only the final label series
#     return df[[obstype + '_final_label' for obstype in checked_obstypes]]

