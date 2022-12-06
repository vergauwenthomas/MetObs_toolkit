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
from .data_import import coarsen_time_resolution
from .landcover_functions import geotiff_point_extraction
from .geometry_functions import find_largest_extent
from .plotting_functions import spatial_plot, timeseries_plot, timeseries_comp_plot, qc_stats_pie
from .qc_checks import gross_value_check, persistance_check, missing_timestamp_check, duplicate_timestamp_check
from .qc_checks import step_check, internal_consistency_check
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

        self.df = pd.DataFrame() #Dataset with observations and labels per observations
        #TODO: a static metadf excludes moving observations :s
        self.metadf = pd.DataFrame() #Dataset with metadata (static)
        self.data_template = pd.DataFrame() #dataframe containing all information on the description and mapping
        
    
        
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
            return Station(name=stationname, df=sta_df, 
                           meta_series=sta_meta_series,
                           data_template=self.data_template)
        except KeyError:
            logger.warning(f'{stationname} not found in the dataset.')
            print(f'{stationname} not found in the dataset.')
            return None
        
        
    
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
            logger.error('The dataset is empty!')
        else: 
            starttimestr = datetime.strftime(min(self.df.index.get_level_values(level='datetime')),
                                             Settings.print_fmt_datetime)
            endtimestr = datetime.strftime(max(self.df.index.get_level_values(level='datetime')), Settings.print_fmt_datetime)
            
            stations_available = list(self.df.index.get_level_values(level='name').unique())
        
            print(f'Observations found for period: {starttimestr} --> {endtimestr}')
            logger.debug(f'Observations found for period: {starttimestr} --> {endtimestr}')
            print(f'Following stations are in dataset: {stations_available}')
            logger.debug(f'Following stations are in dataset: {stations_available}')
        
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
    
    
    def write_to_csv(self, filename=None, add_final_labels=True):
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
            self.add_final_qc_labels()
         
       
        #Get observations and metadata columns in the right order
        logger.debug('Merging data and metadata')
        
        
        #make column ordering
        df_columns = observation_types.copy() #observations
        df_columns.extend(location_info) #metadata
        df_columns.extend([col for col in self.df.columns if col.endswith('_label')]) #add qc labels
        df_columns.insert(0, 'datetime') # timestamp as first column
        
    
        
        #unstack observations and merge with metadf
        df = self.df.reset_index()
        metadf = self.metadf.reset_index()
        writedf = df.merge(metadf, how='left', on='name')
        
        #sort and subset columns
        writedf = writedf[df_columns]
        
                
        #find observation type that are not present
        ignore_obstypes = [col for col in observation_types if writedf[col].isnull().all()]
        
        writedf = writedf.drop(columns=ignore_obstypes)
        
        logger.debug(f'Skip quality labels for obstypes: {ignore_obstypes}.')
        

       
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
        writedf.to_csv(path_or_buf=filepath,
                       sep=';',
                       na_rep='NaN',
                       index=True)        
        
    

    
    # =============================================================================
    #     Quality control
    # =============================================================================
    
    def apply_quality_control(self, obstype='temp',
                              gross_value=True, persistance=True,
                              step=True, internal_consistency=True, ignore_val=np.nan):
        """
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
        #TODO save init observations
        
       
        
        if gross_value:
            print('Applying the gross-value-check on all stations.')
            logger.info('Applying gross value check on the full dataset')
            
       
            checked_obs, qc_flags = gross_value_check(input_series = self.df[obstype],
                                                obstype=obstype,
                                                ignore_val=ignore_val)
            #update the dataset
            self.df[obstype] = checked_obs
            label_column_name = qc_flags.name
            self.df[label_column_name] = qc_flags
            
        if persistance:
            print('Applying the persistance-check on all stations.')
            logger.info('Applying persistance check on the full dataset')
           
            checked_obs, qc_flags = persistance_check(input_series=self.df[obstype],
                                                      obstype=obstype,
                                                      ignore_val=ignore_val)

            #update the dataset
            self.df[obstype] = checked_obs
            label_column_name = qc_flags.name
            self.df[label_column_name] = qc_flags
            
            
        if step:
            print('Applying the step-check on all stations.')
            logger.info('Applying step-check on the full dataset')
           
            checked_obs, qc_flags = step_check(input_series=self.df[obstype],
                                                      obstype=obstype,
                                                      ignore_val=ignore_val)

            #update the dataset
            self.df[obstype] = checked_obs
            label_column_name = qc_flags.name
            self.df[label_column_name] = qc_flags
            
        if internal_consistency:
            print('Applying the internal-concsistency-check on all stations.')
            logger.info('Applying step-check on the full dataset')
           
            checked_obs, qc_flags = internal_consistency_check(input_series=self.df[obstype],
                                                               humidity_series=self.df['humidity'],
                                                               obstype=obstype,
                                                               ignore_val=ignore_val)

            #update the dataset
            self.df[obstype] = checked_obs
            label_column_name = qc_flags.name
            self.df[label_column_name] = qc_flags
    
    def get_qc_stats(self, obstype='temp', stationnames=None, make_plot=True):
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
        
        
        if isinstance(stationnames, type(None)):
            df = self.df
            title=f'Frequency for {obstype}-qc-checks on all stations.'
        else: 
            title=f'Frequency for {obstype}-qc-checks on {stationnames}.'
            if isinstance(stationnames, str):
                stationnames = [stationnames]
                
            df = self.df.loc[self.df.index.get_level_values(level='name').isin(stationnames)]
            
        #Do not include the 'final_qc_label' in the statis
        if obstype + '_final_label' in df.columns:
            logger.debug(f'The {obstype}_final_label, is ingored for the QC-stats. ')
            df = df.loc[:, df.columns != obstype+'_final_label']
            
        #stats on datset level
        dataset_qc_stats = get_qc_effectiveness_stats(df =df,
                                                      obstype=obstype,
                                                      observation_types = observation_types,
                                                      qc_labels=Settings.qc_observation_labels)
    
        if make_plot:
            qc_stats_pie(qc_stats=dataset_qc_stats,
                         figsize=Settings.plot_settings['qc_stats']['figsize'],
                         title=title)
               
        
        return dataset_qc_stats
    
    def add_final_qc_labels(self):
        """
        Add a final qualty label per quality-checked observation type to
        the dataset.df.

        Returns
        -------
        None.

        """
        #add final quality labels
    
        final_labels = final_qc_label_maker(df = self.df,
                                            label_to_numeric_mapper = Settings.qc_numeric_label_mapper)
        if isinstance(final_labels, type(pd.Series())):
            final_labels = final_labels.to_frame()
            
        for final_label_column in final_labels.columns:
            self.df[final_label_column] = final_labels[final_label_column]
            
        
        
    # =============================================================================
    #     importing data        
    # =============================================================================
            
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
        
        
        if coarsen_timeres:
            logger.info(f'Coarsen timeresolution to {Settings.target_time_res} using the {Settings.resample_method}-method.')
            df = coarsen_time_resolution(df=df,
                                          freq=Settings.target_time_res,
                                          method=Settings.resample_method)
            
        
        #convert dataframe to multiindex (datetime - name)
        df = df.set_index(['name', df.index])
        
        
       
        self.update_dataset_by_df(df)
        
    
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
        
        if coarsen_timeres:
            df = coarsen_time_resolution(df=df,
                                          freq=Settings.target_time_res,
                                          method=Settings.resample_method)
            
        #convert dataframe to multiindex (datetime - name)
        df = df.set_index(['name', df.index])
        df = df.sort_index()
        
        #If an ID has changed or not present in the metadatafile, the stationname and metadata is Nan
        #These observations will be removed
        unknown_obs = df[df.index.get_level_values('name').isnull()]
        if not unknown_obs.empty:
            logger.warning('There is an unknown station in the dataset (probaply due to an ID that is not present in the metadata file). This will be removed.')
            df = df[~df.index.get_level_values('name').isnull()]
        
        self.update_dataset_by_df(df)
    
    
    def update_dataset_by_df(self, dataframe):
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
       
        #Check import
        self.df = duplicate_timestamp_check(df=self.df)
        self.df = missing_timestamp_check(df=self.df)
        print(self.df)
        
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




def final_qc_label_maker(df, label_to_numeric_mapper):
    """
    This function creates a final label based on de individual qc labels. If all labels
    are ok, the final label is ok. Else the final label will be that of the individual qc-label
    which rejected the obseration.
    
    This functions converts labels to numeric values, algebra to get final label, and inversly 
    convert to labels. This is faster than looping over the rows.

    Parameters
    ----------
    qc_df : pandas.DataFrame 
        the specific qc_label_df with the datetimeindex, the first column the observations,
        and labels for each QC check per column.
    label_to_numeric_mapper : dict
        The dictionary that maps qc-labels to numeric values (for speedup).

    Returns
    -------
    final_labels : pd.Series
        A series with the final labels and the same index as the qc_df.

    """ 
        
    #invert numeric mapper
    inv_label_to_num = {v: k for k, v in label_to_numeric_mapper.items()}
    
    
    #extract label columns
    qc_labels_columns = [col for col in df.columns if not col in observation_types]
    #Extra savety
    qc_labels_columns = [col for col in qc_labels_columns if col.endswith('_label')]
    
    
    # find the observation types on which QC is applied
    checked_obstypes = [obstype for obstype in observation_types if any([qc_column.startswith(obstype+'_') for qc_column in qc_labels_columns])]
    
    
    #generete final label per obstype
    for obstype in checked_obstypes:
        logger.debug(f'Generating final QC labels for {obstype}.')
        #Get qc column namse specific for this obstype
        specific_columns = [col for col in qc_labels_columns if col.startswith(obstype+'_')]
        #add qc labels that are applicable on all obstypes
        if 'missing_timestamp_label' in qc_labels_columns:
            specific_columns.append('missing_timestamp_label')
        
        
        #get labels dataframe
        qc_df = df[specific_columns]
        num_qc_df = pd.DataFrame()
        num_qc_df = qc_df.applymap(label_to_numeric_mapper.get )
       
    
        df[obstype+'_final_label'] = num_qc_df.sum(axis=1, skipna=True).map(inv_label_to_num)
    
    #return only the final label series
    return df[[obstype + '_final_label' for obstype in checked_obstypes]]


