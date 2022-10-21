#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""
from collections.abc import Iterable

import pandas as pd
import geopandas as gpd
import numpy as np
# import matplotlib.pyplot as plt
# import geoplot as gplt
# import mapclassify as mc
from datetime import datetime

from .settings import Settings
from .data_import import import_data_from_csv, import_data_from_database, template_to_package_space, import_metadata_from_csv
from .data_import import coarsen_time_resolution
from .landcover_functions import geotiff_point_extraction
from .geometry_functions import find_largest_extent
from .plotting_functions import spatial_plot, timeseries_plot, timeseries_comp_plot
# from .qc_checks import duplicate_timestamp, gross_value_check


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




# =============================================================================
# station class
# =============================================================================
class Station:
    def __init__(self, station_name, network_name):
        self.network = network_name
        self.name = station_name
        
        #Meta data without processing
        self.lat = pd.Series()
        self.lon = pd.Series()
        self.call_name = None #ex. Antwerpen Zoo
        self.location = None #ex. Antwerpen 
        
        #Observations
        self.temp = pd.Series()
        self.radiation_temp = pd.Series() 
        
        self.humidity = pd.Series()
        
        self.precip = pd.Series()
        self.precip_sum = pd.Series()
        
        self.wind_speed = pd.Series()
        self.wind_gust = pd.Series()
        self.wind_direction = pd.Series()
        
        self.pressure = pd.Series()
        self.pressure_at_sea_level = pd.Series()
        
        #physiographic data
        self.lcz = None
        
        #Units and descriptions
        self.units = {'temp': None,
                      'radiation_temp': None,
                      'humidity': None,
                      'precip': None,
                      'precip_sum': None,
                      'wind_speed': None,
                      'wind_gust': None,
                      'wind_direction': None,
                      'pressure': None,
                      'pressure_at_sea_level': None}
        self.obs_description = {'temp': None,
                          'radiation_temp': None,
                          'humidity': None,
                          'precip': None,
                          'precip_sum': None,
                          'wind_speed': None,
                          'wind_gust': None,
                          'wind_direction': None,
                          'pressure': None,
                          'pressure_at_sea_level': None}
        
        #attributes that can be filled with info from other functions
        self.qc_labels_df =  {'temp': pd.DataFrame(),
                              'radiation_temp': pd.DataFrame(),
                              'humidity': pd.DataFrame(),
                              'precip': pd.DataFrame(),
                              'precip_sum': pd.DataFrame(),
                              'wind_speed': pd.DataFrame(),
                              'wind_gust': pd.DataFrame(),
                              'wind_direction': pd.DataFrame(),
                              'pressure': pd.DataFrame(),
                              'pressure_at_sea_level': pd.DataFrame()}
        
        #TODO: remove the qc_info and combine all info in the qc_labels_df
        self.qc_info = {} #will be filled by qc checks

    def show(self):
        starttimestr = datetime.strftime(min(self.df().index), Settings.print_fmt_datetime)
        endtimestr = datetime.strftime(max(self.df().index), Settings.print_fmt_datetime)
        
        
        print(' ------- Static info -------')
        print('Stationname: ', self.name)
        print('Network: ', self.network)
        print('Call name: ', self.call_name)
        print('Location: ', self.location)
        print('latitude: ', self.lat)
        print('longtitude: ',self.lon)
       
        print(' ------- Physiography info -------')
        print('LCZ: ', self.lcz)
        print(' ')
        print(' ------- Observations info -------')
    
        print('Observations found for period: ', starttimestr, ' --> ', endtimestr)

        
    def df(self):
        """
        Convert all observations of the station to a pandas dataframe.

        Returns
        -------
        pandas.DataFrame
            A Dataframe containing all observations with a datetime index.

        """
        return pd.DataFrame([self.temp,
                             self.radiation_temp, 
                            self.humidity,
                            self.precip,
                            self.precip_sum,          
                            self.wind_speed,
                            self.wind_gust,
                            self.wind_direction,                   
                            self.pressure,
                            self.pressure_at_sea_level]).transpose()


        
        
    def get_lcz(self):
        
        geo_templates = Settings.geo_datasets_templates
        lcz_file = Settings.geo_lcz_file
        
        if isinstance(lcz_file, type(None)):
            print('No lcz tif location in the settings. Update settings: ')
            print('settings_obj.update_settings(geotiff_lcz_file="...."')
            return
        
        lcz_templates = [geo_templ for geo_templ in geo_templates if geo_templ['usage']=='LCZ']
        
        assert len(lcz_templates)==1, 'More (or no) lcz template found!'
        
        lcz_template = lcz_templates[0]
        
        human_mapper = {num: lcz_template['covers'][num]['cover_name'] 
                        for num in lcz_template['covers'].keys()}
        

        #Check if coordinates are available
        if np.isnan(self.lat.iloc[0]) | np.isnan(self.lon.iloc[0]):
            self.lcz = 'Location unknown'
            return 'Location unknown'
        
        #TODO: lat and lons are time depending, now use first lat, lon

        lcz = geotiff_point_extraction(lat=self.lat.iloc[0],
                                       lon=self.lon.iloc[0],
                                       geotiff_location=lcz_file,
                                       geotiff_crs=lcz_template['epsg'],
                                       class_to_human_mapper=human_mapper)

        self.lcz = lcz
        return lcz

    def make_plot(self, variable='temp', title=None):
        """
        Make a timeseries plot of one attribute.

        Parameters
        ----------
        variable : str, optional
            Name of attribute to plot. Must be one of [temp, radiation_temp, humidity, precip, wind_speed wind_gust, wind_direction, pressure, pressure_at_sea_level].
            The default is 'temp'.
        **kwargs : 
            named-arguments that are passed to matplolib.pyplot.plot()

        Returns
        -------
        ax : AxesSubplot
            AxesSubplot is returned so layer can be added to it.

        """
        
        #Load default plot settings
        default_settings=Settings.plot_settings['time_series']
      
        
        #Make title
        if isinstance(title, type(None)):
            title = self.name + ': ' + self.obs_description[variable]
    
        
        #make figure
        ax = timeseries_plot(dtseries=getattr(self, variable),
                             title=title,
                             xlabel='',
                             ylabel=self.units[variable],
                             figsize=default_settings['figsize'],
                             )
                        
        
        
        return ax
        

# =============================================================================
# Dataset class
# =============================================================================

class Dataset:
    def __init__(self):
        self._stationlist = []
        self.df = pd.DataFrame()
        
        self.data_template = {}
        
    
    def get_station(self, stationname):
        """
        Extract a station object from the dataset.

        Parameters
        ----------
        stationname : String
            Name of the station, example 'vlinder16'

        Returns
        -------
        station_obj : vlinder_toolkit.Station
            

        """
        for station_obj in self._stationlist:
            if stationname == station_obj.name:
                return station_obj
        
        print(stationname, ' not found in the dataset!')
    
    def get_geodataframe(self):
        gdf = gpd.GeoDataFrame(self.df,
                               geometry=gpd.points_from_xy(self.df['lon'],
                                                           self.df['lat']))
        return gdf
    
    def show(self):
        if self.df.empty:
            print("This dataset is empty!")
        else: 
            starttimestr = datetime.strftime(min(self.df.index), Settings.print_fmt_datetime)
            endtimestr = datetime.strftime(max(self.df.index), Settings.print_fmt_datetime)
            
            stations_available = list(self.df.name.unique())
        
            print('Observations found for period: ', starttimestr, ' --> ', endtimestr)
            print('Following stations are in dataset: ', stations_available)
        
        
    def make_plot(self, stationnames, variable='temp',
                                   starttime=None, endtime=None,
                                   title=None, legend=True):
        
        default_settings=Settings.plot_settings['time_series']
        
        plotdf = self.df[self.df['name'].isin(stationnames)]
        
        if not isinstance(starttime, type(None)):
            plotdf = plotdf.loc[(plotdf >= starttime)]
        if not isinstance(endtime, type(None)):
            plotdf = plotdf.loc[(plotdf <= endtime)]
        
        
        relevant_columns = ['name']
        relevant_columns.append(variable)
        plotdf = plotdf[relevant_columns]
        
        plotdf = pd.pivot(plotdf,
                          columns='name',
                          values=variable)
        
        
        if isinstance(title, type(None)):
            title=Settings.display_name_mapper[variable] + ' for stations: ' + str(stationnames)
        
        ax = timeseries_comp_plot(plotdf=plotdf,
                                  title=title,
                                  xlabel='',
                                  ylabel=Settings.display_name_mapper[variable],
                                  figsize=default_settings['figsize'])
        
        return ax
        
        
        
        
        
    def make_geo_plot(self, variable='temp', title=None, timeinstance=None, legend=True,
                      vmin=None, vmax=None):
        
        #Load default plot settings
        default_settings=Settings.plot_settings['spatial_geo']
        
        #get first timeinstance of the dataset if not given
        if isinstance(timeinstance, type(None)):
            timeinstance=self.df.index.min()
        
        #subset to timeinstance
        subdf = self.df.loc[timeinstance]
        #create geodf
        gdf = gpd.GeoDataFrame(subdf,
                               geometry=gpd.points_from_xy(subdf['lon'], subdf['lat']))
        
        gdf = gdf[[variable, 'geometry']]
        
    
        

        #make color scheme for field
        if variable in categorical_fields:
            if variable == 'lcz':
                #use all available LCZ categories
                use_quantiles=False
            else:
                use_quantiles=True
        else:
            use_quantiles=False
     
        
        #if observations extend is contained by default exten, use default else use obs extend
        use_extent=find_largest_extent(geodf=gdf,
                                       extentlist=default_settings['extent'])
        
        
        #Style attributes
        if isinstance(title, type(None)):
            if variable in static_fields:
                title = Settings.display_name_mapper[variable]
            else:
                dtstring = datetime.strftime(timeinstance, default_settings['fmt'])
                title = Settings.display_name_mapper[variable] + ' at ' + dtstring
        
        ax = spatial_plot(gdf=gdf,
                          variable=variable,
                          legend=legend,
                          proj_type=default_settings['proj'],
                          use_quantiles=use_quantiles,
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
            
    # =============================================================================
    #     importing data        
    # =============================================================================
            
    def import_data_from_file(self, network='vlinder', coarsen_timeres=False):
        print('Settings input data file: ', Settings.input_data_file)
        
        
        # Read observations into pandas dataframe
        df, template = import_data_from_csv(input_file = Settings.input_data_file,
                                  file_csv_template=Settings.input_csv_template)
        
        
        if isinstance(Settings.input_metadata_file, type(None)):
            print('WARNING: No metadata file is defined. Add your settings object.')
        else:
            meta_df = import_metadata_from_csv(input_file=Settings.input_metadata_file,
                                               file_csv_template=Settings.input_metadata_template)
            
            #merge additional metadata to observations
            meta_cols = [colname for colname in meta_df.columns if not colname.startswith('_')]
            additional_meta_cols = list(set(meta_cols).difference(df.columns))
            if bool(additional_meta_cols):
                additional_meta_cols.append('name') #merging on name
                df_index = df.index #merge deletes datetime index somehow? so add it back on the merged df
                df = df.merge(right=meta_df[additional_meta_cols],
                              how='left', 
                              on='name')
                df.index = df_index
        
        
        #update dataset object
        self.data_template = template
        
        if coarsen_timeres:
            df = coarsen_time_resolution(df=df,
                                         freq=Settings.target_time_res,
                                         method=Settings.resample_method)
            
        
        self.update_dataset_by_df(df)
        
    
    def import_data_from_database(self,
                              start_datetime=None,
                              end_datetime=None,
                              coarsen_timeres=False):
        if isinstance(start_datetime, type(None)):
            start_datetime=datetime.date.today() - datetime.timedelta(days=1)
        if isinstance(end_datetime, type(None)):
            end_datetime=datetime.date.today()
       
            
        # Read observations into pandas dataframe
        df = import_data_from_database(Settings,
                                       start_datetime=start_datetime,
                                       end_datetime=end_datetime)
        
        
        #Make data template
        self.data_template = template_to_package_space(Settings.vlinder_db_obs_template)
        
        if coarsen_timeres:
            df = coarsen_time_resolution(df=df,
                                         freq=Settings.target_time_res,
                                         method=Settings.resample_method)
        
        self.update_dataset_by_df(df)
    
    
    
            
    def update_dataset_by_df(self, dataframe):
        """
        Update the dataset object and all it attributes by a dataframe.

        Parameters
        ----------
        dataframe : pandas.DataFrame
        A dataframe that has an datetimeindex and following columns: 'name, temp, radiation_temp, humidity, ...'
            

        Returns
        -------
        None.

        """
       
        
        #reset dataset attributes
        self.df = dataframe
        self._stationlist = [] 
        
      
        
        observation_types = ['temp', 'radiation_temp', 'humidity', 'precip',
                             'precip_sum', 'wind_speed', 'wind_gust', 'wind_direction',
                             'pressure', 'pressure_at_sea_level']
        
        # Create a list of station objects
        for stationname in dataframe.name.unique():
            #extract observations
            station_obs = dataframe[dataframe['name'] == stationname].sort_index()
            
            if station_obs.empty:
                print('skip stationname: ', stationname)
                continue
            
            #find network
            if 'network' in station_obs.columns:
                network = station_obs['network'].iloc[0]
            else:
                if 'linder' in stationname:
                    network='vlinder'
                elif 'occa' in stationname:
                    network='mocca'
                else:
                    network='Unknown'
            
            
            
            #initialise station object
            station_obj = Station(station_name=stationname, 
                                  network_name=network)
            
            
            
            for obstype in observation_types:
                #fill attributes of station object
                setattr(station_obj, obstype, station_obs[obstype])
            
                #Fill QC dataframes with observations
                station_obj.qc_labels_df[obstype] = pd.DataFrame(data = {'observations': station_obs[obstype]})
                station_obj.qc_labels_df[obstype]['status'] = 'ok'

            #Apply IO checks
            
            #check for missing timestamps
            checked_df, missing_dt_list = missing_timestamp_check(station_obj)
            if bool(missing_dt_list):
                for obstype in checked_df.columns:
                    #update observations with missing obs as nan's
                    setattr(station_obj, obstype, checked_df[obstype]) 
                    #update QC dataframes
                    station_obj.qc_labels_df[obstype] = pd.DataFrame(data = {'observations': checked_df[obstype]})
                    station_obj.qc_labels_df[obstype]['status'] = ['ok' if dt not in missing_dt_list else 'missing timestamp' for dt in checked_df.index ]
                    
                    
                    
            
            
            
            #check if meta data is available
            if 'lat' in station_obs.columns:
                station_obj.lat = station_obs['lat']
                check_for_nan(station_obj.lat, 'latitude', stationname)
            if 'lon' in station_obs.columns:
                station_obj.lon = station_obs['lon']
                check_for_nan(station_obj.lon, 'longtitude', stationname)
            if 'call_name' in station_obs.columns:
                station_obj.call_name = station_obs['call_name'].iloc[0]
                check_for_nan(station_obj.call_name, 'call_name', stationname)
            if 'location' in station_obs.columns:
                station_obj.location = station_obs['location'].iloc[0]
                check_for_nan(station_obj.location, 'location', stationname)
                
            # Get physiography data if possible
            if not isinstance(Settings.geo_lcz_file, type(None)):
                _ = station_obj.get_lcz()
                
            
            #Update units and description dicts of the station using the used template
            for obs_field in station_obj.units.keys():
                try:
                    station_obj.units[obs_field] = self.data_template[obs_field]['units']
                    station_obj.obs_description[obs_field] = self.data_template[obs_field]['description']
                except:
                   continue 
            
            
            #update stationlist
            self._stationlist.append(station_obj)
            
        
        #Update dataset df with information created on station level
        
        #add LCZ to dataset df
        if not isinstance(Settings.geo_lcz_file, type(None)):
            lcz_dict = {station.name: station.lcz for station in self._stationlist}
            self.df['lcz'] = self.df['name'].map(lcz_dict)
            
            
          
def check_for_nan(value, fieldname, stationname):
    if isinstance(value, float):
        if np.isnan(value):
            print('Nan found for ', fieldname, ' in ', stationname, '!!')
    elif isinstance(value, type(pd.Series())):
        if value.isnull().sum() > 0:
            n_nans = value.isnull().sum()
            print(n_nans, "Nan's found in ", fieldname, '-iterable in ', stationname, '!!')
        

        
def missing_timestamp_check(station):
    """
    Looking for missing timestaps by assuming an observation frequency. The assumed frequency is the most occuring frequency.
    If missing observations are detected, the observations dataframe is extended by these missing timestamps with Nan's as filling values.

    Parameters
    ----------
    station : Station object
        The station you whant to apply this check on.

    Returns
    -------
    df : pandas.DataFrame()
        The observations dataframe (same as Station.df()).
    missing_datetimes : list of datetimes
        The list of the missing timestamps.

    """     
   
    
    df = station.df()
    
    #extrac observed frequencies
    likely_freq = df.index.to_series().diff().value_counts().idxmax()
    
    
    missing_datetimeindices = pd.date_range(start = df.index.min(),
                         end = df.index.max(),
                         freq=likely_freq).difference(df.index)
    
    missing_df = pd.DataFrame(data=np.nan,
                              index=missing_datetimeindices,
                              columns=df.columns)
    
    df = df.append(missing_df)
    
    
    df = df.sort_index()
    
    
    return df, missing_datetimeindices.to_list()
    
    
