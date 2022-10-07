#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

from .settings import Settings
from .data_import import import_data_from_csv, import_data_from_database, template_to_package_space
from .data_import import coarsen_time_resolution
# from .qc_checks import duplicate_timestamp, gross_value_check

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

    # def apply_qc(self, obstype='temp', ignore_val=np.nan):
        
    #     self = duplicate_timestamp(self, obstype)
        
    #     self = gross_value_check(self, obstype, ignore_val)
        
        
        
        


    def make_plot(self, variable='temp', ax=None, **kwargs):
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
        
        data = getattr(self, variable)
        
        if isinstance(ax, type(None)):
            fig, ax = plt.subplots() 
        
        
        
        
        ax=data.plot(ax=ax, **kwargs)
        
        #Add text labels
        ax.set_title(self.name + ': ' + self.obs_description[variable])
        
        ax.set_xlabel('')
        ax.set_ylabel(self.units[variable])
        
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
        
    def show(self):
        if self.df.empty:
            print("This dataset is empty!")
        else: 
            starttimestr = datetime.strftime(min(self.df.index), Settings.print_fmt_datetime)
            endtimestr = datetime.strftime(max(self.df.index), Settings.print_fmt_datetime)
            
            stations_available = list(self.df.name.unique())
            
            print('Observations found for period: ', starttimestr, ' --> ', endtimestr)
            print('Following stations are in dataset: ', stations_available)
            
            
    # =============================================================================
    #     importing data        
    # =============================================================================
            
    def import_data_from_file(self, network='vlinder', coarsen_timeres=False):
        print('Settings input file: ', Settings.input_file)
        # Read observations into pandas dataframe
        df, template = import_data_from_csv(input_file = Settings.input_file,
                                  file_csv_template=Settings.input_csv_template)
        
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
        print('Importing data from ', Settings.input_file)
            
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
        
        print(dataframe.columns)
        
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
            if 'linder' in stationname:
                network='vlinder'
            elif 'occa' in stationname:
                network='mocca'
            else:
                network='Unknonw'
            
            
            #initialise station object
            station_obj = Station(station_name=stationname, 
                                  network_name=network)
            
            
            
            for obstype in observation_types:
                #fill attributes of station object
                setattr(station_obj, obstype, station_obs[obstype])
            
                #Fill QC dataframes with observations
                station_obj.qc_labels_df[obstype] = pd.DataFrame(data = {'observations': station_obs[obstype]})
                station_obj.qc_labels_df[obstype]['status'] = 'ok'

            
            
            #check if meta data is available
            if 'lat' in station_obs.columns:
                station_obj.lat = station_obs['lat']
            if 'lon' in station_obs.columns:
                station_obj.lon = station_obs['lon']
            if 'call_name' in station_obs.columns:
                station_obj.call_name = station_obs['call_name'].iloc[0]
            if 'location' in station_obs.columns:
                station_obj.location = station_obs['location'].iloc[0]
            
            #Update units and description dicts of the station using the used template
            for obs_field in station_obj.units.keys():
                try:
                    station_obj.units[obs_field] = self.data_template[obs_field]['units']
                    station_obj.obs_description[obs_field] = self.data_template[obs_field]['description']
                except:
                   continue 
            
            
            #update stationlist
            self._stationlist.append(station_obj)
            
            
          
            
        
        
    
    
