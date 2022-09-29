#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
The class object for a Vlinder/mocca station
@author: thoverga
"""

import pandas as pd
# import matplotlib.pyplot as plt
from datetime import datetime

from .IO import import_data_from_csv 
from .physical_info import datetime_settings, description_dict, unit_dict


# =============================================================================
# station class
# =============================================================================
class Station:
    def __init__(self, station_name, network_name):
        self.network = network_name
        self.name = station_name
        
        #Meta data without processing
        self.lat = []
        self.lon = []
        self.location = None #ex. Boerenkreek
        
        
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


    def plot(self, variable='temp', **kwargs):
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
        
        ax=data.plot(**kwargs)
        
        #Add text labels
        ax.set_title(self.name + ': ' + description_dict[variable])
        
        ax.set_xlabel('')
        ax.set_ylabel(unit_dict[variable])
        
        return ax
        

# =============================================================================
# Dataset class
# =============================================================================

class Dataset:
    def __init__(self):
        self._stationlist = []
        self.df = pd.DataFrame()
        
    
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
            starttimestr = datetime.strftime(min(self.df.index), datetime_settings['string_representation_format'])
            endtimestr = datetime.strftime(max(self.df.index), datetime_settings['string_representation_format'])
            
            stations_available = list(self.df.name.unique())
            
            print('Observations found for period: ', starttimestr, ' --> ', endtimestr)
            print('Following stations are in dataset: ', stations_available)
            
            
            
            
    def import_data_from_file(self, Settings, network='vlinder'):
        
        # Read observations into pandas dataframe
        df = import_data_from_csv(Settings)
        
        #update dataset object
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
        
        # Create a list of station objects
        for stationname in dataframe.name.unique():
            
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
            #extract observations
            station_obs = dataframe[dataframe['name'] == stationname].sort_index()
            
            #fill attributes of station object
            station_obj.temp = station_obs['temp']
            station_obj.radiation_temp = station_obs['radiation_temp']
            
            station_obj.humidity = station_obs['humidity']
            
            station_obj.precip = station_obs['precip']
            station_obj.precip_sum = station_obs['precip_sum']
            
            station_obj.wind_speed = station_obs['wind_speed']
            station_obj.wind_gust = station_obs['wind_gust']
            station_obj.wind_direction = station_obs['wind_direction']
            
            station_obj.pressure = station_obs['pressure']
            station_obj.pressure_at_sea_level = station_obs['pressure_at_sea_level']
            
            #update stationlist
            self._stationlist.append(station_obj)
        
        
        
    
    
