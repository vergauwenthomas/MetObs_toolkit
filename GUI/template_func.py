#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:00:45 2023

@author: thoverga
"""


Obs_map_values = {
    
     'temp':{'units': ['K', 'Celcius'], 'description':'2mT', 'dtype':'object'},
     'radiation_temp':{'units': ['K', 'Celcius'], 'description':'Blackglobe temperature', 'dtype':'float64'},
     'humidity':{'units': ['%'], 'description':'Relative humidity', 'dtype':'float64'},
     'precip':{'units': ['l/hour x m²'], 'description':'Precipitation intensity', 'dtype':'float64'},     
     'precip_sum':{'units': ['l/m²', 'Celcius'], 'description':'Precipitation cumulated', 'dtype':'float64'},
     'wind_speed':{'units': ['m/s'], 'description':'Average wind speed', 'dtype':'float64'},
     'wind_gust':{'units': ['m/s'], 'description':'Wind gust', 'dtype':'float64'},
     'wind_direction':{'units': ['°'], 'description':'Wind direction (from north CW)', 'dtype':'float64'},
     'pressure':{'units': ['Pa'], 'description':'Air pressure', 'dtype':'float64'},
     'pressure_at_sea_level':{'units': ['Pa'], 'description':'Pressure at sealevel', 'dtype':'float64'}
     
    }

Dt_map_values = {
    'datetime': {'format':'%Y-%m-%d %H:%M:%S', 'dtype': object},
    '_date':{'format':'%Y-%m-%d', 'dtype': object},
    '_time':{'format':'%H:%M:%S', 'dtype': object},
    }

Meta_map_values = {
    'name': {'description': 'Station name/ID', 'dtype':'object'},
    'lat':{'description': 'Latitude', 'dtype':'float64'},
    'lon':{'description': 'Longitude', 'dtype':'float64'},
    'location':{'description': 'Location identifier', 'dtype':'object'},
    'call_name':{'description': 'Pseudo name of station', 'dtype':'object'},
    'network':{'description': 'Name of the network', 'dtype':'object'}
    }



def set_templ_vals(main, data_columns, metadata_columns):
    enable_all_boxes(main)
    csv_columns = data_columns
    csv_columns.insert(0, 'NO MAPPING')
    
    
    # Get all obs boxes for colum mapping
    obs_boxes = [main.temp_col_CB, main.radtemp_col_CB, main.hum_col_CB,
                 main.pre_col_CB, main.pre_s_col_CB, main.wind_col_CB,
                 main.gust_col_CB, main.dir_col_CB, main.p_col_CB,
                 main.psea_col_CB, main.datetime_col_CB, main.date_col_CB, 
                 main.time_col_CB]
    
    for box in obs_boxes:
        box.setEnabled(True)
        box.addItems(csv_columns)
        
    metadata_boxes = [main.name_col_CB, main.lat_col_CB, main.lon_col_CB,
                      main.loc_col_CB, main.call_col_CB, main.network_col_CB]
    
    # combine all possible columns:
    csv_columns.extend(metadata_columns)
    for box in metadata_boxes:
        box.setEnabled(True)
        box.addItems(csv_columns)
    
    
    # =============================================================================
    # Observation types       
    # =============================================================================
    
    
    
    # Fill with defaults
    main.temp_units_CB.addItems(Obs_map_values['temp']['units'])
    main.temp_desc_T.setText(Obs_map_values['temp']['description'])
    
    main.radtemp_units_CB.addItems(Obs_map_values['radiation_temp']['units'])
    main.radtemp_desc_T.setText(Obs_map_values['radiation_temp']['description'])
    
    main.hum_units_CB.addItems(Obs_map_values['humidity']['units'])
    main.hum_desc_T.setText(Obs_map_values['humidity']['description'])

    main.pre_units_CB.addItems(Obs_map_values['precip']['units'])
    main.pre_desc_T.setText(Obs_map_values['precip']['description'])

    main.pre_s_units_CB.addItems(Obs_map_values['precip_sum']['units'])
    main.pre_s_desc_T.setText(Obs_map_values['precip_sum']['description'])
    
    main.wind_units_CB.addItems(Obs_map_values['wind_speed']['units'])
    main.wind_desc_T.setText(Obs_map_values['wind_speed']['description'])
    
    main.gust_units_CB.addItems(Obs_map_values['wind_gust']['units'])
    main.gust_desc_T.setText(Obs_map_values['wind_gust']['description'])
    
    main.dir_units_CB.addItems(Obs_map_values['wind_direction']['units'])
    main.dir_desc_T.setText(Obs_map_values['wind_direction']['description'])
    
    main.p_units_CB.addItems(Obs_map_values['pressure']['units'])
    main.p_desc_T.setText(Obs_map_values['pressure']['description'])
    
    main.psea_units_CB.addItems(Obs_map_values['pressure_at_sea_level']['units'])
    main.psea_desc_T.setText(Obs_map_values['pressure_at_sea_level']['description'])
    
    # =============================================================================
    #     Datetime
    # =============================================================================
     
    # Fill with defaults
    main.datetime_fmt_T.setText(Dt_map_values['datetime']['format'])
    main.date_fmt_T.setText(Dt_map_values['_date']['format'])
    main.time_fmt_T.setText(Dt_map_values['_time']['format'])
    
    
def enable_all_boxes(main):
    main.temp_col_CB.setEnabled(True)
    main.temp_units_CB.setEnabled(True)
    
    main.radtemp_units_CB.setEnabled(True)
    main.radtemp_desc_T.setEnabled(True)
    main.hum_units_CB.setEnabled(True)
    main.hum_desc_T.setEnabled(True)

    main.pre_units_CB.setEnabled(True)
    main.pre_desc_T.setEnabled(True)

    main.pre_s_units_CB.setEnabled(True)
    main.pre_s_desc_T.setEnabled(True)
    
    main.wind_units_CB.setEnabled(True)
    main.wind_desc_T.setEnabled(True)
    
    main.gust_units_CB.setEnabled(True)
    main.gust_desc_T.setEnabled(True)
    
    main.dir_units_CB.setEnabled(True)
    main.dir_desc_T.setEnabled(True)
    
    main.p_units_CB.setEnabled(True)
    main.p_desc_T.setEnabled(True)
    
    main.psea_units_CB.setEnabled(True)
    main.psea_desc_T.setEnabled(True)
   
    main.build_B.setEnabled(True)
    
def make_template_build(main):
    print('in builder')
    
    def get_obs_map(map_CB, unit_CB, desc_T):
        returndict = {
            'map_column': desc_T.getText()
            }
    test = get_obs_map(main.temp_col_CB, main.temp_units_CB, main.temp_desc_T)
    print(test)
    
    