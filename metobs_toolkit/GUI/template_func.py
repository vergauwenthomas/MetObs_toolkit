#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:00:45 2023

@author: thoverga
"""

import metobs_toolkit.GUI.path_handler as path_handler
import os
from pathlib import Path
import numpy as np
import pandas as pd

# =============================================================================
# default settings
# =============================================================================

not_mapable='NO MAPPING'

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

# =============================================================================
# Functions
# =============================================================================

def set_templ_vals(main, data_columns, metadata_columns):
    enable_all_boxes(main)
    csv_columns = data_columns
    csv_columns.insert(0, not_mapable)


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

    # 1. Read in the mapping as dictionaries
    obsmapper = {} #for observationtypes
    dtmapper = {} #for datetimes
    metamapper = {} #for metadata
    def get_obs_map(map_CB, unit_CB, desc_T):
        mapcolname = str(map_CB.currentText())
        if mapcolname != not_mapable:
            returndict = {
                'map_column': mapcolname,
                'map_unit': str(unit_CB.currentText()),
                'map_desc': str(desc_T.text())
                }
        else:
            returndict = {
                'map_column': np.nan,
                'map_unit': np.nan,
                'map_desc': np.nan
                }

        return returndict

    # add obsmapper to obsmapperdict
    obsmapper['temp'] = get_obs_map(main.temp_col_CB,
                                    main.temp_units_CB,
                                    main.temp_desc_T)
    obsmapper['radiation_temp'] = get_obs_map(main.radtemp_col_CB,
                                              main.radtemp_units_CB,
                                              main.radtemp_desc_T)

    obsmapper['humidity'] = get_obs_map(main.hum_col_CB,
                                        main.hum_units_CB,
                                        main.hum_desc_T)
    obsmapper['precip'] = get_obs_map(main.pre_col_CB,
                                      main.pre_units_CB,
                                      main.pre_desc_T)
    obsmapper['precip_sum'] = get_obs_map(main.pre_s_col_CB,
                                          main.pre_s_units_CB,
                                          main.pre_s_desc_T)
    obsmapper['wind_speed'] = get_obs_map(main.wind_col_CB,
                                          main.wind_units_CB,
                                          main.wind_desc_T)
    obsmapper['wind_gust'] = get_obs_map(main.gust_col_CB,
                                         main.gust_units_CB,
                                         main.gust_desc_T)
    obsmapper['wind_direction'] = get_obs_map(main.dir_col_CB,
                                              main.dir_units_CB,
                                              main.dir_desc_T)
    obsmapper['pressure'] = get_obs_map(main.p_col_CB,
                                        main.p_units_CB,
                                        main.p_desc_T)
    obsmapper['pressure_at_sea_level'] = get_obs_map(main.psea_col_CB,
                                                     main.psea_units_CB,
                                                     main.psea_desc_T)

    # add datatime to mapper
    def get_dt_map(map_CB, fmt_T):
        mapcolname = str(map_CB.currentText())
        if mapcolname != not_mapable:
            returndict = {
                'map_column': mapcolname,
                'map_fmt': str(fmt_T.text())
                }
        else:
            returndict = {
                'map_column': np.nan,
                'map_fmt': np.nan
                }

        return returndict

    dtmapper['datetime'] = get_dt_map(main.datetime_col_CB,
                                      main.datetime_fmt_T)
    dtmapper['_date'] = get_dt_map(main.date_col_CB,
                                  main.date_fmt_T)
    dtmapper['_time'] = get_dt_map(main.time_col_CB,
                                  main.time_fmt_T)
    # add metadata to mapper
    def get_meta_map(map_CB):
        mapcolname = str(map_CB.currentText())
        if mapcolname != not_mapable:
            returndict = {'map_column': mapcolname}
        else:
            returndict = {'map_column': np.nan}
        return returndict


    metamapper['name'] = get_meta_map(main.name_col_CB)
    metamapper['lat'] = get_meta_map(main.lat_col_CB)
    metamapper['lon'] = get_meta_map(main.lon_col_CB)
    metamapper['location'] = get_meta_map(main.loc_col_CB)
    metamapper['call_name'] = get_meta_map(main.call_col_CB)
    metamapper['network'] = get_meta_map(main.network_col_CB)

    # 2. Convert to a dataframe
    def to_dataframe(mapper, dtype):
        df = pd.DataFrame(mapper).transpose()
        df.index.name = 'varname'
        df = df.dropna(axis = 0, how = 'all')
        df = df.reset_index()
        df = df.rename(columns={'map_column': 'template column name',
                                'map_unit': 'units',
                                'map_desc': 'description',
                                'map_fmt': 'format'})
        df['dtype'] = dtype
        return df

    mapdf = pd.concat([to_dataframe(dtmapper,'object'),
                       to_dataframe(metamapper,'object'),
                       to_dataframe(obsmapper,'float64')])
    mapdf = mapdf.reset_index(drop=True)

    # 3. check if mapping is valid
    #TODO

    # 4. write mappingtemplate to csv
    tmp_csv_file = os.path.join(path_handler.TMP_dir,'template.csv')

    mapdf.to_csv(path_or_buf=tmp_csv_file, sep=',', index=False)

    return mapdf


def get_all_templates():
    """
    Returns a dict with keys the filenames of all available templates and keys the full paths.
    :return: DESCRIPTION
    :rtype: dict

    """
    template_dict = {}

    # default templates
    template_dict['default_template.csv'] = path_handler.tlk_default_template

    # all templates in cache
    filenames, filepaths = path_handler.list_csv_filenames(path_handler.CACHE_dir)
    template_dict.update(dict(zip(filenames, filepaths)))

    return template_dict









