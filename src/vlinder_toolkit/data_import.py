#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:24:06 2022

@author: thoverga
"""
import sys
# import json
# import datetime
import pandas as pd

import mysql.connector
from mysql.connector import errorcode


def template_to_package_space(specific_template):
    
   returndict = {val['varname']: {'orig_name': key} for key, val in specific_template.items()}
   for key, item in returndict.items():
       orig_dict = dict(specific_template[item['orig_name']])
       orig_dict.pop('varname')
       returndict[key].update(orig_dict)
   return returndict


    
def compress_dict(nested_dict, valuesname):
    """
    This function unnests a nested dictionary for a specific valuename that is a key in the nested dict. 

    Parameters
    ----------
    
    nested_dict : dict 
        Nested dictionary
    
    valuesname : str 
        Nested dict Key-name of nested dict.

    Returns
    -------
    returndict : DICT
        A dictionarry where the keys are kept that have the valuesname as a nesteddict key, 
        and values are the values of the values of the valuesname. 
        {[key-nested_dict-if-exists]: nested_dict[key-nested_dict-if-exists][valuesname]}

    """
    returndict = {}
    for key, item in nested_dict.items():
        if valuesname in item:
            returndict[key] = item[valuesname]
    return returndict    
    






# def get_template_from_df_columns(columns, Settings):
    
#     for test_template in Settings.template_list:
#         columnnames = compress_dict(test_template, 'varname').keys()
        
#         boollist = [col_name_csv in columnnames for col_name_csv in columns]
       
#         if not (False in list(set(boollist))): #all columnnames are found in template
#             return test_template



def import_data_from_csv(input_file, file_csv_template ):
    
    assert not isinstance(input_file, type(None)), "Specify input file in the settings!"    
    df = pd.read_csv(input_file, sep=';')
    
    
    assert not df.empty, "Dataset is empty!"
    
    
    

    # import template
    if isinstance(file_csv_template, type(None)):
        # templ =get_template_from_df_columns(df.columns)
        sys.exit("Not implemented yet !!!")
    else:
        templ = file_csv_template





    # rename columns to toolkit attriute names
    df = df.rename(columns=compress_dict(templ, 'varname'))
    
    
    #COnvert template to package-space
    template =template_to_package_space(templ)
    
    
    #format columns
    df = df.astype(dtype=compress_dict(template, 'dtype'))
    
    #create datetime column
    datetime_fmt = template['_date']['fmt'] + ' ' + template['_time']['fmt']
    df['datetime'] =pd.to_datetime(df['_date'] +' ' + df['_time'],
                                    format=datetime_fmt) 
    #TODO implement timezone settings
    
    
    #Set datetime index
    df = df.set_index('datetime', drop=True, verify_integrity=False)
    
    
    #drop 'date' and 'time' columns
    df = df.drop(columns=['_date', '_time'])
    
    
    #Keep only columns as defined in the template
    for column in df.columns:
        if not (column in template.keys()):
            df = df.drop(columns=[column])
    
    # add template to the return
    
    return df, template

#%%
def import_data_from_database(Settings,
                              start_datetime,
                              end_datetime):
   
    # =============================================================================
    # Make connection to database
    # =============================================================================
    
    #Make connection with database (needs ugent VPN active)
    
    try:
        connection = mysql.connector.connect(host=Settings.db_host,
                                      database=Settings.db_database,
                                      user=Settings.db_user,
                                      password=Settings.db_passw,
                                      connection_timeout=1)
    except mysql.connector.Error as err:
        if err.errno == errorcode.ER_ACCESS_DENIED_ERROR:
            print("Something is wrong with your user name or password!")
            print('Make shure the following envrionment variables are defind:')
            print('    VLINDER_DB_USER_NAME')
            print('    VLINDER_DB_USER_PASW')
            print('or update the Settings.db_user and Settings.db_passw')
            sys.exit()
        elif err.errno == 2003:
            print("Can't connect to ", Settings.db_host, ' host. Make shure your Ugent VPN is on!')
            sys.exit()
        
    
    
    # =============================================================================
    # Read all meta data from database    
    # =============================================================================
    
    
    metadata_Query = "select * from " + Settings.db_meta_table
    cursor = connection.cursor()
    cursor.execute(metadata_Query)
    metadata = cursor.fetchall()
    metadata = pd.DataFrame(metadata)
    # metadata_columns = list(cursor.column_names)
    
    metadata.columns = list(cursor.column_names)
    
    #subset relevent columns
    metadata = metadata[list(Settings.vlinder_db_meta_template.keys())]
    
    #rename columns to standards
    metadata = metadata.rename(columns=compress_dict(Settings.vlinder_db_meta_template, 'varname'))
    
    #COnvert template to package-space
    template = template_to_package_space(Settings.vlinder_db_meta_template)
    
    
    #format columns
    metadata = metadata.astype(dtype=compress_dict(template, 'dtype'))
    
    
    # =============================================================================
    # Read observations data 
    # =============================================================================
    
    
    
    assert start_datetime < end_datetime, "start_datetime is not earlier thand end_datetime!"
    
    observation_types = ['all']
    
    
    #observation types to strig
    if observation_types[0] == 'all':
        obs_type_query_str = '*'    
    else: #TODO
        print( "NOT IMPLEMENTED YET")
        obs_type_query_str = '*' 
    
    #format datetime
    
    datetime_db_info = [item for item in Settings.vlinder_db_obs_template.values() if item['varname'] == 'datetime'][0]
    
    startstring = start_datetime.strftime(format=datetime_db_info['fmt']) #datetime to string
    endstring = end_datetime.strftime(format=datetime_db_info['fmt']) #datetime to string
    _inverted_template = template_to_package_space(Settings.vlinder_db_obs_template)
    datetime_column_name = _inverted_template['datetime']['orig_name']
    
    
    #select all stations
    obsdata_Query=str(r"SELECT ") + obs_type_query_str + ' ' + \
            str(r"FROM ") + Settings.db_obs_table + str(' ') + \
            str(r"WHERE ") + datetime_column_name +  str(r">='") + startstring  + str(r"' AND ") + \
            datetime_column_name +  str(r"<='") + endstring  + str(r"'  ") + \
            str(r"ORDER BY ") + datetime_column_name
    
    
    
    
    
    
    print(obsdata_Query)
    
    
    
    cursor.execute(obsdata_Query)
    obsdata = cursor.fetchall()
    obsdata = pd.DataFrame(obsdata)
    
    
    
    obsdata.columns = list(cursor.column_names)
    
    #subset relevent columns
    obsdata = obsdata[list(Settings.vlinder_db_obs_template.keys())]
    
    #format columns
    obsdata = obsdata.astype(dtype=compress_dict(Settings.vlinder_db_obs_template, 'dtype'))
    
    #rename columns to standards
    obsdata = obsdata.rename(columns=compress_dict(Settings.vlinder_db_obs_template, 'varname'))
    
    
    connection.close()
    
    
    # =============================================================================
    # merge Observatios and metadata
    # =============================================================================
    
    combdata = obsdata.merge(metadata,
                             how='left',
                             on='id')
    combdata = combdata.drop(columns=['id'])
    combdata['datetime'] =pd.to_datetime(combdata['datetime'],
                                    format=datetime_db_info['fmt']) 
    #TODO implement timezone settings
    
    
    #Set datetime index
    combdata = combdata.set_index('datetime', drop=True, verify_integrity=False)
    
    return combdata

