#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct  6 13:25:02 2022
@author: thoverga
"""

#%%

# import metobs_toolkit

import os
import sys
from pathlib import Path


lib_folder = Path(__file__).resolve().parents[1]
sys.path.insert(0,str(lib_folder))


import metobs_toolkit



#%%




#%%



# #%% WIDE
# # wide
# file = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide.csv"
# # file_template ="/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_template.csv"
# file_metadata = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_metadata.csv"


# # long
# file = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/vlinderdata_small.csv"
# # file_template ="/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_template.csv"
# file_metadata = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/vlinder_metadata.csv"



# #%%
# import pandas as pd
# import numpy as np
# from datetime import datetime
# from metobs_toolkit.data_import import  _read_csv_file
# from metobs_toolkit import observation_types


# def col_option_input(columns):
#     mapper={}
#     i=1
#     for col in columns:
#         print(f'  {i}. {col}')
#         mapper[i] = col
#         i+=1

#     print(f'  x. -- not valid --')

#     if i <= 3:
#         repr_str = '('
#         for i in np.arange(1, i):
#             repr_str += str(i)+', '
#         # remove last comma
#         repr_str = repr_str[:-2] + ') : '
#         num_ans = (input(f'{repr_str}'))
#     else:
#         num_ans = (input(f'(1 - {i-1}) : '))

#     if num_ans == 'x':
#         print(' ... This setting is not provided! ...')
#         return None

#     print(f' ... {mapper[int(num_ans)]} selected ... \n')


#     return mapper[int(num_ans)]



# template_dict = {}
# print('This prompt will help to build a template for your data and metadata. Answer the prompt and hit Enter. \n \n')

# # =============================================================================
# # Map data file
# # =============================================================================

# print(' *******      Data File   ***********')

# datafilepath = input('Give the full path to your data file : ')
# datafilepath = file #debug
# print(' ... opening the data file ...')
# data = _read_csv_file(datafilepath)
# columnnames = data.columns.to_list()


# format_dict = {1: 'Long format', 2: 'Wide format', 3: 'Single station format'}

# format_option = int(input('How is your dataset structured : \n \
#     1. Long format (station observations are stacked as rows) \n \
#     2. Wide format (columns represent different stations) \n \
#     3. Single station (columns represent observation(s) of one station) \n \
# (1, 2 or 3) : '))

# print(f' \n... oke, {format_dict[format_option]} selected ...\n')

# #Datatime mapping
# datetime_option = int(input(' How are the timestamps represented in your dataset: \n \
#     1. In a single column (ex: 2023/06/07 16:12:30) \n \
#     2. By a column with dates, and another column with times \n \
#   (1 or 2) : '))

# if datetime_option ==1:

#     # Datetime mapping
#     template_dict['datetime'] = {}
#     print('\n Which is your timestamp columnname: ')
#     template_dict['datetime']['orig_name'] = col_option_input(columnnames)
#     columnnames.remove(template_dict['datetime']['orig_name'])

#     template_dict['datetime']['format'] = input('Type your datetime format (ex. %Y-%m-%d %H:%M:%S) : ')

# else:
#     # Date mapping
#     template_dict['_date'] = {}
#     print('Which column represents the DATES : ')
#     template_dict['_date']['orig_name'] = col_option_input(columnnames)
#     columnnames.remove(template_dict['_date']['orig_name'])
#     template_dict['_date']['format'] = input('Type your date format (ex. %Y-%m-%d) : ')
#     print(' \n')

#     # Time mapping
#     template_dict['_time'] = {}
#     print('Which column represents the TIMES : ')
#     template_dict['_time']['orig_name'] = col_option_input(columnnames)
#     columnnames.remove(template_dict['_time']['orig_name'])

#     template_dict['_time']['format'] = input('Type your time format (ex. %H:%M:%S) : ')

# # Obstype mapping in long format:
# obstype_desc = {
#     'name': 'name (name of the stations represented by strings)',
#     'temp': "temp (temperature)",
#     'radiation_temp' : "radiation_temp (radiation temperature)",
#     'humidity' : "humidity (humidity)",
#     'precip' : "precip (precipitation intensity)",
#     'precip_sum' : "precip_sum (precipitation cumulated)",
#     'wind_speed' : "wind_speed (wind speed)",
#     'wind_gust': "wind_gust (wind gust)",
#     'wind_direction' : "wind_direction (wind direction in degrees)" ,
#     'pressure' : "pressure (measured pressure)",
#     'pressure_at_sea_level' : "pressure_at_sea_level (altitude corrected pressure)"}
# inv_obstype_desc = {val: key for key, val in obstype_desc.items()}

# if (format_option == 1) | (format_option == 3):
#     # long format
#     print('What do the following columns represent: \n')
#     obstype_options=list(obstype_desc.values())
#     for col in columnnames:

#         contin = input(f'  add column {col} to the template? (y/n) ')
#         if contin != 'y':
#             continue

#         print(f'\n {col} : ')

#         desc_return = col_option_input(obstype_options)
#         if desc_return is None:
#             continue #when enter x
#         obstype = inv_obstype_desc[desc_return]


#         if obstype == 'temp':
#             _unit_num = input(' In Celcius (1), or Kelvin (2) : ')
#             units = {1:'Celcius', 2: 'Kelvin'}[int(_unit_num)]

#             print(units)
#         elif obstype == 'name':
#             template_dict['name'] = {'orig_name': col}
#             continue
#         else:
#             units = input(' What are the units : ')


#         description = input('Some more details on the observation : ')

#         # update template
#         template_dict[obstype] = {
#             'orig_name' : col,
#             'units': units,
#             'description': description
#             }
#         obstype_options.remove(obstype_desc[obstype])

# if format_option == 2:
#     print('Does these columns represent stations: ')
#     for col in columnnames:
#         print(f'  col')
#     cont = input('(y/n) ')
#     if cont!='y':

#         print('\n In a Wide-format, REMOVE THE COLUMNS that do not represent differnt satations, before proceding! \n')


# # debug
# print(f'\n   {template_dict} \n')
# # debug


# # =============================================================================
# # Map metadatafile
# # =============================================================================

# print('\n \n *******      Meta Data File   ***********')
# metatemplate_dict = {}
# meta_avail = input('The metadata contains metadata for each station, where each column represent a metadata type. \n Do you have a metadatafile? (y/n) ')
# if meta_avail == 'y':
#     metadatafilepath = input('Give the full path to your metadata file : ')
#     metadatafilepath = file_metadata #debug
#     print(' ... opening the metadata file ...')
#     metadata = _read_csv_file(metadatafilepath)
#     metacolumnnames = metadata.columns.to_list()


#     meta_desc = {
#         'name': 'name (the column with the stationnames, must be unique for each station)',
#         'lat':	'lat (the latitudes of the stations as a numeric values)',
#         'lon':	'lon (The longtitudes of the stations as a numeric values)',
#         'location': 'location (the city/region of the stations) (OPTIONAL)',
#         'call_name': 'call_name (an informal name of the stations) (OPTIONAL)',
#         'network': 'network (the name of the network the stations belong to) (OPTIONAL)',
#         }
#     inv_meta_desc = {val: key for key, val in meta_desc.items()}

#     print('What do the following columns represent: \n')
#     meta_options=list(meta_desc.values())
#     for col in metacolumnnames:
#         # debug
#         print(f'\n   {metatemplate_dict} \n')
#         # debug

#         contin = input(f'  add {col} to the template? (y/n) ')
#         if contin != 'y':
#             continue

#         print(f'\n {col} : ')
#         desc_return =col_option_input(meta_options)
#         if desc_return is None:
#             continue #when enter x
#         metatype = inv_meta_desc[desc_return]

#         # check if the name column is equalt in the data template to avoid creating
#         # two templates
#         if metatype == 'name':
#             if 'name' in template_dict:
#                 if not col == template_dict['name']['orig_name']:
#                     print(f'WARNING, the "name" column in the datafile is different than in the metadatafile! \
# Rename in your metadatafile : {col} ---> {template_dict["name"]["orig_name"]}')
#                     cont = input ('Renaming done? (y/n) ')
#                     if cont != 'y':
#                         sys.exit(f'Please rename {col} ---> {template_dict["name"]["orig_name"]} in your metadata file.')

#         metatemplate_dict[metatype] = {'orig_name': col}
#         meta_options.remove(meta_desc[metatype])


# print('\n   ... Oke, that is all the info for the mapping. Now i will do some basic tests to see if the mapping works.')


# #%%

# # template_dict = {'_date': {'orig_name': 'Datum', 'format': 'feumef'}, '_time': {'orig_name': 'Tijd (UTC)', 'format': 'moijmoij'}, 'name': {'orig_name': 'Vlinder'}}
# # metatemplate_dict = {'name': {'orig_name': 'Vlinder'}, 'lat': {'orig_name': 'lat'}, 'lon': {'orig_name': 'lon'}}

# # datafilepath = file #debug
# # print(' ... opening the data file ...')
# # data = _read_csv_file(datafilepath)

# # metadatafilepath = file_metadata #debug
# # print(' ... opening the metadata file ...')
# # metadata = _read_csv_file(metadatafilepath)

# # format_option=1
# #%%



# # =============================================================================
# # Apply tests
# # =============================================================================

# # apply tests the first row
# data_test = data.iloc[0].to_dict()
# metadata_test = metadata.iloc[0].to_dict()


# # test if name is in the metadat
# print (' *  ... checking metadata columns ... ')
# if ((bool(metatemplate_dict)) & (not 'name' in metatemplate_dict)):
#     print(f'Error! There is no metadata column containing the station names in the template! Add this column to the metadatafile of the template.')
#     sys.exit('Template invalid, see last message. ')



# # test if all stationnames are present in the metadata
# print (' *  ... checking compatible station names ... ')
# if ((bool(metatemplate_dict)) & (format_option == 1) & ('name' in template_dict)):
#     stanames_data = data[template_dict['name']['orig_name']].unique()
#     stanames_metadata = metadata[metatemplate_dict['name']['orig_name']].unique()

#     unmapped = [sta for sta in stanames_data if not sta in stanames_metadata]
#     if bool(unmapped):
#         print(f'Warning! The following stations are found in the data, but not in the metadata: {unmapped}')


# if ((bool(metatemplate_dict)) & (format_option == 2)):
#     stanames_data = data[template_dict['name']['orig_name']].unique().to_list()
#     stanames_metadata = metadata[metatemplate_dict['name']['orig_name']].unique().to_list()

#     unmapped = [sta for sta in stanames_data if not sta in stanames_metadata]
#     print(f'Warning! The following stations are found in the data, but not in the metadata: {unmapped}')


# # test if a stationname column is available in a long format
# print (' *  ... checking data columns ... ')
# if ((format_option == 1) & (not 'name' in template_dict)):
#     print(' \n WARNING: There is no information which column in the data file represents the names of the stations. The toolkit will assume that the observations are from ONE station! \n')

# # check if a least one mapped observation type exist
# if (format_option != 2):
#     present_obs = [key for key in template_dict.keys() if key in observation_types]
#     if not bool(present_obs):
#         print('ERROR! There is no observation type included in the template! Add at least one observation type when mapping the data file.')
#         sys.exit('Template invalid, see last message. ')





# # test datetime format
# print (' *  ... checking timestamps formats ... ')
# if 'datetime' in template_dict:
#     escape = False
#     while not escape:
#         test_dt = data_test[template_dict['datetime']['orig_name']]
#         try:
#             _ = datetime.strptime(test_dt,
#                                   template_dict['datetime']['format'])
#             print ('   ... testing datetime format is ...  OK!')
#             escape=True
#         except:
#             print(f'ERROR: the {template_dict["datetime"]["format"]} does not work for {test_dt}')
#             template_dict['datetime']['format'] = input('\n Try new timestamp format (ex. %Y-%m-%d %H:%M:%S) : ')

# if '_date' in template_dict:
#     escape = False
#     while not escape:
#         test_dt = data_test[template_dict['_date']['orig_name']]
#         try:
#             _ = datetime.strptime(test_dt,
#                                   template_dict['_date']['format'])
#             print ('   ... testing date format is OK!')
#             escape=True
#         except:
#             print(f'ERROR: the {template_dict["_date"]["format"]} does not work for {test_dt}')
#             template_dict['_date']['format'] = input('\n Try new date format (ex. %Y-%m-%d) : ')
# if '_time' in template_dict:
#     escape = False
#     while not escape:
#         test_dt = data_test[template_dict['_time']['orig_name']]
#         try:
#             _ = datetime.strptime(test_dt,
#                                   template_dict['_time']['format'])
#             print ('   ... testing time format is OK!')
#             escape=True
#         except:
#             print(f'ERROR: the {template_dict["_time"]["format"]} does not work for {test_dt}')
#             template_dict['_time']['format'] = input('\n Try new time format (ex. %H:%M:%S) : ')


# # check if all data columns are mapped
# print (' *  ... checking for unmapped data columns ... ')
# if (format_option == 1) | (format_option == 3):
#     present_columns = list(data_test.keys())
#     mapped_cols = [val['orig_name'] for val in template_dict.values()]
#     for col in present_columns:
#         if not col in mapped_cols:
#             print(f' Warning! {col} in the datafile is not present in the template, and thus it will not be used.')



# # check if all metadata columns are mapped
# print (' *  ... checking for unmapped metadata columns ... ')
# present_columns = list(metadata_test.keys())
# mapped_cols = [val['orig_name'] for val in metatemplate_dict.values()]
# for col in present_columns:
#     if not col in mapped_cols:
#         print(f' Warning! {col} in the metadatafile is not present in the template, and thus it will not be used.')

# # =============================================================================
# # Saving the template
# # =============================================================================

# print('\n ------ Saving the template ----- \n')
# save_dir = input('Give a directory where to save the template : ')
# template_dict.update(metatemplate_dict) #this is why name in data and metadata should have the same mapping !!


# # Convert to dataframe


# df = pd.DataFrame(template_dict).transpose()
# df.index.name = 'varname'
# df = df.rename(columns={'orig_name' : 'template column name'})
# df = df.reset_index()

# # write to csv
# filepath = os.path.join(save_dir, 'template.csv')
# df.to_csv(filepath, na_rep = '', index=False)
# print(f' DONE! The template is writen here: {filepath}')




# #%%
# # Convert to dataframe
# temp = {'_date': {'orig_name': 'Datum', 'format': '%Y-%m-%d'}, '_time': {'orig_name': 'Tijd (UTC)', 'format': '%H:%M:%S'}, 'temp': {'orig_name': 'Temperatuur', 'units': 'Celcius', 'description': '2mT passive'}, 'humidity': {'orig_name': 'Vochtigheid', 'units': 'percentage', 'description': 'at 2m'}, 'name': {'orig_name': 'Vlinder'}, 'lat': {'orig_name': 'lat'}, 'lon': {'orig_name': 'lon'}}

# df = pd.DataFrame(temp).transpose()
# df.index.name = 'varname'
# df = df.rename(columns={'orig_name' : 'template column name'})
# df = df.reset_index()

# # write to csv
















#%%

# # Make an empty dataset
# dataset = metobs_toolkit.Dataset()

# # Add the demo data files to the dataset settings
# dataset.update_settings(input_data_file = file,
#                         input_metadata_file = file_metadata,
#                         data_template_file = file_template,
#                         metadata_template_file = file_template, # Contains also the metadata mapping
#                         )

# # Load the data from the demo data files
# dataset.import_data_from_file(long_format=False,
#                               obstype='temp',
#                               obstype_dtype='float',
#                               obstype_description='oijmojoijomm',
#                               obstype_unit='Celcius')

# print(dataset.data_template)
# # dataset.coarsen_time_resolution()
# #%% Long

# file = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/vlinderdata_small.csv"
# # file_template ="/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/tests/test_data/debug_wide_template.csv"
# file_metadata = "/home/thoverga/Documents/VLINDER_github/MetObs_toolkit/static_data/vlinder_metadata.csv"

# # Make an empty dataset
# dataset_long = metobs_toolkit.Dataset()

# # Add the demo data files to the dataset settings
# dataset_long.update_settings(input_data_file = file,
#                         input_metadata_file = file_metadata,
#                         # data_template_file = file_template,
#                         # metadata_template_file = file_template, # Contains also the metadata mapping
#                         )

# # Load the data from the demo data files
# dataset_long.import_data_from_file(long_format=True)

#%%
# dataset.make_plot()
