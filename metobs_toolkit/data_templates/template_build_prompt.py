#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:25:45 2023

@author: thoverga
"""

import os, sys
import pandas as pd
import numpy as np
from datetime import datetime
import pytz

from metobs_toolkit.data_import import  _read_csv_file
from metobs_toolkit import observation_types



def col_option_input(columns):
    mapper={}
    i=1
    for col in columns:
        print(f'  {i}. {col}')
        mapper[i] = col
        i+=1

    print(f'  x. -- not valid --')
    valid_input = False
    while valid_input == False:
        if i <= 3:
            repr_str = '('
            for i in np.arange(1, i):
                repr_str += str(i)+', '
            # remove last comma
            repr_str = repr_str[:-2] + ') : '
            num_ans = (input(f'{repr_str}'))
        else:
            num_ans = (input(f'(1 - {i-1}) : '))

        if num_ans == 'x':
            print(' ... This setting is not provided! ...')
            return None

        try:
            _ = mapper[int(num_ans)]
            valid_input = True
        except:
            valid_input = False
            print(f'{num_ans} is not a valid input.')


    print(f' ... {mapper[int(num_ans)]} selected ... \n')


    return mapper[int(num_ans)]








def build_template_prompt():
    template_dict = {}
    print('This prompt will help to build a template for your data and metadata. Answer the prompt and hit Enter. \n \n')

    # =============================================================================
    # Map data file
    # =============================================================================

    print(' *******      Data File   ***********')

    datafilepath = input('Give the full path to your data file : ')
    print(' ... opening the data file ...')
    data = _read_csv_file(datafilepath)
    columnnames = data.columns.to_list()


    format_dict = {1: 'Long format', 2: 'Wide format', 3: 'Single station format'}

    format_option = int(input('How is your dataset structured : \n \
        1. Long format (station observations are stacked as rows) \n \
        2. Wide format (columns represent different stations) \n \
        3. Single station (columns represent observation(s) of one station) \n \
    (1, 2 or 3) : '))

    print(f' \n... oke, {format_dict[format_option]} selected ...\n')

    #Datatime mapping
    datetime_option = int(input(' How are the timestamps represented in your dataset: \n \
        1. In a single column (ex: 2023/06/07 16:12:30) \n \
        2. By a column with dates, and another column with times \n \
      (1 or 2) : '))

    if datetime_option ==1:

        # Datetime mapping
        template_dict['datetime'] = {}
        print('\n Which is your timestamp columnname: ')
        template_dict['datetime']['orig_name'] = col_option_input(columnnames)
        columnnames.remove(template_dict['datetime']['orig_name'])

        template_dict['datetime']['format'] = input('Type your datetime format (ex. %Y-%m-%d %H:%M:%S) : ')

    else:
        # Date mapping
        template_dict['_date'] = {}
        print('Which column represents the DATES : ')
        template_dict['_date']['orig_name'] = col_option_input(columnnames)
        columnnames.remove(template_dict['_date']['orig_name'])
        template_dict['_date']['format'] = input('Type your date format (ex. %Y-%m-%d) : ')
        print(' \n')

        # Time mapping
        template_dict['_time'] = {}
        print('Which column represents the TIMES : ')
        template_dict['_time']['orig_name'] = col_option_input(columnnames)
        columnnames.remove(template_dict['_time']['orig_name'])

        template_dict['_time']['format'] = input('Type your time format (ex. %H:%M:%S) : ')

    # Obstype mapping in long format:
    obstype_desc = {
        'name': 'name (name of the stations represented by strings)',
        'temp': "temp (temperature)",
        'radiation_temp' : "radiation_temp (radiation temperature)",
        'humidity' : "humidity (humidity)",
        'precip' : "precip (precipitation intensity)",
        'precip_sum' : "precip_sum (precipitation cumulated)",
        'wind_speed' : "wind_speed (wind speed)",
        'wind_gust': "wind_gust (wind gust)",
        'wind_direction' : "wind_direction (wind direction in degrees)" ,
        'pressure' : "pressure (measured pressure)",
        'pressure_at_sea_level' : "pressure_at_sea_level (altitude corrected pressure)"}
    inv_obstype_desc = {val: key for key, val in obstype_desc.items()}

    if (format_option == 1) | (format_option == 3):
        # long format
        print('What do the following columns represent: \n')
        obstype_options=list(obstype_desc.values())
        for col in columnnames:

            contin = input(f'  add column {col} to the template? (y/n) ')
            if contin != 'y':
                continue

            print(f'\n {col} : ')

            desc_return = col_option_input(obstype_options)
            if desc_return is None:
                continue #when enter x
            obstype = inv_obstype_desc[desc_return]


            if obstype == 'temp':
                _unit_num = input(' In Celcius (1), or Kelvin (2) : ')
                units = {1:'Celcius', 2: 'Kelvin'}[int(_unit_num)]

                print(units)
            elif obstype == 'name':
                template_dict['name'] = {'orig_name': col}
                continue
            else:
                units = input(' What are the units : ')


            description = input('Some more details on the observation : ')

            # update template
            template_dict[obstype] = {
                'orig_name' : col,
                'units': units,
                'description': description
                }
            obstype_options.remove(obstype_desc[obstype])

    if format_option == 2:
        print('Does these columns represent stations: ')
        for col in columnnames:
            print(f'  col')
        cont = input('(y/n) ')
        if cont!='y':

            print('\n In a Wide-format, REMOVE THE COLUMNS that do not represent differnt satations, before proceding! \n')

        print('What observation type does you data represent : ')
        desc_return = col_option_input(obstype_options)
        if desc_return is None:
            print('This is not an option, select an observation type.')
            sys.exit('invalid obstype for wide dataset, see last message. ')
        wide_obstype = inv_obstype_desc[desc_return]



    # =============================================================================
    # Map metadatafile
    # =============================================================================

    print('\n \n *******      Meta Data File   ***********')
    metatemplate_dict = {}
    meta_avail = input('The metadata contains metadata for each station, where each column represent a metadata type. \n Do you have a metadatafile? (y/n) ')
    if meta_avail == 'y':
        metadatafilepath = input('Give the full path to your metadata file : ')
        print(' ... opening the metadata file ...')
        metadata = _read_csv_file(metadatafilepath)
        metacolumnnames = metadata.columns.to_list()


        meta_desc = {
            'name': 'name (the column with the stationnames, must be unique for each station)',
            'lat':	'lat (the latitudes of the stations as a numeric values)',
            'lon':	'lon (The longtitudes of the stations as a numeric values)',
            'location': 'location (the city/region of the stations) (OPTIONAL)',
            'call_name': 'call_name (an informal name of the stations) (OPTIONAL)',
            'network': 'network (the name of the network the stations belong to) (OPTIONAL)',
            }
        inv_meta_desc = {val: key for key, val in meta_desc.items()}

        print('What do the following columns represent: \n')
        meta_options=list(meta_desc.values())
        for col in metacolumnnames:

            contin = input(f'  add {col} to the template? (y/n) ')
            if contin != 'y':
                continue

            print(f'\n {col} : ')
            desc_return =col_option_input(meta_options)
            if desc_return is None:
                continue #when enter x
            metatype = inv_meta_desc[desc_return]

            # check if the name column is equalt in the data template to avoid creating
            # two templates
            if metatype == 'name':
                if 'name' in template_dict:
                    if not col == template_dict['name']['orig_name']:
                        print(f'WARNING, the "name" column in the datafile is different than in the metadatafile! \
    Rename in your metadatafile : {col} ---> {template_dict["name"]["orig_name"]}')
                        cont = input ('Renaming done? (y/n) ')
                        if cont != 'y':
                            sys.exit(f'Please rename {col} ---> {template_dict["name"]["orig_name"]} in your metadata file.')

            metatemplate_dict[metatype] = {'orig_name': col}
            meta_options.remove(meta_desc[metatype])


    print('\n   ... Oke, that is all the info for the mapping. Now i will do some basic tests to see if the mapping works.')


    # =============================================================================
    # Apply tests
    # =============================================================================


    #  ------- tests on data ---------
    # apply tests the first row
    data_test = data.iloc[0].to_dict()

    # test if a stationname column is available in a long format
    print (' *  ... checking data columns ... ')
    if ((format_option == 1) & (not 'name' in template_dict)):
        print(' \n WARNING: There is no information which column in the data file represents the names of the stations. The toolkit will assume that the observations are from ONE station! \n')
        format_option = 3

    # check if a least one mapped observation type exist
    if (format_option != 2):
        present_obs = [key for key in template_dict.keys() if key in observation_types]
        if not bool(present_obs):
            print('ERROR! There is no observation type included in the template! Add at least one observation type when mapping the data file.')
            sys.exit('Template invalid, see last message. ')





    # test datetime format
    print (' *  ... checking timestamps formats ... ')
    if 'datetime' in template_dict:
        escape = False
        while not escape:
            test_dt = data_test[template_dict['datetime']['orig_name']]
            try:
                _ = datetime.strptime(test_dt,
                                      template_dict['datetime']['format'])
                print ('   ... testing datetime format is ...  OK!')
                escape=True
            except:
                print(f'ERROR: the {template_dict["datetime"]["format"]} does not work for {test_dt}')
                template_dict['datetime']['format'] = input('\n Try new timestamp format (ex. %Y-%m-%d %H:%M:%S) : ')

    if '_date' in template_dict:
        escape = False
        while not escape:
            test_dt = data_test[template_dict['_date']['orig_name']]
            try:
                _ = datetime.strptime(test_dt,
                                      template_dict['_date']['format'])
                print ('   ... testing date format is OK!')
                escape=True
            except:
                print(f'ERROR: the {template_dict["_date"]["format"]} does not work for {test_dt}')
                template_dict['_date']['format'] = input('\n Try new date format (ex. %Y-%m-%d) : ')
    if '_time' in template_dict:
        escape = False
        while not escape:
            test_dt = data_test[template_dict['_time']['orig_name']]
            try:
                _ = datetime.strptime(test_dt,
                                      template_dict['_time']['format'])
                print ('   ... testing time format is OK!')
                escape=True
            except:
                print(f'ERROR: the {template_dict["_time"]["format"]} does not work for {test_dt}')
                template_dict['_time']['format'] = input('\n Try new time format (ex. %H:%M:%S) : ')


    # check if all data columns are mapped
    print (' *  ... checking for unmapped data columns ... ')
    if (format_option == 1) | (format_option == 3):
        present_columns = list(data_test.keys())
        mapped_cols = [val['orig_name'] for val in template_dict.values()]
        for col in present_columns:
            if not col in mapped_cols:
                print(f' Warning! {col} in the datafile is not present in the template, and thus it will not be used.')



    # -------- tests on metadata ----------
    if bool(metatemplate_dict):
        # apply tests the first row
        metadata_test = metadata.iloc[0].to_dict()


        # test if name is in the metadat
        print (' *  ... checking metadata columns ... ')
        if ((not 'name' in metatemplate_dict)):
            print(f'Error! There is no metadata column containing the station names in the template! Add this column to the metadatafile of the template.')
            sys.exit('Template invalid, see last message. ')



        # test if all stationnames are present in the metadata
        print (' *  ... checking compatible station names ... ')
        if ((format_option == 1) & ('name' in template_dict)):
            stanames_data = data[template_dict['name']['orig_name']].unique()
            stanames_metadata = metadata[metatemplate_dict['name']['orig_name']].unique()

            unmapped = [sta for sta in stanames_data if not sta in stanames_metadata]
            if bool(unmapped):
                print(f'Warning! The following stations are found in the data, but not in the metadata: {unmapped}')


        if ((format_option == 2)):
            stanames_data = data[template_dict['name']['orig_name']].unique().to_list()
            stanames_metadata = metadata[metatemplate_dict['name']['orig_name']].unique().to_list()

            unmapped = [sta for sta in stanames_data if not sta in stanames_metadata]
            print(f'Warning! The following stations are found in the data, but not in the metadata: {unmapped}')






        # check if all metadata columns are mapped
        print (' *  ... checking for unmapped metadata columns ... ')
        present_columns = list(metadata_test.keys())
        mapped_cols = [val['orig_name'] for val in metatemplate_dict.values()]
        for col in present_columns:
            if not col in mapped_cols:
                print(f' Warning! {col} in the metadatafile is not present in the template, and thus it will not be used.')




    # =============================================================================
    # Saving the template
    # =============================================================================

    print('\n ------ Saving the template ----- \n')
    save_dir = input('Give a directory where to save the template (as template.csv) : ')

    template_dict.update(metatemplate_dict) #this is why name in data and metadata should have the same mapping !!


    # Convert to dataframe


    df = pd.DataFrame(template_dict).transpose()
    df.index.name = 'varname'
    df = df.rename(columns={'orig_name' : 'template column name'})
    df = df.reset_index()

    # write to csv
    templatefilepath = os.path.join(save_dir, 'template.csv')
    df.to_csv(templatefilepath, na_rep = '', index=False)
    print(f' DONE! The template is writen here: {templatefilepath}')



    # =============================================================================
    # Tips for the user
    # =============================================================================
    apply_tips = input("Do you want some help creating your Dataset? (y/n) : ")
    if apply_tips == 'y':

        print('\n ------ How to use the template ----- ')

        print('(Some questions will be asked that are case-specific) \n')
        tzchange = input(' Are the timestamps in UTC? (y/n) : ')
        tz_update = False
        if tzchange != 'y':
            print('Select a timezone: ')
            tzstring = col_option_input(pytz.all_timezones)
            tz_update = True

        output_change = input(' Do you plan to save images to a direcory? (y/n) : ')
        output_update = False
        if output_change == 'y':
            output_folder = input(' Give the path of your output direcory : ')
            output_update = True

        staname_update = False
        if ((format_option == 3) & (not 'name' in template_dict)):#single station with no name information
            staname = input(' What is the name of your station : ')
            staname_update=True

        gaps_change = input(' Do you want to use the default gaps defenition? (y/n) : ')
        gaps_update = False
        if gaps_change != 'y':
            gapsize = int(input(' What is the minimum number of consecutive missing records to define as a gap? (default=40) : '))
            gaps_update = True



        print('\n#1. Define the paths to your files: \n')
        print(f'data_file = "{datafilepath}"')
        if bool(metatemplate_dict):
            print(f'meta_data_file = "{metadatafilepath}"')

        print(f'data_template = "{templatefilepath}"')
        if bool(metatemplate_dict):
            print(f'meta_data_template = "{templatefilepath}"')

        print('\n#2. initiate a dataset: \n')
        print('your_dataset = metobs_toolkit.Dataset()')

        print('\n#3. Update the paths to your files: \n')
        print('your_dataset.update_settings(')
        print('    input_data_file = data_file,')
        print('    data_template_file = data_template,')
        if bool(metatemplate_dict):
            print('    input_metadata_file = meta_data_file,')
            print('    metadata_template_file = meta_data_template,')
        if output_update:
            print(f'    output_folder = "{output_folder}",')
        print('    )')

        # extra case specific options
        if ((tz_update) | (staname_update) | (gaps_update)):
            print('\n#3B. Update specific settings (optional): \n')

        if tz_update:
            print(f'your_dataset.update_timezone(timezonestr = "{tzstring}")')
        if staname_update:
            print(f'your_dataset.update_default_name(default_name = "{staname}")')
        if gaps_update:
            print(f'your_dataset.update_qc_settings(gapsize_in_records = {gapsize})')


        print('\n#4. Import your data : \n')

        print('your_dataset.import_data_from_file(')
        if format_option == 2:
            print('    long_format=False,')
            print(f'    obstype={wide_obstype},')
        else:
            print('    long_format=True,')

        print('    freq_estimation_method=None, #None for default(="highest"), "highest" or "median"')
        print('    freq_estimation_simplify=None, #None for default(=True), True or False')
        print('    freq_estimation_simplify_error=None, #None for default(="2T"), or timedelta string.')

        print('    )')





    return df