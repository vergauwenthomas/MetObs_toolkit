#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:25:45 2023

@author: thoverga
"""

import os
import sys
import pandas as pd
import numpy as np
import copy
from datetime import datetime
import pytz


from metobs_toolkit.template import _get_empty_templ_dict, _pwrite_templdict_to_json


from metobs_toolkit.obstypes import Obstype, tlk_obstypes
from metobs_toolkit.data_import import _read_csv_to_df


def add_new_obstype():

    print("\n --- Adding a new observation type --- \n")

    # get obsname
    name_ok = False
    while not name_ok:
        obsname = str(input("Give the name of your observation type: "))
        if obsname in tlk_obstypes.keys():
            print(
                f"!! {obsname} is already a knonw observation type. This cannot be added."
            )
        else:
            name_ok = True

    # get std unit
    std_unit = str(
        input(
            "Give the standard unit (how the toolkit should store/present the data): "
        )
    )

    # Get input data unit
    is_std_unit = yes_no_ques(f" Are the {obsname} values in your data in {std_unit}")
    if is_std_unit:
        cur_unit = std_unit
        unit_conv = {std_unit: ["x"]}
    else:
        cur_unit = str(input("Give the unit your data is in: "))
        print(
            f"Give the expression on how to convert {cur_unit} values to {std_unit}. "
        )
        print("  * Example: Kelvin (= new unit) to °C :  x - 273.15 ")
        print(
            "  * Example: Farenheit to °C : x-32.0; x/1.8    (executed left to right)"
        )

        conv_str = str(input(" : "))
        # cleanup and make list if needend
        conv_str = list(conv_str.replace(" ", "").split(";"))

        unit_conv = {cur_unit: conv_str}
    # Description
    description = str(
        input(f"Give a detailed description of the {obsname} type (optional): ")
    )

    # Aliases and coversions

    # Do not add this in the prompt, the prompt should not check the more advanced
    # settigns. If the prompt could cover 95% of all user needs, that would be great.
    # The others should help themself with the documentation to create aliases
    # and conversions

    unit_aliases = {}

    # create obstype:
    new_obstype = Obstype(
        obsname=obsname,
        std_unit=std_unit,
        description=description,
        unit_aliases=unit_aliases,
        unit_conversions=unit_conv,
    )
    return new_obstype, cur_unit


def get_unit(obstype):

    available_units = obstype.get_all_units()
    available_units.append("ADD A NEW UNIT")

    print(f"\n Select the unit your {obstype.name} data is in:  \n")
    conv_str = None
    unit = col_option_input(available_units)
    if unit == "ADD A NEW UNIT":
        unit = str(input("Give the unit your data is in: "))
        print(
            f"Give the expression on how to convert {unit} values to {obstype.get_standard_unit()}. "
        )
        print("  * Example: Kelvin (= new unit) to °C :  x - 273.15 ")
        print(
            "  * Example: Farenheit to °C : x-32.0; x/1.8    (executed left to right)"
        )

        conv_str = str(input(" : "))
        # cleanup and make list if needend
        conv_str = list(conv_str.replace(" ", "").split(";"))

    return unit, conv_str


def col_option_input(columns):
    """Convert options to numerics and ask for input."""
    mapper = {}
    i = 1
    for col in columns:
        if col == "Unnamed: 0":
            print(f"  {i}. {col} (--> this is the index of your csv file)")
        else:
            print(f"  {i}. {col}")
        mapper[i] = col
        i += 1

    print("  x. -- not valid --")
    valid_input = False
    while valid_input is False:
        if i <= 3:
            repr_str = "("
            for i in np.arange(1, i):
                repr_str += str(i) + ", "
            # remove last comma
            repr_str = repr_str[:-2] + ") : "
            num_ans = input(f"{repr_str}")
        else:
            num_ans = input(f"(1 - {i-1}) : ")

        if num_ans == "x":
            print(" ... This setting is not provided! ...")
            return None

        try:
            _ = mapper[int(num_ans)]
            valid_input = True
        except KeyError:
            valid_input = False
            print(f"{num_ans} is not a valid input.")

    print(f" ... {mapper[int(num_ans)]} selected ... \n")
    return mapper[int(num_ans)]


def yes_no_ques(text):
    """Get yes/no input."""
    valid_input = False

    while valid_input is False:
        prompt = input(f" {text}. (y/n) : ")

        if (prompt == "y") | (prompt == "Y"):
            valid_input = True
            return True
        elif (prompt == "n") | (prompt == "N"):
            valid_input = True
            return False
        else:
            print(f" {prompt} is not y or n, give a suitable answer.")


def usr_input_dir(text):
    """Prompt directory path.

    question and check if the answer is a directory, return the path
    if it is a directory, repeat else.
    """
    is_dir = False
    while is_dir is False:
        inp_dir = input(f"{text} : ")
        if os.path.isdir(inp_dir):
            is_dir = True
        else:
            print(f"{inp_dir} is not a directory, try again.")
    return inp_dir


def usr_input_file(text):
    """Prompt file path.

    Prompt question and check if the answer is a file, return the path if it
    exists, repeat else.
    """
    is_file = False
    while is_file is False:
        inp_file = input(f"{text} : ")
        if os.path.isfile(inp_file):
            is_file = True
        else:
            print(f"{inp_file} is not found, try again.")
    return inp_file


def build_template_prompt(debug=False):
    """Launch the prompt to help make a template."""
    tmpl_dict = _get_empty_templ_dict()
    tmpl_dict["data_related"]["obstype_mapping"] = []

    known_obstypes = copy.copy(tlk_obstypes)
    new_units = {}

    print(
        "This prompt will help to build a template for your data and metadata. Answer the prompt and hit Enter. \n \n"
    )

    print(" *******      File locations   *********** \n")

    datafilepath = usr_input_file("Give the full path to your data file")
    meta_avail = yes_no_ques("Do you have a file with the metadata?")
    if meta_avail:
        metadatafilepath = usr_input_file("Give the full path to your metadata file")

    # =============================================================================
    # Map data file
    # =============================================================================

    print("\n\n *******      Data File   ***********")

    # datafilepath = usr_input_file('Give the full path to your data file')
    print(" ... opening the data file ...")
    data = _read_csv_to_df(datafilepath, {"nrows": 10})
    columnnames = data.columns.to_list()

    format_dict = {
        "Long format (station observations are stacked as rows)": 1,
        "Wide format (columns represent different stations)": 2,
        "Single station format (columns represent observation(s) of one station)": 3,
    }

    print("How is your dataset structured : \n")
    format_option = col_option_input(format_dict.keys())
    print(f" \n... oke, {format_option} selected ...\n")
    format_option = format_dict[format_option]
    if debug:
        print(f"format numeric option: {format_option}")
    if format_option == 1:
        tmpl_dict["data_related"]["structure"] = "long"
        # options_dict["data_structure"] = "long"
    if format_option == 2:
        tmpl_dict["data_related"]["structure"] = "wide"
        # options_dict["data_structure"] = "wide"
    if format_option == 3:
        tmpl_dict["data_related"]["structure"] = "single_station"
        # options_dict["data_structure"] = "single_station"

    # Datatime mapping
    dt_dict = {
        "In a single column (ex: 2023/06/07 16:12:30)": 1,
        "By a column with dates, and another column with times": 2,
    }
    print("How are the timestamps present in your data file : \n")
    datetime_option = col_option_input(dt_dict.keys())
    datetime_option = dt_dict[datetime_option]

    if datetime_option == 1:

        # Datetime mapping

        print("\n Which is your timestamp columnname: ")

        datetimecolumn = col_option_input(columnnames)
        tmpl_dict["data_related"]["timestamp"]["datetime_column"] = datetimecolumn
        columnnames.remove(datetimecolumn)

        example = data[datetimecolumn].iloc[0]
        tmpl_dict["data_related"]["timestamp"]["datetime_fmt"] = input(
            f"Type your datetime format (ex. %Y-%m-%d %H:%M:%S), (your first timestamp: {example}) : "
        )

    else:
        # Date mapping
        print("Which column represents the DATES : ")
        datecolumn = col_option_input(columnnames)
        tmpl_dict["data_related"]["timestamp"]["date_column"] = datecolumn
        columnnames.remove(datecolumn)

        example = data[datecolumn].iloc[0]
        tmpl_dict["data_related"]["timestamp"]["date_fmt"] = input(
            f"Type your date format (ex. %Y-%m-%d), (your first timestamp: {example}) : "
        )

        print(" \n")

        # Time mapping

        print("Which column represents the TIMES : ")
        timecolumn = col_option_input(columnnames)
        tmpl_dict["data_related"]["timestamp"]["time_column"] = timecolumn
        columnnames.remove(timecolumn)

        example = data[timecolumn].iloc[0]
        tmpl_dict["data_related"]["timestamp"]["time_fmt"] = input(
            f"Type your time format (ex. %H:%M:%S), (your first timestamp: {example}) : "
        )

    # Obstype mapping in long format:
    obstype_desc = {"name": "name (name of the stations represented by strings)"}
    obstype_desc.update(
        {
            ob.name: f"{ob.name} : {ob.get_description()}"
            for ob in known_obstypes.values()
        }
    )
    obstype_desc.update(
        {
            "ADD NEW OBSERVATION TYPE": "add a new observation type if it is not present in this list."
        }
    )

    inv_obstype_desc = {val: key for key, val in obstype_desc.items()}
    obstype_options = list(obstype_desc.values())

    if (format_option == 1) | (format_option == 3):
        # long format
        print("What do the following columns represent: \n")

        for col in columnnames:
            if col == "Unnamed: 0":
                contin = yes_no_ques(
                    f"\n add column {col} (: probably this is the index of the csv file) to the template?"
                )
            else:
                contin = yes_no_ques(f"\n add column {col} to the template?")

            if contin is False:
                continue

            print(f"\n {col} : ")

            desc_return = col_option_input(obstype_options)
            if desc_return is None:
                continue  # when enter x

            # 1) add a new obstype
            if inv_obstype_desc[desc_return] == "ADD NEW OBSERVATION TYPE":
                new_obstype, cur_unit = add_new_obstype()

                known_obstypes[new_obstype.name] = new_obstype  # add to knonw obstypes
                obstype = new_obstype.name
                units = cur_unit
                description = new_obstype.get_description()

            # 2) name column is mapped
            elif inv_obstype_desc[desc_return] == "name":
                tmpl_dict["data_related"]["name_column"] = col
                obstype_options.remove(
                    "name (name of the stations represented by strings)"
                )
                continue

            # 3) existing obstype
            else:
                obstype = inv_obstype_desc[desc_return]

                # add unit
                units, conv_str = get_unit(known_obstypes[obstype])
                if conv_str is not None:
                    # add new units to the dict
                    new_units[obstype] = {"unit": units, "conv": conv_str}

                description = input("Some more details on the observation (optional): ")
                obstype_options.remove(obstype_desc[obstype])

            # update template

            obsdict = {
                "tlk_obstype": obstype,
                "columnname": col,
                "unit": str(units),
                "description": str(description),
            }

            tmpl_dict["data_related"]["obstype_mapping"].append(obsdict)

    if format_option == 2:
        print("\n Does these columns represent stations: ")
        for col in columnnames:
            print(f"  {col} ")

        cont = yes_no_ques("")
        if cont is False:
            print(
                "\n In a Wide-format, REMOVE THE COLUMNS that do not represent different satations, before proceding! \n"
            )
        else:
            stationnames = columnnames

        print("\n What observation type does you data represent : ")
        obstype_options.remove(obstype_desc["name"])
        desc_return = col_option_input(obstype_options)
        if desc_return is None:
            print("This is not an option, select an observation type.")
            sys.exit("invalid obstype for wide dataset, see last message. ")
        wide_obstype = inv_obstype_desc[desc_return]

        # 1) add a new obstype
        if wide_obstype == "ADD NEW OBSERVATION TYPE":
            new_obstype, cur_unit = add_new_obstype()
            wide_obstype = new_obstype.name
            known_obstypes[new_obstype.name] = new_obstype  # add to knonw obstypes
            units = cur_unit
            description = new_obstype.get_description()

        # 2) Knonw obstype
        else:
            # add unit
            units, conv_str = get_unit(known_obstypes[wide_obstype])
            if conv_str is not None:
                # add new units to the dict
                new_units[wide_obstype] = {"unit": units, "conv": conv_str}

            description = input("Some more details on the observation (optional): ")

        # update template

        obsdict = {
            "tlk_obstype": wide_obstype,
            "columnname": None,
            "unit": str(units),
            "description": str(description),
        }

        tmpl_dict["data_related"]["obstype_mapping"].append(obsdict)

    if debug:
        print(f"template_dict: {tmpl_dict}")

    # =============================================================================
    # Map metadatafile
    # =============================================================================

    print("\n \n *******      Meta Data   ***********")

    metatemplate_dict = {}

    if meta_avail:
        print(" ... opening the metadata file ...")
        metadata = _read_csv_to_df(metadatafilepath, {"nrows": 10})
        metacolumnnames = metadata.columns.to_list()

        # map the required columns (name)

        # if multiple stations are in the dataset, this column is required
        if format_option != 3:
            print("Which column does represent the NAMES of the stations?")
            name_column = col_option_input(metacolumnnames)
            tmpl_dict["metadata_related"]["name_column"] = name_column
            metacolumnnames.remove(name_column)

        # if the data is a single station, this column is ignored
        else:
            staname = input("\n What is the name of your station : ")
            tmpl_dict["single_station_name"] = staname

        # map columns that are used by the toolit (lat, lon)
        with_coords = yes_no_ques(
            "\n are there coordinates (latitude, longitude) columns in the metadata?"
        )
        if with_coords:
            print("Which column does represent the LATITUDES?")
            lat_column = col_option_input(metacolumnnames)
            tmpl_dict["metadata_related"]["lat_column"] = lat_column
            metacolumnnames.remove(lat_column)

            print("Which column does represent the LONGITUDES?")
            lon_column = col_option_input(metacolumnnames)
            tmpl_dict["metadata_related"]["lon_column"] = lon_column
            metacolumnnames.remove(lon_column)

        # Which other (not used by the toolkit) to add.
        if len(metacolumnnames) > 0:
            add_cols = yes_no_ques(
                f"\n Do you want to include/use remaining columns in the metadatafile? \n ({str(metacolumnnames)})"
            )
            if add_cols:
                for col in metacolumnnames:
                    add_bool = yes_no_ques(f"\n Add {col} in the metada?")
                    if add_bool:
                        tmpl_dict["metadata_related"]["columns_to_include"].append(
                            str(col)
                        )

    #     meta_desc = {
    #         "name": "name (the column with the stationnames, must be unique for each station)",
    #         "lat": "lat (the latitudes of the stations as a numeric values)",
    #         "lon": "lon (The longtitudes of the stations as a numeric values)",
    #         "location": "location (the city/region of the stations) (OPTIONAL)",
    #         "call_name": "call_name (an informal name of the stations) (OPTIONAL)",
    #         "network": "network (the name of the network the stations belong to) (OPTIONAL)",
    #     }
    #     inv_meta_desc = {val: key for key, val in meta_desc.items()}

    #     print("What do the following columns represent: \n")
    #     meta_options = list(meta_desc.values())
    #     for col in metacolumnnames:
    #         if col == 'Unnamed: 0':
    #             contin = yes_no_ques(f"\n add column {col} (: probably this is the index of the csv file) to the template?")
    #         else:
    #             contin = yes_no_ques(f"add {col} to the template?")
    #         if contin is False:
    #             continue

    #         print(f"\n {col} : ")
    #         desc_return = col_option_input(meta_options)
    #         if desc_return is None:
    #             continue  # when enter x
    #         metatype = inv_meta_desc[desc_return]

    #         # check if the name column is equalt in the data template to avoid creating
    #         # two templates
    #         if metatype == "name":
    #             if "name" in template_dict:
    #                 if not col == template_dict["name"]["orig_name"]:
    #                     print(
    #                         f'WARNING, the "name" column in the datafile is different than in the metadatafile! \
    # Rename in your metadatafile : {col} ---> {template_dict["name"]["orig_name"]}'
    #                     )
    #                     cont = yes_no_ques("Renaming done?")
    #                     if cont is False:
    #                         sys.exit(
    #                             f'Please rename {col} ---> {template_dict["name"]["orig_name"]} in your metadata file.'
    #                         )

    #         metatemplate_dict[metatype] = {"orig_name": col}
    #         meta_options.remove(meta_desc[metatype])
    # if debug:
    #     print(f"metatemplate_dict : {metatemplate_dict}")

    # =============================================================================
    # Apply tests
    # =============================================================================
    # print("\n \n *******      Testing template compatibility   ***********")
    # print(
    #     "\n   ... Oke, that is all the info for the mapping. Now i will do some basic tests to see if the mapping works."
    # )

    # #  ------- tests on data ---------
    # # apply tests the first row
    # data_test = data.iloc[0].to_dict()

    # # test if a stationname column is available in a long format
    # print(" *  ... checking data columns ... ")
    # if (format_option == 1) & (not "name" in template_dict):
    #     print(
    #         " \n WARNING: There is no information which column in the data file represents the names of the stations. The toolkit will assume that the observations are from ONE station! \n"
    #     )
    #     format_option = 3

    # # check if a least one mapped observation type exist
    # if format_option != 2:
    #     present_obs = [
    #         key for key in template_dict.keys() if key in known_obstypes.keys()
    #     ]
    #     if not bool(present_obs):
    #         print(
    #             "ERROR! There is no observation type included in the template! Add at least one observation type when mapping the data file."
    #         )
    #         sys.exit("Template invalid, see last message. ")

    # # test datetime format
    # print(" *  ... checking timestamps formats ... ")
    # if "datetime" in template_dict:
    #     escape = False
    #     while not escape:
    #         test_dt = data_test[template_dict["datetime"]["orig_name"]]
    #         try:
    #             _ = datetime.strptime(test_dt, template_dict["datetime"]["format"])
    #             print("   ... testing datetime format is ...  OK!")
    #             escape = True
    #         except:
    #             print(
    #                 f'ERROR: the {template_dict["datetime"]["format"]} does not work for {test_dt}'
    #             )
    #             template_dict["datetime"]["format"] = input(
    #                 "\n Try new timestamp format (ex. %Y-%m-%d %H:%M:%S) : "
    #             )

    # if "_date" in template_dict:
    #     escape = False
    #     while not escape:
    #         test_dt = data_test[template_dict["_date"]["orig_name"]]
    #         try:
    #             _ = datetime.strptime(test_dt, template_dict["_date"]["format"])
    #             print("   ... testing date format is OK!")
    #             escape = True
    #         except:
    #             print(
    #                 f'ERROR: the {template_dict["_date"]["format"]} does not work for {test_dt}'
    #             )
    #             template_dict["_date"]["format"] = input(
    #                 "\n Try new date format (ex. %Y-%m-%d) : "
    #             )
    # if "_time" in template_dict:
    #     escape = False
    #     while not escape:
    #         test_dt = data_test[template_dict["_time"]["orig_name"]]
    #         try:
    #             _ = datetime.strptime(test_dt, template_dict["_time"]["format"])
    #             print("   ... testing time format is OK!")
    #             escape = True
    #         except:
    #             print(
    #                 f'ERROR: the {template_dict["_time"]["format"]} does not work for {test_dt}'
    #             )
    #             template_dict["_time"]["format"] = input(
    #                 "\n Try new time format (ex. %H:%M:%S) : "
    #             )

    # # check if all data columns are mapped
    # print(" *  ... checking for unmapped data columns ... ")
    # if (format_option == 1) | (format_option == 3):
    #     present_columns = list(data_test.keys())
    #     mapped_cols = [val["orig_name"] for val in template_dict.values()]
    #     for col in present_columns:
    #         if not col in mapped_cols:
    #             print(
    #                 f" Warning! {col} in the datafile is not present in the template, and thus it will not be used."
    #             )

    # # -------- tests on metadata ----------
    # if bool(metatemplate_dict):
    #     # apply tests the first row
    #     metadata_test = metadata.iloc[0].to_dict()

    #     # test if name is in the metadat in a long format
    #     print(" *  ... checking metadata columns ... ")
    #     if (not "name" in metatemplate_dict) & ((format_option in [1, 2])):
    #         print(
    #             f"Error! There is no metadata column containing the station names in the template! Add this column to the metadatafile of the template."
    #         )
    #         sys.exit("Template invalid, see last message. ")

    #     print(" *  ... checking metadata name duplicates... ")
    #     if format_option in [1, 2]:
    #         stanames_metadata = metadata[metatemplate_dict["name"]["orig_name"]]
    #         if stanames_metadata.duplicated().any():
    #             dubs = stanames_metadata[stanames_metadata.duplicated()]
    #             print(
    #                 f"Error! There are duplicated names present in the metadatafile {dubs}. Remove the duplicates manually."
    #             )
    #             sys.exit("Template invalid, see last message. ")

    #     # test if all stationnames are present in the metadata
    #     print(" *  ... checking compatible station names ... ")
    #     if (format_option == 1) & ("name" in template_dict):
    #         stanames_data = data[template_dict["name"]["orig_name"]].unique()
    #         stanames_metadata = metadata[
    #             metatemplate_dict["name"]["orig_name"]
    #         ].unique()

    #         unmapped = [sta for sta in stanames_data if not sta in stanames_metadata]
    #         if bool(unmapped):
    #             print(
    #                 f"Warning! The following stations are found in the data, but not in the metadata: {unmapped}"
    #             )

    #     if format_option == 2:
    #         # 1. no duplicates in stationnames
    #         if not (len(stationnames) == len(set(stationnames))):
    #             print(
    #                 f"Error! Duplicated station names found in the columns of the dataset: {stationnames}"
    #             )
    #             sys.exit("Template invalid, see last message. ")

    #         # 2. check if all stationname in the data are defined in the metadata,
    #         # If there are no mapped stationnames give error, else give warning

    #         stanames_metadata = metadata[
    #             metatemplate_dict["name"]["orig_name"]
    #         ].to_list()
    #         unmapped = [
    #             staname for staname in stationnames if not staname in stanames_metadata
    #         ]

    #         if len(unmapped) == len(stationnames):
    #             print(
    #                 f"Error! None of the stationnames in the dataset ({stationnames}), are found in the metadataset ({stanames_metadata})."
    #             )
    #             sys.exit("Template invalid, see last message. ")

    #         if len(unmapped) < len(stationnames):
    #             print(f" unmapped: {unmapped}")
    #             print(f" stationnames: {stationnames}")
    #             print(f" stationnames metadta: {stanames_metadata}")
    #             print(
    #                 f"Warning! The following stations are present in the data but not in the metadata: {unmapped}"
    #             )

    #     # check if all metadata columns are mapped
    #     print(" *  ... checking for unmapped metadata columns ... ")
    #     present_columns = list(metadata_test.keys())
    #     mapped_cols = [val["orig_name"] for val in metatemplate_dict.values()]
    #     for col in present_columns:
    #         if not col in mapped_cols:
    #             print(
    #                 f" Warning! {col} in the metadatafile is not present in the template, and thus it will not be used."
    #             )

    # # make sure the stationname is unique in single station datafile
    # if format_option == 3:
    #     print(" *  ... checking if stationname is unique ... ")
    #     if bool(metatemplate_dict):
    #         if "name" in metatemplate_dict:
    #             names = metadata[metatemplate_dict["name"]["orig_name"]].unique()
    #             if len(names) > 1:
    #                 print(
    #                     f"Error! multiple station names found in the {metatemplate_dict['name']['orig_name']} metadata column."
    #                 )
    #                 sys.exit("Template invalid, see last message. ")
    #     else:
    #         if "name" in template_dict:
    #             names = data[template_dict["name"]["orig_name"]].unique()
    #             if len(names) > 1:
    #                 print(
    #                     f"Error! multiple station names found in the {template_dict['name']['orig_name']} data column."
    #                 )
    #                 sys.exit("Template invalid, see last message. ")

    # =============================================================================
    #     Some extra options
    # =============================================================================

    # template_dict.update(
    #     metatemplate_dict
    # )  # this is why name in data and metadata should have the same mapping !!

    print("\n \n *******      Extra options    ***********")
    # if (tmpl_dict["data_related"]["structure"] == "single_station") & (
    #     tmpl_dict["metadata_related"]["name_column"] is None
    # ):

    #     # single station with no name information
    #     staname = input("\n What is the name of your station : ")
    #     tmpl_dict["single_station_name"] = staname

    tzchange = yes_no_ques("\n Are the timestamps in UTC?")
    if tzchange is False:
        print("\n Select a timezone: ")
        tzstring = col_option_input(pytz.all_timezones)
        tmpl_dict["data_related"]["timestamp"]["timezone"] = tzstring
    else:
        tmpl_dict["data_related"]["timestamp"]["timezone"] = "UTC"

    # =============================================================================
    # Saving the template
    # =============================================================================

    print("\n ------ Saving the template ----- \n")
    save_dir = usr_input_dir(
        "Give a directory where to save the template (as template.json)"
    )

    # write to csv
    templatefilepath = os.path.join(save_dir, "template.json")
    _pwrite_templdict_to_json(templdict=tmpl_dict, trgfile=templatefilepath)

    print(f" DONE! The template is written here: {templatefilepath}")

    # =============================================================================
    # Tips for the user
    # =============================================================================

    apply_tips = yes_no_ques("Do you want some help creating your Dataset?")
    if apply_tips is True:

        print("\n ------ How to use the template ----- ")

        print("(Some questions will be asked that are case-specific) \n")

        output_change = yes_no_ques("Do you plan to save images to a direcory?")
        output_update = False
        if output_change is True:
            output_folder = input(" Give the path of your output direcory : ")
            output_update = True

        gaps_change = yes_no_ques("Do you want to use the default gaps defenition?")
        gaps_update = False
        if gaps_change is False:
            gapsize = int(
                input(
                    " What is the minimum number of consecutive missing records to define as a gap? (default=40) : "
                )
            )
            gaps_update = True

        print("\n\n ========= RUN THIS CODE ========= \n\n")

        print("\n#1. Define the paths to your files: \n")
        print(f'data_file = r"{datafilepath}"')
        if meta_avail:
            print(f'meta_data_file = r"{metadatafilepath}"')

        print(f'template = r"{templatefilepath}"')

        print("\n#2. initiate a dataset: \n")
        print("your_dataset = metobs_toolkit.Dataset()")

        print("\n#3. Update the paths to your files: \n")
        print("your_dataset.update_settings(")
        print("    input_data_file = data_file,")
        if meta_avail:
            print("    input_metadata_file = meta_data_file,")
        print("    template_file = template,")
        if output_update:
            print(f'    output_folder = "{output_folder}",')
        print("    )")

        # extra case specific options
        if gaps_update:
            print("\n#3B. Update specific settings (optional): \n")

        if gaps_update:
            print(f"your_dataset.update_qc_settings(gapsize_in_records = {gapsize})")

        # add new obstypes if needed
        to_add_obstypes = [
            newobsname
            for newobsname in known_obstypes.keys()
            if newobsname not in tlk_obstypes.keys()
        ]
        if bool(to_add_obstypes):
            print(
                "\n# Define non-standard observation types, and add them to the dataset: \n"
            )
            for newob in to_add_obstypes:
                new_obstype = known_obstypes[newob]
                print("new_obstype = metobs_toolkit.Obstype(")
                print(f'                 obsname="{new_obstype.name}",')
                print(f'                 std_unit="{new_obstype.get_standard_unit()}",')
                print(
                    f'                 description="{new_obstype.get_description()}",'
                )
                print(f"                 unit_aliases={new_obstype.units_aliases},")
                print(f"                 unit_conversions={new_obstype.conv_table})")
                print("\n\n #add the new obstype to your dataset. \n")
                print("your_dataset.add_new_observationtype(Obstype=new_obstype)")
                print("\n\n")

        # add new units if needed

        if bool(new_units):
            print(
                "\n# Define non-standard units, and add them to the corresponding units: \n"
            )
            for obstype, unit_info in new_units.items():
                print("your_dataset.add_new_unit(")
                print(f'                         obstype="{obstype}",')
                print(f'                         new_unit="{unit_info["unit"]}",')
                print(
                    f'                         conversion_expression={unit_info["conv"]})'
                )

                print("\n\n")

        print("\n#4. Import your data : \n")

        print("your_dataset.import_data_from_file()")

    return
