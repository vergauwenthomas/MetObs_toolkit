#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 24 10:25:45 2023

@author: thoverga
"""

import os
import sys
import copy
import logging

import pandas as pd
import numpy as np
import pytz

# from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.backend_collection.dev_collection import get_function_defaults
from metobs_toolkit.template import _get_empty_templ_dict
from metobs_toolkit.dataset import Dataset
from metobs_toolkit.obstypes import Obstype, tlk_obstypes, MetObsUnitUnknown
from metobs_toolkit.io_collection.filereaders import find_suitable_reader
from metobs_toolkit.io_collection.filewriters import write_dict_to_json

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def test_unit(unitstring: str) -> bool:
    """
    Test if a unit can be interpreted by pint.

    Parameters
    ----------
    unitstring : str
        The unit string to test.

    Returns
    -------
    bool
        True if the unit is valid, False otherwise.
    """
    try:
        Obstype(obsname="dummy", std_unit=unitstring, description="dummy")
        return True
    except MetObsUnitUnknown:
        print(
            f"!! {unitstring} is not a known Pint-unit. \
See https://github.com/hgrecco/pint/blob/master/pint/default_en.txt \
for a complete list.\n (Note that you can use prefixes (kilo, hecto, ... ) \
and expressions (meter/second, liter/second**2, ...)"
        )
        return False


@log_entry
def add_new_obstype() -> tuple:
    """
    Interactively add a new observation type.

    Returns
    -------
    tuple
        Tuple containing the new Obstype object and the current unit string.
    """
    print("\n --- Adding a new observation type --- \n")

    # get obsname
    name_ok = False
    while not name_ok:
        obsname = str(input("Give the name of your observation type: "))
        if obsname in tlk_obstypes.keys():
            print(
                f"!! {obsname} is already a known observation type. This \
cannot be added."
            )
        else:
            name_ok = True

    # get std unit
    print(
        "(Note: Units are handled by the pint-package. See a list of all \
available units here: \
https://github.com/hgrecco/pint/blob/master/pint/default_en.txt) "
    )
    stdunit_ok = False
    while not stdunit_ok:
        std_unit = str(
            input(
                "Give the standard unit (how the toolkit should \
store/present the data): "
            )
        )
        # test if the unit is valid
        stdunit_ok = test_unit(std_unit)

    # Get input data unit
    cur_unit = get_unit(
        trgobstype=Obstype(obsname=obsname, std_unit=std_unit, description="_dummy")
    )

    # Description
    description = str(
        input(f"Give a detailed description of the {obsname} type (optional): ")
    )
    # create obstype:
    new_obstype = Obstype(
        obsname=obsname,
        std_unit=std_unit,
        description=description,
    )
    return new_obstype, cur_unit


@log_entry
def get_unit(trgobstype: Obstype) -> str:
    """
    Prompt the user to specify the unit of the data for a given observation
    type.

    Parameters
    ----------
    trgobstype : Obstype
        The target observation type.

    Returns
    -------
    str
        The unit string for the data.
    """
    is_std_unit = yes_no_ques(
        f" Are the {trgobstype} values in your data in {trgobstype.std_unit}?"
    )
    if is_std_unit:
        cur_unit = trgobstype.std_unit
        return cur_unit

    else:
        compatible_units = trgobstype.get_compatible_units()
        curunit_ok = False
        while not curunit_ok:
            print(
                f"The following units are compatible with \
{trgobstype}: \n {compatible_units}"
            )
            cur_unit = str(
                input(
                    "Give the unit your data is in (you can add prefixes \
like kilo/hecto/etc if you need): "
                )
            )

            # test unit
            try:
                trgobstype.original_unit = cur_unit  # checks compatibility
                curunit_ok = True
            except MetObsUnitUnknown:
                print(
                    f"!! {cur_unit} is not compatible with {trgobstype}. \
Provide a compatible unit."
                )
                pass
        return cur_unit


@log_entry
def col_option_input(columns) -> str:
    """
    prompt a list of options.

    Parameters
    ----------
    columns : iterable
        Iterable of column names or options.

    Returns
    -------
    str
        The selected column name or option. Returns None if 'x' is selected.
    """
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
            for j in np.arange(1, i):
                repr_str += str(j) + ", "
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
        except (KeyError, ValueError):
            valid_input = False
            print(f"{num_ans} is not a valid input.")

    print(f" ... {mapper[int(num_ans)]} selected ... \n")
    return mapper[int(num_ans)]


@log_entry
def yes_no_ques(text: str) -> bool:
    """
    Get yes/no input from the user.

    Parameters
    ----------
    text : str
        The question to prompt.

    Returns
    -------
    bool
        True if 'y' is selected, False if 'n' is selected.
    """
    valid_input = False

    while valid_input is False:
        prompt = input(f" {text}. (y/n) : ")

        if (prompt == "y") or (prompt == "Y"):
            valid_input = True
            return True
        elif (prompt == "n") or (prompt == "N"):
            valid_input = True
            return False
        else:
            print(f" {prompt} is not y or n, give a suitable answer.")


@log_entry
def usr_input_dir(text: str) -> str:
    """
    Prompt for a directory path and check if it exists.

    Parameters
    ----------
    text : str
        The prompt text.

    Returns
    -------
    str
        The valid directory path.
    """
    is_dir = False
    while is_dir is False:
        inp_dir = input(f"{text} : ")
        if os.path.isdir(inp_dir):
            is_dir = True
        else:
            print(f"{inp_dir} is not a directory, try again.")
    return inp_dir


@log_entry
def usr_input_file(text: str) -> str:
    """
    Prompt for a file path and check if it exists.

    Parameters
    ----------
    text : str
        The prompt text.

    Returns
    -------
    str
        The valid file path.
    """
    is_file = False
    while is_file is False:
        inp_file = input(f"{text} : ")
        if os.path.isfile(inp_file):
            is_file = True
        else:
            print(f"{inp_file} is not found, try again.")
    return inp_file


@log_entry
def build_template_prompt() -> None:
    """
    Launch an interactive prompt to construct a template.json file.

    When called, an interactive prompt will start. Answer the questions, and
    hit Enter to continue. At the end of the prompt, you can specify a
    location where to save the template.json file.

    Returns
    -------
    None

    Notes
    -----
    It is good practice to rename the template.json file to specify the
    corresponding data file(s).

    At the end, the prompt asks if you need further assistance. If you do, the
    prompt will print out code that you can copy and run to create a
    `Dataset()`.

    Warning
    -------
    In previous versions (<=v0.2.1) the template file was a CSV. Thus you have
    to create the template again to be compatible with this version of the
      toolkit.
    """
    tmpl_dict = _get_empty_templ_dict()
    tmpl_dict["data_related"]["obstype_mapping"] = []

    known_obstypes = copy.copy(tlk_obstypes)

    print(
        "This prompt will help to build a template for your data and metadata. \
Answer the prompt and hit Enter. \n \n"
    )

    print(" *******      File locations   *********** \n")

    data_avail = yes_no_ques("Do you have a file with the OBSERVATIONAL DATA?")
    if data_avail:
        datafilepath = usr_input_file("Give the full path to your data file")

    meta_avail = yes_no_ques("Do you have a file with the METADATA?")
    if meta_avail:
        metadatafilepath = usr_input_file("Give the full path to your metadata file")

    if (not data_avail) and (not meta_avail):
        sys.exit("No template can be built without a data or metadata file.")

    # ==========================================================================
    # Map data file
    # ==========================================================================
    if data_avail:
        print("\n\n *******      Data File   ***********")

        print(" ... opening the data file ...")
        datareader = find_suitable_reader(filepath=datafilepath, is_url=False)

        data = datareader.read(nrows=10)
        columnnames = data.columns.to_list()

        format_dict = {
            "Long format (station observations are stacked as rows)": 1,
            "Wide format (columns represent different stations)": 2,
            "Single station format (columns represent \
observation(s) of one station)": 3,
        }

        print("How is your dataset structured : \n")
        format_option = col_option_input(format_dict.keys())
        print(f" \n... oke, {format_option} selected ...\n")
        format_option = format_dict[format_option]

        if format_option == 1:
            tmpl_dict["data_related"]["structure"] = "long"
        if format_option == 2:
            tmpl_dict["data_related"]["structure"] = "wide"
        if format_option == 3:
            tmpl_dict["data_related"]["structure"] = "single_station"

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

            fmt_is_ok = False
            while not fmt_is_ok:
                example = data[datetimecolumn].iloc[0]
                datetimefmt = input(
                    f"Type your datetime format (ex. %Y-%m-%d %H:%M:%S), \
(your first timestamp: {example}) : "
                )
                # Test datetime format
                try:
                    _ = pd.to_datetime(data[datetimecolumn], format=datetimefmt)
                    fmt_is_ok = True
                except ValueError:
                    print(
                        f" !! {datetimefmt} is not a suitable format for your \
{datetimecolumn}-column, check your data and input a suitable format."
                    )
                    pass

            tmpl_dict["data_related"]["timestamp"]["datetime_fmt"] = datetimefmt
        else:
            # Date mapping
            print("Which column represents the DATES : ")
            datecolumn = col_option_input(columnnames)
            tmpl_dict["data_related"]["timestamp"]["date_column"] = datecolumn
            columnnames.remove(datecolumn)

            fmt_is_ok = False
            while not fmt_is_ok:
                example = data[datecolumn].iloc[0]
                datefmt = input(
                    f"Type your date format (ex. %Y-%m-%d), (your first timestamp: {example}) : "
                )
                # test format
                try:
                    _ = pd.to_datetime(data[datecolumn], format=datefmt)
                    fmt_is_ok = True
                except ValueError:
                    print(
                        f" !! {datefmt} is not a suitable format for your {datecolumn}-column, check your data and input a suitable format."
                    )
                    pass
            tmpl_dict["data_related"]["timestamp"]["date_fmt"] = datefmt

            print(" \n")

            # Time mapping
            print("Which column represents the TIMES : ")
            timecolumn = col_option_input(columnnames)
            tmpl_dict["data_related"]["timestamp"]["time_column"] = timecolumn
            columnnames.remove(timecolumn)

            fmt_is_ok = False
            while not fmt_is_ok:
                example = data[timecolumn].iloc[0]
                timefmt = input(
                    f"Type your time format (ex. %H:%M:%S), (your first timestamp: {example}) : "
                )

                # test format
                try:
                    _ = pd.to_datetime(data[timecolumn], format=timefmt)
                    fmt_is_ok = True
                except ValueError:
                    print(
                        f" !! {timefmt} is not a suitable format for your {timecolumn}-column, check your data and input a suitable format."
                    )
                    pass

            tmpl_dict["data_related"]["timestamp"]["time_fmt"] = timefmt
        # Set the timezone
        tzchange = yes_no_ques("\n Are the timestamps in UTC?")
        if tzchange is False:
            print("\n Select a timezone: ")
            tzstring = col_option_input(pytz.all_timezones)
            tmpl_dict["data_related"]["timestamp"]["timezone"] = tzstring
        else:
            tmpl_dict["data_related"]["timestamp"]["timezone"] = "UTC"

        # Obstype mapping in long format:
        obstype_desc = {"name": "name (name of the stations represented by strings)"}
        obstype_desc.update(
            {ob.name: f"{ob.name} : {ob.description}" for ob in known_obstypes.values()}
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

                    known_obstypes[new_obstype.name] = (
                        new_obstype  # add to known obstypes
                    )
                    obstype = new_obstype.name
                    units = cur_unit
                    description = new_obstype.description

                # 2) name column is mapped
                elif inv_obstype_desc[desc_return] == "name":
                    tmpl_dict["data_related"]["name_column"] = col
                    obstype_options.remove(
                        "name (name of the stations represented by strings)"
                    )
                    continue

                # 3) existing obstype
                else:
                    knownobstype = known_obstypes[inv_obstype_desc[desc_return]]
                    obstype = knownobstype.name
                    # get unit
                    units = get_unit(knownobstype)

                    description = input(
                        "Some more details on the observation (optional): "
                    )
                    obstype_options.remove(obstype_desc[knownobstype.name])

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
                    "\n In a Wide-format, REMOVE THE COLUMNS that do not represent different stations, before proceeding! \n"
                )
            else:
                pass

            print("\n What observation type does you data represent : ")
            obstype_options.remove(obstype_desc["name"])
            desc_return = col_option_input(obstype_options)
            if desc_return is None:
                print("This is not an option, select an observation type.")
                sys.exit("invalid obstype for wide dataset, see last message. ")
            wide_obstype_name = inv_obstype_desc[desc_return]

            # 1) add a new obstype
            if wide_obstype_name == "ADD NEW OBSERVATION TYPE":
                new_obstype, cur_unit = add_new_obstype()
                wide_obstype_name = new_obstype.name
                known_obstypes[new_obstype.name] = new_obstype  # add to knonw obstypes
                units = cur_unit
                description = new_obstype.description

            # 2) Known obstype
            else:
                # add unit
                units = get_unit(known_obstypes[wide_obstype_name])
                description = input("Some more details on the observation (optional): ")

            # update template

            obsdict = {
                "tlk_obstype": wide_obstype_name,
                "columnname": None,
                "unit": str(units),
                "description": str(description),
            }

            tmpl_dict["data_related"]["obstype_mapping"].append(obsdict)

    # =============================================================================
    # Map metadatafile
    # =============================================================================

    print("\n \n *******      Meta Data   ***********")

    if meta_avail:
        print(" ... opening the metadata file ...")
        metadatareader = find_suitable_reader(filepath=metadatafilepath, is_url=False)

        metadata = metadatareader.read(nrows=10)
        metacolumnnames = metadata.columns.to_list()

        # map the required columns (name)

        if data_avail:
            if format_option != 3:
                # if multiple stations are in the dataset, this column is required
                print("Which column does represent the NAMES of the stations?")
                name_column = col_option_input(metacolumnnames)
                tmpl_dict["metadata_related"]["name_column"] = name_column
                metacolumnnames.remove(name_column)

            # if the data is a single station, this column is ignored
            else:
                staname = input("\n What is the name of your station : ")
                tmpl_dict["single_station_name"] = staname
        else:
            # no data is available, but still needs the 'name' column
            print("Which column does represent the NAMES of the stations?")
            name_column = col_option_input(metacolumnnames)
            tmpl_dict["metadata_related"]["name_column"] = name_column
            metacolumnnames.remove(name_column)

        # map columns that required for the template (lat, lon, altitude)
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

        with_altitude = yes_no_ques(
            "\n Is there a columns in the metadata representing the ALTITUDE (in meter) of the stations?"
        )
        if with_altitude:
            print("Which column does represent the ALTITUDE?")
            altitude_column = col_option_input(metacolumnnames)
            tmpl_dict["metadata_related"]["altitude_column"] = altitude_column
            metacolumnnames.remove(altitude_column)

        # Which other (not used by the toolkit) to add.
        if len(metacolumnnames) > 0:
            add_cols = yes_no_ques(
                f"\n Do you want to include one or more of the remaining columns in the metadatafile? \n ({str(metacolumnnames)})"
            )
            if add_cols:
                for col in metacolumnnames:
                    add_bool = yes_no_ques(f"\n Add {col} in the metada?")
                    if add_bool:
                        tmpl_dict["metadata_related"]["columns_to_include"].append(
                            str(col)
                        )

    # =============================================================================
    # Saving the template
    # =============================================================================

    print("\n ------ Saving the template ----- \n")
    save_dir = usr_input_dir(
        "Give a directory where to save the template (as template.json)"
    )

    # write to csv
    templatefilepath = os.path.join(save_dir, "template.json")
    write_dict_to_json(dictionary=tmpl_dict, trgfile=templatefilepath)
    print(f" DONE! The template is written here: {templatefilepath}")

    # =============================================================================
    # Tips for the user
    # =============================================================================

    apply_tips = yes_no_ques("Do you want some help creating your Dataset?")
    if apply_tips is True:
        print("\n\n ========= RUN THIS CODE ========= \n\n")

        print("\n# Define the paths to your files: ")
        if data_avail:
            print(f'data_file = r"{datafilepath}"')
        if meta_avail:
            print(f'meta_data_file = r"{metadatafilepath}"')

        print(f'template = r"{templatefilepath}"')

        print("\n# Initiate a dataset:")
        print("your_dataset = metobs_toolkit.Dataset()")

        # add new obstypes if needed
        to_add_obstypes = [
            newobsname
            for newobsname in known_obstypes.keys()
            if newobsname not in tlk_obstypes.keys()
        ]
        if bool(to_add_obstypes):
            print(
                "\n# Define non-standard observation types, and add them to the dataset: "
            )
            for newob in to_add_obstypes:
                new_obstype = known_obstypes[newob]
                print("new_obstype = metobs_toolkit.Obstype(")
                print(f'                 obsname="{new_obstype.name}",')
                print(f'                 std_unit="{new_obstype.std_unit}",')
                print(f'                 description="{new_obstype.description}",')
                print("                 )")
                print("\n\n #add the new obstype to your dataset. \n")
                print("your_dataset.add_new_observationtype(obstype=new_obstype)")
                print("\n\n")

        # add new units if needed

        print("\n# Import your data :")
        if data_avail:
            print("your_dataset.import_data_from_file(")
            print("    input_data_file = data_file,")
            if meta_avail:
                print("    input_metadata_file = meta_data_file,")

            print("    template_file = template,")
            print("    #The following arguments are filled with default values.")
            # add defualt arguments
            skiplist = ["input_data_file", "input_metadata_file", "template_file"]
            kwargdict = get_function_defaults(Dataset.import_data_from_file)
            for kwarg, kwargvalue in kwargdict.items():
                if kwarg in skiplist:
                    continue
                if isinstance(kwargvalue, str):
                    print(f'    {kwarg} = "{kwargvalue}",')
                else:
                    print(f"    {kwarg} = {kwargvalue},")

            print(")")

        else:
            print("your_dataset.import_only_metadata_from_file(")
            print("    input_metadata_file = meta_data_file,")
            print("    template_file = template,")
            print(")")

    return
