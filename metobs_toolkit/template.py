#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 09:48:24 2024

@author: thoverga
"""

import os
import sys
import logging


import pandas as pd
from pytz import all_timezones

from metobs_toolkit.data_import import _read_csv_to_df


class Template:
    """Contains all info and methods to work with a template."""

    def __init__(self):
        # to renmae the columns
        self.obscolumnmap = {}  # toolkitname --> datacolumnname
        self.namemap = {"name": None}  # name --> name column in data AND metadata

        # obstypes details
        self.obsdetails = {}  # obsname(tlk) --> {unit: , description: ..}

        # Special always required
        self.dataformat = "long"  # long or wide or single_station?
        self.timestampinfo = {
            "datetimecolumn": None,
            "time_column": None,
            "date_column": None,
            "fmt": None,
        }

        # Extra options (not required)
        self.tz = None

        # Not activaly used attributes
        self.filepath = None

    def get_info(self):
        key_len = 7
        print("------ Template columns map ---------")
        for key, val in self.obscolumnmap.items():
            print(f" * {key.ljust(key_len)} <---> {val.ljust(key_len)}")

    # =============================================================================
    # Setters
    # =============================================================================
    def _set_filepath(self, filepath):
        self.filepath = str(filepath)

    def _set_dataformat(self, datafmt):
        if str(datafmt).lower() == "long":
            self.dataformat = "long"
            return
        if str(datafmt).lower() == "wide":
            self.dataformat = "wide"
            return
        if str(datafmt).lower() == "single_station":
            self.dataformat = "single_station"
            return
        sys.exit(
            f"The data format specified in the template ({datafmt}) is neither long, wide or single_station"
        )

    def _set_obs_info(self, obstempldf):
        for _idx, row in obstempldf.iterrows():
            self.obscolumnmap[str(row["varname"])] = str(row["template column name"])

            self.obsdetails[str(row["varname"])] = {
                "unit": str(row["units"]),
                "description": str(row["description"]),
            }

    def _set_name(self, namecolumn):
        self.namemap["name"] = str(namecolumn)

    def _set_tz(self, tzstring):
        if str(tzstring) in all_timezones:
            self.tz = str(tzstring)
            return
        sys.exit(
            f"{tzstring} is not a valid timezone string. See pytz.all_timezones for all valid timezone strings."
        )

    def _set_timestampinfo(self, templdf):

        if "datetime" in templdf["varname"].values:
            # get columnname of the timestamps
            self.timestampinfo["datetimecolumn"] = str(
                templdf.loc[
                    templdf["varname"] == "datetime", "template column name"
                ].squeeze()
            )
            # get the datetime format
            self.timestampinfo["fmt"] = str(
                templdf.loc[templdf["varname"] == "datetime", "format"].squeeze()
            )
            return
        if "_date" in templdf["varname"].values:
            if "_time" in templdf["varname"].values:
                # get columnname of the timestamps
                self.timestampinfo["date_column"] = str(
                    templdf.loc[
                        templdf["varname"] == "_date", "template column name"
                    ].squeeze()
                )
                self.timestampinfo["time_column"] = str(
                    templdf.loc[
                        templdf["varname"] == "_time", "template column name"
                    ].squeeze()
                )

                # get the datetime format
                date_fmt = str(
                    templdf.loc[templdf["varname"] == "_date", "format"].squeeze()
                )
                time_fmt = str(
                    templdf.loc[templdf["varname"] == "_time", "format"].squeeze()
                )
                self.timestampinfo["fmt"] = (
                    f"{date_fmt} {time_fmt}"  # concat with space inbetween
                )
                return
            else:
                sys.exit(
                    "No time column is mapped by the template (_date is mapped, but not _time)."
                )

        else:
            sys.exit(
                "No datetime columns are mapped by the template. Neither datetime or a _date and _time combination is found."
            )

    # =============================================================================
    # Getters
    # =============================================================================


def read_csv_template(file, known_obstypes, data_long_format=True):
    """
    Import a template from a csv file.

    Format options will be stored in a seperate dictionary. (Because these
    do not relate to any of the data columns.)

    Parameters
    ----------
    file : str
        Path to the csv template file.
    known_obstypes : list
        A list of known observation types. These consist of the default
        obstypes and the ones added by the user.
    data_long_format : bool, optional
        If True, this format structure has priority over the format structure
        in the template file. The default is True.

    Returns
    -------
    template : dict
        The template related to the data/metadata columns.
    opt_kwargs : dict
        Options and settings present in the template.

    """
    template = Template()

    templdf = _read_csv_to_df(filepath=file, kwargsdict={})
    # Drop emty rows
    templdf = templdf.dropna(axis="index", how="all")

    # Extracting general settings
    assert (
        "options" in templdf.columns
    ), 'The "options" column is not present in the template.'
    assert (
        "options_values" in templdf.columns
    ), 'The "options" column is not present in the template.'

    optionsdf = templdf[["options", "options_values"]].dropna(axis="index", how="all")
    options = dict(zip(optionsdf["options"], optionsdf["options_values"]))

    # Updatet template attributes
    template._set_filepath(file)

    assert (
        "data_structure" in options.keys()
    ), 'the "data_structure" is a required option that must be in the "options" column of the template.'
    template._set_dataformat(datafmt=options["data_structure"])

    # Get timestamps info (how they are represented in the data)
    template._set_timestampinfo(templdf=templdf)

    if "timezone" in options:
        template._set_tz(tzstring=options["timezone"])

    assert (
        "name" in templdf["varname"].values
    ), '"name" is required in the "varname" column of the template.'
    template._set_name(
        namecolumn=templdf.loc[
            templdf["varname"] == "name", "template column name"
        ].squeeze()
    )

    obstempldf = templdf.loc[templdf["varname"].isin(known_obstypes)]
    template._set_obs_info(obstempldf=obstempldf)

    return template


# def extract_options_from_template(templ, known_obstypes):
#     """Filter out options settings from the template dataframe.

#     Parameters
#     ----------
#     templ : pandas.DataFrame()
#         Template in a dataframe structure
#     known_obstypes : list
#         A list of known observation types. These consist of the default
#         obstypes and the ones added by the user.

#     Returns
#     -------
#     new_templ : pandas.DataFrame()
#         The template dataframe with optioncolumns removed.
#     opt_kwargs : dict
#         Options and settings present in the template dataframe.

#     """
#     opt_kwargs = {}
#     if "options" in templ.columns:
#         if "options_values" in templ.columns:
#             opt = templ[["options", "options_values"]]
#             # drop nan columns
#             opt = opt[opt["options"].notna()]
#             # convert to dict
#             opt = opt.set_index("options")["options_values"].to_dict()

#             # check options if valid
#             possible_options = {
#                 "data_structure": ["long", "wide", "single_station"],
#                 "stationname": "_any_",
#                 "obstype": known_obstypes,
#                 "obstype_unit": "_any_",
#                 "obstype_description": "_any_",
#                 "timezone": all_timezones,
#             }
#             for key, val in opt.items():
#                 key, val = str(key), str(val)
#                 if key not in possible_options:
#                     sys.exit(
#                         f"{key} is not a known option in the template. These are the possible options: {list(possible_options.keys())}"
#                     )

#                 if possible_options[key] == "_any_":
#                     pass  # value can be any string

#                 else:
#                     if val not in possible_options[key]:
#                         sys.exit(
#                             f"{val} is not a possible value for {key}. These values are possible for {key}: {possible_options[key]}"
#                         )

#                 # overload to kwargs:

#                 if key == "data_structure":
#                     if val == "long":
#                         opt_kwargs["long_format"] = True
#                     elif val == "wide":
#                         opt_kwargs["long_format"] = False
#                     else:
#                         # single station
#                         opt_kwargs["long_format"] = True
#                 if key == "stationname":
#                     if not opt["data_structure"] == "single_station":
#                         logger.warning(
#                             f'{val} as {key} in the template options will be ignored because the datastructure is not "single_station" (but {opt["data_structure"]})'
#                         )
#                     else:
#                         opt_kwargs["single"] = val
#                 if key == "obstype":
#                     opt_kwargs["obstype"] = val
#                 if key == "obstype_unit":
#                     opt_kwargs["obstype_unit"] = val
#                 if key == "obstype_description":
#                     opt_kwargs["obstype_description"] = val
#                 if key == "timezone":
#                     opt_kwargs["timezone"] = val

#         else:
#             sys.exit(
#                 '"options" column found in the template, but no "options_values" found!'
#             )

#     # remove the options from the template
#     new_templ = templ.drop(columns=["options", "options_values"], errors="ignore")
#     return new_templ, opt_kwargs
