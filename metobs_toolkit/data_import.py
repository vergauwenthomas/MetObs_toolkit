#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 16:24:06 2022

@author: thoverga
"""
import sys
import warnings
import logging
import pandas as pd

# import mysql.connector
# from mysql.connector import errorcode
from pytz import all_timezones

logger = logging.getLogger(__name__)


# =============================================================================
# Helpers
# =============================================================================


def _remove_keys_from_dict(dictionary, keys):
    for key in keys:
        dictionary.pop(key, None)
    return dictionary


def template_to_package_space(specific_template):
    """Invert template dictionary."""
    returndict = {
        val["varname"]: {"orig_name": key} for key, val in specific_template.items()
    }
    for key, item in returndict.items():
        orig_dict = dict(specific_template[item["orig_name"]])
        orig_dict.pop("varname")
        returndict[key].update(orig_dict)
    return returndict


def find_compatible_templatefor(df_columns, template_list):
    """Test if template is compatible with dataaframe columns."""
    for templ in template_list:
        found = all(keys in list(df_columns) for keys in templ.keys())
        if found:
            logger.info("Compatible template found. ")
            return templ
    sys.exit("No compatible teplate found!")


def compress_dict(nested_dict, valuesname):
    """Unnest dictionary info for valuename.

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


def _read_csv_to_df(filepath, kwargsdict):
    assert not isinstance(filepath, type(None)), f"No filepath is specified: {filepath}"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if bool(kwargsdict):
            df = pd.read_csv(filepath_or_buffer=filepath, **kwargsdict)
        else:
            common_seperators = [None, ";", ",", "    ", "."]
            for sep in common_seperators:
                df = pd.read_csv(filepath, sep=sep)
                assert not df.empty, f"{filepath} is empty!"

                if len(df.columns) > 1:
                    break

    assert (
        len(df.columns) > 1
    ), f"Only one column detected from import using these seperators: {common_seperators}. See if csv template is correct."

    return df


# =============================================================================
# Template
# =============================================================================


def check_template_compatibility(template, df_columns, filetype):
    """Log template compatiblity with dataframe columns.

    Parameters
    ----------
    template : dict
        Template dictionary.
    df_columns : pd.index
        Dataframe columns to map.
    filetype : str
        'data', 'metadata' or other description of the dataframe.

    Returns
    -------
    None.

    """
    # ignore datetime because this is already mapped
    present_cols = [col for col in df_columns if col != "datetime"]
    assumed_cols = [key for key in template.keys() if key != "datetime"]

    # in mapper but not in df
    unmapped_assumed = [
        templ_var for templ_var in assumed_cols if templ_var not in present_cols
    ]

    if len(unmapped_assumed) > 0:
        logger.info(
            f"The following columns are not present in the {filetype},\
 and cannot be mapped: {unmapped_assumed}"
        )

    # in df but not in mapper
    unmapped_appearing = [col for col in present_cols if col not in assumed_cols]
    if len(unmapped_appearing) > 0:
        logger.info(
            f"The following columns in the {filetype} cannot be mapped with the template: {unmapped_appearing}."
        )

    # check if at least one column is mapped
    if len(list(set(present_cols) - set(assumed_cols))) == len(present_cols):
        sys.exit(
            f"Fatal: The given template: {assumed_cols} does not match with any of the {filetype} columns: {present_cols}."
        )


def extract_options_from_template(templ, known_obstypes):
    """Filter out options settings from the template dataframe.

    Parameters
    ----------
    templ : pandas.DataFrame()
        Template in a dataframe structure
    known_obstypes : list
        A list of known observation types. These consist of the default
        obstypes and the ones added by the user.

    Returns
    -------
    new_templ : pandas.DataFrame()
        The template dataframe with optioncolumns removed.
    opt_kwargs : dict
        Options and settings present in the template dataframe.

    """
    opt_kwargs = {}
    if "options" in templ.columns:
        if "options_values" in templ.columns:
            opt = templ[["options", "options_values"]]
            # drop nan columns
            opt = opt[opt["options"].notna()]
            # convert to dict
            opt = opt.set_index("options")["options_values"].to_dict()

            # check options if valid
            possible_options = {
                "data_structure": ["long", "wide", "single_station"],
                "stationname": "_any_",
                "obstype": known_obstypes,
                "obstype_unit": "_any_",
                "obstype_description": "_any_",
                "timezone": all_timezones,
            }
            for key, val in opt.items():
                key, val = str(key), str(val)
                if key not in possible_options:
                    sys.exit(
                        f"{key} is not a known option in the template. These are the possible options: {list(possible_options.keys())}"
                    )

                if possible_options[key] == "_any_":
                    pass  # value can be any string

                else:
                    if val not in possible_options[key]:
                        sys.exit(
                            f"{val} is not a possible value for {key}. These values are possible for {key}: {possible_options[key]}"
                        )

                # overload to kwargs:

                if key == "data_structure":
                    if val == "long":
                        opt_kwargs["long_format"] = True
                    elif val == "wide":
                        opt_kwargs["long_format"] = False
                    else:
                        # single station
                        opt_kwargs["long_format"] = True
                if key == "stationname":
                    if not opt["data_structure"] == "single_station":
                        logger.warning(
                            f'{val} as {key} in the template options will be ignored because the datastructure is not "single_station" (but {opt["data_structure"]})'
                        )
                    else:
                        opt_kwargs["single"] = val
                if key == "obstype":
                    opt_kwargs["obstype"] = val
                if key == "obstype_unit":
                    opt_kwargs["obstype_unit"] = val
                if key == "obstype_description":
                    opt_kwargs["obstype_description"] = val
                if key == "timezone":
                    opt_kwargs["timezone"] = val

        else:
            sys.exit(
                '"options" column found in the template, but no "options_values" found!'
            )

    # remove the options from the template
    new_templ = templ.drop(columns=["options", "options_values"], errors="ignore")
    return new_templ, opt_kwargs


def read_csv_template(file, known_obstypes, data_long_format=True):
    """Import a template from a csv file.

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
    templ = _read_csv_to_df(filepath=file, kwargsdict={})

    # Extract structure options from template
    templ, opt_kwargs = extract_options_from_template(templ, known_obstypes)

    # Drop emty rows
    templ = templ.dropna(axis="index", how="all")

    if "long_format" in opt_kwargs.keys():
        data_long_format = opt_kwargs["long_format"]

    if data_long_format:
        # Drop variables that are not present in templ
        templ = templ[templ["template column name"].notna()]

    # templates have nested dict structure where the keys are the column names in the csv file, and the
    # values contain the mapping information to the toolkit classes and names.

    # create dictionary from templframe
    templ = templ.set_index("template column name")

    # create a dict from the dataframe, remove Nan value row wise
    template = {}
    for idx, row in templ.iterrows():
        template[idx] = row[~row.isnull()].to_dict()

    return template, opt_kwargs


# =============================================================================
# Metadata
# =============================================================================


def import_metadata_from_csv(input_file, template, kwargs_metadata_read):
    """Import metadata as a dataframe.

    Parameters
    ----------
    input_file : str
        Path to the metadata (csv) file.
    template : dict
        Template dictionary.
    kwargs_metadata_read : dict
        Extra user-specific kwargs to pass to the pd.read_csv() function.

    Returns
    -------
    df : pandas.DataFrame()
        The metadata in a pandas dataframe with columnnames in the toolkit
        standards.

    """
    assert not isinstance(input_file, type(None)), "Specify input file in the settings!"

    df = _read_csv_to_df(input_file, kwargs_metadata_read)

    # validate template
    # template = read_csv_template(template_file)
    check_template_compatibility(template, df.columns, filetype="metadata")

    # rename columns to toolkit attriute names
    column_mapper = {val["orig_name"]: key for key, val in template.items()}
    df = df.rename(columns=column_mapper)

    return df


# =============================================================================
# Data
# =============================================================================


def wide_to_long(df, template, obstype):
    """Convert a wide dataframe to a long format.

    Convert a wide dataframe that represents obstype-observations to a long
    dataframe (=standard toolkit structure).

    Parameters
    ----------
    df : pandas.DataFrame()
        Wide dataframe.
    template : dict
        The dictionary to update the 'name' key on.
    obstype : str
        A MetObs obstype.

    Returns
    -------
    longdf : pandas.DataFrame
        Long dataframe.
    template : dict
        Updateted template dictionary.

    """
    # the df is assumed to have one datetime column, and the others represent
    # stations with their obstype values

    stationnames = df.columns.to_list()
    stationnames.remove("datetime")

    longdf = pd.melt(
        df,
        id_vars=["datetime"],
        value_vars=stationnames,
        var_name="name",
        value_name=obstype,
    )

    # # update template
    # template[obstype] = template["_wide_dummy"]
    # del template["_wide_dummy"]

    # add name to the template
    template["name"] = {"varname": "name", "dtype": "object"}

    return longdf, template


def import_data_from_csv(
    input_file,
    template,
    long_format,
    obstype,
    obstype_units,
    obstype_description,
    known_obstypes,
    kwargs_data_read,
):
    """Import data as a dataframe.

    Parameters
    ----------
    input_file : str
        Path to the data (csv) file.
    template : dict
        template dictionary.
    long_format : bool
        If True, a long format is assumed else wide.
    obstype : str
        If format is wide, this is the observationtype.
    obstype_units : str
       If format is wide, this is the observation unit.
    obstype_description : str
        If format is wide, this is the observation description.
    known_obstypes : list
        A list of known observation types. These consist of the default
        obstypes and the ones added by the user.
    kwargs_data_read : dict
        Kwargs passed to the pd.read_csv() function.

    Returns
    -------
    df : pandas.DataFrame()
        A long dataframe containing the observations.
    invtemplate : dict
        Template in toolkit space.

    """
    # 1. Read data into df
    df = _read_csv_to_df(filepath=input_file, kwargsdict=kwargs_data_read)

    # 2. Read template
    invtemplate = template_to_package_space(template)

    # 3. Make datetime column (needed for wide to long conversion)
    if "datetime" in invtemplate.keys():

        df = df.rename(columns={invtemplate["datetime"]["orig_name"]: "datetime"})
        df["datetime"] = pd.to_datetime(
            df["datetime"], format=invtemplate["datetime"]["format"]
        )

        inv_temp_remove_keys = ["datetime"]
        temp_remove_keys = [invtemplate["datetime"]["orig_name"]]
    elif ("_date" in invtemplate.keys()) & ("_time" in invtemplate.keys()):

        datetime_fmt = (
            invtemplate["_date"]["format"] + " " + invtemplate["_time"]["format"]
        )
        df["datetime"] = pd.to_datetime(
            df[invtemplate["_date"]["orig_name"]]
            + " "
            + df[invtemplate["_time"]["orig_name"]],
            format=datetime_fmt,
        )
        df = df.drop(
            columns=[
                invtemplate["_date"]["orig_name"],
                invtemplate["_time"]["orig_name"],
            ]
        )

        inv_temp_remove_keys = ["_time", "_date"]
        temp_remove_keys = [
            invtemplate["_date"]["orig_name"],
            invtemplate["_time"]["orig_name"],
        ]
    else:
        sys.exit(
            "Impossible to map the dataset to a datetime column, verify your template please."
        )

    # 3.b Remove the datetime keys from the template

    invtemplate = _remove_keys_from_dict(invtemplate, inv_temp_remove_keys)
    template = _remove_keys_from_dict(template, temp_remove_keys)

    # 4. convert wide data to long if needed
    if not long_format:
        template[obstype] = {}
        invtemplate[obstype] = {}
        template[obstype]["varname"] = obstype
        invtemplate[obstype]["orig_name"] = obstype  # use default as orig name
        if obstype_units is not None:
            template[obstype]["units"] = obstype_units
            invtemplate[obstype]["units"] = obstype_units
        if obstype_description is not None:
            template[obstype]["description"] = obstype_description
            invtemplate[obstype]["description"] = obstype_description

        df, template = wide_to_long(df, template, obstype)

    # 5. check compatibility
    check_template_compatibility(template, df.columns, filetype="data")

    # 6. map to default name space
    df = df.rename(columns=compress_dict(template, "varname"))

    # 7. Keep only columns as defined in the template
    cols_to_keep = list(invtemplate.keys())
    cols_to_keep.append("datetime")
    cols_to_keep.append("name")
    cols_to_keep = list(set(cols_to_keep))
    df = df.loc[:, df.columns.isin(cols_to_keep)]

    # 8. Set index
    df = df.reset_index()
    df = df.drop(columns=["index"], errors="ignore")
    df = df.set_index("datetime")

    # 8. map to numeric dtypes
    for col in df.columns:
        if col in known_obstypes:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        if col in ["lon", "lat"]:
            df[col] = pd.to_numeric(df[col], errors="coerce")

    # add template to the return
    return df, invtemplate
