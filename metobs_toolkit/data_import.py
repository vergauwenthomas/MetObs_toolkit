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
from pytz import all_timezones
from metobs_toolkit.template import _create_datetime_column

logger = logging.getLogger(__name__)


# =============================================================================
# Helpers
# =============================================================================


def _read_csv_to_df(filepath, kwargsdict):
    assert not isinstance(filepath, type(None)), f"No filepath is specified: {filepath}"
    common_seperators = [None, ";", ",", "    ", "."]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")

        if bool(kwargsdict):
            if "sep" not in kwargsdict:
                for sep in common_seperators:
                    df = pd.read_csv(filepath, sep=sep, **kwargsdict)
                    assert not df.empty, f"{filepath} is empty!"
                    if len(df.columns) > 1:
                        break
            else:
                df = pd.read_csv(filepath_or_buffer=filepath, **kwargsdict)
        else:
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

    # 1. Read metadata into df
    df = _read_csv_to_df(filepath=input_file, kwargsdict=kwargs_metadata_read)

    # 2. Check template on blacklist columnnames and compatibility
    template._metadata_template_compatibility_test(df.columns)

    blacklist_remapper = template._apply_blacklist(columns=df.columns, on_data=False)
    df.rename(columns=blacklist_remapper, inplace=True)

    # 3. map the name column
    # If single stations, AND name column is None (not mapped), then add name column with default name
    if (template._is_data_single_station()) & (
        template.metadata_namemap["name"] is None
    ):
        df["name"] = template.single_station_name
    else:
        df.rename(columns=template._get_metadata_name_map(), inplace=True)

    # 4. map all observation column names
    metacolmap = template._get_metadata_column_map()
    df.rename(columns=metacolmap, inplace=True)

    # 5. subset to relevant columns (name + datetime + all_obstypes_in_tlk_space)
    relev_columns = ["name"]
    relev_columns.extend(list(metacolmap.values()))
    try:
        df = df[list(set(relev_columns))]
    except KeyError as e:
        raise MetobsDataImportError(
            "The template refers to columns not present in the metadata file."
        )

    # make sure the names are strings
    df["name"] = df["name"].astype(str)

    # 8. map to numeric dtypes
    for col in df.columns:
        if col in ["lat", "lon"]:  # all other columns are observations
            df[col] = pd.to_numeric(df[col], errors="coerce")

    return df


# =============================================================================
# Data
# =============================================================================


def wide_to_long(df, obstypename):
    """Convert a wide dataframe to a long format.

    Convert a wide dataframe that represents obstype-observations to a long
    dataframe (=standard toolkit structure).

    Parameters
    ----------
    df : pandas.DataFrame()
        Wide dataframe with a "datetime" column.
    obstypename : str
        A MetObs obstype name.

    Returns
    -------
    longdf : pandas.DataFrame
        Long dataframe.
    template : dict
        Updated template dictionary.

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
        value_name=obstypename,
    )
    return longdf


def import_data_from_csv(
    input_file,
    template,
    known_obstypes,
    kwargs_data_read,
):
    """Import data as a dataframe.

    Parameters
    ----------
    input_file : str
        Path to the data (csv) file.
    template : metobs_toolkit.Template
        template for mapping the data.
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
    logger.debug(
        f"Reading the data records from {input_file}, with kwargs: {kwargs_data_read}"
    )
    df = _read_csv_to_df(filepath=input_file, kwargsdict=kwargs_data_read)

    # Test blacklist columns and apply compatibility check of template and data
    blacklist_remapper = template._apply_blacklist(columns=df.columns, on_data=True)
    df.rename(columns=blacklist_remapper, inplace=True)

    template._data_template_compatibility_test(datacolumns=df.columns)

    # 2. Create a singel "datetime" column (needed for wide-long transfer)
    df = _create_datetime_column(df=df, template=template)

    # 2B if the data is singel station, add a name column if is is not present
    if template._is_data_single_station():
        # check if there is a column indicating the name of the station that is mapped
        assumed_name_col = list(template._get_data_name_map().keys())[0]
        if assumed_name_col is None:
            df["_dummy_name_column"] = template._get_single_station_default_name()
            # add it to the template
            template._set_dataname("_dummy_name_column")

    # 2C if the data is wide, convert it into a long format
    elif not template._is_data_long():
        df = wide_to_long(df=df, obstypename=template._get_wide_obstype())

    # 3. map the name column
    df.rename(columns=template._get_data_name_map(), inplace=True)
    # make sure the names are strings
    df["name"] = df["name"].astype(str)
    # 4. map all observation column names
    df.rename(columns=template._get_obs_column_map(), inplace=True)

    # 5. subset to relevant columns (name + datetime + all_obstypes_in_tlk_space)
    try:
        df = df[template._get_all_mapped_data_cols_in_tlk_space()]
    except Exception as e:
        raise MetobsDataImportError(
            "Template-Datafile mismatch (template refers to non-existing column in the dataset)."
        )

    # 6.. Set index
    df.set_index("datetime", inplace=True)

    # 7. create timezone-aware datetime index
    df.index = df.index.tz_localize(tz=template._get_tz())
    logger.debug(f"df head: \n {df.head()}")
    return df


# =============================================================================
# Exceptions
# =============================================================================


class MetobsDataImportError(Exception):
    """Exception raised for errors on importing Data."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
