#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 28 08:12:02 2022

@author: thoverga
"""
import os, sys
from pathlib import Path
import pandas as pd


default_template_file = os.path.join(
    str(Path(__file__).parent), "template_defaults", "default_template.csv"
)


# =============================================================================
# templates
# =============================================================================


# templates have nested dict structure where the keys are the column names in the csv file, and the
# values contain the mapping information to the toolkit classes and names.


def read_csv_template(file, data_long_format=True, obstype=None):
    common_seperators = [";", ",", "    ", "."]
    for sep in common_seperators:
        templ = pd.read_csv(file, sep=sep)
        assert not templ.empty, "Template is empty!"

        if len(templ.columns) > 1:
            break

    # Drop emty rows
    templ = templ.dropna(axis="index", how="all")

    if data_long_format:
        # Drop variables that are not present in templ
        templ = templ[templ["template column name"].notna()]

    else:
        # Do not do this for wide dataframes since the present obstype
        # is not specified in template oclumn name, but the defenition and datatype do.
        templ.loc[templ["varname"] == obstype, "template column name"] = "_wide_dummy"


    # create dictionary from templframe
    templ = templ.set_index("template column name")

    # create a dict from the dataframe, remove Nan value row wise
    template = {}
    for idx, row in templ.iterrows():
        template[idx] = row[~row.isnull()].to_dict()

    return template
