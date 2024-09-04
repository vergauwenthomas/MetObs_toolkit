#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:37:33 2024

@author: thoverga
"""

import pandas as pd

print(f"pandas version: {pd.__version__}")

import metobs_toolkit

print("toolkit version: ", metobs_toolkit.__version__)


# %%

demofile = metobs_toolkit.demo_datafile

print(demofile)

df = pd.read_csv(demofile)  # naive import

print(f"naive df read: {df}")


df = pd.read_csv(demofile, sep=";", engine="python")  # python engine

print(f"naive df read: {df}")


# %%


def _read_csv_to_df(filepath, kwargsdict):
    assert not isinstance(filepath, type(None)), f"No filepath is specified: {filepath}"
    common_seperators = [None, ";", ",", "    ", "."]
    # with warnings.catch_warnings():
    #     warnings.simplefilter("ignore")

    if bool(kwargsdict):
        if "sep" not in kwargsdict:
            for sep in common_seperators:
                print(f" seperator: {sep}")
                df = pd.read_csv(filepath, sep=sep, **kwargsdict)
                print(f"df :{df}")
                assert not df.empty, f"{filepath} is empty!"
                if len(df.columns) > 1:
                    break
        else:
            df = pd.read_csv(filepath_or_buffer=filepath, **kwargsdict)
    else:
        for sep in common_seperators:
            print(f" seperator: {sep}")
            df = pd.read_csv(filepath, sep=sep)
            print(f"df :{df}")
            assert not df.empty, f"{filepath} is empty!"

            if len(df.columns) > 1:
                break
    print(df)  # REPORT ME IF YOU SEE ME
    assert (
        len(df.columns) > 1
    ), f"Only one column detected from import using these seperators: {common_seperators}. See if csv template is correct."

    return df


_read_csv_to_df(demofile, kwargsdict={})
