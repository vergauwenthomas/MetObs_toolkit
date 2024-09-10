#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:02:08 2024

@author: thoverga
"""

import os
import sys
from pathlib import Path


def _set_path_to_toolkit():
    lib_folder = Path(__file__).resolve().parents[1]
    if sys.path[0] != lib_folder:
        sys.path.insert(0, str(lib_folder))

    import metobs_toolkit

    print(f" Using metobs: v{metobs_toolkit.__version__}")


def _set_pandas_display_options():
    import pandas as pd

    # pd.options.display.float_format = '{:,.2f}'.format
    pd.set_option("display.width", 50)
    pd.set_option("display.max_columns", 8)
    pd.set_option("display.max_colwidth", 25)
    pd.set_option("display.max_rows", 10)
    pd.set_option("display.max_dir_items", 30)
    pd.set_option("display.precision", 2)
    pd.set_option("display.expand_frame_repr", False)


def setup_and_run_doctest():
    import doctest

    # set path to the this metobstoolkit
    _set_path_to_toolkit()
    # Set display options
    _set_pandas_display_options()

    # Get all files currently present in cwd
    cur_files = os.listdir(os.getcwd())

    # Run doctest
    doctest.testmod(
        verbose=False,  # print only failures
        report=True,
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
    )

    # delete the files that are generated (and stored in cwd) by the doctest
    files_after_check = os.listdir(os.getcwd())
    created = list(set(files_after_check) - set(cur_files))
    for f in created:
        if not f.endswith(".py"):  # better be sure
            print(f"Deleting doctest output artifact: {f}")
            os.remove(f)
