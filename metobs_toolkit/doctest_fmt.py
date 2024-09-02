#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:02:08 2024

@author: thoverga
"""


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

    # Set display options
    _set_pandas_display_options()

    # Run doctest
    doctest.testmod(
        verbose=False,  # print only failures
        report=True,
        optionflags=doctest.ELLIPSIS | doctest.NORMALIZE_WHITESPACE,
    )
