#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Module with functions for writing csv files.

@author: thoverga
"""

import os
import logging

logger = logging.getLogger(__name__)


# =============================================================================
# IO simple helpers
# =============================================================================
def _does_trg_file_exist(trg):
    """Check if a target file/folder already exists. Return bool."""
    return os.path.isfile(trg)


def _remove_file(trg):
    logger.info(f"Removing {trg}")
    os.remove(trg)
    return


def write_df_to_csv(df, trgfile, to_csv_kwargs={}):
    # write to csv in output folder
    logger.info(f"write to file: {trgfile}")
    print(f"write to file: {trgfile}")
    df.to_csv(path_or_buf=trgfile, **to_csv_kwargs)


class MetobsWritingError(Exception):
    """Exception raised for errors in the template."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
