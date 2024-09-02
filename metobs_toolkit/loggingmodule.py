#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  2 14:23:30 2024

@author: thoverga
"""

import os
import logging
from datetime import datetime


def add_FileHandler(trglogfile, setlvl="DEBUG", clearlog=True):

    if clearlog:
        if os.path.isfile(trglogfile):
            os.remove(trglogfile)

    # log_path = os.path.join(str(Path(__file__).parent.parent.parent), "logfile.log")
    # # Create the Handler for logging data to a file - will be hereted for children
    file_handler = logging.FileHandler(filename=trglogfile)
    file_handler.setLevel(setlvl.upper())  # set handler level on debug
    # # Create a Formatter for formatting the log messages
    file_logger_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(file_logger_formatter)
    # Add the Handler to the Logger
    rootlog = logging.getLogger("<metobs_toolkit>")
    rootlog.addHandler(file_handler)

    rootlog.info(f"FileHandler set at {datetime.now()}")


def add_StreamHandler(setlvl="DEBUG"):

    # add Streamhandler
    streamhandler = logging.StreamHandler()
    streamhandler.setLevel(setlvl.upper())
    stream_logger_formatter = logging.Formatter(
        "LOG:: %(name)s - %(levelname)s - %(message)s"
    )
    streamhandler.setFormatter(stream_logger_formatter)
    # Add the Handler to the Logger
    rootlog = logging.getLogger("<metobs_toolkit>")
    rootlog.addHandler(streamhandler)

    rootlog.info(f"StreamHandler set at {datetime.now()}")


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
