#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import logging
from pathlib import Path


# Create the Logger
loggers = logging.getLogger(__name__)  # logger name is <metobs-toolkit>
loggers.setLevel(logging.DEBUG)
loggers.handlers.clear()  # clear all handlers

# File handler
log_path = os.path.join(str(Path(__file__).parent.parent.parent), "logfile.log")
# # Create the Handler for logging data to a file - will be hereted for children
file_handler = logging.FileHandler(filename=log_path)
file_handler.setLevel(logging.DEBUG)
# # Create a Formatter for formatting the log messages
file_logger_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
file_handler.setFormatter(file_logger_formatter)
# Add the Handler to the Logger
loggers.addHandler(file_handler)

# add Streamhandler
streamhandler = logging.StreamHandler()
streamhandler.setLevel(logging.DEBUG)
stream_logger_formatter = logging.Formatter(
    "LOG:: %(name)s - %(levelname)s - %(message)s"
)
streamhandler.setFormatter(stream_logger_formatter)
loggers.addHandler(streamhandler)


loggers.info("Logger initiated")


BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(BASE_PATH)

# demo files
demo_datafile = os.path.join(
    BASE_PATH, "metobs_toolkit", "datafiles", "demo_datafile.csv"
)
demo_metadatafile = os.path.join(
    BASE_PATH, "metobs_toolkit", "datafiles", "demo_metadatafile.csv"
)
demo_template = os.path.join(
    BASE_PATH, "metobs_toolkit", "datafiles", "demo_template.json"
)


# =============================================================================
# Import classes and function to be used by the user
# =============================================================================

# import the Dataset core + extensions
# (Do not change order!!)
from metobs_toolkit.dataset_core import Dataset
from metobs_toolkit.dataset_settings_updater import Dataset
from metobs_toolkit.dataset_visuals import Dataset
from metobs_toolkit.dataset_gap_handling import Dataset
from metobs_toolkit.dataset_qc_handling import Dataset


# User accesable classes
from metobs_toolkit.station import Station  # after all Dataset extensions !!
from metobs_toolkit.modeldata import Modeldata
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.obstype_modeldata import ModelObstype, ModelObstype_Vectorfield
from metobs_toolkit.analysis import Analysis


# Special functions that can be directly called by te user
from metobs_toolkit.template_build_prompt import build_template_prompt
from metobs_toolkit.landcover_functions import connect_to_gee


# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.2.3"
