#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import logging
from pathlib import Path


# Create the Logger
loggers = logging.getLogger(__name__)  # logger name is <metobs-toolkit>
loggers.setLevel(logging.DEBUG)


# Adding Handlers

# File handler
log_path = os.path.join(str(Path(__file__).parent.parent.parent), "logfile.log")
# # Create the Handler for logging data to a file - will be hereted for children
file_handler = logging.FileHandler(filename=log_path)
file_handler.setLevel(logging.DEBUG)
# # Create a Formatter for formatting the log messages
logger_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
file_handler.setFormatter(logger_formatter)
# Add the Handler to the Logger
loggers.addHandler(file_handler)


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
    BASE_PATH, "metobs_toolkit", "datafiles", "demo_templatefile.csv"
)


# =============================================================================
# Import classes and function to be used by the user
# =============================================================================

from metobs_toolkit.dataset import Dataset
from metobs_toolkit.station import Station
from metobs_toolkit.modeldata import Modeldata
from metobs_toolkit.obstypes import Obstype
from metobs_toolkit.obstype_modeldata import ModelObstype, ModelObstype_Vectorfield


# Special functions that can be directly called by te user
from metobs_toolkit.data_templates.template_build_prompt import build_template_prompt
from metobs_toolkit.landcover_functions import connect_to_gee

# =============================================================================
# Import extenders
# =============================================================================
from metobs_toolkit.dataset_settings_updater import Dataset

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.2.0"
