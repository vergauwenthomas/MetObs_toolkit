#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# __init__.py

import logging
from pathlib import Path
import os, sys

# Create the Logger
loggers = logging.getLogger(__name__)  # logger name is <vlinder-toolkit>
loggers.setLevel(logging.DEBUG)

log_path = os.path.join(str(Path(__file__).parent.parent.parent), "logfile.log")

# # Create the Handler for logging data to a file - will be hereted for children
logger_handler = logging.FileHandler(filename=log_path)
logger_handler.setLevel(logging.DEBUG)

# # Create a Formatter for formatting the log messages
logger_formatter = logging.Formatter("%(name)s - %(levelname)s - %(message)s")
logger_handler.setFormatter(logger_formatter)

# Add the Handler to the Logger
loggers.addHandler(logger_handler)
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
#  Static variables to be reached by users
# =============================================================================
observation_types = [
    "temp",
    "radiation_temp",
    "humidity",
    "precip",
    "precip_sum",
    "wind_speed",
    "wind_gust",
    "wind_direction",
    "pressure",
    "pressure_at_sea_level",
]


# =============================================================================
# Import classes and function to be used by the user
# =============================================================================

from metobs_toolkit.dataset import Dataset
from metobs_toolkit.station import Station
from metobs_toolkit.modeldata import Modeldata

# import GUI
from metobs_toolkit.gui_launcher import launch_gui

# =============================================================================
# Import extenders
# =============================================================================
from metobs_toolkit.dataset_settings_updater import Dataset

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!
__version__ = "0.1.1a0"

