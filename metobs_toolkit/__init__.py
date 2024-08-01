#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import logging
from pathlib import Path


# =============================================================================
# Setup logs
# =============================================================================
# Create the Root logger
rootlog = logging.getLogger(__name__)  # logger name is <metobs-toolkit>
rootlog.setLevel(logging.DEBUG)  # set rootlogger on debug
rootlog.handlers.clear()  # clear all handlers

from metobs_toolkit.logging import add_FileHandler, add_StreamHandler

add_FileHandler(
    trglogfile=os.path.join(str(Path(__file__).parent.parent.parent), "logfile.log")
)
rootlog.info("Logger initiated")


# =============================================================================
# Specify demo data paths
# =============================================================================

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

# import the Dataset core

from metobs_toolkit.dataset_core import Dataset

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
