#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import os
import sys
import logging
from datetime import datetime
from pathlib import Path


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

# from metobs_toolkit.dataset_core import Dataset
from metobs_toolkit.Newdataset import Dataset


# User accesable classes
from metobs_toolkit.oldstation import Station  # after all Dataset extensions !!

# from metobs_toolkit.modeldata import Modeldata
from metobs_toolkit.obstypes import Obstype, ModelObstype, ModelObstype_Vectorfield

# from metobs_toolkit.analysis import Analysis
from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
)
from metobs_toolkit.gap import (
    Gap,
)  # No direct usecase, but needed for creation of artificial gaps (+ doc api)


from metobs_toolkit.geedatasetmanagers import default_datasets as default_GEE_datasets


# Special functions that can be directly called by te user
from metobs_toolkit.Newdataset import import_dataset_from_pkl
from metobs_toolkit.template_build_prompt import build_template_prompt
from metobs_toolkit.gee_api import connect_to_gee

# =============================================================================
# Version
# =============================================================================

# DO not change this manually!

__version__ = "0.4.0a"


# =============================================================================
# Setup logs
# =============================================================================


# Create the Root logger
rootlog = logging.getLogger(__name__)  # logger name is <metobs-toolkit>
rootlog.setLevel(logging.DEBUG)  # set rootlogger on debug
rootlog.handlers.clear()  # clear all handlers

from metobs_toolkit.loggingmodule import add_FileHandler, add_StreamHandler

add_FileHandler(
    trglogfile=os.path.join(str(Path(__file__).parent.parent.parent), "logfile.log")
)
rootlog.info("Logger initiated")
