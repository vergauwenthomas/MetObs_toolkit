#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# flake8: noqa: F401

import os
import sys
import logging

# =============================================================================
# Import classes and function to be used by the user
# =============================================================================


from metobs_toolkit.dataset import Dataset
from metobs_toolkit.obstypes import Obstype, ModelObstype, ModelObstype_Vectorfield

from metobs_toolkit.analysis import Analysis
from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
)

from metobs_toolkit.sensordata import SensorData  # give user acces
from metobs_toolkit.modeltimeseries import ModelTimeSeries  # give user acces

from metobs_toolkit.geedatasetmanagers import (
    default_datasets as default_GEE_datasets,
)

# Special functions that can be directly called by te user
from metobs_toolkit.dataset import import_dataset_from_pkl
from metobs_toolkit.template_build_prompt import build_template_prompt
from metobs_toolkit.gee_api import connect_to_gee

from metobs_toolkit.backend_collection.loggingmodule import (
    add_FileHandler,
    add_StreamHandler,
)


# =============================================================================
# Specify demo data paths
# =============================================================================

BASE_PATH = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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
# Setup logs
# =============================================================================

# Create the Root logger
rootlog = logging.getLogger(__name__)  # logger name is <metobs-toolkit>
rootlog.setLevel(logging.DEBUG)  # set rootlogger on debug
rootlog.handlers.clear()  # clear all handlers
# rootlog.setLevel(logging.DEBUG)  # set rootlogger on debug


# Set the default handler
console_handler = logging.StreamHandler()
console_handler.setLevel("WARNING")
formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
console_handler.setFormatter(formatter)
rootlog.addHandler(console_handler)


# =============================================================================
# Version
# =============================================================================

from metobs_toolkit.settings_collection.version import __version__
