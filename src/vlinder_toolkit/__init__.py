#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# __init__.py


#Import classes
from .stations import Station, Dataset
from .settings import Settings


#Import QC checks
from . import qc_checks


#import templates and settingsfiles
# from .data_templates import csv_templates, db_templates
# from .settings_files import server_login.json


__version__ = "version debug"


#import functions
# from .IO import import_data 


# =============================================================================
# 
# =============================================================================
