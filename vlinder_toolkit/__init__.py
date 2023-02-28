#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# __init__.py

import logging
from pathlib import Path
import os
# Create the Logger
loggers = logging.getLogger(__name__) #logger name is <vlinder-toolkit>
loggers.setLevel(logging.DEBUG)

log_path = os.path.join(str(Path(__file__).parent.parent.parent), 'logfile.log')

# # Create the Handler for logging data to a file - will be hereted for children
logger_handler = logging.FileHandler(filename=log_path)
logger_handler.setLevel(logging.DEBUG)

# # Create a Formatter for formatting the log messages
logger_formatter = logging.Formatter('%(name)s - %(levelname)s - %(message)s')
logger_handler.setFormatter(logger_formatter)

# Add the Handler to the Logger
loggers.addHandler(logger_handler)
loggers.info('Logger initiated') 






#Import classes
from .dataset import Dataset
from .settings import Settings


#Import QC checks
from . import qc_checks




__version__ = "version debug"



