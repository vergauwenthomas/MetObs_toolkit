#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:28:20 2024

@author: thoverga
"""

import copy
import logging
import pandas as pd


from metobs_toolkit.settings import Settings
from metobs_toolkit.obstypes import tlk_obstypes


logger = logging.getLogger(__name__)


class _DatasetBase(object):
    """Abstract class for holding data and metadata."""

    def __init__(self):
        """Construct all the necessary attributes for Dataset object."""
        logger.info("Initialise dataset")

        # Dataset with 'good' observations
        self.df = pd.DataFrame()

        # Dataset with metadata (static)
        self.metadf = pd.DataFrame()

        # dictionary storing present observationtypes
        self.obstypes = copy.copy(tlk_obstypes)  # init with all tlk obstypes
        self.settings = copy.deepcopy(Settings())

    # =============================================================================
    #     attribute setters
    # =============================================================================
    def _set_df(self, df):
        self.df = df

    def _set_metadf(self, metadf):
        self.metadf = metadf

    def _set_obstypes(self, obstypes):
        self.obstypes = obstypes

    def _set_settings(self, settings):
        self.settings = settings
