#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 26 09:31:51 2023

@author: thoverga
"""

from metobs_toolkit import dataset

class Station(dataset.Dataset):
    def __init__(self, name, df, outliersdf, gaps, missing_obs, gapfilldf,
                 metadf, data_template, settings):
        self.name = name
        self.df = df
        self.outliersdf = outliersdf
        self.gaps = gaps
        self.missing_obs = missing_obs
        self.gapfilldf = gapfilldf
        self.metadf = metadf
        self.data_template = data_template
        self.settings=settings


        self._istype = 'Station'

