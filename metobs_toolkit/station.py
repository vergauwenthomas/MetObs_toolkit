#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Station class that inherits all methods of the Dataset
class.

A Station holds all observations of one station.
"""

from metobs_toolkit import dataset


class Station(dataset.Dataset):
    def __init__(
        self,
        name,
        df,
        outliersdf,
        gaps,
        missing_obs,
        gapfilldf,
        missing_fill_df,
        metadf,
        data_template,
        settings,
        _qc_checked_obstypes,
        _applied_qc,
    ):
        self.name = name
        self.df = df
        self.outliersdf = outliersdf
        self.gaps = gaps
        self.missing_obs = missing_obs
        self.gapfilldf = gapfilldf
        self.missing_fill_df = missing_fill_df
        self.metadf = metadf
        self.data_template = data_template
        self.settings = settings
        self._qc_checked_obstypes = _qc_checked_obstypes
        self._applied_qc = _applied_qc

        self._istype = "Station"
