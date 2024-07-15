#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Station class that inherits all methods of the Dataset
class.

A Station holds all observations of one station.
"""
import pandas as pd
from metobs_toolkit import dataset


class Station(dataset.Dataset):
    """A class holding all information of one station. Inherit all from Dataset."""

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
        obstypes,
        template,
        settings,
        _qc_checked_obstypes,
        _applied_qc,
    ):
        """Initiate the Station object."""
        dataset.Dataset.__init__(self)

        # Set data attributes (using abstract dataset class)
        self._set_df(df)
        self._set_metadf(metadf)
        self._set_obstypes(obstypes)
        self._set_settings(settings)

        # Set data attribute
        self.gaps = gaps
        self.missing_obs = missing_obs
        self.gapfilldf = gapfilldf
        self.missing_fill_df = missing_fill_df
        self.template = template
        self._qc_checked_obstypes = _qc_checked_obstypes
        self._applied_qc = _applied_qc

        # Station specific
        self.name = name
        self._istype = "Station"

        # setup of object
        self.setup_metadata_dtyes()

    def setup_metadata_dtyes(self):
        """Make sure the dtypes are not lost when subsetting."""
        numeric_columns = ["lat", "lon"]
        timedelta_columns = ["assumed_import_frequency", "dataset_resolution"]

        for col in numeric_columns:
            if col in self.metadf.columns:
                self.metadf[col] = pd.to_numeric(self.metadf[col])

        for col in timedelta_columns:
            if col in self.metadf.columns:
                self.metadf[col] = pd.to_timedelta(self.metadf[col])
