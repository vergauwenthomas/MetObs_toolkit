#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 26 10:14:38 2024

@author: thoverga
"""

import sys
import logging
import numpy as np
import pandas as pd
import ee

from metobs_toolkit.landcover_functions import (
    extract_pointvalues,
    extract_buffer_frequencies,
)

logger = logging.getLogger(__name__)


class GeeExtractor:
    """Class holding methods for extracting data from Google Earth Engine."""

    def __init__(self):
        logger.info("Initiate GeeExtractor instance")
        self.location = None  # GEE location
        self.usage = None  # Human readable application domain
        self.scale = None  # scale of the gee dataset ( = highest resolution)
        self.value_type = "numeric"  # numeric or categorical
        self.is_image = False
        self.is_imagecollection = True
        self.dynamical = True  # time evolution? To be used for timeseries

        self.missing_value_label = np.nan
        # Optional
        self.band_of_use = (
            None  # band to use for imagecollections (or None if no band available)
        )
        self.credentials = None

        self.categorical_map = (
            {}
        )  # dictionary to map categories {gee_Value: target_str}
        self.categorical_aggregation = (
            {}
        )  # dictionary to aggregate categories to {target_str : [list of gee_values]}

    def __str__(self):
        if self.location is None:
            return f"GeeExtractor with unknown target gee location."
        else:
            return f"GeeExtractor with for {self.usage} use: \n \
  location : {self.location}"

    def __repr__(self):
        return self.__str__()

    # =============================================================================
    # Extraction methods
    # =============================================================================

    def connect_to_gee(self, **kwargs):
        """Authenticate to GEE if needed."""
        if not ee.data._credentials:  # check if ee connection is initialized
            ee.Authenticate(**kwargs)
            ee.Initialize()
        return

    def extract_static_point_values(self, metadf, aggregate=False):

        logger.info(f"Extracting static point values from {self.location}")
        # check if all attributes are set and correct
        if self.location is None:
            sys.exit(f"No values could be extracted from {self}.")

        if self.scale is None:
            sys.exit(
                f"Scale is {self.scale}, no values could be extracted from {self}."
            )
        if self.dynamical:
            sys.exit(f"Could not extract static values of a dynamic geedataset: {self}")

        # check if metadf is valid
        if not _is_metadf_valid(metadf):
            sys.exit()

        filtered_metadf = metadf[metadf[["lat", "lon"]].notnull().all(1)]

        # authenticate to GEE
        self.connect_to_gee()

        valuesdf = extract_pointvalues(
            metadf=filtered_metadf,
            scale=self.scale,
            trg_gee_loc=self.location,
            band_of_use=self.band_of_use,
            is_imagecollection=self.is_imagecollection,
            is_image=self.is_image,
        )
        # Format to a series
        combdf = metadf.merge(
            valuesdf, how="left", left_index=True, right_on="feature_idx"
        )
        combdf = combdf.set_index(metadf.index)
        values_series = combdf[self.band_of_use]

        # Map to human labels if needed
        if aggregate:
            if self.value_type != "categorical":
                sys.exit(
                    f"Aggregating covers is only valid for gee datasets with value_type = categorical."
                )
            if not bool(self.categorical_aggregation):
                sys.exit(
                    "Cannot aggregate categorical values, since no categorical_aggregation is defined."
                )
            aggmap_to_human = _create_agg_map_dict(self.categorical_aggregation)
            values_series = values_series.map(aggmap_to_human)
        elif bool(self.categorical_map):
            values_series = values_series.map(self.categorical_map)

        # Replace missing values
        values_series = values_series.fillna(self.missing_value_label)

        return values_series

    def extract_static_buffer_frequencies(self, metadf, buffer, aggregate=False):

        logger.info(f"Extracting static buffer frequencies from {self.location}")
        # check if all attributes are set and correct
        if self.location is None:
            sys.exit(f"No values could be extracted from {self}.")

        if self.scale is None:
            sys.exit(
                f"Scale is {self.scale}, no values could be extracted from {self}."
            )
        if self.dynamical:
            sys.exit(f"Could not extract static values of a dynamic geedataset: {self}")
        if self.value_type != "categorical":
            sys.exit(
                f"Extarcting buffer frequencies is only valid for Categorical dataset, not {self.value_type}."
            )

        # check if metadf is valid
        if not _is_metadf_valid(metadf):
            sys.exit()

        filtered_metadf = metadf[metadf[["lat", "lon"]].notnull().all(1)]

        # authenticate to GEE
        self.connect_to_gee()

        valuesdf = extract_buffer_frequencies(
            metadf=filtered_metadf,
            scale=self.scale,
            trg_gee_loc=self.location,
            band_of_use=self.band_of_use,
            is_imagecollection=self.is_imagecollection,
            is_image=self.is_image,
            bufferradius=buffer,
        )

        # Format to a Dataframe
        combdf = metadf.index.to_frame().merge(
            valuesdf, how="left", left_index=True, right_on="feature_idx"
        )
        combdf = combdf.set_index(metadf.index)
        valuesdf = combdf.drop(
            columns=["feature_idx", metadf.index.name], errors="ignore"
        )

        # Aggregate if needed
        if aggregate:
            if self.value_type != "categorical":
                sys.exit(
                    f"Aggregating covers is only valid for gee datasets with value_type = categorical."
                )
            if not bool(self.categorical_aggregation):
                sys.exit(
                    "Cannot aggregate categorical values, since no categorical_aggregation is defined."
                )
            aggdf = pd.DataFrame(index=valuesdf.index)
            for aggcat, coverlist in self.categorical_aggregation.items():
                present_agg_classes = [
                    str(num) for num in coverlist if str(num) in valuesdf.columns
                ]
                aggdf[aggcat] = valuesdf[present_agg_classes].sum(axis=1, skipna=False)

            valuesdf = aggdf
        # Map to human labels if needed
        elif bool(self.categorical_map):
            # typecast values to str
            renamer = {str(key): val for key, val in self.categorical_map.items()}
            valuesdf = valuesdf.rename(columns=renamer)

        # TODO: aggregate if needed

        # Replace missing values
        valuesdf = valuesdf.fillna(self.missing_value_label)
        return valuesdf

    # =============================================================================
    # default targets
    # =============================================================================
    def activate_LCZ(self):
        self.location = "RUB/RUBCLIM/LCZ/global_lcz_map/v1"
        self.usage = "LCZ"
        self.scale = 100
        self.value_type = "categorical"
        self.is_image = False
        self.is_imagecollection = True
        self.dynamical = False

        self.missing_value_label = "Unknown"
        # Optional
        self.band_of_use = "LCZ_Filter"
        self.credentials = "Demuzere M.; Kittner J.; Martilli A.; Mills, G.; Moede, C.; Stewart, I.D.; van Vliet, J.; Bechtel, B. A global map of local climate zones to support earth system modelling and urban-scale environmental science. Earth System Science Data 2022, 14 Volume 8: 3835-3873. doi:10.5194/essd-14-3835-2022"
        self.categorical_map = {
            1: "Compact highrise",  # mapvalue: (color, human class)
            2: "Compact midrise",
            3: "Compact lowrise",
            4: "Open highrise",
            5: "Open midrise",
            6: "Open lowrise",
            7: "Lightweight lowrise",
            8: "Large lowrise",
            9: "Sparsely built",
            10: "Heavy industry",
            11: "Dense Trees (LCZ A)",
            12: "Scattered Trees (LCZ B)",
            13: "Bush, scrub (LCZ C)",
            14: "Low plants (LCZ D)",
            15: "Bare rock or paved (LCZ E)",
            16: "Bare soil or sand (LCZ F)",
            17: "Water (LCZ G)",
        }

        self.categorical_aggregation = {}

    def activate_worldcover(self):
        self.location = "ESA/WorldCover/v200"
        self.usage = "landcover"
        self.scale = 10
        self.value_type = "categorical"
        self.is_image = False
        self.is_imagecollection = True
        self.dynamical = False

        self.missing_value_label = np.nan
        # Optional
        self.band_of_use = "Map"
        self.credentials = "https://spdx.org/licenses/CC-BY-4.0.html"
        self.categorical_map = {
            10: "Tree cover",  # mapvalue: (color, human class)
            20: "Shrubland",
            30: "Grassland",
            40: "Cropland",
            50: "Built-up",
            60: "Bare / sparse vegetation",
            70: "Snow and ice",
            80: "Permanent water bodies",
            90: "Herbaceous wetland",
            95: "Mangroves",
            100: "Moss and lichen",
        }

        self.categorical_aggregation = {
            "water": [70, 80, 90, 95],
            "pervious": [10, 20, 30, 40, 60, 100],
            "impervious": [50],
        }

    def activate_DEM(self):

        self.location = "CGIAR/SRTM90_V4"
        self.usage = "elevation"
        self.scale = 100
        self.value_type = "numeric"
        self.is_image = True
        self.is_imagecollection = False
        self.dynamical = False

        self.missing_value_label = np.nan
        # Optional
        self.band_of_use = "elevation"
        self.credentials = "SRTM Digital Elevation Data Version 4"


# =============================================================================
# Helpers
# =============================================================================
def _create_agg_map_dict(aggdict):
    invdict = {}
    for key, val in aggdict.items():
        for cov in val:
            invdict[cov] = key

    return invdict


def _is_metadf_valid(metadf):
    assert not metadf.empty, "Metadf is an empty dataframe."
    assert (
        "lat" in metadf.columns
    ), f"lat column not found in columns of metadf {metadf.columns}"
    assert (
        "lon" in metadf.columns
    ), f"lon column not found in columns of metadf {metadf.columns}"
    assert (
        metadf.index.is_unique
    ), f"The index of the metadf is not unique: {metadf.index}"

    # TODO check if lat and lon columns are numerical
    return True
