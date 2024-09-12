#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Modeldata class and all its methods.

A Modeldata holds all timeseries coming from a model and methods to use them.
"""
import os
import copy
import sys
import pickle
import pandas as pd
import logging
import datetime as datetimemodule
from time import sleep
import ee


from metobs_toolkit.df_helpers import (
    init_multiindexdf,
    conv_tz_multiidxdf,
    xs_save,
    multiindexdf_datetime_subsetting,
)

import metobs_toolkit.gee_api as gee_api
from metobs_toolkit.obstype_modeldata import default_era5_obstypes

from metobs_toolkit.plotting_functions import (
    model_timeseries_plot,
    timeseries_plot,
    folium_map,
    add_title_to_folium_map,
)
from metobs_toolkit.gee_api import connect_to_gee

# from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.obstype_modeldata import (
    # model_obstypes,
    ModelObstype,
    ModelObstype_Vectorfield,
)

logger = logging.getLogger(__name__)

# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================


class _GeeModelData:
    """Parent class for working with a GEE modeldataset."""

    def __init__(
        self,
        name,
        location,
        value_type,
        scale,
        is_static,
        is_image,
        credentials,
        is_mosaic=False,
    ):
        """
        Create a GeeModelData abstract instance.

        Parameters
        ----------
        name : str
            The user-defined name for referring to this GEE dataset.
        location : str
            The location of the dataset on GEE. Navigate the GEE datasets online,
            to get the location of a GEE dataset.
        value_type : "numeric" or "categorical"
            Specify how to interpret the values of the GEE dataset.
        scale : int
            The Scale (See GEE doc) of the dataset to extract values of.
        is_static : bool
            If True, the GEE dataset is static and has no time-evolution component.
            Else the GEE dataset has a time component.
        is_image : bool
            If True, the GEE dataset is opened as ee.Image(), else ee.ImageCollection().
        credentials : str
            Credentials of the GEE dataset.
        is_mosaic : bool, optional
            If True, ee.mosaic() is applied on the GEE dataset. The default is False.


        Returns
        -------
        None.

        """

        self.metadf = pd.DataFrame()  # will hold static data

        self.name = str(name)
        self.location = str(location)
        # self.band_of_use = str(band_of_use)

        if str(value_type) not in ["categorical", "numeric"]:
            raise MetobsModelDataError(
                f'value_type: {value_type} is not "categorical" or "numeric"'
            )
        self.value_type = str(value_type)

        self.scale = int(scale)
        self.is_static = bool(is_static)  # Time dependence or not
        self.is_image = bool(is_image)  # ee.image or ee.imagecollection
        self._is_mosaic = bool(
            is_mosaic
        )  # if true, .mosaic() is applied (stitching the maps into one big map)

        self.credentials = str(credentials)

    # =============================================================================
    # Checks
    # =============================================================================
    def _check_metadf_validity(self, metadf):
        """
        Check if a metadf is valid (coordinates and structure-wise). If
        it is not valid, an error is raised.

        Parameters
        ----------
        metadf : pandas.DataFrame()
            The metadata as a (geo)pandas dataframe.

        Returns
        -------
        None.

        """
        if metadf.empty:
            raise MetobsModelDataError(f"There is no metadata provided for {self}.")
        if metadf.index.name != "name":
            raise MetobsModelDataError(
                f"Wrong index name for setting {metadf} to {self}"
            )
        if "lat" not in metadf.columns:
            raise MetobsModelDataError(f'No "lat" column in the metadf of {self}.')
        if metadf["lat"].isnull().all():
            raise MetobsModelDataError(
                f'All values of the "lat" column in the metadf are Nan.'
            )
        if "lon" not in metadf.columns:
            raise MetobsModelDataError(f'No "lon" column in the metadf of {self}.')
        if metadf["lon"].isnull().all():
            raise MetobsModelDataError(
                f'All values of the "lon" column in the metadf are Nan.'
            )

    # =============================================================================
    # Setters
    # =============================================================================

    def set_metadf(self, metadf):
        """Setter for the metadf."""
        self._check_metadf_validity(metadf)
        self.metadf = metadf[["lon", "lat", "geometry"]]

    # =============================================================================
    # Getters
    # =============================================================================
    def _get_all_gee_bandnames(self):
        """Return a list of all the bandnames of the GEE dataset."""
        if self.is_image:
            return list(ee.Image(self.location).bandNames().getInfo())
        else:
            return list(ee.ImageCollection(self.location).first().bandNames().getInfo())

    def _get_base_details(self):
        """Print out basis details of the GEE dataset."""
        print(f"------ Details --------- \n")
        print(f" * name: {self.name}")
        print(f" * location: {self.location}")
        # print(f' * band of use: {self.band_of_use}')
        print(f" * value_type: {self.value_type}")
        print(f" * scale: {self.scale}")
        print(f" * is_static: {self.is_static}")
        print(f" * is_image: {self.is_image}")
        print(f" * is_mosaic: {self._is_mosaic}")
        print(f" * credentials: {self.credentials}")
        return

    def _clear_metadata(self):
        self.metadf = pd.DataFrame()  # will hold static data


class GeeStaticModelData(_GeeModelData):
    """Class for working with static GEE modeldatasets."""

    def __init__(
        self,
        name,
        location,
        band_of_use,
        value_type,
        scale,
        is_image,
        is_mosaic=False,
        credentials="",
        class_map={},
        agg_scheme={},
        col_scheme={},
    ):
        """
        Create a GeeStaticModelData instance representing a GEE dataset without
        a time dimensions.

        Parameters
        ----------
        name : str
            The user-defined name for referring to this GEE dataset.
        location : str
            The location of the dataset on GEE. Navigate the GEE datasets online,
            to get the location of a GEE dataset.
        band_of_use : str
            The name of the band to use. Navigate the GEE datasets online, and
            click on bands to get a table and description of them.
        value_type : "numeric" or "categorical"
            Specify how to interpret the values of the GEE dataset.
        scale : int
            The Scale (See GEE doc) of the dataset to extract values of. This
            can be found on the GEE dataset information page.
        is_image : bool
            If True, the GEE dataset is opened as ee.Image(), else
            ee.ImageCollection(). This can be found on the GEE dataset
            information page.
        is_mosaic : bool, optional
            If True, ee.mosaic() is applied on the GEE dataset. The default is False.
        credentials : str, optional
            Credentials of the GEE dataset. The default is "".
        class_map : dict, optional
            If value_type is categorical, then the class_map defines how the
            numeric values are mapped to 'human-categories'. The keys are the
            numeric values, the values are the human-labels. The default is {}.
        agg_scheme : dict, optional
            If value_types is categorical, then the agg scheme defines custom-
            made classes, which are aggregates of the present classes. The
            keys are the names of the custom-classes, the values are lists, with
            the corresponding numeric values. The default is {}.
        col_scheme : dict, optional
            if value_types is categorical, the col_sheme defines the colors used
            for each class. The keys are the numeric values, the values are
            the colors (in hex form). The default is {}.

        Returns
        -------
        None.

        See Also
        -----------
        GeeDynamicModelData: Gee Modeldata dataset for time-evolving data.
        GeeStaticModelData.set_metadf: Set metadata (station locations).
        GeeStaticModelData.get_info: Print out detailed info method.
        GeeStaticModelData.extract_static_point_data: Extract point values.
        GeeStaticModelData.extract_static_buffer_frac_data: Extract buffer fractions.
        GeeStaticModelData.make_gee_plot: Make an interactive spatial gee plot.

        Note
        -------
        In general, specifying a scale smaller than the true scale of the GEE
        dataset has no impact on the results (but can affect the computation time).

        Examples
        --------
        As an example, we will create a GeeStaticModelData instance representing
        representing the Copernicus Corine landcover dataset. This Dataset is
        available as a GEE dataset: https://developers.google.com/earth-engine/datasets/catalog/COPERNICUS_CORINE_V20_100m#bands.

        >>> import metobs_toolkit
        >>> corine = metobs_toolkit.GeeStaticModelData(
        ...             name="corine",
        ...             location="COPERNICUS/CORINE/V20/100m", #See GEE
        ...             band_of_use="landcover", #See GEE
        ...             value_type="categorical",
        ...             scale=100,
        ...             is_image=False,
        ...             is_mosaic=True,  # test this
        ...             credentials="Copernicus team:  https://land.copernicus.eu/pan-european/corine-land-cover/clc2018?tab=metadata.",
        ...             # Class definitions: this is given for categorical maps, see GEE.
        ...             class_map={
        ...                 111: "Artificial surfaces > Urban fabric > Continuous urban fabric",
        ...                 112: "Artificial surfaces > Urban fabric > Discontinuous urban fabric",
        ...                 121: "Artificial surfaces > Industrial, commercial, and transport units > Industrial or commercial units",
        ...                 122: "Artificial surfaces > Industrial, commercial, and transport units > Road and rail networks and associated land",
        ...                 123: "Artificial surfaces > Industrial, commercial, and transport units > Port areas",
        ...                 124: "Artificial surfaces > Industrial, commercial, and transport units > Airports",
        ...                 131: "Artificial surfaces > Mine, dump, and construction sites > Mineral extraction sites",
        ...                 132: "Artificial surfaces > Mine, dump, and construction sites > Dump sites",
        ...                 133: "Artificial surfaces > Mine, dump, and construction sites > Construction sites",
        ...                 141: "Artificial surfaces > Artificial, non-Agricultural vegetated areas > Green urban areas",
        ...                 142: "Artificial surfaces > Artificial, non-Agricultural vegetated areas > Sport and leisure facilities",
        ...                 211: "Agricultural areas > Arable land > Non-irrigated arable land",
        ...                 212: "Agricultural areas > Arable land > Permanently irrigated land",
        ...                 213: "Agricultural areas > Arable land > Rice fields",
        ...                 221: "Agricultural areas > Permanent crops > Vineyards",
        ...                 222: "Agricultural areas > Permanent crops > Fruit trees and berry plantations",
        ...                 223: "Agricultural areas > Permanent crops > Olive groves",
        ...                 231: "Agricultural areas > Pastures > Pastures",
        ...                 241: "Agricultural areas > Heterogeneous Agricultural areas > Annual crops associated with permanent crops",
        ...                 242: "Agricultural areas > Heterogeneous Agricultural areas > Complex cultivation patterns",
        ...                 243: "Agricultural areas > Heterogeneous Agricultural areas > Land principally occupied by Agriculture, with significant areas of natural vegetation",
        ...                 244: "Agricultural areas > Heterogeneous Agricultural areas > Agro-forestry areas",
        ...                 311: "Forest and semi natural areas > Forests > Broad-leaved Forest",
        ...                 312: "Forest and semi natural areas > Forests > Coniferous Forest",
        ...                 313: "Forest and semi natural areas > Forests > Mixed Forest",
        ...                 321: "Forest and semi natural areas > Scrub and/or herbaceous vegetation associations > Natural grasslands",
        ...                 322: "Forest and semi natural areas > Scrub and/or herbaceous vegetation associations > Moors and heathland",
        ...                 323: "Forest and semi natural areas > Scrub and/or herbaceous vegetation associations > Sclerophyllous vegetation",
        ...                 324: "Forest and semi natural areas > Scrub and/or herbaceous vegetation associations > Transitional woodland-shrub",
        ...                 331: "Forest and semi natural areas > Open spaces with little or no vegetation > Beaches, dunes, sands",
        ...                 332: "Forest and semi natural areas > Open spaces with little or no vegetation > Bare rocks",
        ...                 333: "Forest and semi natural areas > Open spaces with little or no vegetation > Sparsely vegetated areas",
        ...                 334: "Forest and semi natural areas > Open spaces with little or no vegetation > Burnt areas",
        ...                 335: "Forest and semi natural areas > Open spaces with little or no vegetation > Glaciers and perpetual snow",
        ...                 411: "Wetlands > Inland wetlands > Inland marshes",
        ...                 412: "Wetlands > Inland wetlands > Peat bogs",
        ...                 421: "Wetlands > Maritime wetlands > Salt marshes",
        ...                 422: "Wetlands > Maritime wetlands > Salines",
        ...                 423: "Wetlands > Maritime wetlands > Intertidal flats",
        ...                 511: "Water bodies > Inland waters > Water courses",
        ...                 512: "Water bodies > Inland waters > Water bodies",
        ...                 521: "Water bodies > Marine waters > Coastal lagoons",
        ...                 522: "Water bodies > Marine waters > Estuaries",
        ...                 523: "Water bodies > Marine waters > Sea and ocean",
        ...             },
        ...             # Aggregation scheme: leave blank, or create your own aggregation scheme.
        ...             agg_scheme={
        ...                 "Artificial_surface": [111,112,121,122,123,124,131,132,133,141,142],
        ...                 "Agricultural_surface": [211,212,213,221,222,223,231,241,242,243,244],
        ...                 "Forest_and_natural": [311,312,313,321,322,323,324,331,332,333,334,335],
        ...                 "Wetlands": [411,412,421,422,423],
        ...                 "Water": [511,512,521,522,523],
        ...             },
        ...             # color scheme: often an example is given, see GEE
        ...             col_scheme={
        ...                 111: "#e6004d",
        ...                 112: "#ff0000",
        ...                 121: "#cc4df2",
        ...                 122: "#cc0000",
        ...                 123: "#e6cccc",
        ...                 124: "#e6cce6",
        ...                 131: "#a600cc",
        ...                 132: "#a64dcc",
        ...                 133: "#ff4dff",
        ...                 141: "#ffa6ff",
        ...                 142: "#ffe6ff",
        ...                 211: "#ffffa8",
        ...                 212: "#ffff00",
        ...                 213: "#e6e600",
        ...                 221: "#e68000",
        ...                 222: "#f2a64d",
        ...                 223: "#e6a600",
        ...                 231: "#e6e64d",
        ...                 241: "#ffe6a6",
        ...                 242: "#ffe64d",
        ...                 243: "#e6cc4d",
        ...                 244: "#f2cca6",
        ...                 311: "#80ff00",
        ...                 312: "#00a600",
        ...                 313: "#4dff00",
        ...                 321: "#ccf24d",
        ...                 322: "#a6ff80",
        ...                 323: "#a6e64d",
        ...                 324: "#a6f200",
        ...                 331: "#e6e6e6",
        ...                 332: "#cccccc",
        ...                 333: "#ccffcc",
        ...                 334: "#000000",
        ...                 335: "#a6e6cc",
        ...                 411: "#a6a6ff",
        ...                 412: "#4d4dff",
        ...                 421: "#ccccff",
        ...                 422: "#e6e6ff",
        ...                 423: "#a6a6e6",
        ...                 511: "#00ccf2",
        ...                 512: "#80f2e6",
        ...                 521: "#00ffa6",
        ...                 522: "#a6ffe6",
        ...                 523: "#e6f2ff"}
        ...             )
        >>> corine
        GeeStaticModelData instance of corine  (no metadata has been set)

        If you want to use this map with your `Dataset`, add them to the known
        Modeldatasets.

        >>> dataset = metobs_toolkit.Dataset()
        >>> dataset.add_new_geemodeldata(Modeldata=corine)
        >>> dataset.gee_datasets
        {'lcz': GeeStaticModelData instance of lcz  (no metadata has been set) , 'altitude': GeeStaticModelData instance of altitude  (no metadata has been set) , 'worldcover': GeeStaticModelData instance of worldcover  (no metadata has been set) , 'ERA5-land': Empty GeeDynamicModelData instance of ERA5-land , 'corine': GeeStaticModelData instance of corine  (no metadata has been set) }
        """
        super().__init__(
            name=name,
            location=location,
            # band_of_use=band_of_use,
            value_type=value_type,
            scale=scale,
            is_static=True,
            is_image=is_image,
            is_mosaic=is_mosaic,
            credentials=credentials,
        )

        self.class_map = class_map
        self.agg_scheme = agg_scheme
        self.col_scheme = col_scheme
        self.band_of_use = band_of_use

        self.__name__ = "GeeStaticModelData"

    # ========================================================================
    # Special methods
    # ========================================================================
    def __str__(self):

        if self.metadf.empty:
            return (
                f"{self.__name__} instance of {self.name}  (no metadata has been set) "
            )
        else:
            return f"{self.__name__} instance of {self.name}  (known metadata)"

    def __repr__(self):
        """Print overview information of the modeldata."""
        return self.__str__()

    def _clear_data(self):
        """Clear all attributes that hold meta and extracted data."""
        self._clear_metadata()

    # =========================================================================
    # Setters
    # =========================================================================

    # =========================================================================
    # Getters
    # =========================================================================
    def get_info(self):
        """
        Prints out detailed information of the GeeStaticModelData.

        Returns
        -------
        None.

        Examples
        --------
        As an example, we will use the LCZ map, which is a default `GeeStaticModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.gee_datasets
        {'lcz': GeeStaticModelData instance of lcz  (no metadata has been set) , 'altitude': GeeStaticModelData instance of altitude  (no metadata has been set) , 'worldcover': GeeStaticModelData instance of worldcover  (no metadata has been set) , 'ERA5-land': Empty GeeDynamicModelData instance of ERA5-land }

        These are the defaults. We select the lcz one:

        >>> lcz_model = dataset.gee_datasets['lcz']
        >>> lcz_model
        GeeStaticModelData instance of lcz  (no metadata has been set)

        To print out detailed information use the `GeeStaticModelData.get_info()`
        method.

        >>> lcz_model.get_info()
        GeeStaticModelData instance of lcz  (no metadata has been set)
        ------ Details ---------
        <BLANKLINE>
         * name: lcz
         * location: RUB/RUBCLIM/LCZ/global_lcz_map/latest
         * value_type: categorical
         * scale: 100
         * is_static: True
         * is_image: False
         * is_mosaic: True
         * credentials: Demuzere M.; Kittner J.; Martilli A.; Mills, G.; Moede, C.; Stewart, I.D.; van Vliet, J.; Bechtel, B. A global map of local climate zones to support earth system modelling and urban-scale environmental science. Earth System Science Data 2022, 14 Volume 8: 3835-3873. doi:10.5194/essd-14-3835-2022
         -- Band --
         LCZ_Filter
        <BLANKLINE>
         -- classification --
         {1: 'Compact highrise', 2: 'Compact midrise', 3: 'Compact lowrise', 4: 'Open highrise', 5: 'Open midrise', 6: 'Open lowrise', 7: 'Lightweight lowrise', 8: 'Large lowrise', 9: 'Sparsely built', 10: 'Heavy industry', 11: 'Dense Trees (LCZ A)', 12: 'Scattered Trees (LCZ B)', 13: 'Bush, scrub (LCZ C)', 14: 'Low plants (LCZ D)', 15: 'Bare rock or paved (LCZ E)', 16: 'Bare soil or sand (LCZ F)', 17: 'Water (LCZ G)'}
        <BLANKLINE>
         -- aggregation --
         {}
        <BLANKLINE>
         -- colors --
         {1: '#8c0000', 2: '#d10000', 3: '#ff0000', 4: '#bf4d00', 5: '#ff6600', 6: '#ff9955', 7: '#faee05', 8: '#bcbcbc', 9: '#ffccaa', 10: '#555555', 11: '#006a00', 12: '#00aa00', 13: '#648525', 14: '#b9db79', 15: '#000000', 16: '#fbf7ae', 17: '#6a6aff'}
        <BLANKLINE>
        """
        print(str(self))
        self._get_base_details()
        print(f" -- Band -- \n {self.band_of_use} \n")
        print(f" -- classification -- \n {self.class_map} \n")
        print(f" -- aggregation -- \n {self.agg_scheme} \n")
        print(f" -- colors -- \n {self.col_scheme} \n")

        return

    # =========================================================================
    # Functionality
    # =========================================================================

    def extract_static_point_data(self):
        """
        Extract point values at the locations in the metadata.

        First, a connection with the gee server is made. Then the coordinates
        of the metadata and the details of the GEE dataset are sent to the
        GEE server. There the point values are extracted and are sent back.


        Returns
        -------
        df : pandas.DataFrame
            A dataframe with the stationnames as index, and one column with
            values (mapped to human-labels if the dataset is categorical and
            a cat_map is defined).

        See Also
        -----------
        metobs_toolkit.connect_to_gee: Establish connection with GEE services.
        GeeStaticModelData.set_metadf: Set metadata (station locations).
        GeeStaticModelData.get_info: Print out detailed info method.
        GeeStaticModelData.extract_static_buffer_frac_data: Extract buffer fractions.


        Note
        -------
        Make sure that the metadata is set. Use the
        `GeeStaticModelData.set_metadata()` for this.

        Note
        -------
        Make sure that you are authenticated to the GEE services.
        Run `metobs_toolkit.connect_to_gee()`, to do this. See the Documentation
        for more details.

        Examples
        --------
        As an example, we will use the LCZ map, which is a default `GeeStaticModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset

        In order to extract point data, we need metadata (i.g. coordinates
        of stations). Therefore we will import the demo dataset so that we
        can use the metadata to extract values.

        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)
        >>> dataset.metadf
                     lat   lon               school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.98  3.82                UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.02  3.71                UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.32  4.95          Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.34  4.93          Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.05  3.68         Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        ...          ...   ...                  ...                       ...                ...                       ...                       ...
        vlinder24  51.17  3.57        OLV ten Doorn  POINT (3.57206 51.16702)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder25  51.15  3.71    Einstein Atheneum  POINT (3.70861 51.15472)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder26  51.16  5.00          Sint Dimpna  POINT (4.99765 51.16176)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder27  51.06  3.73  Sec. Kunstinstituut   POINT (3.72807 51.0581)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder28  51.04  3.77             GO! Ath.  POINT (3.76974 51.03529)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        <BLANKLINE>
        [28 rows x 7 columns]

        Now we add the metadata to the Modeldataset.

        >>> lcz_model = dataset.gee_datasets['lcz']
        >>> lcz_model
        GeeStaticModelData instance of lcz  (no metadata has been set)
        >>> lcz_model.set_metadf(dataset.metadf) #overload the metadf
        >>> lcz_model
        GeeStaticModelData instance of lcz  (known metadata)

        To extract point-values, use the
        `GeeStaticModelData.extract_static_point_data()` method. First, make
        sure you are authenticated with the GEE services.

        >>> metobs_toolkit.connect_to_gee() #only required once per session
        >>> lcz_df = lcz_model.extract_static_point_data()
        >>> lcz_df
                                   lcz
        name
        vlinder01   Low plants (LCZ D)
        vlinder02        Large lowrise
        vlinder03         Open midrise
        vlinder04       Sparsely built
        vlinder05        Water (LCZ G)
        ...                        ...
        vlinder24  Dense Trees (LCZ A)
        vlinder25        Water (LCZ G)
        vlinder26         Open midrise
        vlinder27      Compact midrise
        vlinder28         Open lowrise
        <BLANKLINE>
        [28 rows x 1 columns]

        """
        if self.metadf.empty:
            raise MetobsModelDataError(
                f"No metadata is present for the GeeStaticModelData. No extraction possible."
            )

        # test if coordiantes are available
        self._check_metadf_validity(self.metadf)

        # =============================================================================
        # df to featurecollection
        # =============================================================================

        ee_fc = gee_api._df_to_features_point_collection(self.metadf)

        # =============================================================================
        # extract raster values
        # =============================================================================

        raster = gee_api.get_ee_obj(
            self
        )  # raster is image/imagecollection filtered to bands

        # extract points
        results = raster.sampleRegions(
            collection=ee_fc, scale=self.scale  # feature collection here
        ).getInfo()

        # extract properties
        if not bool(results["features"]):
            # no data retrieved
            logger.warning(
                f"Something went wrong, gee did not return any data: {results}"
            )
            logger.info(
                f"(Could it be that (one) these coordinates are not on the map: {self.metadf}?)"
            )
            return pd.DataFrame()

        # =============================================================================
        # to dataframe
        # =============================================================================

        properties = [x["properties"] for x in results["features"]]
        df = pd.DataFrame(properties)

        # map to human space if categorical
        if bool(self.class_map):
            df[self.band_of_use] = df[self.band_of_use].map(self.class_map)

        # rename to values to toolkit space
        df = df.rename(columns={self.band_of_use: self.name})

        # #format index
        df = df.set_index(["name"])

        return df

    def extract_static_buffer_frac_data(self, bufferradius, agg_bool=False):
        """Extract cover frequencies in circular buffers at stations.

        This methods will create (circular) buffers, specified by the radius,
        on each station. All the gridpoint in the radius are sampled and
        occurence frequency is computed for each landcover category in the
        buffer.

        If an aggregation scheme is known (see the .agg_scheme attribute), then
        the frequencies can be computed for the aggregated classes instead of the
        dataset classes.

        Parameters
        ----------
        bufferradius : int
            The radius (in meters) of the buffer.
        agg_bool : bool, optional
            If True, and if an aggregation scheme is known, then the frequencies
            are computed on the aggregated classes. Else it is computed for the
            GEE dataset classes. The default is False.

        Returns
        -------
        pandas.DataFrame
            A dataframe with the stationnames (from the metadf) as index. The
            columns are the class names. An additional "buffer_radius" column
            is added with the numeric value of the radius.

        See Also
        -----------
        metobs_toolkit.connect_to_gee: Establish connection with GEE services.
        GeeStaticModelData.set_metadf: Set metadata (station locations).
        GeeStaticModelData.get_info: Print out detailed info method.
        GeeStaticModelData.extract_static_point_data: Extract point values.

        Note
        -------
        Make sure that the metadata is set. Use the
        `GeeStaticModelData.set_metadata()` for this.

        Note
        -------
        Make sure that you are authenticated to the GEE services.
        Run `metobs_toolkit.connect_to_gee()`, to do this. See the Documentation
        for more details.

        Examples
        --------
        As an example, we will use the ESA worldcovermap, which is a default `GeeStaticModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset

        In order to extract point data, we need metadata (i.g. coordinates
        of stations). Therefore we will import the demo dataset, so that we
        can use the metadata to extract values.

        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)
        >>> dataset.metadf
                     lat   lon               school                  geometry dataset_resolution                  dt_start                    dt_end
        name
        vlinder01  50.98  3.82                UGent  POINT (3.81576 50.98044)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder02  51.02  3.71                UGent   POINT (3.7097 51.02238)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder03  51.32  4.95          Heilig Graf  POINT (4.95211 51.32458)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder04  51.34  4.93          Heilig Graf  POINT (4.93473 51.33552)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder05  51.05  3.68         Sint-Barbara  POINT (3.67518 51.05266)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        ...          ...   ...                  ...                       ...                ...                       ...                       ...
        vlinder24  51.17  3.57        OLV ten Doorn  POINT (3.57206 51.16702)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder25  51.15  3.71    Einstein Atheneum  POINT (3.70861 51.15472)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder26  51.16  5.00          Sint Dimpna  POINT (4.99765 51.16176)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder27  51.06  3.73  Sec. Kunstinstituut   POINT (3.72807 51.0581)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        vlinder28  51.04  3.77             GO! Ath.  POINT (3.76974 51.03529)    0 days 00:05:00 2022-09-01 00:00:00+00:00 2022-09-15 23:55:00+00:00
        <BLANKLINE>
        [28 rows x 7 columns]

        Now we add the metadata to the Modeldataset.

        >>> worldcover_model = dataset.gee_datasets['worldcover']
        >>> worldcover_model
        GeeStaticModelData instance of worldcover  (no metadata has been set)
        >>> worldcover_model.set_metadf(dataset.metadf) #overload the metadf
        >>> worldcover_model
        GeeStaticModelData instance of worldcover  (known metadata)

        For more details on the Modeldataset you can use the
        `GeeStaticModelData.get_info()` method.

        >>> worldcover_model.get_info()
        GeeStaticModelData instance of worldcover  (known metadata)
        ------ Details ---------
        <BLANKLINE>
         * name: worldcover
         * location: ESA/WorldCover/v200
         * value_type: categorical
         * scale: 10
         * is_static: True
         * is_image: False
         * is_mosaic: True
         * credentials: https://spdx.org/licenses/CC-BY-4.0.html
         -- Band --
         Map
        <BLANKLINE>
         -- classification --
         {10: 'Tree cover', 20: 'Shrubland', 30: 'Grassland', 40: 'Cropland', 50: 'Built-up', 60: 'Bare / sparse vegetation', 70: 'Snow and ice', 80: 'Permanent water bodies', 90: 'Herbaceous wetland', 95: 'Mangroves', 100: 'Moss and lichen'}
        <BLANKLINE>
         -- aggregation --
         {'water': [70, 80, 90, 95], 'pervious': [10, 20, 30, 40, 60, 100], 'impervious': [50]}
        <BLANKLINE>
         -- colors --
         {10: '006400', 20: 'ffbb22', 30: 'ffff4c', 40: 'f096ff', 50: 'fa0000', 60: 'b4b4b4', 70: 'f0f0f0', 80: '0064c8', 90: '0096a0', 95: '00cf75', 100: 'fae6a0'}
        <BLANKLINE>

        To extract buffer frequencies, use the
        `GeeStaticModelData.extract_static_buffer_frac_data()` method. First, make
        sure you are authenticated with the GEE services.

        >>> metobs_toolkit.connect_to_gee() #only required once per session

        You can choose a buffer radius (in meters) to compute cover fractions in.

        >>> cover_frac = worldcover_model.extract_static_buffer_frac_data(
        ...                             bufferradius=150)
        >>> cover_frac
                                  Grassland  Cropland  Built-up  Tree cover  Permanent water bodies  Bare / sparse vegetation  Herbaceous wetland
        name      buffer_radius
        vlinder01 150                 0.29      0.70  8.10e-03        0.00                    0.00                       0.0                 0.0
        vlinder02 150                 0.12      0.00  5.29e-01        0.35                    0.00                       0.0                 0.0
        vlinder03 150                 0.05      0.00  7.57e-01        0.19                    0.00                       0.0                 0.0
        vlinder04 150                 0.77      0.08  1.08e-01        0.05                    0.00                       0.0                 0.0
        vlinder05 150                 0.15      0.00  3.06e-01        0.15                    0.39                       0.0                 0.0
        ...                            ...       ...       ...         ...                     ...                       ...                 ...
        vlinder24 150                 0.03      0.00  8.73e-02        0.88                    0.00                       0.0                 0.0
        vlinder25 150                 0.07      0.01  5.24e-02        0.02                    0.84                       0.0                 0.0
        vlinder26 150                 0.03      0.00  7.84e-01        0.19                    0.00                       0.0                 0.0
        vlinder27 150                 0.01      0.00  9.39e-01        0.05                    0.00                       0.0                 0.0
        vlinder28 150                 0.15      0.00  3.88e-01        0.47                    0.00                       0.0                 0.0
        <BLANKLINE>
        [28 rows x 7 columns]


        As could be seen in the output of `get_info()`, there is an aggregation
        scheme present. We can apply that scheme to get aggregated fractions.

        >>> cover_frac = worldcover_model.extract_static_buffer_frac_data(
        ...                             bufferradius=150,
        ...                             agg_bool=True)
        >>> cover_frac
                                 water  pervious  impervious
        name      buffer_radius
        vlinder01 150             0.00      0.99    8.10e-03
        vlinder02 150             0.00      0.47    5.29e-01
        vlinder03 150             0.00      0.24    7.57e-01
        vlinder04 150             0.00      0.89    1.08e-01
        vlinder05 150             0.39      0.30    3.06e-01
        ...                        ...       ...         ...
        vlinder24 150             0.00      0.91    8.73e-02
        vlinder25 150             0.84      0.11    5.24e-02
        vlinder26 150             0.00      0.22    7.84e-01
        vlinder27 150             0.00      0.06    9.39e-01
        vlinder28 150             0.00      0.61    3.88e-01
        <BLANKLINE>
        [28 rows x 3 columns]

        """

        if self.metadf.empty:
            raise MetobsModelDataError(
                f"No metadata is present for the GeeStaticModelData. No extraction possible."
            )

        # test if coordiantes are available
        self._check_metadf_validity(self.metadf)

        # =============================================================================
        # df to featurecollection
        # =============================================================================
        ee_fc = gee_api._df_to_features_buffer_collection(
            self.metadf, int(bufferradius)
        )

        # =============================================================================
        # extract raster frequencies
        # =============================================================================

        def rasterExtraction(image):
            feature = image.reduceRegions(
                reducer=ee.Reducer.frequencyHistogram(),
                collection=ee_fc,  # feature collection here
                scale=self.scale,  # Cell size of raster
            )
            return feature

        raster = gee_api.get_ee_obj(self, force_mosaic=False)  # dataset
        results = raster.map(rasterExtraction).flatten().getInfo()
        # =============================================================================
        # to dataframe
        # =============================================================================

        freqs = {
            staprop["properties"]["name"]: staprop["properties"]["histogram"]
            for staprop in results["features"]
        }
        freqsdf = pd.DataFrame(freqs)

        # format frequency df
        freqsdf = freqsdf.transpose().fillna(0)
        freqsdf.index.name = "name"

        # normalize freqs
        freqsdf = freqsdf.div(freqsdf.sum(axis=1), axis=0)

        # =============================================================================
        # Aggregation and classification
        # =============================================================================
        if agg_bool:
            # create agg map
            agg_df = pd.DataFrame()
            for agg_name, agg_classes in self.agg_scheme.items():
                present_agg_classes = [
                    str(num) for num in agg_classes if str(num) in freqsdf.columns
                ]
                agg_df[agg_name] = freqsdf[present_agg_classes].sum(axis=1)

            freqsdf = agg_df

        # # map to human space if categorical
        elif bool(self.class_map):
            freqsdf = freqsdf.rename(
                columns={str(key): str(val) for key, val in self.class_map.items()}
            )
        else:
            pass

        # =============================================================================
        # Format dataframe
        # =============================================================================
        freqsdf["buffer_radius"] = bufferradius
        freqsdf = freqsdf.reset_index().set_index(["name", "buffer_radius"])

        return freqsdf

    def make_gee_plot(
        self, outputfolder=None, filename=None, vmin=None, vmax=None, overwrite=False
    ):
        """Make an interactive spatial plot of the GEE dataset and the stations.

        This method will create an interactive plot of the GEE dataset. If
        metadata is present, it will be displayed as markers on the map.

        The interactive map can be saved as an HTML file, by specifying the
        target path.


        Parameters
        ----------
        outputfolder : str or None, optional
            Path to the folder to save the HTML file. If None, the map will
            not be saved as an HTML file. The default is None.
        filename : str or None, optional
            The filename for the HTML file. If a filename is given, if it does
            not end with ".html", the prefix is added. If None, the map will not
            be saved as an HTML file. The default is None.
        vmin : num or None, optional
            If the dataset is not categorical, vmin is the minimum values
            assigned to the colormap. If None, vmin is computed by extracting
            the values at the locations of the stations. If no metadata is
            available, and vmin is None then vmin is set to 0. The default is None.
        vmax : num or None, optional
            If the dataset is not categorical, vmax is the minimum values
            assigned to the colormap. If None, vmax is computed by extracting
            the values at the locations of the stations. If no metadata is
            available, and vmax is None then vmax is set to 1. The default is None.
        overwrite : bool, optional
            If True, and if the target file exists, then it will be overwritten.
            Else, an error will be raised when the target file already exists.
            The default is False.

        Returns
        -------
        MAP : geemap.foliummap.Map
            The interactive map of the GeeStaticModelData.

        See Also
        -----------
        metobs_toolkit.connect_to_gee: Establish connection with GEE services.
        GeeStaticModelData.set_metadf: Set metadata (station locations).
        GeeStaticModelData.get_info: Print out detailed info method.

        Warning
        ---------
        To display the interactive map a graphical interactive backend is
        required, which could be missing. You can recognice this when no
        map is displayed, but the python console prints out a message similar
        to `<geemap.foliumap.Map at 0x7ff7586b8d90>`.

        In that case, you can specify a `outputfolder` and `outputfile`, save the map as a HTML file, and
        open in with a browser.

        Examples
        --------
        As an example we will make a plot of the LCZ map, which is a default `GeeStaticModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> lcz_model = dataset.gee_datasets['lcz']

        If you want your stations present in the map, then you must add the
        metadf to the Modeldata. This step is not required.

        >>> #we will use the demo metadata
        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)

        >>> lcz_model.set_metadf(dataset.metadf) #overload the metadf
        >>> lcz_model
        GeeStaticModelData instance of lcz  (known metadata)

        We will save the map as a (HTML) file. You can specify where so save it,
        for this example we will store it in the current working directory
        (`os.getcwd()`)

        >>> import os
        >>> map = lcz_model.make_gee_plot(
        ...                  outputfolder = os.getcwd(),
        ...                  filename = 'LCZ_map.html',
        ...                  overwrite=True)


        """

        # check if outputfile exists
        if (outputfolder is not None) | (filename is not None):
            save = True

            if not os.path.isdir(outputfolder):
                raise MetobsModelDataError(f"{outputfolder} is not a directory!")

            # check file extension in the filename:
            if filename[-5:] != ".html":
                filename += ".html"

            full_path = os.path.join(outputfolder, filename)

            # check if file exists
            if os.path.isfile(full_path):
                if overwrite:
                    logger.info(f"Overwrite the file at {full_path}.")
                    os.remove(full_path)
                else:
                    raise MetobsModelDataError(f"{full_path} is already a file!")
        else:
            save = False

        # Connect to Gee
        connect_to_gee()

        # Get image
        im = gee_api.get_ee_obj(self)

        # make empty map (FOLIUM BACKEND !! )
        MAP = folium_map()

        # show stations
        if self.metadf.empty:
            pass

        else:
            # Extract point values (for hover info and vmin/vmax in cont case)
            metadf = copy.deepcopy(self.metadf)
            df = self.extract_static_point_data()
            metadf = metadf.merge(df, how="left", left_index=True, right_index=True)

            _add_stations_to_folium_map(
                Map=MAP, metadf=metadf, display_cols=["name", self.name]
            )
            # fix center
            centroid = self.metadf.dissolve().centroid
            MAP.setCenter(lon=centroid.x.item(), lat=centroid.y.item(), zoom=8)

        if bool(self.class_map):
            # categorical color scheme
            vmin = min(self.class_map.keys())
            vmax = max(self.class_map.keys())
            if bool(self.col_scheme):
                var_visualization = {
                    "bands": [self.band_of_use],
                    "min": vmin,
                    "max": vmax,
                    "palette": list(self.col_scheme.values()),
                }

                # legend dict is human-label : color
                MAP.add_legend(
                    title="NLCD Land Cover Classification",
                    legend_dict={
                        self.class_map[numval]: self.col_scheme[numval]
                        for numval in self.class_map.keys()
                    },
                )

            else:
                var_visualization = {
                    "bands": [self.band_of_use],
                    "min": vmin,
                    "max": vmax,
                }

        else:
            if self.metadf.empty:
                if vmin == None:
                    vmin = 0.0
                if vmax == None:
                    vmax = 1.0
            else:
                # continuous color scheme
                obsmin = df[self.name].min()
                obsmax = df[self.name].max()
                if vmin == None:
                    vmin = obsmin - ((obsmax - obsmin) * 0.15)
                if vmax == None:
                    vmax = obsmax + ((obsmax - obsmin) * 0.15)

            var_visualization = {
                "bands": [self.band_of_use],
                "min": vmin,
                "max": vmax,
                "palette": [
                    "000080",
                    "0000d9",
                    "4000ff",
                    "8000ff",
                    "0080ff",
                    "00ffff",
                    "00ff80",
                    "80ff00",
                    "daff00",
                    "ffff00",
                    "fff500",
                    "ffda00",
                    "ffb000",
                    "ffa400",
                    "ff4f00",
                    "ff2500",
                    "ff0a00",
                    "ff00ff",
                ],
            }

            MAP.add_colorbar_branca(
                vis_params=var_visualization,
                index=None,
                label="",
                categorical=False,
                step=None,
                background_color=None,
            )

        # add layer
        MAP.add_layer(
            ee_object=im,
            vis_params=var_visualization,
            name=self.name,
        )

        # add layer control
        MAP.addLayerControl()

        # Save
        if save:
            logger.info(f"Saving {self.name} gee plot at: {full_path}")
            MAP.save(full_path)

        return MAP


class GeeDynamicModelData(_GeeModelData):
    """Class for working with Dynamic GEE modeldatasets."""

    def __init__(
        self,
        name,
        location,
        value_type,
        scale,
        time_res,
        modelobstypes,
        is_image=False,
        is_mosaic=False,
        credentials="",
        class_map={},
        agg_scheme={},
        col_scheme={},
    ):
        """
        Create a GeeDynamicModelData instance representing a GEE dataset with
        a time dimension.

        Parameters
        ----------
        name : str
            The user-defined name for referring to this GEE dataset.
        location : str
            The location of the dataset on GEE. Navigate the GEE datasets online,
            to get the location of a GEE dataset.
        value_type : "numeric" or "categorical"
            Specify how to interpret the values of the GEE dataset.
        scale : int
            The Scale (See GEE doc) of the dataset to extract values of. This
            can be found on the GEE dataset information page.
        time_res : str
            The time resolution of the dataset is represented as a timedelta string.
            Common resolutions are '1h' or '3h'. This can be found on the GEE
            dataset information page.
        modelobstypes : list of ModelObstype and ModelObstype_Vectorfield
            The ModelObstype's (scalar and vector) defined for this GEE dataset. These ModelObstype
            define how Obstypes are linked to bands of the GEE dataset.
        is_image : bool
            If True, the GEE dataset is opened as ee.Image(), else
            ee.ImageCollection(). This can be found on the GEE dataset
            information page.
        is_mosaic : bool, optional
            If True, ee.mosaic() is applied on the GEE dataset. The default is False.
        credentials : str, optional
            Credentials of the GEE dataset. The default is "".
        class_map : dict, optional
            If value_type is categorical, than the class_map defines how the
            numeric values are mapped to 'human-categories'. The keys are the
            numeric values, the values are the human-labels. The default is {}.
        agg_scheme : dict, optional
            If value_types is categorical, then the agg scheme defines custom-
            made classes, which are aggregates of the present classes. The
            keys are the names of the custom-classes, the values are lists, with
            the corresponding numeric values. The default is {}.
        col_scheme : dict, optional
            if value_types is categorical, the col_sheme defines the colors used
            for each class. The keys are the numeric values, the values are
            the colors (in hex form). The default is {}.

        Returns
        -------
        None.

        See Also
        -----------
        GeeStaticModelData : Gee Modeldata dataset without time dimension.
        metobs_toolkit.ModelObstype: A Obstype for Modeldata (linked to a band).
        metobs_toolkit.ModelObstype_Vectorfield: A Vectorfield version of ModelObstype.
        GeeDynamicModelData.set_metadf: Set metadata (station locations).
        GeeDynamicModelData.get_info: Print out detailed info method.
        GeeDynamicModelData.extract_timeseries_data: Extract data from GEE.
        GeeDynamicModelData.save_modeldata: Save Modeldata as pickle.
        metobs_toolkit.import_modeldata_from_pkl: Import Modeldata from pickle.


        Note
        -------
        In general, specifying a scale smaller than the true scale of the GEE
        dataset has no impact on the results (but can affect the computation time).

        """
        super().__init__(
            name=name,
            location=location,
            value_type=value_type,
            scale=scale,
            is_static=False,
            is_image=is_image,
            is_mosaic=is_mosaic,
            credentials=credentials,
        )
        # name - datetime index, and (tlk-renamed) bandnames as columns (=wide since perfect freq assumed)
        self.modeldf = pd.DataFrame()  # will hold time-related records (timeseries)

        self.modelobstypes = {}
        for obs in modelobstypes:
            if not (
                (isinstance(obs, ModelObstype))
                | (isinstance(obs, ModelObstype_Vectorfield))
            ):
                raise MetobsModelDataError(
                    f"{obs} is not an instance of ModelObstype or ModelObstype_Vectorfield but of type {type(obs)}."
                )
            self.modelobstypes[obs.name] = obs

        self.time_res = str(time_res)

        self.__name__ = "GeeDynamicModelData"

    def __str__(self):
        if self.modeldf.empty:
            return f"Empty {self.__name__} instance of {self.name} "
        else:
            return f"{self.__name__} instance of {self.name} with modeldata "

    def __repr__(self):
        """Print overview information of the modeldata."""
        return self.__str__()

    def _clear_data(self):
        """Clear all attributes that holds meta and extracted data."""
        self._clear_metadata()
        self.modeldf = pd.DataFrame()

    # =============================================================================
    # Setters
    # =============================================================================

    def _set_modeldf(self, modeldf):
        """Set the modeldf attribute"""

        if not list(modeldf.index.names) == ["name", "datetime"]:
            raise MetobsModelDataError()(
                f"A dataframe is being set as .modeldf with wrong modeldf.index.names: {modeldf.index.names}"
            )
        self.modeldf = modeldf

    def add_modelobstype(self, modelobstype):
        """Add a new ModelObstype to the GeeDynamicModelData.


        Parameters
        ----------
        modelobstype : ModelObstype or ModelObstype_Vectorfield
            The new modelobstype to add.

        Returns
        -------
        None.

        See Also
        -----------
        metobs_toolkit.ModelObstype: A Obstype for Modeldata (linked to a band).
        metobs_toolkit.ModelObstype_Vectorfield: A Vectorfield version of ModelObstype.

        Examples
        --------
        As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        The ERA5 GeeDyanmicModelData is equipped with a few Modelobstypes.
        >>> era5.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

        Now we create a ModelObstype, which is a regular Obstype, but with
        extra information on how this is linked to a band in the Gee dataset. So
        we start by looking op the gee dataset online and inspect the present
        bands: https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY#bands

        As example, we add the 'surface_solar_radiation_downwards' band. Therefore
        we need to create an Obstype representing solar radiation (downward) first.

        >>> solar_rad = metobs_toolkit.Obstype(obsname='solar_radiation',
        ...                                    std_unit='J/m²',
        ...                                    description= 'Downward solar flux (perpendicular to earth surface)',
        ...                                    unit_aliases={},
        ...                                    unit_conversions={})
        >>> solar_rad
        Obstype instance of solar_radiation

        Now we make a ModelObstype from it. We know that the values of solar_rad
        are scalar, so we create a `ModelObstype` instance (for vectorfields,
        use the `ModelObstype_Vectorfield` class.)

        >>> solar_rad_in_era5 = metobs_toolkit.ModelObstype(obstype=solar_rad,
        ...                                                 model_unit='J/m²',
        ...                                                 model_band='surface_solar_radiation_downwards')
        >>> solar_rad_in_era5
        ModelObstype instance of solar_radiation (linked to band: surface_solar_radiation_downwards)

        At last, we can add the solar radiation to the Geedataset by using
        the `GeeDynamicModelData.add_modelobstype()` method.

        >>> era5.add_modelobstype(modelobstype=solar_rad_in_era5)

        >>> era5.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m), 'solar_radiation': ModelObstype instance of solar_radiation (linked to band: surface_solar_radiation_downwards)}

        """
        if not (
            (isinstance(modelobstype, ModelObstype))
            | (isinstance(modelobstype, ModelObstype_Vectorfield))
        ):
            raise MetobsModelDataError(
                f"{modelobstype} is not a ModelObstype of ModelObstype_Vectorfield"
            )
        if modelobstype.name in self.modelobstypes.keys():
            # Check equlity
            if modelobstype == self.modelobstypes[modelobstype.name]:
                return  # modelobstye is already present
            else:
                raise MetobsModelDataError(
                    f"There is already a known ModelObstype with {modelobstype.name} as a name: {self.modelobstypes[modelobstype.name]}"
                )
        else:
            self.modelobstypes[modelobstype.name] = modelobstype

    # =============================================================================
    # Convertors/formatters
    # =============================================================================
    def _convert_units(self):
        """Convert the units of the modeldf attr to toolkit space, and set the
        modeldf again.

        """

        modeldf = self.modeldf
        for obs in self.modelobstypes.values():
            if obs.name in modeldf.columns:
                modeldf[obs.name] = obs.convert_to_standard_units(
                    input_data=modeldf[obs.name], input_unit=obs.get_modelunit()
                )
        self._set_modeldf(modeldf)
        return

    def _format_gee_df_structure(self, geedf):
        """Format a dataframe (constructed directly from gee), to a modeldf (
        name, datetime index and modelobstype names as columns). Then set the
        modeldf attribute.
        """

        # format datetime
        geedf["datetime"] = pd.to_datetime(geedf["datetime"], format="%Y%m%d%H%M%S")
        # (assume all gee dataset are in UTC)
        geedf["datetime"] = geedf["datetime"].dt.tz_localize("UTC")

        # 2. Format dataframe
        # format index
        geedf = geedf.set_index(["name", "datetime"])
        geedf = geedf.sort_index()

        # 3. Look for vector components, add the columns in the geedf and add the
        # scalar equivalent obstypes
        new_obs = []
        for obs in self.modelobstypes.values():
            if isinstance(obs, ModelObstype_Vectorfield):
                if (obs.get_modelband_u() in geedf.columns) & (
                    obs.get_modelband_v() in geedf.columns
                ):
                    # 1 add aplitude column + obstype
                    amp_series, amp_obstype = obs.compute_amplitude(df=geedf)
                    geedf[amp_obstype.name] = amp_series
                    new_obs.append(amp_obstype)

                    # 2. add direction column + obstype
                    dir_series, dir_obstype = obs._compute_angle(df=geedf)
                    geedf[dir_obstype.name] = dir_series
                    new_obs.append(dir_obstype)

        for obs in new_obs:
            self.add_modelobstype(obs)

        # rename to values to toolkit space (only for scalar fields)
        scalar_obs = [
            obs for obs in self.modelobstypes.values() if isinstance(obs, ModelObstype)
        ]
        band_mapper = {obs.get_modelband(): obs.name for obs in scalar_obs}
        geedf = geedf.rename(columns=band_mapper)

        modeldf = geedf
        self._set_modeldf(modeldf)

    def _subset_to_obstypes(self, trg_obstypes):
        """Subset the modeldf (atttr) to list of target obstypes."""
        keep_columns = []

        for obs in trg_obstypes:
            if isinstance(obs, ModelObstype):
                keep_columns.append(obs.name)
            elif isinstance(obs, ModelObstype_Vectorfield):
                keep_columns.append(obs._amp_obs_name)
                keep_columns.append(obs._dir_obs_name)
            else:
                raise MetobsModelDataError(
                    f"{obs} is not a ModelObstype or ModelObstype_Vectorfield."
                )

        self.modeldf = self.modeldf[keep_columns]

    # =============================================================================
    # Getters
    # =============================================================================

    def _get_bandnames(self, trg_obstypes):
        """Get a list of all known target band names."""

        trg_bands = []
        for obs in trg_obstypes:
            if isinstance(obs, ModelObstype):
                trg_bands.append(obs.get_modelband())
            elif isinstance(obs, ModelObstype_Vectorfield):
                trg_bands.append(obs.get_modelband_u())
                trg_bands.append(obs.get_modelband_v())
            else:
                raise MetobsModelDataError(
                    f"{obs} is not an instance of ModelObstype or ModelObstype_Vectorfield."
                )
        return trg_bands

    def _get_time_res(self):
        return str(self.time_res)

    def _get_unknown_modelcolumns(self):
        """Get a list of all columns in the modeldf, that are unknown obstypes."""
        return list(set(self.modeldf.columns) - set(self.modelobstypes.keys()))

    def get_info(self):
        """Print out detailed information about the GeeDynamicModelData.

        Returns
        -------
        None.

        See Also
        -----------
        GeeDynamicModelData.set_metadf: Set metadata (station locations).
        GeeDynamicModelData.add_modelobstype: Add a new ModelObstype.
        GeeDynamicModelData.extract_timeseries_data: Extract data from GEE.

        Examples
        --------
        As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()`

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        To get a detailed overview, use the `GeeDynamicModelData.get_info()`
        method, to print out the info.

        >>> era5.get_info()
        Empty GeeDynamicModelData instance of ERA5-land
        ------ Details ---------
        <BLANKLINE>
         * name: ERA5-land
         * location: ECMWF/ERA5_LAND/HOURLY
         * value_type: numeric
         * scale: 2500
         * is_static: False
         * is_image: False
         * is_mosaic: False
         * credentials:
         * time res: 1h
        <BLANKLINE>
         -- Known Modelobstypes --
        <BLANKLINE>
         * temp : ModelObstype instance of temp (linked to band: temperature_2m)
            (conversion: Kelvin --> Celsius)
         * pressure : ModelObstype instance of pressure (linked to band: surface_pressure)
            (conversion: pa --> pa)
         * wind : ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)
            vectorfield that will be converted to:
              * wind_speed
              * wind_direction
            (conversion: m/s --> m/s)
        <BLANKLINE>
         -- Metadata --
        <BLANKLINE>
        No metadf is set.
        <BLANKLINE>
         -- Modeldata --
        <BLANKLINE>
        No model data is set.

        """
        print(str(self))
        self._get_base_details()
        print(f" * time res: {self.time_res}")
        print("\n -- Known Modelobstypes -- \n")
        for obs in self.modelobstypes.values():
            print(f" * {obs.name} : {obs}")
            if isinstance(obs, ModelObstype_Vectorfield):
                print("    vectorfield that will be converted to: ")
                print(f"      * {obs._amp_obs_name}")
                print(f"      * {obs._dir_obs_name}")
            print(
                f"    (conversion: {obs.get_modelunit()} --> {obs.get_standard_unit()})"
            )

        print("\n -- Metadata -- \n")
        if self.metadf.empty:
            print("No metadf is set.")
        else:
            print(f"{self.metadf}")

        print("\n -- Modeldata -- \n")

        if self.modeldf.empty:
            print("No model data is set.")
        else:
            print(f"{self.modeldf}")

        if bool(self._get_unknown_modelcolumns()):
            print(
                "\n (The following data columns (bandnames) are present, without a corresponding Modelobstype: "
            )
            print(f"{self._get_unknown_modelcolumns()}")

    def make_gee_plot(
        self,
        timeinstance,
        modelobstype="temp",
        outputfolder=None,
        filename=None,
        vmin=None,
        vmax=None,
        overwrite=False,
    ):
        """Make an interactive spatial plot of the GEE dataset and the stations.

        This method will create an interactive plot of the GEE dataset at an
        instance in time. If metadata is present, it will be displayed as
        markers on the map.

        The interactive map can be saved as an HTML file, by specifying the
        target path.


        Parameters
        ----------
        timeinstance : datetime.datetime or pandas.Timestamp
            The timeinstance to plot the GEE dataset. This timestamp is
            rounded down with the time resolution (.time_res). The timeinstance,
            is interpreted as UTC.
        modelobstype : str, optional
            The name of the ModelObstype to plot. The modelobstype name must be
            known. The default is "temp".
        outputfolder : str or None, optional
            Path to the folder to save the HTML file. If None, the map will
            not be saved as an HTML file. The default is None.
        filename : str or None, optional
            The filename for the HTML file. If a filename is given, if it does
            not end with ".html", the prefix is added. If None, the map will not
            be saved as an HTML file. The default is None.
        vmin : num or None, optional
            vmin is the minimum value assigned to the colormap. If None, and
            metadata is set, vmin is computed by computing the minimum
            modelvalue in a boundbox defined by the locations of the stations.
            If no metadata is available, and vmin is None then vmin is set to
            0. The default is None.
        vmax : num or None, optional
            vmax is the minimum value assigned to the colormap. If None, and
            metadata is set, vmax is computed by computing the minimum
            modelvalue in a boundbox defined by the locations of the stations.
            If no metadata is available, and vmax is None then vmax is set to
            1. The default is None.
        overwrite : bool, optional
            If True, and if the target file exists, then it will be overwritten.
            Else, an error will be raised when the target file already exists.
            The default is False.

        Returns
        -------
        MAP : geemap.foliummap.Map
            The interactive map of the GeeStaticModelData.

        See Also
        -----------
        metobs_toolkit.connect_to_gee: Establish connection with GEE services.
        GeeDynamicModelData.set_metadf: Set metadata (station locations).
        GeeDynamicModelData.get_info: Print out detailed info method.
        GeeDynamicModelData.make_plot: Make a timeseries plot.

        Warning
        ---------
        To display the interactive map a graphical interactive backend is
        required, which could be missing. You can recognice this when no
        map is displayed, but the Python console prints out a message similar
        to `<geemap.foliumap.Map at 0x7ff7586b8d90>`.

        In that case, you can specify an `outputfolder` and `outputfile`, save the map as an HTML file, and
        open it with a browser.

        Note
        ------
        Be aware that it is not possible to plot a ModelObstype_Vectorfield.

        Examples
        --------
        As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()` and we will plot the 2m-temperature field.

        >>> import metobs_toolkit
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        For illustration we will add metadata (station locations) to the era5
        Modeldata, so that the locations appear on the map. We use the demo
        dataset's metadata for this.

        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)
        >>> era5.set_metadf(dataset.metadf)

        We can no use the `GeeDynamicModelData.make_gee_plot()` method to make
        an interactive spatial plot of the era5 Modeldata (and since we added
        metadata, the locations of the stations are added as an extra layer).

        We will save the output as a (HTML) file and store it in the
        current working directory (`os.getcwd`) as illustration.

        We specify a time instance that is present in the dataset (see gee
        dataset info online). The time instance is rounded down, respecting
        the time resolution of the GEE dataset.

        >>> import os
        >>> import datetime
        >>> dt = datetime.datetime(2006,11,18, 20, 15)
        >>> str(dt)
        '2006-11-18 20:15:00'

        >>> era5.make_gee_plot(
        ...     timeinstance=dt, #will be rounded down to 18/11/2006 20:00:00
        ...     modelobstype="temp", #which modelobstype to plot
        ...     outputfolder=os.getcwd(),
        ...     filename=f'era5_temp_at_{dt}.html',
        ...     vmin=None, #because metadata is present, a vmin will be computed
        ...     vmax=None, #because metadata is present, a vmax will be computed
        ...     overwrite=True,)
        <geemap.foliumap.Map object at ...



        """

        # Format datetime
        if isinstance(timeinstance, datetimemodule.datetime):
            dt = pd.Timestamp(timeinstance)
        elif isinstance(timeinstance, pd.Timestamp):
            dt = timeinstance
        else:
            raise MetobsModelDataError(
                f"{timeinstance} is not in a datetime format (datetime.datetime or pandas.Timestamp)."
            )
        # check timezone and make tz-awer
        if dt.tz is None:
            # tz naive timestamp
            dt = dt.tz_localize(tz="UTC")
        else:
            # tz aware timestamp --> convert to tz of records
            dt = dt.tz_convert(tz="UTC")

        # floor to resolution
        dt = dt.floor(self.time_res)

        # Check modelobstype
        if modelobstype not in self.modelobstypes:
            raise MetobsModelDataError(
                f"{modelobstype} is not a known modelobstype ({self.modelobstypes})"
            )
        modelobstype = self.modelobstypes[modelobstype]

        if not isinstance(modelobstype, ModelObstype):
            raise MetobsModelDataError(f"{modelobstype} is not a ModelObstype.")

        # check if outputfile exists
        if (outputfolder is not None) | (filename is not None):
            save = True

            if not os.path.isdir(outputfolder):
                raise MetobsModelDataError(f"{outputfolder} is not a directory!")

            # check file extension in the filename:
            if filename[-5:] != ".html":
                filename += ".html"

            full_path = os.path.join(outputfolder, filename)

            # check if file exists
            if os.path.isfile(full_path):
                if overwrite:
                    logger.info(f"Overwrite the file at {full_path}.")
                    os.remove(full_path)
                else:
                    raise MetobsModelDataError(f"{full_path} is already a file!")
        else:
            save = False

        # Connect to Gee
        connect_to_gee()

        # Get image
        im = gee_api.get_ee_obj(self, target_bands=[modelobstype.get_modelband()])
        # filter to timestamp
        im = im.filterDate(
            start=dt.isoformat(), end=(dt + pd.Timedelta(self.time_res)).isoformat()
        )
        im = im.first()

        # make empty map (FOLIUM BACKEND !! )
        MAP = folium_map()

        # show stations
        if self.metadf.empty:
            vmin = 0.0
            vmax = 1.0

        else:
            _add_stations_to_folium_map(
                Map=MAP, metadf=self.metadf, display_cols=["name"]
            )
            # fix center
            centroid = self.metadf.dissolve().centroid
            MAP.setCenter(lon=centroid.x.item(), lat=centroid.y.item(), zoom=8)

            # Create GEE boundbox of the ROI

            metadf = self.metadf.to_crs("epsg:4326")
            (xmin, ymin, xmax, ymax) = metadf.total_bounds

            if (vmin is None) | (vmax is None):
                roi = ee.Geometry.BBox(west=xmin, south=ymin, east=xmax, north=ymax)
                roi_min = im.reduceRegion(
                    ee.Reducer.min(), roi, scale=self.scale
                ).getInfo()[modelobstype.get_modelband()]
                roi_max = im.reduceRegion(
                    ee.Reducer.max(), roi, scale=self.scale
                ).getInfo()[modelobstype.get_modelband()]
            if vmin is None:
                vmin = roi_min - ((roi_max - roi_min) * 0.15)
                if vmin == vmax:
                    vmin = vmax - 1.0
            if vmax is None:
                vmax = roi_max + ((roi_max - roi_min) * 0.15)
                if vmax == vmin:
                    vmax = vmin + 1.0

        var_visualization = {
            "bands": [modelobstype.get_modelband()],
            "min": vmin,
            "max": vmax,
            "palette": [
                "000080",
                "0000d9",
                "4000ff",
                "8000ff",
                "0080ff",
                "00ffff",
                "00ff80",
                "80ff00",
                "daff00",
                "ffff00",
                "fff500",
                "ffda00",
                "ffb000",
                "ffa400",
                "ff4f00",
                "ff2500",
                "ff0a00",
                "ff00ff",
            ],
        }

        MAP.add_colorbar_branca(
            vis_params=var_visualization,
            index=None,
            label="",
            categorical=False,
            step=None,
            background_color=None,
        )

        # add layer
        MAP.add_layer(
            ee_object=im,
            vis_params=var_visualization,
            name=f"{self.name}: {modelobstype.get_modelband()}",
        )

        # add title
        title = f"{self.name} GEE plot of {modelobstype.get_modelband()} (in {modelobstype.get_modelunit()}) at {dt}."
        MAP = add_title_to_folium_map(title=title, Map=MAP)

        # add layer control
        MAP.addLayerControl()

        # Save
        if save:
            logger.info(f"Saving {self.name} gee plot at: {full_path}")
            MAP.save(full_path)

        return MAP

    def make_plot(
        self,
        obstype_model="temp",
        Dataset=None,
        obstype_dataset=None,
        stationnames=None,
        starttime=None,
        endtime=None,
        title=None,
        show_outliers=True,
        show_filled=True,
        legend=True,
        _ax=None,  # needed for GUI, not recommended use
    ):
        """Plot timeseries of the modeldata.

        This function creates a timeseries plot for the Modeldata. When a
        metobs_toolkit.Dataset is provided, it is plotted in the same figure.

        The line colors represent the timeseries for different locations.



        Parameters
        ----------
        obstype_model : string, optional
            The name of the ModelObstype to plot. The default is 'temp'.
        Dataset : metobs_toolkit.Dataset, optional
            A Dataset instance with observations plotted in the same figure.
            Observations are represented by solid lines and modeldata by dashed
            lines. The default is None.
        obstype_dataset : string, optional
            Fieldname of the Dataset to visualize. Only relevant when a dataset
            is provided. If None, obsype_dataset = obstype_model. The default
            is None.
        stationnames : list, optional
            A list with stationnames to include in the timeseries. If None is
            given, all the stations are used, defaults to None.
        starttime : datetime.datetime, optional
             Specify the start datetime for the plot. If None is given it will
             use the start datetime of the dataset. The default to None.
        endtime : datetime.datetime, optional
             Specify the end datetime for the plot. If None is given it will
             use the end datetime of the dataset. The default to None.
        title : string, optional
             Title of the figure, if None a default title is generated. The
             default is None.
        show_outliers : bool, optional
             If true the observations labeled as outliers will be included in
             the plot. Only relevant when a dataset is provided. The default
             is True.
        show_filled : bool, optional
             If true the filled values for gaps and missing observations will
             be included in the plot. Only relevant when a dataset is provided.
             The default is True.
        legend : bool, optional
             If True, a legend is added to the plot. The default is True.


        Returns
        -------
        axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

        See Also
        -----------
        metobs_toolkit.connect_to_gee: Establish connection with GEE services.
        GeeDynamicModelData.set_metadf: Set metadata (station locations).
        GeeDynamicModelData.extract_timeseries_data: Extract data from GEE.
        GeeDynamicModelData.get_info: Print out detailed info method.
        GeeDynamicModelData.make_gee_plot: Make interactive spatial GEE plot.

        Examples
        --------

        .. plot::
            :context: close-figs

            As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
            present in a `metobs_toolkit.Dataset()`. We will use the demo dataset,
            extract era5 timeseries data for the stations in the demo dataset, and
            plot the temperature.

            >>> import metobs_toolkit
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                input_data_file=metobs_toolkit.demo_datafile,
            ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                template_file=metobs_toolkit.demo_template)
            >>> era5 = dataset.gee_datasets['ERA5-land']
            >>> era5.set_metadf(dataset.metadf)
            >>> era5
            Empty GeeDynamicModelData instance of ERA5-land

            Now we will extract temperature timeseries from ERA5. We already have
            a ModelObstype representing the temperature present in the GeeDynamicModelData:

            >>> era5.modelobstypes
            {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

            We define a (short) period and extract ERA5 timeseries.

            >>> import datetime
            >>> tstart = datetime.datetime(2022,9,4)
            >>> tend = datetime.datetime(2022,9,5)
            >>> era5.extract_timeseries_data(
            ...                            startdt_utc=tstart,
            ...                            enddt_utc=tend,
            ...                            obstypes=['temp'])

            If the data request is small, the timeseries are present. (If the
            data request is larger, a CSV file is written to your google drive.
            Download that file, and use the `GeeDynamicModelData.set_modeldata_from_csv()`
            method.)

            >>> era5.modeldf
                                                  temp
            name      datetime
            vlinder01 2022-09-04 00:00:00+00:00  17.71
                      2022-09-04 01:00:00+00:00  17.56
                      2022-09-04 02:00:00+00:00  17.24
                      2022-09-04 03:00:00+00:00  16.75
                      2022-09-04 04:00:00+00:00  16.33
            ...                                    ...
            vlinder28 2022-09-04 20:00:00+00:00  21.40
                      2022-09-04 21:00:00+00:00  20.91
                      2022-09-04 22:00:00+00:00  20.51
                      2022-09-04 23:00:00+00:00  20.28
                      2022-09-05 00:00:00+00:00  20.13
            <BLANKLINE>
            [700 rows x 1 columns]


            (As you can see, the timeseries are automatically converted to the
             toolkit standards of the Obstype: °C)

            Now we can plot the timeseries. We can plot the observations in the
            same figure as well.

            >>> era5.make_plot(
            ...     obstype_model="temp", #plot temperature timeseries of the model (=era5)
            ...     Dataset=dataset,
            ...     obstype_dataset='temp') #plot temperature observations
            <Axes: title={'center': 'ERA5-land and temp observations.'}, ylabel='temp ...

        """
        logger.info(f"Make {obstype_model}-timeseries plot of model data")

        # Basic test
        if obstype_model not in self.modeldf.columns:
            raise MetobsModelDataError(
                f"{obstype_model} is not foud in the modeldata df (columns = {self.modeldf.columns})."
            )

        if self.modeldf.empty:
            logger.warning("The modeldata is empty.")
            return
        if obstype_dataset is None:
            obstype_dataset = obstype_model

        if Dataset is not None:
            if obstype_dataset not in Dataset._get_present_obstypes():
                logger.warning(f"{obstype_dataset} is not foud in the Dataframe df.")
                return

        model_df = self.modeldf

        # ------ filter model ------------

        # Filter on obstype
        model_df = model_df[[obstype_model]]

        # Subset on stationnames
        if stationnames is not None:
            model_df = model_df[
                model_df.index.get_level_values("name").isin(stationnames)
            ]

        # Subset on start and endtime
        model_df = multiindexdf_datetime_subsetting(model_df, starttime, endtime)

        #  -------- Filter dataset (if available) -----------
        if Dataset is not None:
            # check if there is data
            Dataset._data_is_required_check()
            # combine all dataframes
            mergedf = Dataset.get_full_status_df(return_as_wide=False)

            # subset to obstype
            mergedf = xs_save(mergedf, obstype_dataset, level="obstype")

            # Subset on stationnames
            if stationnames is not None:
                mergedf = mergedf[
                    mergedf.index.get_level_values("name").isin(stationnames)
                ]

            # Subset on start and endtime
            mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

        # Generate ylabel
        y_label = self.modelobstypes[obstype_model]._get_plot_y_label()

        # Generate title
        title = f"{self.name}"
        if Dataset is not None:
            title = (
                f"{title} and {Dataset.obstypes[obstype_dataset].name} observations."
            )

        # make plot of the observations
        if Dataset is not None:
            # make plot of the observations
            _ax, col_map = timeseries_plot(
                mergedf=mergedf,
                title=title,
                ylabel=y_label,
                colorby="name",
                show_legend=legend,
                show_outliers=show_outliers,
                show_filled=show_filled,
                settings=Dataset.settings,
                _ax=_ax,
            )

            # Make plot of the model on the previous axes
            ax, col_map = model_timeseries_plot(
                df=model_df,
                obstype=obstype_model,
                title=title,
                ylabel=y_label,
                show_primary_legend=False,
                add_second_legend=True,
                _ax=_ax,
                colorby_name_colordict=col_map,
            )

        else:
            # Make plot of model on empty axes
            ax, _colmap = model_timeseries_plot(
                df=model_df,
                obstype=obstype_model,
                title=title,
                ylabel=y_label,
                show_primary_legend=legend,
                add_second_legend=False,
                _ax=_ax,
            )

        return ax

    # =============================================================================
    # Functionallity
    # =============================================================================

    def extract_timeseries_data(
        self,
        startdt_utc,
        enddt_utc,
        obstypes=["temp"],
        get_all_bands=False,
        drive_filename=None,
        drive_folder="gee_timeseries_data",
        force_direct_transfer=False,
        force_to_drive=False,
    ):
        """Extract timeseries data and set the modeldf.

        The timeseries are extracted at the location of the stations. Nan's are
        used when there is no data at that location (out of data-mask).

        There are two possibilities for extracting the timeseries:
            * If the data request is not too big, then the timeseries are
              direct imported. No additional steps are required. The Modeldata
              is converted to Metobs-toolkit standards (name and units), and
              amplitude/direction fields are computed for requested vectorfields.

            * If the data request is to big, GEE services will write a datafile
              (CSV) directly to your Google Drive. Wait until the file appears
              on your Drive, then download it. You can import the data by using
              the `GeeDynamicModelData.set_modeldata_from_csv()` method.

        Parameters
        ----------
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstypes : list of strings, optional
            A list of ModelObstype names to extract modeldata for. These obstypes
            must be known. The default is ['temp']
        get_all_bands : bool, optional
            If True, all values (over all bands) are extracted. If the band is
            linked to a ModelObstye, then the name of the modelObstype is used
            instead of the band name. If True, the obstypes argument is ignored.
            The default is False.
        drive_filename : str or None, optional
            If given, the data will be saved as this filename on your Google Drive.
            This argument will only take effect when the data is written to
            Google Drive. If None, a custom filename is created. The default is
            None.
        drive_folder: str
            The name of the folder on your Google Drive to save the drive_filename
            in. If the folder, does not exists it will be created instead. This
            argument will only take effect when the data is written to Google
            Drive.
        force_direct_transfer: bool, optional
            If True, the data is demanded as a direct transfer (no file writen
            to google drive). If the request is to large, an GEE error is raised.
            The default is False.
        force_to_drive: bool, optional
            If True, The gee data is writen to a file on your drive. Direct
            transfer of data is prohibited. The default is False.

        Returns
        -------
        None.

        See Also
        -----------
        metobs_toolkit.connect_to_gee: Establish connection with GEE services.
        GeeDynamicModelData.set_metadf: Set metadata (station locations).
        GeeDynamicModelData.get_info: Print out detailed info method.
        GeeDynamicModelData.make_plot: Make a timeseries plot.

        Note
        -------
        Make sure that the metadata is set. Use the
        `GeeDynamicModelData.set_metadata()` for this.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your Google Drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        Examples
        --------
        As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()`. We will extract era5 timeseries
        data for the stations in the demo dataset.

        >>> import metobs_toolkit
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5.set_metadf(dataset.metadf)
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        Now we will extract temperature and wind timeseries from ERA5. We already have
        a ModelObstype representing the temperature present in the GeeDynamicModelData:

        >>> era5.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

        Note that the "wind" Modelobstype is a ModelObstype_Vectorfield!

        >>> era5.modelobstypes['wind']
        ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)

        When "wind" is requested for ERA5, the toolkit will request the u and v
        components of the wind, download them, convert the units and compute
        the amplitude and direction (scalar) values. At last, two new ModelObstypes
        are (automatically) created for the amplitude and direction values.


        We define a (short) period and extract ERA5 timeseries.

        >>> import datetime
        >>> tstart = datetime.datetime(2022,9,4)
        >>> tend = datetime.datetime(2022,9,5)
        >>> era5.extract_timeseries_data(
        ...                            startdt_utc=tstart,
        ...                            enddt_utc=tend,
        ...                            obstypes=['temp', 'wind'])

        If the data request is small, the timeseries are present. (If the
        data request is larger, a CSV file is writen to your google drive.
        Download that file, and use the `GeeDynamicModelData.set_modeldata_from_csv()`
        method.)

        >>> era5.modeldf
                                              temp  wind_speed  wind_direction
        name      datetime
        vlinder01 2022-09-04 00:00:00+00:00  17.71        2.44          359.39
                  2022-09-04 01:00:00+00:00  17.56        2.47           11.91
                  2022-09-04 02:00:00+00:00  17.24        2.28           19.15
                  2022-09-04 03:00:00+00:00  16.75        2.22           22.71
                  2022-09-04 04:00:00+00:00  16.33        2.22           22.72
        ...                                    ...         ...             ...
        vlinder28 2022-09-04 20:00:00+00:00  21.40        2.37          251.99
                  2022-09-04 21:00:00+00:00  20.91        2.47          270.47
                  2022-09-04 22:00:00+00:00  20.51        2.50          298.30
                  2022-09-04 23:00:00+00:00  20.28        2.76          322.41
                  2022-09-05 00:00:00+00:00  20.13        2.97          331.08
        <BLANKLINE>
        [700 rows x 3 columns]



        As you can see, the timeseries are automatically converted to the
        toolkit standards of the Obstype.

        We can also see that the amplitude and direction of the windfield,
        are added as ModelObstypes.

        >>> era5.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m), 'wind_speed': ModelObstype instance of wind_speed (linked to band: wind_speed), 'wind_direction': ModelObstype instance of wind_direction (linked to band: wind_direction)}


        """
        # ====================================================================
        # Test input
        # ====================================================================
        self._check_metadf_validity(self.metadf)

        if (force_direct_transfer) & (force_to_drive):
            raise MetobsModelDataError(
                "Both force_direct_transfer and force_to_drive could not be True at the same time."
            )
        # Check obstypes
        if isinstance(obstypes, str):
            obstypes = [obstypes]  # convert to list

        for obstype in obstypes:
            # is obstype mapped?
            if obstype not in self.modelobstypes.keys():
                raise MetobsModelDataError(
                    f"{obstype} is an unknown modelobservation type of the {self}."
                )
        # convert to model obstypes
        obstypes = [self.modelobstypes[obs] for obs in obstypes]

        # convert timestamps to pd.Timestamps with utc as timezone if unaware
        startdt_utc = pd.Timestamp(startdt_utc)
        if startdt_utc.tz is None:
            startdt_utc = startdt_utc.tz_localize(tz="UTC")
        enddt_utc = pd.Timestamp(enddt_utc)
        if enddt_utc.tz is None:
            enddt_utc = enddt_utc.tz_localize(tz="UTC")

        # ====================================================================
        # GEE api extraction
        # ====================================================================

        # Connect to Gee
        connect_to_gee()

        # Get bandnames
        if get_all_bands:
            bandnames = self._get_all_gee_bandnames()
        else:
            bandnames = self._get_bandnames(trg_obstypes=obstypes)

        logger.info(f"{bandnames} are extracted from {self}.")

        _est_data_size = gee_api._estimate_data_size(
            metadf=self.metadf,
            startdt=startdt_utc,
            enddt=enddt_utc,
            time_res=self.time_res,
            n_bands=len(bandnames),
        )

        if force_direct_transfer:
            use_drive = False
        elif force_to_drive:
            use_drive = True

        elif _est_data_size > 4900:
            print(
                "THE DATA AMOUNT IS TOO LARGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
            )
            logger.info(
                "THE DATA AMOUNT IS TOO LARGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
            )

            use_drive = True
        else:
            use_drive = False
        # =============================================================================
        # df to featurecollection
        # =============================================================================

        ee_fc = gee_api._df_to_features_point_collection(self.metadf)

        # =============================================================================
        # extract raster values
        # =============================================================================

        def rasterExtraction(image):
            feature = image.sampleRegions(
                collection=ee_fc,  # feature collection here
                scale=self.scale,  # Cell size of raster
            )
            return feature

        # Because the daterange is maxdate exclusive, add the time resolution to the enddt
        enddt_utc = enddt_utc + pd.Timedelta(self.time_res)

        raster = gee_api.get_ee_obj(self, target_bands=bandnames)  # dataset
        results = (
            raster.filter(
                ee.Filter.date(
                    gee_api._datetime_to_gee_datetime(startdt_utc),
                    gee_api._datetime_to_gee_datetime(enddt_utc),
                )
            )
            .map(gee_api._addDate)
            .map(rasterExtraction)
            .flatten()
        )

        if not use_drive:
            results = results.getInfo()

            # extract properties
            properties = [x["properties"] for x in results["features"]]
            df = pd.DataFrame(properties)

            if df.empty:
                sys.exit("ERROR: the returned timeseries from GEE are empty.")

            self._format_gee_df_structure(df)  # format to wide structure + update attr

            # Subset
            if not get_all_bands:
                self._subset_to_obstypes(trg_obstypes=obstypes)

            # convert units + update attr
            self._convert_units()

        else:
            if drive_filename is None:
                _filename = f"{self.name}_timeseries_data"
            else:
                _filename = str(drive_filename)
            print(
                f"The timeseries will be written to your Drive in {drive_folder}/{_filename} "
            )
            logger.info(
                f"The timeseries will be written to your Drive in {drive_folder}/{_filename} "
            )

            data_columns = ["datetime", "name"]
            data_columns.extend(bandnames)

            task = ee.batch.Export.table.toDrive(
                collection=results,
                description="MetObs_extracting_timeseries",
                folder=drive_folder,
                fileNamePrefix=_filename,
                fileFormat="CSV",
                selectors=data_columns,
            )

            task.start()
            logger.info("The google server is handling your request ...")
            sleep(3)
            finished = False
            while finished is False:
                if task.status()["state"] == "READY":
                    logger.info("Awaitening execution ...")
                    sleep(4)
                elif task.status()["state"] == "RUNNING":
                    logger.info("Running ...")
                    sleep(4)
                else:
                    logger.info("finished")
                    finished = True

            doc_folder_id = task.status()["destination_uris"][0]
            print(
                "The data is transfered! Open the following link in your browser: \n\n"
            )
            print(f"{doc_folder_id} \n\n")
            print(
                "To upload the data to the model, use the Modeldata.set_model_from_csv() method"
            )

            return

    def _modeldf_to_long(self):
        """Convert the modeldf attribute to a long format (name, datetime, obstype) index."""
        longdf = self.modeldf.stack(future_stack=True)
        longdf = (
            longdf.reset_index()
            .rename(columns={"level_2": "obstype", 0: "value"})
            .set_index(["name", "datetime", "obstype"])
        )
        return longdf

    def _interpolate_modeldata(self, to_multiidx, method="time"):
        """Interpolate modeldata in time.

        Interpolate the modeldata timeseries, to a given name-datetime
        multiindex.

        The modeldata will be converted to the timezone of the multiindex.

        If no interpolation can be done, Nan values are used.

        Parameters
        ----------
        to_multiidx : pandas.MultiIndex
            A name - datetime (tz-aware) multiindex to interpolate the
            modeldata timeseries to.
        method : str, optional
            The interpolation method to use. Is passed to
            pandas.DataFrame.interpolate() method. The default is 'time'.

        Returns
        -------
        returndf : pandas.DataFrame
            A dataframe with to_multiidx as an index.
            The values are the interpolated values.

        """
        returndf = init_multiindexdf()

        recordsdf = init_multiindexdf()
        recordsdf.index = to_multiidx
        # iterate over stations check to avoid extrapolation is done per stations
        for sta in recordsdf.index.get_level_values("name").unique():
            sta_recordsdf = xs_save(recordsdf, sta, level="name", drop_level=False)
            sta_moddf = xs_save(self.modeldf, sta, level="name", drop_level=False)

            if sta_moddf.empty:
                logger.warning(f"There are not modeldata records for {sta}!")
                # empyt sta_moddg --> not model data --> empyt return
                mergedf = sta_recordsdf  # overload the records index
                mergedf = mergedf.reindex(
                    columns=list(self.modeldf.columns)
                )  # add all present obstypes as nan columns

            else:
                # convert modeldata to timezone of observations
                sta_moddf = conv_tz_multiidxdf(
                    df=sta_moddf,
                    timezone=sta_recordsdf.index.get_level_values("datetime").tz,
                )

                # check if modeldata is will not be extrapolated !
                if min(sta_recordsdf.index.get_level_values("datetime")) < min(
                    sta_moddf.index.get_level_values("datetime")
                ):
                    logger.warning("Modeldata will be extrapolated")
                if max(sta_recordsdf.index.get_level_values("datetime")) > max(
                    sta_moddf.index.get_level_values("datetime")
                ):
                    logger.warning("Modeldata will be extrapolated")

                # combine model and records
                mergedf = sta_recordsdf.merge(
                    sta_moddf, how="outer", left_index=True, right_index=True
                )

                # reset index for time interpolation
                mergedf = mergedf.reset_index().set_index("datetime").sort_index()

                # interpolate missing modeldata
                mergedf = mergedf.drop(columns=["name"])
                mergedf.interpolate(method=method, limit_area="inside", inplace=True)
                mergedf["name"] = sta

                # convert back to multiindex
                mergedf = (
                    mergedf.reset_index().set_index(["name", "datetime"]).sort_index()
                )
                # filter only records
                mergedf = mergedf.loc[sta_recordsdf.index]

            returndf = pd.concat([returndf, mergedf])
        return returndf

    def save_modeldata(
        self,
        outputfolder,
        filename="saved_modeldata.pkl",
        overwrite=False,
    ):
        """Save a Modeldata instance to a (pickle) file.

        Parameters
        ----------
        outputfolder : str
            The path to the folder to save the file.
        filename : str, optional
            The name of the output file. The default is 'saved_modeldata.pkl'.
        overwrite : bool, optional
            If True, the target file will be overwritten if it should exist.
            The default is False.

        See Also
        -----------
        GeeDynamicModelData.set_modeldata_from_csv: Set modeldata from a csv datafile.
        metobs_toolkit.import_modeldata_from_pkl: Import modeldata from a pkl file.

        Returns
        -------
        None.

        Examples
        --------
        As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()`. We will use the demo dataset,
        extract era5 timeseries data for the stations in the demo dataset, and
        save the GeeDynamicModelData.

        >>> import metobs_toolkit
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5.set_metadf(dataset.metadf)
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        Now we will extract temperature timeseries from ERA5. We already have
        a ModelObstype representing the temperature present in the GeeDynamicModelData:

        >>> era5.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

        We define a (short) period and extract ERA5 timeseries.

        >>> import datetime
        >>> tstart = datetime.datetime(2022,9,4)
        >>> tend = datetime.datetime(2022,9,5)
        >>> era5.extract_timeseries_data(
        ...                            startdt_utc=tstart,
        ...                            enddt_utc=tend,
        ...                            obstypes=['temp'])

        If the data request is small, the timeseries are present. (If the
        data request is larger, a CSV file is writen to your google drive.
        Download that file, and use the `GeeDynamicModelData.set_modeldata_from_csv()`
        method.)

        >>> era5.get_info()
        GeeDynamicModelData instance of ERA5-land with modeldata
        ------ Details ---------
        <BLANKLINE>
         * name: ERA5-land
         * location: ECMWF/ERA5_LAND/HOURLY
         * value_type: numeric
         * scale: 2500
         * is_static: False
         * is_image: False
         * is_mosaic: False
         * credentials:
         * time res: 1h
        <BLANKLINE>
         -- Known Modelobstypes --
        <BLANKLINE>
         * temp : ModelObstype instance of temp (linked to band: temperature_2m)
            (conversion: Kelvin --> Celsius)
         * pressure : ModelObstype instance of pressure (linked to band: surface_pressure)
            (conversion: pa --> pa)
         * wind : ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)
            vectorfield that will be converted to:
              * wind_speed
              * wind_direction
            (conversion: m/s --> m/s)
        <BLANKLINE>
         -- Metadata --
        <BLANKLINE>
                    lon    lat                  geometry
        name
        vlinder01  3.82  50.98  POINT (3.81576 50.98044)
        vlinder02  3.71  51.02   POINT (3.7097 51.02238)
        vlinder03  4.95  51.32  POINT (4.95211 51.32458)
        vlinder04  4.93  51.34  POINT (4.93473 51.33552)
        vlinder05  3.68  51.05  POINT (3.67518 51.05266)
        ...         ...    ...                       ...
        vlinder24  3.57  51.17  POINT (3.57206 51.16702)
        vlinder25  3.71  51.15  POINT (3.70861 51.15472)
        vlinder26  5.00  51.16  POINT (4.99765 51.16176)
        vlinder27  3.73  51.06   POINT (3.72807 51.0581)
        vlinder28  3.77  51.04  POINT (3.76974 51.03529)
        <BLANKLINE>
        [28 rows x 3 columns]
        <BLANKLINE>
         -- Modeldata --
        <BLANKLINE>
                                              temp
        name      datetime
        vlinder01 2022-09-04 00:00:00+00:00  17.71
                  2022-09-04 01:00:00+00:00  17.56
                  2022-09-04 02:00:00+00:00  17.24
                  2022-09-04 03:00:00+00:00  16.75
                  2022-09-04 04:00:00+00:00  16.33
        ...                                    ...
        vlinder28 2022-09-04 20:00:00+00:00  21.40
                  2022-09-04 21:00:00+00:00  20.91
                  2022-09-04 22:00:00+00:00  20.51
                  2022-09-04 23:00:00+00:00  20.28
                  2022-09-05 00:00:00+00:00  20.13
        <BLANKLINE>
        [700 rows x 1 columns]

        We will save the era5 GeeDynamicModelData now as a (pkl) file. As an
        example, we will store it in the current working directory (`os.getcwd()`)

        >>> import os
        >>> era5.save_modeldata(
        ...                outputfolder=os.getcwd(),
        ...                filename="era5_modeldata.pkl",
        ...                overwrite=True)
        Modeldata saved in ...

        If we later want to open the saved era5 Modeldata, we use the
        `GeeDynamicModelData.import_modeldata_from_pkl()` method.

        >>> your_era5 = metobs_toolkit.import_modeldata_from_pkl(
        ...                     folder_path=os.getcwd(),
        ...                     filename="era5_modeldata.pkl")

        """
        # check if outputfolder exists
        if not os.path.isdir(outputfolder):
            raise MetobsModelDataError(f"{outputfolder} is not a directory!")

        # check file extension in the filename:
        if filename[-4:] != ".pkl":
            filename += ".pkl"

        full_path = os.path.join(outputfolder, filename)

        # check if file exists
        if os.path.isfile(full_path):
            if overwrite:
                logger.info(f"Overwrite the file at {full_path}.")
                os.remove(full_path)
            else:
                raise MetobsModelDataError(f"{full_path} is already a file!")

        with open(full_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

        print(f"Modeldata saved in {full_path}")
        logger.info(f"Modeldata saved in {full_path}")

    def set_modeldata_from_csv(self, csvpath):
        """Import timeseries data that is stored in a CSV file.

        The name of the gee dataset the timeseries are coming from must be the
        same as the .name attribute of the Modeldata.


        The timeseries will be formatted and converted to standard toolkit
        units. If band names of known ModelObstype_Vectorfields are detecte,
        the corresponding amplitude and direction fields are computed aswell.

        Parameters
        ----------
        csvpath : str
            Path of the CSV file containing the modeldata timeseries. (This is
            the path to the file that you have downloaded from your Google
            Drive.)

        Returns
        -------
        None.

        See Also
        -----------
        GeeDynamicModelData.save_modeldata: Save modeldata as a pkl file.
        metobs_toolkit.import_modeldata_from_pkl: Import modeldata from a pkl file.

        Examples
        --------
        As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
        present in a `metobs_toolkit.Dataset()`. We will use the demo dataset,
        and extract era5 timeseries data for the stations in the demo dataset.

        >>> import metobs_toolkit
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                input_data_file=metobs_toolkit.demo_datafile,
        ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                template_file=metobs_toolkit.demo_template)
        >>> era5 = dataset.gee_datasets['ERA5-land']
        >>> era5.set_metadf(dataset.metadf)
        >>> era5
        Empty GeeDynamicModelData instance of ERA5-land

        Now we will extract temperature timeseries from ERA5. We already have
        a ModelObstype representing the temperature present in the GeeDynamicModelData:

        >>> era5.modelobstypes
        {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

        We define a period, and extract ERA5 timeseries.

        >>> import datetime
        >>> tstart = datetime.datetime(2022,9,4)
        >>> tend = datetime.datetime(2022,9,14)
        >>> era5.extract_timeseries_data(
        ...                            startdt_utc=tstart,
        ...                            enddt_utc=tend,
        ...                            obstypes=['temp'],
        ...                            drive_filename="your_era5_temperature_data.csv",
        ...                             drive_folder="gee_timeseries_data")
        THE DATA AMOUNT IS TOO LARGE FOR INTERACTIVE SESSION, ...

        If the data request is small, the timeseries are present. If the
        data request is larger, a CSV file is writen to your google drive.
        Download that file, and use the `GeeDynamicModelData.set_modeldata_from_csv()`
        method.

        Let's assume that you downloaded the CSV file ("your_era5_temperature_data.csv")
        and saved it locally on your computer.

        >>> era5.set_modeldata_from_csv(csvpath="<Path to the donwloaded csv file>") # doctest: +SKIP

        """

        # 1. Read csv and set timezone
        df = pd.read_csv(csvpath, sep=",")
        # format to wide structure
        self._format_gee_df_structure(df)
        # convert units
        self._convert_units()


def import_modeldata_from_pkl(folder_path, filename="saved_modeldata.pkl"):
    """Import a GeeDynamicModelData instance from a (pickle) file.

    Parameters
    ----------
    folder_path : str
        The path to the folder where the GeeDynamicModelData pickle file is
        located.
    filename : str, optional
        The name of the output file. The default is 'saved_modeldata.pkl'.

    Returns
    -------
    metobs_toolkit.GeeDynamicModelData
        The modeldata instance.

    See Also
    -----------
    GeeDynamicModelData: Gee Modeldata dataset for time-evolving data.
    GeeDynamicModelData.save_modeldata: Save modeldata as a pkl file.
    GeeDynamicModelData.set_modeldata_from_csv: Set modeldata from a csv datafile.

    Examples
    --------
    As an example, we will use the ERA5-Land dataset, which is a default `GeeDynamicModelData`
    present in a `metobs_toolkit.Dataset()`. We will use the demo dataset,
    extract era5 timeseries data for the stations in the demo dataset, and
    save the GeeDynamicModelData. Then we wil import the data from a pkl file.

    >>> import metobs_toolkit
    >>> #Create your Dataset
    >>> dataset = metobs_toolkit.Dataset() #empty Dataset
    >>> dataset.import_data_from_file(
    ...                input_data_file=metobs_toolkit.demo_datafile,
    ...                input_metadata_file=metobs_toolkit.demo_metadatafile,
    ...                template_file=metobs_toolkit.demo_template)
    >>> era5 = dataset.gee_datasets['ERA5-land']
    >>> era5.set_metadf(dataset.metadf)

    Now we will extract temperature timeseries from ERA5. We already have
    a ModelObstype representing the temperature present in the GeeDynamicModelData:

    >>> era5.modelobstypes
    {'temp': ModelObstype instance of temp (linked to band: temperature_2m), 'pressure': ModelObstype instance of pressure (linked to band: surface_pressure), 'wind': ModelObstype_Vectorfield instance of wind (linked to bands: u_component_of_wind_10m and v_component_of_wind_10m)}

    We define a (short) period and extract ERA5 timeseries.

    >>> import datetime
    >>> tstart = datetime.datetime(2022,9,4)
    >>> tend = datetime.datetime(2022,9,5)
    >>> era5.extract_timeseries_data(
    ...                            startdt_utc=tstart,
    ...                            enddt_utc=tend,
    ...                            obstypes=['temp'])

    If the datarequest is small, the timeseries are present. (If the
    datarequest is larger, a csv file is writen to your google drive.
    Download that file, and use the `GeeDynamicModelData.set_modeldata_from_csv()`
    method.)


    We will save the era5 GeeDynamicModelData now as a (pkl) file. As an
    example we will store it in the current working directory (`os.getcwd()`)

    >>> import os
    >>> era5.save_modeldata(
    ...                outputfolder=os.getcwd(),
    ...                filename="era5_modeldata.pkl",
    ...                overwrite=True)
    Modeldata saved in ...

    If we later want to open the saved era5 Modeldata, we use the
    `GeeDynamicModelData.import_modeldata_from_pkl()` method.

    >>> your_era5 = metobs_toolkit.import_modeldata_from_pkl(
    ...                     folder_path=os.getcwd(),
    ...                     filename="era5_modeldata.pkl")

    """
    # check if folder_path is known and exists
    if not os.path.isdir(folder_path):
        raise MetobsModelDataError(f"{folder_path} is not a directory!")

    full_path = os.path.join(folder_path, filename)

    # check if file exists
    if not os.path.isfile(full_path):
        raise MetobsModelDataError(f"{full_path} does not exist.")

    with open(full_path, "rb") as inp:
        modeldata = pickle.load(inp)

    return modeldata


def _add_stations_to_folium_map(Map, metadf, display_cols=["name"]):
    """Add stations as markers to the folium map."""

    import folium

    metadf = metadf.reset_index()
    metadf["geometry"] = metadf["geometry"].to_crs("epsg:4326")
    for _, row in metadf.iterrows():
        point = row["geometry"]
        popuptext = ""
        for disp in display_cols:
            popuptext = f"{popuptext}{row[disp]}\n"
        folium.Marker(
            location=[point.y, point.x], fill_color="#43d9de", popup=popuptext, radius=8
        ).add_to(Map)

    return Map


# =============================================================================
# Errors
# =============================================================================


class MetobsModelDataError(Exception):
    """Exception raised for errors with ModelData."""

    pass


# =============================================================================
# Define default datasets
# =============================================================================

global_lcz_map = GeeStaticModelData(
    name="lcz",
    location="RUB/RUBCLIM/LCZ/global_lcz_map/latest",
    band_of_use="LCZ_Filter",
    value_type="categorical",
    scale=100,
    is_image=False,
    is_mosaic=True,
    credentials="Demuzere M.; Kittner J.; Martilli A.; Mills, G.; Moede, C.; Stewart, I.D.; van Vliet, J.; Bechtel, B. A global map of local climate zones to support earth system modelling and urban-scale environmental science. Earth System Science Data 2022, 14 Volume 8: 3835-3873. doi:10.5194/essd-14-3835-2022",
    class_map={
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
    },
    col_scheme={
        1: "#8c0000",
        2: "#d10000",
        3: "#ff0000",
        4: "#bf4d00",
        5: "#ff6600",
        6: "#ff9955",
        7: "#faee05",
        8: "#bcbcbc",
        9: "#ffccaa",
        10: "#555555",
        11: "#006a00",
        12: "#00aa00",
        13: "#648525",
        14: "#b9db79",
        15: "#000000",
        16: "#fbf7ae",
        17: "#6a6aff",
    },
)


global_dem = GeeStaticModelData(
    name="altitude",
    location="CGIAR/SRTM90_V4",
    band_of_use="elevation",
    value_type="numeric",
    scale=100,
    is_image=True,
    is_mosaic=False,
    credentials="SRTM Digital Elevation Data Version 4",
)

global_worldcover = GeeStaticModelData(
    name="worldcover",
    location="ESA/WorldCover/v200",
    band_of_use="Map",
    value_type="categorical",
    scale=10,
    is_image=False,
    is_mosaic=True,  # test this
    credentials="https://spdx.org/licenses/CC-BY-4.0.html",
    class_map={
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
    },
    agg_scheme={
        "water": [70, 80, 90, 95],
        "pervious": [10, 20, 30, 40, 60, 100],
        "impervious": [50],
    },
    col_scheme={
        10: "006400",
        20: "ffbb22",
        30: "ffff4c",
        40: "f096ff",
        50: "fa0000",
        60: "b4b4b4",
        70: "f0f0f0",
        80: "0064c8",
        90: "0096a0",
        95: "00cf75",
        100: "fae6a0",
    },
)

era5_land = GeeDynamicModelData(
    name="ERA5-land",
    location="ECMWF/ERA5_LAND/HOURLY",
    value_type="numeric",
    scale=2500,
    is_image=False,
    is_mosaic=False,  # test this
    credentials="",
    time_res="1h",
    modelobstypes=copy.deepcopy(default_era5_obstypes),
)


default_datasets = {  # Static datasets
    global_lcz_map.name: global_lcz_map,
    global_dem.name: global_dem,
    global_worldcover.name: global_worldcover,
    # Dynamic datasets
    era5_land.name: era5_land,
}

# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
