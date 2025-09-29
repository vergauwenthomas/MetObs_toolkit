#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Modeldata class and all its methods.

A Modeldata holds all timeseries coming from a model and methods to use them.
"""

# Standard library imports
from pathlib import Path
import copy
import sys
from typing import Union

import logging
from time import sleep

# Third-party imports
import pandas as pd
import numpy as np
import ee

# Local imports
import metobs_toolkit.gee_api as gee_api
from metobs_toolkit.backend_collection.errorclasses import MetObsModelDataError
import metobs_toolkit.backend_collection.printing_collection as printing
from metobs_toolkit.obstypes import default_era5_obstypes
from metobs_toolkit.plot_collection import (
    folium_map,
    add_title_to_folium_map,
    add_stations_to_folium_map,
)

# Fallback legend imports
try:
    from branca.element import MacroElement, Template
except ImportError:
    MacroElement = None
    Template = None
from metobs_toolkit.gee_api import connect_to_gee

from metobs_toolkit.obstypes import (
    ModelObstype,
    ModelObstype_Vectorfield,
)

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================


class _GEEDatasetManager:
    """Parent class for working with a GEE modeldataset.

    This class is abstract and holds methods and attributes that are applicable to all GEE Datasets.
    """

    def __init__(
        self,
        name: str,
        location: str,
        value_type: str,
        scale: int,
        is_static: bool,
        is_image: bool,
        credentials: str,
        is_mosaic: bool = False,
    ):
        """
        Create a GeeModelData abstract instance.

        Parameters
        ----------
        name : str
            The user-defined name for referring to this GEE dataset.
        location : str
            The location of the dataset on GEE.
        value_type : str
            Specify how to interpret the values of the GEE dataset ("numeric" or "categorical").
        scale : int
            The scale of the dataset to extract values of.
        is_static : bool
            If True, the GEE dataset is static and has no time-evolution component.
        is_image : bool
            If True, the GEE dataset is opened as ee.Image(), else ee.ImageCollection().
        credentials : str
            Credentials of the GEE dataset.
        is_mosaic : bool, optional
            If True, ee.mosaic() is applied on the GEE dataset. Default is False.

        Returns
        -------
        None
        """
        self.name = str(name)
        self.location = str(location)

        if str(value_type) not in ["categorical", "numeric"]:
            raise MetObsModelDataError(
                f'value_type: {value_type} is not "categorical" or "numeric"'
            )
        self.value_type = str(value_type)

        self.scale = int(scale)
        self.is_static = bool(is_static)
        self.is_image = bool(is_image)
        self._is_mosaic = bool(is_mosaic)
        self.credentials = str(credentials)

    # ------------------------------------------
    #    Specials
    # ------------------------------------------

    def __str__(self):
        return f"{self.__name__} representation of {self.name} "

    def __repr__(self):
        """Return string representation of the object."""
        return f"{type(self).__name__}(name={self.name}, location={self.location})"

    # =============================================================================
    # Checks
    # =============================================================================
    def _check_metadf_validity(self, metadf) -> pd.DataFrame:
        """
        Check if a metadf is valid (coordinates and structure-wise). If
        it is not valid, an error is raised.

        It filters out the rows in the metadf for which there are
        no 'lat' or 'lon' values.

        Parameters
        ----------
        metadf : pandas.DataFrame
            The metadata as a (geo)pandas dataframe.

        Returns
        -------
        pandas.Dataframe:
            The formatted metadf.
        """
        if metadf.empty:
            raise MetObsModelDataError(f"There is no metadata provided for {self}.")
        if metadf.index.name != "name":
            raise MetObsModelDataError(
                f"Wrong index name for setting {metadf} to {self}"
            )
        if "lat" not in metadf.columns:
            raise MetObsModelDataError(f'No "lat" column in the metadf of {self}.')
        if metadf["lat"].isnull().all():
            raise MetObsModelDataError(
                'All values of the "lat" column in the metadf are Nan.'
            )
        if "lon" not in metadf.columns:
            raise MetObsModelDataError(f'No "lon" column in the metadf of {self}.')
        if metadf["lon"].isnull().all():
            raise MetObsModelDataError(
                'All values of the "lon" column in the metadf are Nan.'
            )

        # at this point, it can happen that some stations do not have
        # a coordinate and others does. Drop these rows in the metadf,
        # since a call to Nan coordinates results in an error.

        if metadf[["lat", "lon"]].isnull().any(axis=1).any():
            missing_coords = metadf[
                metadf[["lat", "lon"]].isnull().any(axis=1)
            ].index.tolist()
            logger.warning(
                f"The following stations have missing coordinates, no data will be extracted: {missing_coords}"
            )
        metadf = metadf.dropna(subset=["lat", "lon"])
        return metadf

    def _get_all_gee_bandnames(self) -> list:
        """
        Return a list of all the bandnames of the GEE dataset.

        Returns
        -------
        list
            List of band names.
        """
        if self.is_image:
            return list(ee.Image(self.location).bandNames().getInfo())
        else:
            return list(ee.ImageCollection(self.location).first().bandNames().getInfo())

    def _get_base_details(self) -> str:
        """
        Print out basic details of the GEE dataset.

        Returns
        -------
        str
            String with dataset details.
        """
        retstr = ""
        retstr += printing.print_fmt_section("GEE Dataset details")
        retstr += printing.print_fmt_line(f"name: {self.name}")
        retstr += printing.print_fmt_line(f"location: {self.location}")
        retstr += printing.print_fmt_line(f"value_type: {self.value_type}")
        retstr += printing.print_fmt_line(f"scale: {self.scale}")
        retstr += printing.print_fmt_line(f"is_static: {self.is_static}")
        retstr += printing.print_fmt_line(f"is_image: {self.is_image}")
        retstr += printing.print_fmt_line(f"is_mosaic: {self._is_mosaic}")
        retstr += printing.print_fmt_line(f"credentials: {self.credentials}")
        return retstr


class GEEStaticDatasetManager(_GEEDatasetManager):
    """Class for working with static GEE modeldatasets."""

    def __init__(
        self,
        name: str,
        location: str,
        band_of_use: str,
        value_type: str,
        scale: int,
        is_image: bool,
        is_mosaic: bool = False,
        credentials: str = "",
        class_map: dict = {},
        agg_scheme: dict = {},
        col_scheme: dict = {},
    ):
        """
        Create a GeeStaticDataset instance representing a GEE dataset without a time dimension.

        Parameters
        ----------
        name : str
            The user-defined name for referring to this GEE dataset.
        location : str
            The location of the dataset on GEE.
        band_of_use : str
            The name of the band to use.
        value_type : str
            Specify how to interpret the values of the GEE dataset.
        scale : int
            The scale of the dataset to extract values of.
        is_image : bool
            If True, the GEE dataset is opened as ee.Image(), else ee.ImageCollection().
        is_mosaic : bool, optional
            If True, ee.mosaic() is applied on the GEE dataset. Default is False.
        credentials : str, optional
            Credentials of the GEE dataset. Default is "".
        class_map : dict, optional
            Mapping of numeric values to human-labels for categorical datasets. Default is {}.
        agg_scheme : dict, optional
            Aggregation scheme for custom classes. Default is {}.
        col_scheme : dict, optional
            Color scheme for classes. Default is {}.

        Returns
        -------
        None
        """
        super().__init__(
            name=name,
            location=location,
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

        self.__name__ = "GeeStaticDatasetManager"

    @log_entry
    def get_info(self, printout: bool = True) -> None:
        """
        Print out detailed information of the GeeStaticDataset.

        Parameters
        ----------
        printout : bool, optional
            If True, prints the information. If False, returns the string. Default is True.

        Returns
        -------
        None or str
        """

        retstr = ""
        retstr += printing.print_fmt_title("General info of GEEStaticDataset")

        retstr += self._get_base_details()
        retstr += printing.print_fmt_line(f"target band: {self.band_of_use}")

        retstr += printing.print_fmt_line("classification: ")
        retstr += printing.print_fmt_dict(self.class_map, identlvl=2)

        retstr += printing.print_fmt_line("aggregation: ")
        retstr += printing.print_fmt_dict(self.agg_scheme, identlvl=2)

        retstr += printing.print_fmt_line("colors: ")
        retstr += printing.print_fmt_dict(self.col_scheme, identlvl=2)

        if printout:
            print(retstr)
        else:
            return retstr

    @log_entry
    def extract_static_point_data(self, metadf: pd.DataFrame) -> pd.DataFrame:
        """
        Extract point values at the locations in the metadata.

        Parameters
        ----------
        metadf : pandas.DataFrame
            Metadata dataframe with station locations.

        Returns
        -------
        pandas.DataFrame
            Dataframe with station names as index and one column with values.
        """
        logger.debug(
            f"Entering GEEStaticDatasetManager.extract_static_point_data for {self}"
        )
        metadf = self._check_metadf_validity(metadf)

        ee_fc = gee_api._df_to_features_point_collection(metadf)

        raster = gee_api.get_ee_obj(self)

        results = raster.sampleRegions(collection=ee_fc, scale=self.scale).getInfo()

        if not bool(results["features"]):
            logger.warning(
                f"Something went wrong, gee did not return any data: {results}"
            )
            logger.info(
                f"(Could it be that (one) these coordinates are not on the map: {metadf}?)"
            )
            return pd.DataFrame()

        properties = [x["properties"] for x in results["features"]]
        df = pd.DataFrame(properties)

        if bool(self.class_map):
            df[self.band_of_use] = df[self.band_of_use].map(self.class_map)

        df = df.rename(columns={self.band_of_use: self.name})
        df = df.set_index(["name"])

        # Make sure all the index elements from metadf are in the df
        # Note: if gee does not return a value, that station is not present in the
        # df. So add them as Nan's
        df = df.reindex(metadf.index)

        return df

    @log_entry
    def extract_static_buffer_frac_data(
        self, metadf: pd.DataFrame, bufferradius: int, agg_bool: bool = False
    ) -> pd.DataFrame:
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
        metadf : pandas.DataFrame
            Metadata dataframe with station locations.
        bufferradius : int
            The radius (in meters) of the buffer.
        agg_bool : bool, optional
            If True, frequencies are computed on aggregated classes. Default is False.

        Returns
        -------
        pandas.DataFrame
            Dataframe with station names as index and class frequencies as columns.
        """
        logger.debug(
            f"Entering GEEStaticDatasetManager.extract_static_buffer_frac_data for {self}"
        )
        if metadf.empty:
            raise MetObsModelDataError(
                "No metadata is present for the GeeStaticDataset. No extraction possible."
            )

        metadf = self._check_metadf_validity(metadf)

        ee_fc = gee_api._df_to_features_buffer_collection(metadf, int(bufferradius))

        @log_entry
        def rasterExtraction(image):
            feature = image.reduceRegions(
                reducer=ee.Reducer.frequencyHistogram(),
                collection=ee_fc,
                scale=self.scale,
            )
            return feature

        raster = gee_api.get_ee_obj(self, force_mosaic=False)
        results = raster.map(rasterExtraction).flatten().getInfo()

        freqs = {
            staprop["properties"]["name"]: staprop["properties"]["histogram"]
            for staprop in results["features"]
        }
        freqsdf = pd.DataFrame(freqs)

        freqsdf = freqsdf.transpose().fillna(0)
        freqsdf.index.name = "name"

        freqsdf = freqsdf.div(freqsdf.sum(axis=1), axis=0)

        if agg_bool:
            agg_df = pd.DataFrame()
            for agg_name, agg_classes in self.agg_scheme.items():
                present_agg_classes = [
                    str(num) for num in agg_classes if str(num) in freqsdf.columns
                ]
                agg_df[agg_name] = freqsdf[present_agg_classes].sum(axis=1)

            freqsdf = agg_df

        elif bool(self.class_map):
            freqsdf = freqsdf.rename(
                columns={str(key): str(val) for key, val in self.class_map.items()}
            )
        else:
            pass

        freqsdf["buffer_radius"] = bufferradius
        freqsdf = freqsdf.reset_index().set_index(["name", "buffer_radius"])

        return freqsdf

    @log_entry
    def make_gee_plot(
        self,
        metadf: pd.DataFrame,
        save: bool = False,
        outputfolder: str = None,
        filename: str = None,
        vmin: Union[float, int, None] = None,
        vmax: Union[float, int, None] = None,
        overwrite: bool = False,
    ):
        """Make an interactive spatial plot of the GEE dataset and the stations.

        This method will create an interactive plot of the GEE dataset. If
        metadata is present, it will be displayed as markers on the map.

        The interactive map can be saved as an HTML file, by specifying the
        target path.


        Parameters
        ----------
        metadf : pandas.DataFrame
            Metadata dataframe with station locations.
        save : bool, optional
            If True, saves the map as an HTML file. Default is False.
        outputfolder : str or None, optional
            Path to the folder to save the HTML file. Default is None.
        filename : str or None, optional
            The filename for the HTML file. Default is None.
        vmin : numeric or None, optional
            Minimum value for colormap. Default is None.
        vmax : numeric or None, optional
            Maximum value for colormap. Default is None.
        overwrite : bool, optional
            If True, overwrites existing file. Default is False.

        Returns
        -------
        geemap.foliumap.Map
            The interactive map of the GeeStaticDataset.
        """
        if save:
            if outputfolder is None:
                raise MetObsModelDataError(
                    "If save is True, then outputfolder must be specified."
                )
            if filename is None:
                raise MetObsModelDataError(
                    "If save is True, then filename must be specified."
                )
            if filename[-5:] != ".html":
                filename += ".html"

            target_path = Path(outputfolder).joinpath(filename)
            if target_path.exists():
                if overwrite:
                    logger.info(f"Overwrite the file at {target_path}.")
                    target_path.unlink()
                else:
                    raise MetObsModelDataError(
                        f"{target_path} is already a file and overwrite is set to False!"
                    )

        connect_to_gee()

        im = gee_api.get_ee_obj(self)

        MAP = folium_map()

        if metadf.empty:
            pass

        else:
            if self.name in metadf.columns:
                logger.debug(
                    f"{self.name} is already found in the metadf, no point extraction is needed."
                )
            else:
                logger.debug(
                    f"{self.name} is extracted as point values (to present in markers)."
                )
                df = self.extract_static_point_data(metadf=metadf)
                metadf = metadf.merge(df, how="left", left_index=True, right_index=True)

            add_stations_to_folium_map(
                Map=MAP, metadf=metadf, display_cols=["name", self.name]
            )

            centroid = metadf.to_crs("EPSG:3857").dissolve().centroid.to_crs(metadf.crs)
            MAP.setCenter(lon=centroid.x.item(), lat=centroid.y.item(), zoom=8)

        if bool(self.class_map):
            vmin = min(self.class_map.keys())
            vmax = max(self.class_map.keys())
            if bool(self.col_scheme):
                var_visualization = {
                    "bands": [self.band_of_use],
                    "min": vmin,
                    "max": vmax,
                    "palette": list(self.col_scheme.values()),
                }

                try:
                    MAP.add_legend(
                        title="NLCD Land Cover Classification",
                        legend_dict={
                            self.class_map[numval]: self.col_scheme[numval]
                            for numval in self.class_map.keys()
                        },
                    )
                except FileNotFoundError:
                    # Fallback if geemap template is missing (Python 3.12+ / packaging issue)
                    logger.warning(
                        "geemap legend template not found, using fallback legend"
                    )
                    _add_fallback_legend(
                        MAP,
                        title="NLCD Land Cover Classification",
                        class_map=self.class_map,
                        col_scheme=self.col_scheme,
                    )

            else:
                var_visualization = {
                    "bands": [self.band_of_use],
                    "min": vmin,
                    "max": vmax,
                }

        else:
            if metadf.empty:
                if vmin is None:
                    vmin = 0.0
                if vmax is None:
                    vmax = 1.0
            else:
                obsmin = np.nanmin(metadf[self.name])
                obsmax = np.nanmax(metadf[self.name])
                if np.isnan(obsmin):
                    obsmin = 0.0
                if np.isnan(obsmax):
                    obsmax = 1.0

                if vmin is None:
                    vmin = obsmin - ((obsmax - obsmin) * 0.15)
                if vmax is None:
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

        MAP.add_layer(
            ee_object=im,
            vis_params=var_visualization,
            name=self.name,
        )

        MAP.addLayerControl()

        if save:
            logger.info(f"Saving {self.name} gee plot at: {target_path}")
            MAP.save(target_path)

        return MAP


class GEEDynamicDatasetManager(_GEEDatasetManager):
    """Class for working with Dynamic GEE modeldatasets."""

    def __init__(
        self,
        name: str,
        location: str,
        value_type: str,
        scale: int,
        time_res: str,
        modelobstypes: list,
        is_image: bool = False,
        is_mosaic: bool = False,
        credentials: str = "",
        # class_map: dict = {},
        # agg_scheme: dict = {},
        # col_scheme: dict = {},
    ):
        """
        Create a GeeDynamicDataset instance representing a GEE dataset with a time dimension.

        Parameters
        ----------
        name : str
            The user-defined name for referring to this GEE dataset.
        location : str
            The location of the dataset on GEE.
        value_type : str
            Specify how to interpret the values of the GEE dataset.
        scale : int
            The scale of the dataset to extract values of.
        time_res : str
            The time resolution of the dataset as a timedelta string.
        modelobstypes : list
            List of ModelObstype and ModelObstype_Vectorfield for this dataset.
        is_image : bool, optional
            If True, the GEE dataset is opened as ee.Image(), else ee.ImageCollection(). Default is False.
        is_mosaic : bool, optional
            If True, ee.mosaic() is applied on the GEE dataset. Default is False.
        credentials : str, optional
            Credentials of the GEE dataset. Default is "".


        Returns
        -------
        None
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

        self.modelobstypes = {}
        for obs in modelobstypes:
            if not (
                (isinstance(obs, ModelObstype))
                | (isinstance(obs, ModelObstype_Vectorfield))
            ):
                raise MetObsModelDataError(
                    f"{obs} is not an instance of ModelObstype or ModelObstype_Vectorfield but of type {type(obs)}."
                )
            self.modelobstypes[obs.name] = obs

        self.time_res = str(time_res)

        self.__name__ = "GeeDynamicDatasetManager"

    @log_entry
    def add_modelobstype(self, modelobstype) -> None:
        """
        Add a new ModelObstype to the GeeDynamicDataset.

        Parameters
        ----------
        modelobstype : ModelObstype or ModelObstype_Vectorfield
            The new modelobstype to add.

        Returns
        -------
        None
        """
        if not (
            (isinstance(modelobstype, ModelObstype))
            | (isinstance(modelobstype, ModelObstype_Vectorfield))
        ):
            raise MetObsModelDataError(
                f"{modelobstype} is not a ModelObstype of ModelObstype_Vectorfield"
            )
        if modelobstype.name in self.modelobstypes.keys():
            if modelobstype == self.modelobstypes[modelobstype.name]:
                return
            else:
                raise MetObsModelDataError(
                    f"There is already a known ModelObstype with {modelobstype.name} as a name: {self.modelobstypes[modelobstype.name]}"
                )
        else:
            self.modelobstypes[modelobstype.name] = modelobstype

    def _format_gee_df_structure(self, geedf: pd.DataFrame) -> pd.DataFrame:
        """
        Format a dataframe (constructed directly from GEE) to a modeldf.

        Parameters
        ----------
        geedf : pandas.DataFrame
            Dataframe from GEE.

        Returns
        -------
        pandas.DataFrame
            Formatted dataframe.
        """
        logger.debug(
            f"Entering GEEDynamicDatasetManager._format_gee_df_structure for {self}"
        )
        geedf["datetime"] = pd.to_datetime(geedf["datetime"], format="%Y%m%d%H%M%S")
        geedf["datetime"] = geedf["datetime"].dt.tz_localize("UTC")

        geedf = geedf.set_index(["name", "datetime"])
        geedf = geedf.sort_index()

        new_obs = []
        for obs in self.modelobstypes.values():
            if isinstance(obs, ModelObstype_Vectorfield):
                if (obs.model_band_u in geedf.columns) & (
                    obs.model_band_v in geedf.columns
                ):
                    amp_series, amp_obstype = obs.compute_amplitude(df=geedf)
                    geedf[amp_obstype.name] = amp_series
                    new_obs.append(amp_obstype)

                    dir_series, dir_obstype = obs._compute_angle(df=geedf)
                    geedf[dir_obstype.name] = dir_series
                    new_obs.append(dir_obstype)

        for obs in new_obs:
            self.add_modelobstype(obs)

        scalar_obs = [
            obs for obs in self.modelobstypes.values() if isinstance(obs, ModelObstype)
        ]
        band_mapper = {obs.model_band: obs.name for obs in scalar_obs}
        geedf = geedf.rename(columns=band_mapper)

        modeldf = geedf
        return modeldf

    def _subset_to_obstypes(self, df: pd.DataFrame, trg_obstypes: list) -> pd.DataFrame:
        """
        Subset the modeldf to a list of target obstypes.

        Parameters
        ----------
        df : pandas.DataFrame
            Dataframe to subset.
        trg_obstypes : list
            List of target obstypes.

        Returns
        -------
        pandas.DataFrame
            Subsetted dataframe.
        """
        logger.debug(
            f"Entering GEEDynamicDatasetManager._subset_to_obstypes for {self}"
        )
        keep_columns = []

        for obs in trg_obstypes:
            if isinstance(obs, ModelObstype):
                keep_columns.append(obs.name)
            elif isinstance(obs, ModelObstype_Vectorfield):
                keep_columns.append(obs._amp_obs_name)
                keep_columns.append(obs._dir_obs_name)
            else:
                raise MetObsModelDataError(
                    f"{obs} is not a ModelObstype or ModelObstype_Vectorfield."
                )

        return df[keep_columns]

    def _get_bandnames(self, trg_obstypes: list) -> list:
        """
        Get a list of all known target band names.

        Parameters
        ----------
        trg_obstypes : list
            List of target obstypes.

        Returns
        -------
        list
            List of band names.
        """
        trg_bands = []
        for obs in trg_obstypes:
            if isinstance(obs, ModelObstype):
                trg_bands.append(obs.model_band)
            elif isinstance(obs, ModelObstype_Vectorfield):
                trg_bands.append(obs.model_band_u)
                trg_bands.append(obs.model_band_v)
            else:
                raise MetObsModelDataError(
                    f"{obs} is not an instance of ModelObstype or ModelObstype_Vectorfield."
                )
        return trg_bands

    def _get_time_res(self) -> str:
        """Return the time resolution as a string."""
        return str(self.time_res)

    @log_entry
    def get_info(self, printout: bool = True) -> Union[None, str]:
        """
        Print out detailed information about the GeeDynamicDataset.

        Parameters
        ----------
        printout : bool, optional
            If True, prints the information. If False, returns the string. Default is True.

        Returns
        -------
        None or str
        """
        retstr = ""
        retstr += printing.print_fmt_title("General info of GEEDynamicDataset")

        retstr += self._get_base_details()
        retstr += printing.print_fmt_line(f"time res: {self.time_res}")

        retstr += printing.print_fmt_section("Known Modelobstypes")
        for obs in self.modelobstypes.values():
            retstr += printing.print_fmt_line(f"{obs.name} : {obs}")
            if isinstance(obs, ModelObstype_Vectorfield):
                retstr += printing.print_fmt_line(
                    "vectorfield that will be converted to: ", 2
                )
                retstr += printing.print_fmt_line(f"{obs._amp_obs_name}", 3)
                retstr += printing.print_fmt_line(f"{obs._dir_obs_name}", 3)
            retstr += printing.print_fmt_line(
                f"conversion: {obs.model_unit} --> {obs.std_unit}", 2
            )

        if printout:
            print(retstr)
        else:
            return retstr

    @log_entry
    def make_gee_plot(
        self,
        metadf: pd.DataFrame,
        timeinstance: pd.Timestamp,
        modelobstype: str = "temp",
        save: bool = False,
        outputfolder: str = None,
        filename: str = None,
        vmin: Union[float, int, None] = None,
        vmax: Union[float, int, None] = None,
        overwrite: bool = False,
    ):
        """
        Make an interactive spatial plot of the GEE dataset and the stations.

        Parameters
        ----------
        metadf : pandas.DataFrame
            Metadata dataframe with station locations.
        timeinstance : datetime.datetime or pandas.Timestamp
            The time instance to plot the GEE dataset.
        modelobstype : str, optional
            The name of the ModelObstype to plot. Default is "temp".
        save : bool, optional
            If True, saves the map as an HTML file. Default is False.
        outputfolder : str or None, optional
            Path to the folder to save the HTML file. Default is None.
        filename : str or None, optional
            The filename for the HTML file. Default is None.
        vmin : numeric or None, optional
            Minimum value for colormap. Default is None.
        vmax : numeric or None, optional
            Maximum value for colormap. Default is None.
        overwrite : bool, optional
            If True, overwrites existing file. Default is False.

        Returns
        -------
        geemap.foliumap.Map
            The interactive map of the GeeDynamicDataset.
        """
        if save:
            if outputfolder is None:
                raise MetObsModelDataError(
                    "If save is True, then outputfolder must be specified."
                )
            if filename is None:
                raise MetObsModelDataError(
                    "If save is True, then filename must be specified."
                )
            if filename[-5:] != ".html":
                filename += ".html"

            target_path = Path(outputfolder).joinpath(filename)
            if target_path.exists():
                if overwrite:
                    logger.info(f"Overwrite the file at {target_path}.")
                    target_path.unlink()
                else:
                    raise MetObsModelDataError(
                        f"{target_path} is already a file and overwrite is set to False!"
                    )

        if timeinstance.tz is None:
            timeinstance = timeinstance.tz_localize(tz="UTC")
        else:
            timeinstance = timeinstance.tz_convert(tz="UTC")

        timeinstance = timeinstance.floor(self.time_res)

        if modelobstype not in self.modelobstypes:
            raise MetObsModelDataError(
                f"{modelobstype} is not a known modelobstype ({self.modelobstypes})"
            )
        modelobstype = self.modelobstypes[modelobstype]

        if not isinstance(modelobstype, ModelObstype):
            raise MetObsModelDataError(f"{modelobstype} is not a ModelObstype.")

        connect_to_gee()

        # Get image
        im = gee_api.get_ee_obj(self, target_bands=[modelobstype.model_band])
        if gee_api._is_eeobj_empty(im):
            raise gee_api.MetObsGEEDatasetError(
                f"An empty GEE dataset is returned. Please check the validity of {self} and the target bands {[modelobstype.model_band]}"
            )

        # filter to timestamp
        im = im.filterDate(
            start=timeinstance.isoformat(),
            end=(timeinstance + pd.Timedelta(self.time_res)).isoformat(),
        )
        if gee_api._is_eeobj_empty(im):
            raise gee_api.MetObsGEEDatasetError(
                f"An empty GEE dataset detected after filtering to date. Please check if {timeinstance} is in the geedataset {self}."
            )
        im = im.first()

        # make empty map (FOLIUM BACKEND !! )
        MAP = folium_map()

        # show stations
        if metadf.empty:
            vmin = 0.0
            vmax = 1.0

        else:
            MAP = add_stations_to_folium_map(
                Map=MAP, metadf=metadf, display_cols=["name"]
            )
            # fix center
            centroid = metadf.to_crs("EPSG:3857").dissolve().centroid.to_crs(metadf.crs)
            MAP.setCenter(lon=centroid.x.item(), lat=centroid.y.item(), zoom=8)

            metadf = metadf.to_crs("epsg:4326")
            (xmin, ymin, xmax, ymax) = metadf.total_bounds

            if (vmin is None) | (vmax is None):
                roi = ee.Geometry.BBox(west=xmin, south=ymin, east=xmax, north=ymax)
                roi_min = im.reduceRegion(
                    ee.Reducer.min(), roi, scale=self.scale
                ).getInfo()[modelobstype.model_band]
                roi_max = im.reduceRegion(
                    ee.Reducer.max(), roi, scale=self.scale
                ).getInfo()[modelobstype.model_band]
            if vmin is None:
                vmin = roi_min - ((roi_max - roi_min) * 0.15)
                if vmin == vmax:
                    vmin = vmax - 1.0
            if vmax is None:
                vmax = roi_max + ((roi_max - roi_min) * 0.15)
                if vmax == vmin:
                    vmax = vmin + 1.0

        var_visualization = {
            "bands": [modelobstype.model_band],
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
            name=f"{self.name}: {modelobstype.model_band}",
        )

        title = f"{self.name} GEE plot of {modelobstype.model_band} (in {modelobstype.model_unit}) at {timeinstance}."
        MAP = add_title_to_folium_map(title=title, Map=MAP)

        # add layer control
        MAP.addLayerControl()

        # Save
        if save:
            logger.info(f"Saving {self.name} gee plot at: {target_path}")
            MAP.save(target_path)

        return MAP

    @log_entry
    def extract_timeseries_data(
        self,
        metadf: pd.DataFrame,
        startdt_utc,
        enddt_utc,
        obstypes: list = ["temp"],
        get_all_bands: bool = False,
        drive_filename: str = None,
        drive_folder: str = "gee_timeseries_data",
        force_direct_transfer: bool = False,
        force_to_drive: bool = False,
    ):
        """
        Extract timeseries data and set the modeldf.

        Parameters
        ----------
        metadf : pandas.DataFrame
            Metadata dataframe with station locations.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstypes : list of str, optional
            List of ModelObstype names to extract modeldata for. Default is ['temp'].
        get_all_bands : bool, optional
            If True, all values (over all bands) are extracted. Default is False.
        drive_filename : str or None, optional
            Filename for saving data on Google Drive. Default is None.
        drive_folder : str, optional
            Folder on Google Drive to save the file. Default is "gee_timeseries_data".
        force_direct_transfer : bool, optional
            If True, forces direct transfer. Default is False.
        force_to_drive : bool, optional
            If True, forces writing to Google Drive. Default is False.

        Returns
        -------
        pandas.DataFrame or None
            Extracted timeseries dataframe or None if written to Drive.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your Google Drive. In this case, use the
        ``metobs_toolkit.Dataset.import_gee_data_from_file()`` method to import the data.
        """
        logger.debug(
            f"Entering GEEDynamicDatasetManager.extract_timeseries_data for {self}"
        )
        metadf = self._check_metadf_validity(metadf)

        if (force_direct_transfer) & (force_to_drive):
            raise MetObsModelDataError(
                "Both force_direct_transfer and force_to_drive could not be True at the same time."
            )
        # Check obstypes
        if isinstance(obstypes, str):
            obstypes = [obstypes]

        for obstype in obstypes:
            if obstype not in self.modelobstypes.keys():
                raise MetObsModelDataError(
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

        connect_to_gee()

        if get_all_bands:
            bandnames = self._get_all_gee_bandnames()
        else:
            bandnames = self._get_bandnames(trg_obstypes=obstypes)

        logger.info(f"{bandnames} are extracted from {self}.")

        _est_data_size = gee_api._estimate_data_size(
            metadf=metadf,
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
            logger.warning(
                "THE DATA AMOUNT IS TOO LARGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
            )

            use_drive = True
        else:
            use_drive = False

        ee_fc = gee_api._df_to_features_point_collection(metadf)

        @log_entry
        def rasterExtraction(image):
            feature = image.sampleRegions(
                collection=ee_fc,
                scale=self.scale,
            )
            return feature

        # Because the daterange is maxdate exclusive, add the time resolution to the enddt
        enddt_utc = enddt_utc + pd.Timedelta(self.time_res)

        raster = gee_api.get_ee_obj(self, target_bands=bandnames)
        results = (
            raster.filter(
                ee.Filter.date(
                    gee_api.datetime_to_gee_datetime(startdt_utc),
                    gee_api.datetime_to_gee_datetime(enddt_utc),
                )
            )
            .map(gee_api._addDate)
            .map(rasterExtraction)
            .flatten()
        )

        if not use_drive:
            results = results.getInfo()

            properties = [x["properties"] for x in results["features"]]
            df = pd.DataFrame(properties)

            if df.empty:
                sys.exit("ERROR: the returned timeseries from GEE are empty.")

            df = self._format_gee_df_structure(
                df
            )  # format to wide structure + rename columns etc

            if not get_all_bands:
                df = self._subset_to_obstypes(df=df, trg_obstypes=obstypes)

            # convert units + update attr

            return df

        else:
            if drive_filename is None:
                _filename = f"{self.name}_timeseries_data"
            else:
                _filename = str(drive_filename)
            logger.warning(
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
            logger.warning(
                "The data is transfered! Open the following link in your browser: \n\n"
            )
            logger.warning(f"{doc_folder_id} \n\n")
            logger.warning(
                "To upload the data to the model, use the \
Dataset.import_gee_data_from_file() method."
            )

            return


# =============================================================================
# Helper functions
# =============================================================================


def _add_fallback_legend(Map, title: str, class_map: dict, col_scheme: dict):
    """
    Fallback legend implementation when geemap's add_legend template is missing.

    Parameters
    ----------
    Map : geemap.foliumap.Map
        The map object to add the legend to.
    title : str
        The title for the legend.
    class_map : dict
        Mapping of numeric values to class labels.
    col_scheme : dict
        Mapping of numeric values to color codes.

    Returns
    -------
    geemap.foliumap.Map
        The map with the added legend.
    """
    if MacroElement is None or Template is None:
        logger.warning("Cannot create fallback legend: branca not available")
        return Map

    rows = ""
    for key, label in class_map.items():
        color = col_scheme.get(key, "#FFFFFF")
        rows += f"""
        <tr>
          <td style="text-align:left;">
            <i style="background:{color};opacity:0.85;border:1px solid #555;display:inline-block;width:18px;height:18px;"></i>
          </td>
          <td style="padding-left:4px;">{label}</td>
        </tr>"""

    html = f"""
    <div id="metobs-toolkit-legend" style="position: fixed;
         bottom: 20px; left: 20px; z-index: 9999;
         background: white; padding: 10px 12px; border:2px solid #444;
         font-size: 12px; font-family: Arial, sans-serif; max-height:300px; overflow:auto;">
      <b>{title}</b>
      <table style="border:none; margin-top:4px;">{rows}</table>
    </div>
    """

    tpl = Template(html)
    macro = MacroElement()
    macro._template = tpl
    Map.get_root().add_child(macro)
    return Map


# =============================================================================
# Define default datasets
# =============================================================================

global_LCZ_map = GEEStaticDatasetManager(
    name="LCZ",
    location="RUB/RUBCLIM/LCZ/global_lcz_map/latest",
    band_of_use="LCZ_Filter",
    value_type="categorical",
    scale=100,
    is_image=False,
    is_mosaic=True,
    credentials="Demuzere M.; Kittner J.; Martilli A.; Mills, G.; Moede, C.; Stewart, I.D.; van Vliet, J.; Bechtel, B. A global map of local climate zones to support earth system modelling and urban-scale environmental science. Earth System Science Data 2022, 14 Volume 8: 3835-3873. doi:10.5194/essd-14-3835-2022",
    class_map={
        1: "Compact highrise",
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


global_dem = GEEStaticDatasetManager(
    name="altitude",
    location="CGIAR/SRTM90_V4",
    band_of_use="elevation",
    value_type="numeric",
    scale=100,
    is_image=True,
    is_mosaic=False,
    credentials="SRTM Digital Elevation Data Version 4",
)

global_worldcover = GEEStaticDatasetManager(
    name="worldcover",
    location="ESA/WorldCover/v200",
    band_of_use="Map",
    value_type="categorical",
    scale=10,
    is_image=False,
    is_mosaic=True,
    credentials="https://spdx.org/licenses/CC-BY-4.0.html",
    class_map={
        10: "Tree cover",
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

era5_land = GEEDynamicDatasetManager(
    name="ERA5-land",
    location="ECMWF/ERA5_LAND/HOURLY",
    value_type="numeric",
    scale=2500,
    is_image=False,
    is_mosaic=False,
    credentials="",
    time_res="1h",
    modelobstypes=copy.deepcopy(default_era5_obstypes),
)


default_datasets = {
    global_LCZ_map.name: global_LCZ_map,
    global_dem.name: global_dem,
    global_worldcover.name: global_worldcover,
    era5_land.name: era5_land,
}
