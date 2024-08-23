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

# from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.obstypes import Obstype as Obstype_class
from metobs_toolkit.obstype_modeldata import (
    # model_obstypes,
    ModelObstype,
    ModelObstype_Vectorfield,
)

# from metobs_toolkit.obstype_modeldata import compute_amplitude, compute_angle
from metobs_toolkit.settings import Settings

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
            The user-defined name for refering to this GEE dataset.
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
            If True, ee.mosaic() is appied on the GEE dataset. The default is False.


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
        Check if a metadf is valid (coordinates and structure wise). If
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
            The user-defined name for refering to this GEE dataset.
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
            If True, ee.mosaic() is appied on the GEE dataset. The default is False.
        credentials : str, optional
            Credentials of the GEE dataset. The default is "".
        class_map : dict, optional
            If value_type is categorical, than the class_map defines how the
            numeric values are mapped to 'human-categories'. The keys are the
            numeric values, the values are the human-labels. The default is {}.
        agg_scheme : dict, optional
            If value_types is categorical, than the agg scheme defines custom-
            made classes, which are aggregates of the present classes. The
            keys are the names of the custom-classes, the values are lists, with
            the corresponing numeric values. The default is {}.
        col_scheme : dict, optional
            if value_types is categorical, the col_sheme defines the colors used
            for each class. The keys are the numeric values, the values are
            the colors (in hex form). The default is {}.

        Returns
        -------
        None.

        Note
        -------
        On general, specifying a scale smaller than the true scale of the GEE
        dataset has no impact on the results (but can effect the computation time).

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
            return f"Empty {self.__name__} instance of {self.name} "
        else:
            return f"{self.__name__} instance of {self.name} "

    def __repr__(self):
        """Print overview information of the modeldata."""
        return self.__str__()

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

        First a connection with the gee server is made. Then the coordinates
        of the metadata and the details of the GEE dataset are send to the
        gee server. There the point values are extracted and are send back.


        Returns
        -------
        df : pandas.DataFrame
            A dataframe with the stationnames as index, and one column with
            values (mapped to human-labels if the datataset is categorical and
            a cat_map is defined).

        Note
        -------
        Make sure that the metadata is set. Use the
        `GeeStaticModelData.set_metadata()` for this.

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

        raster = gee_api.get_ee_obj(self)  # dataset
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

            add_stations_to_folium_map(
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


def add_stations_to_folium_map(Map, metadf, display_cols=["name"]):
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


class GeeDynamicModelData(GeeModelData):
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

        self.modelobstypes = modelobstypes
        self.time_res = str(time_res)

        self.__name__ = "GeeDynamicModelData"

    def __str__(self):
        if self.modeldf.empty:
            return f"Empty {self.__name__} instance of {self.name} "
        else:
            return f"{self.__name__} instance of {self.name} with modeldata "

    # =============================================================================
    # Setters
    # =============================================================================

    def _set_modeldf(self, modeldf):

        if not list(modeldf.index.names) == ["name", "datetime"]:
            raise MetobsModelDataError()(
                f"A dataframe is being set as .modeldf with wrong modeldf.index.names: {modeldf.index.names}"
            )
        self.modeldf = modeldf

    def add_modelobstype(self, modelobstype):
        # TODO checks
        if not (
            (isinstance(modelobstype, ModelObstype))
            | (isinstance(modelobstype, ModelObstype_Vectorfield))
        ):
            raise MetobsModelDataError(
                f"{modelobstype} is not a ModelObstype of ModelObstype_Vectorfield"
            )
        self.modelobstypes[modelobstype.name] = modelobstype

    # =============================================================================
    # Convertors/formatters
    # =============================================================================
    def _convert_units(self):

        modeldf = self.modeldf
        for obs in self.modelobstypes.values():
            if obs.name in modeldf.columns:
                modeldf[obs.name] = obs.convert_to_standard_units(
                    input_data=modeldf[obs.name], input_unit=obs.get_modelunit()
                )
        self._set_modeldf(modeldf)
        return

    def _format_gee_df_structure(self, geedf):
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
                    dir_series, dir_obstype = obs.compute_angle(df=geedf)
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
        keep_columns = []

        for obs in trg_obstypes:
            if isinstance(obs, ModelObstype):
                keep_columns.append(obs.name)
            elif isinstance(obs, ModelObstype_Vectorfield):
                keep_columns.append(obs._amp_obs_name)
            else:
                raise MetobsModelDataError(
                    f"{obs} is not a ModelObstype or ModelObstype_Vectorfield."
                )

        self.modeldf = self.modeldf[keep_columns]

    # =============================================================================
    # Getters
    # =============================================================================

    def _get_bandnames(self, trg_obstypes):
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

    def _get_unknonw_modelcolumns(self):
        return list(set(self.modeldf.columns) - set(self.modelobstypes.keys()))

    def get_info(self):
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

        if bool(self._get_unknonw_modelcolumns()):
            print(
                "\n (The following data columns (bandnames) are present, without a corresponding Modelobstype: "
            )
            print(f"{self._get_unknonw_modelcolumns()}")

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
        im = gee_api.get_ee_obj(self, target_bands=modelobstype.get_modelband())
        # filter to timestamp
        im = im.filterDate(
            start=dt.isoformat(), end=(dt + pd.Timedelta(self.time_res)).isoformat()
        )
        im = im.first()

        # make empty map (FOLIUM BACKEND !! )
        MAP = folium_map()

        # show stations
        if self.metadf.empty:
            pass

        else:
            add_stations_to_folium_map(
                Map=MAP, metadf=self.metadf, display_cols=["name"]
            )
            # fix center
            centroid = self.metadf.dissolve().centroid
            MAP.setCenter(lon=centroid.x.item(), lat=centroid.y.item(), zoom=8)

            # Create GEE boundbox of the ROI

            metadf = self.metadf.to_crs("epsg:4326")
            (xmin, ymin, xmax, ymax) = metadf.total_bounds

            roi = ee.Geometry.BBox(west=xmin, south=ymin, east=xmax, north=ymax)

        roi_max = im.reduceRegion(ee.Reducer.max(), roi).getInfo()[
            modelobstype.get_modelband()
        ]
        roi_min = im.reduceRegion(ee.Reducer.min(), roi).getInfo()[
            modelobstype.get_modelband()
        ]

        if vmin == None:
            vmin = roi_min - ((roi_max - roi_min) * 0.15)
        if vmax == None:
            vmax = roi_max + ((roi_max - roi_min) * 0.15)

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

        The line colors represent the timesries for different locations.



        Parameters
        ----------
        obstype_model : string, optional
             Fieldname of the Modeldata to visualise. The default is 'temp'.
        Dataset : metobs_toolkit.Dataset, optional
            A Dataset instance with observations plotted in the same figure.
            Observations are represented by solid line and modeldata by dashed
            lines. The default is None.
        obstype_dataset : string, optional
            Fieldname of the Dataset to visualise. Only relevent when a dataset
            is provided. If None, obsype_dataset = obstype_model. The default
            is None.
        stationnames : list, optional
            A list with stationnames to include in the timeseries. If None is
            given, all the stations are used, defaults to None.
        starttime : datetime.datetime, optional
             Specifiy the start datetime for the plot. If None is given it will
             use the start datetime of the dataset, defaults to None.
        endtime : datetime.datetime, optional
             Specifiy the end datetime for the plot. If None is given it will
             use the end datetime of the dataset, defaults to None.
        title : string, optional
             Title of the figure, if None a default title is generated. The
             default is None.
        show_outliers : bool, optional
             If true the observations labeld as outliers will be included in
             the plot. Only relevent when a dataset is provided. The default
             is True.
        show_filled : bool, optional
             If true the filled values for gaps and missing observations will
             be included in the plot. Only relevent when a dataset is provided.
             The default is True.
        legend : bool, optional
             If True, a legend is added to the plot. The default is True.


        Returns
        -------
        axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

        """
        logger.info(f"Make {obstype_model}-timeseries plot of model data")

        # Basic test
        if obstype_model not in self.modeldf.columns:
            logger.warning(
                f"{obstype_model} is not foud in the modeldata df (columns = {self.modeldf.columns})."
            )
            return
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
        y_label = self.modelobstypes[obstype_model].get_plot_y_label()

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
    ):
        """Extract timeseries of a gee dataset.

        The extraction can only be done if the gee dataset bandname (and units)
        corresponding to the obstype is known.

        The units are converted to the toolkit standard units!!

        Parameters
        ----------
        mapname : str
            Mapname of choice of the GEE dataset to extract data from.
        metadf : pandas.DataFrame
            A dataframe with a 'name' index and  'lat', 'lon' columns.
            Timeseries are extracted for these locations.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstypes : str or list of strings, optional
            Toolkit observation type to extract data from. There should be a
            bandname mapped to this obstype for the gee map. Multiple obstypes
            can be given in a list. The default is 'temp'.


        Returns
        -------
        None.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        """
        # ====================================================================
        # Test input
        # ====================================================================
        self._check_metadf_validity(self.metadf)

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

        if _est_data_size > 4000:
            print(
                "THE DATA AMOUT IS TO LAREGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
            )
            logger.info(
                "THE DATA AMOUT IS TO LAREGE FOR INTERACTIVE SESSION, THE DATA WILL BE EXPORTED TO YOUR GOOGLE DRIVE!"
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

        Returns
        -------
        None.

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
        """Import timeseries data that is stored in a csv file.

        The name of the gee dataset the timeseries are coming from must be the
        same as the .name attribute of the Modeldata.


        The timeseries will be formatted and converted to standard toolkit
        units.

        Parameters
        ----------
        csvpath : str
            Path of the csv file containing the modeldata timeseries.

        Returns
        -------
        None.

        """

        # 1. Read csv and set timezone
        df = pd.read_csv(csvpath, sep=",")
        # format to wide structure
        self._format_gee_df_structure(df)
        # convert units
        self._convert_units()


def import_modeldata_from_pkl(folder_path, filename="saved_modeldata.pkl"):
    """Import a modeldata instance from a (pickle) file.

    Parameters
    ----------
    folder_path : str
        The path to the folder where the Modeldata pickle file is located.
    filename : str, optional
        The name of the output file. The default is 'saved_modeldata.pkl'.

    Returns
    -------
    metobs_toolkit.Modeldata
        The modeldata instance.

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


# =============================================================================
# GEE api
# =============================================================================


def connect_to_gee(**kwargs):
    """
    Setup authentication for the use of the GEE Python API.

    For a fresh kernel, without stored credentials, a prompt/browser window
    will appear with further instructions for the authentication.


    Parameters
    ----------
    **kwargs : Kwargs passed to ee.Authenticate()
        Kwargs are only used by the user, for resetting the gee connection. See
        the Note below.

    Returns
    -------
    None.

    Note
    ------
    Upon calling, this function assumes you have a Google developers account,
    and a project with the Google Earth Engine API enabled.
    See the * Using Google Earth Engine * page for more info.

    Note
    ------
    During the Authentication, you will be asked if you want a read-only scope.
    A read-only scope is sufficient when the data is transferred directly to your
    machine (small data transfers), but will not be sufficient when extracting
    large amounts of data (typical for extracting Modeldata). This is because
    modeldata is written directly to your Google Drive, and therefore
    the read-only scope is insufficient.

    Note
    ------
    Due to several reasons, an EEExeption may be thrown. This is
    likely because of an invalid credential file. To fix this, you
    can update your credential file, and specify a specific authentication method.
    We found that the "notebook" authentication method works best for most users.

    Here is an example on how to update the credentials:

    .. code-block:: python

        import metobs_toolkit

        metobs_toolkit.connect_to_gee(force=True, #create new credentials
                                      auth_mode='notebook', # 'notebook', 'localhost', 'gcloud' (requires gcloud installed) or 'colab' (works only in colab)
                                      )

    """

    if bool(kwargs):  # kwargs are always passed by user, so reinitialize
        ee.Authenticate(**kwargs)
        ee.Initialize()
        return

    if not ee.data._credentials:  # check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return


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


# TODO: first check the static gee dataset class
# then adjust the modeldatatype class so that each instance works with
# one dataset (simplification), inherit from obstype

# then fix the dynamic class


# %%


class Modeldata_old:
    """Class holding data and methods for a modeldata-timeseries."""

    def __init__(self, modelname):
        """Initialize modeldata."""
        self.df = init_multiindexdf()
        self.modelname = modelname

        self._settings = Settings()
        self.mapinfo = self._settings.gee["gee_dataset_info"]

        self.df_tz = "UTC"  # the timezone of the datetimes stored in the df

        self.obstypes = model_obstypes  # Dict name: Obstype-instance

    def __str__(self):
        """Print overview information of the modeldata."""
        if self.df.empty:
            return "Empty Modeldata instance."
        n_stations = self.df.index.get_level_values("name").unique().shape[0]
        obstypes = self.df.columns.to_list()
        startdt = self.df.index.get_level_values("datetime").min()
        enddt = self.df.index.get_level_values("datetime").max()
        data_units = [self.obstypes[col].get_standard_unit() for col in self.df.columns]

        return f"Modeldata instance containing: \n \
    * Modelname: {self.modelname} \n \
    * {n_stations} timeseries \n \
    * The following obstypes are available: {obstypes} \n \
    * Data has these units: {data_units} \n \
    * From {startdt} --> {enddt} (with tz={self.df_tz}) \n \n (Data is stored in the .df attribute)"

    def __repr__(self):
        """Print overview information of the modeldata."""
        return self.__str__()

    def get_info(self):
        """Print out detailed information on the Modeldata."""
        print(str(self))

        print("\n ------ Known gee datasets -----------")
        self.list_gee_datasets()

    def add_obstype(self, Obstype, bandname, band_units, band_description=None):
        """Add a new Observation type for the current Modeldata.


        Parameters
        ----------
        Obstype : metobs_toolkit.obstype.Obstype
            The new Obstype to add.
        bandname : str
            The name of the band that represents the obstype.
        band_units : str
            The unit the band is in. This unit must be a knonw-unit in the
            Obstype.
        band_description : str, optional
            A detailed description of the band. The default is None.

        Returns
        -------
        None.

        """
        if not isinstance(Obstype, Obstype_class):
            sys.exit(
                f"{Obstype} is not an instance of metobs_toolkit.obstypes.Obstype."
            )

        obs = Obstype

        # Test if the band unit is a knonw unit
        if not obs.test_if_unit_is_known(band_units):
            sys.exit(
                f"The {bandname} unit: {band_units} is not a knonw unit for {obs.name}"
            )

        # Make the modeldata extension
        equiv_dict = {
            self.modelname: {
                "name": str(bandname),
                "units": str(band_units),
                "band_desc": str(band_description),
            }
        }

        modeldata_obstype = ModelObstype(obstype=obs, model_equivalent_dict=equiv_dict)

        # add Obstype
        self.obstypes[obs.name] = modeldata_obstype
        logger.info(f"{obs.name} added to the known observation types.")

    def _df_as_triple_idx(self):
        """Return the data in the same form as Dataset records.

        The returned dataframe has ['name', 'obstype', 'datetime'] as index,
        and one column: 'value'.

        """
        df = self.df.stack()
        df.index = df.index.set_names("obstype", level=-1)
        df.name = "value"
        df = df.reset_index().set_index(["name", "obstype", "datetime"]).sort_index()
        return df

    def add_gee_dataset(
        self,
        mapname,
        gee_location,
        obstype,
        bandname,
        units,
        scale,
        band_desc=None,
        time_res="1h",
        is_image=False,
        is_numeric=True,
        credentials="",
    ):
        """Add a new gee dataset to the available gee datasets.

        Parameters
        ----------
        mapname : str
            Mapname of choice for the GEE dataset to add.
        gee_location : str
            Location of the gee dataset (like "ECMWF/ERA5_LAND/HOURLY" for ERA5).
        obstype : str
            The observation type name the band corresponds to.
        bandname : str
            Name of the dataset band as stored on the GEE.
        units : str
            The units of the band.
        scale : int
            The scale to represent the dataset in. (This is a GEE concept that
            is similar to the resolution in meters).
        band_desc : str or None, optional
            Add a descrition to of the band. The default is None.
        time_res : timedelta string, optional
            Time reoslution of the dataset, if is_image == False. The default is '1H'.
        is_image : bool, optional
            If True, the dataset is a ee.Image, else it is assumed to be an
            ee.ImageCollection. The default is False.
        is_numeric : bool, optional
            If True, the bandvalues are interpreted as numerical values rather
            than categorical.. The default is True.
        credentials : str, optional
            Extra credentials of the dataset. The default is ''.

        Returns
        -------
        None.

        Note
        -------
        To list all available gee dataset, use the .list_gee_dataset() method.

        Note
        -------
        Currently no unit conversion is perfomed automatically other than K -->
        Celsius. This will be implemented in the future.

        """
        # check if mapname exists
        if mapname in self.mapinfo.keys():
            logger.warning(
                f"{mapname} is found in the list of known gee datasets: {list(self.mapinfo.keys())}, choose a different mapname."
            )
            return

        if is_numeric:
            val_typ = "numeric"
        else:
            val_typ = "categorical"

        # Dataset defenition
        new_info = {
            mapname: {
                "location": f"{gee_location}",
                "usage": "user defined addition",
                "value_type": val_typ,
                "dynamical": not bool(is_image),
                "scale": int(scale),
                "is_image": bool(is_image),
                "is_imagecollection": not bool(is_image),
                "credentials": f"{credentials}",
            }
        }

        if not is_image:
            new_info[mapname]["time_res"] = f"{time_res}"

        # obstype defenition
        # 1. if obstype exists, update the obstype
        if obstype in self.obstypes:
            self.obstypes[obstype].add_new_band(
                mapname=mapname, bandname=bandname, bandunit=units, band_desc=band_desc
            )

        # 2. if obstype does not exist, create the obstype
        else:
            sys.exit(
                f"{obstype} is an unknown obstype. First add this obstype to the Modeldata, and than add a gee dataset."
            )

        self.mapinfo.update(new_info)
        logger.info(
            f"{mapname} is added to the list of available gee dataset with: {new_info}"
        )
        return

    def list_gee_datasets(self):
        """Print out all the available gee datasets.

        Returns
        -------
        None.

        """
        print("The following datasets are found: ")
        for geename, info in self.mapinfo.items():
            print("\n --------------------------------")
            print(f"{geename} : \n")
            # find which observations that are mappd
            mapped_obs = [
                obstype
                for obstype in self.obstypes.values()
                if obstype.has_mapped_band(geename)
            ]
            if len(mapped_obs) == 0:
                print(f" No mapped observation types for {geename}.")
            else:
                for obs in mapped_obs:
                    obs.get_info()
            print("\n INFO: \n")
            print(f"{info}")

    def _conv_to_timezone(self, tzstr):
        """Convert the timezone of the datetime index of the df attribute.

        Parameters
        ----------
        tzstr : str
            TImezonstring from the pytz module.

        Returns
        -------
        None.

        """
        # get tzstr by datetimindex.tz.zone

        df = self.df
        df["datetime_utc"] = df.index.get_level_values("datetime").tz_convert(tzstr)
        df = df.reset_index()
        df = df.drop(columns=["datetime"])
        df = df.rename(columns={"datetime_utc": "datetime"})
        df = df.set_index(["name", "datetime"])
        self.df = df
        self.df_tz = tzstr

    def _convert_units_to_tlk(self, obstype):
        """Convert the model data of one observation to the standard units.

        The data attributes will be updated.

        Parameters
        ----------
        obstype : str
            Observation type to convert to standard units.

        Returns
        -------
        None.

        """
        # chech if data is available
        if self.df.empty:
            logger.warning("No data to set units for.")
            return

        if obstype not in self.obstypes:
            logger.warning(
                f"{obstype} not found as a known observationtype in the Modeldata."
            )
            return

        if isinstance(self.obstypes[obstype], ModelObstype):
            # scalar obstype
            if obstype not in self.df.columns:
                logger.warning(
                    f"{obstype} not found as observationtype in the Modeldata."
                )
                return
        if isinstance(self.obstypes[obstype], ModelObstype_Vectorfield):
            # vector obstype
            if self.obstypes[obstype].get_u_column() not in self.df.columns:
                logger.warning(
                    f"{self.obstypes[obstype].get_u_column()} not found as observationtype in the Modeldata."
                )
                return
            if self.obstypes[obstype].get_v_column() not in self.df.columns:
                logger.warning(
                    f"{self.obstypes[obstype].get_v_column()} not found as observationtype in the Modeldata."
                )
                return

        cur_unit = self.obstypes[obstype].get_current_data_unit()

        if isinstance(self.obstypes[obstype], ModelObstype):
            converted_data = self.obstypes[obstype].convert_to_standard_units(
                input_data=self.df[obstype], input_unit=cur_unit
            )
            # Update the data and the current unit
            self.df[obstype] = converted_data
            self.obstypes[obstype].set_current_data_unit(
                current_data_unit=self.obstypes[obstype].get_standard_unit()
            )
        if isinstance(self.obstypes[obstype], ModelObstype_Vectorfield):
            u_comp_name = self.obstypes[obstype].get_u_column()
            v_comp_name = self.obstypes[obstype].get_v_column()
            u_comp, v_comp = self.obstypes[obstype].convert_to_standard_units(
                input_df=self.df, input_unit=cur_unit
            )

            self.df[u_comp_name] = u_comp
            self.df[v_comp_name] = v_comp
        logger.info(
            f"{obstype} are converted from {cur_unit} --> {self.obstypes[obstype].get_standard_unit()}."
        )

    def _exploid_2d_vector_field(self, obstype):
        """Compute amplitude and direction of 2D vector field components.

        The amplitude and directions are added to the data attribute, and their
        equivalent observationtypes are added to the known ModelObstypes.

        (The vector components are not saved.)
        Parameters
        ----------
        obstype : str
            The name of the observationtype that is a ModelObstype_Vectorfield.

        Returns
        -------
        None.

        """
        # check if the obstype is a vector field
        if not isinstance(self.obstypes[obstype], ModelObstype_Vectorfield):
            logger.warning(
                f"{obstype} is not a 2D vector field, so it can not be exploided."
            )
            return

        # get amplitude of 2D vectors
        logger.info(f"Computing the amplited of the 2D vector field of {obstype}")
        amp_data, amp_obstype = compute_amplitude(
            modelobs_vectorfield=copy.deepcopy(self.obstypes[obstype]), df=self.df
        )

        # # get direction of 2D vectors
        logger.info(f"Computing the direction of the 2D vector field of {obstype}")
        dir_data, dir_obstype = compute_angle(
            modelobs_vectorfield=copy.deepcopy(self.obstypes[obstype]), df=self.df
        )

        #  ------ update the attributes ---------

        # add new columns to the df
        self.df[amp_obstype.name] = amp_data
        self.df[dir_obstype.name] = dir_data

        # remove components from the df (Needed because they are not linked to an obstype)
        self.df = self.df.drop(
            columns=[
                self.obstypes[obstype].get_u_column(),
                self.obstypes[obstype].get_v_column(),
            ]
        )

        # add the aggregated obstypes to the known obsytpes
        self.obstypes[amp_obstype.name] = amp_obstype
        self.obstypes[dir_obstype.name] = dir_obstype

    def get_gee_dataset_data(
        self, mapname, metadf, startdt_utc, enddt_utc, obstypes=["temp"]
    ):
        """Extract timeseries of a gee dataset.

        The extraction can only be done if the gee dataset bandname (and units)
        corresponding to the obstype is known.

        The units are converted to the toolkit standard units!!

        Parameters
        ----------
        mapname : str
            Mapname of choice of the GEE dataset to extract data from.
        metadf : pandas.DataFrame
            A dataframe with a 'name' index and  'lat', 'lon' columns.
            Timeseries are extracted for these locations.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstypes : str or list of strings, optional
            Toolkit observation type to extract data from. There should be a
            bandname mapped to this obstype for the gee map. Multiple obstypes
            can be given in a list. The default is 'temp'.


        Returns
        -------
        None.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        """
        # ====================================================================
        # Test input
        # ====================================================================
        if metadf.empty:
            logger.warning("The metadf is empty!")
            return

        # Subset metadf to stations with coordinates
        no_coord_meta = metadf[metadf[["lat", "lon"]].isna().any(axis=1)]
        if not no_coord_meta.empty:
            logger.warning(
                f"Following stations do not have coordinates, and thus no modeldata extraction is possible: {no_coord_meta.index.to_list()}"
            )
            metadf = metadf[~metadf[["lat", "lon"]].isna().any(axis=1)]

        # is mapinfo available
        if mapname not in self.mapinfo.keys():
            logger.warning(f"{mapname} is not a known gee dataset.")
            return

        geeinfo = self.mapinfo[mapname]

        # does dataset contain time evolution
        if not geeinfo["dynamical"]:
            logger.warning(
                f"{mapname} is a static dataset, this method does not work on static datasets"
            )
            return

        # Check obstypes
        if isinstance(obstypes, str):
            obstypes = [obstypes]  # convert to list

        for obstype in obstypes:
            # is obstype mapped?
            if obstype not in self.obstypes.keys():
                logger.warning(
                    f"{obstype} is an unknown observation type of the modeldata."
                )
                return
            if not self.obstypes[obstype].has_mapped_band(mapname):
                logger.warning(
                    f"{obstype} is not yet mapped to a bandname in the {mapname} dataset."
                )
                return

        # ====================================================================
        # GEE api extraction
        # ====================================================================

        # Connect to Gee
        connect_to_gee()

        # Get bandname mapper ({bandname1: obstypename1, ...})
        band_mapper = {}
        for obstype in obstypes:
            band_mapper.update(self.obstypes[obstype].get_bandname_mapper(mapname))

        logger.info(f"{band_mapper} are extracted from {mapname}.")
        # Get data using GEE
        df = gee_extract_timeseries(
            metadf=metadf,
            band_mapper=band_mapper,
            mapinfo=geeinfo,
            startdt=startdt_utc,
            enddt=enddt_utc,
            latcolname="lat",
            loncolname="lon",
        )

        self.df = df
        self.modelname = mapname

        if not self.df.empty:
            # Set extra attributes
            self.df_tz = "UTC"

            # convert to standard units
            for obstype in obstypes:
                # Set the current unit (at this point, it is the unit as define in the band)
                self.obstypes[obstype].setup_current_data_unit(mapname=mapname)
                self._convert_units_to_tlk(obstype)
                if isinstance(self.obstypes[obstype], ModelObstype_Vectorfield):
                    self._exploid_2d_vector_field(obstype)
        else:
            self._data_stored_at_drive = True

    def get_ERA5_data(self, metadf, startdt_utc, enddt_utc, obstypes="temp"):
        """Extract timeseries of the ERA5_hourly dataset.

        The units are converted to the toolkit standard units.

        (This method is a specific ERA5_hourly wrapper on the
         get_gee_dataset_data() method)

        Parameters
        ----------
        metadf : pandas.DataFrame
            A dataframe with a 'name' index and  'lat', 'lon' columns.
            Timeseries are extracted for these locations.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        obstypes : str or list of str, optional
            Toolkit observation type to extract data from. There should be a
            bandname mapped to this obstype for the gee map. Multiple
            observation types can be extracted if given as a list. The default is
            'temp'.


        Returns
        -------
        None.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        written to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .set_model_from_csv()
        method.

        """
        # Check obstypes
        if isinstance(obstypes, str):
            obstypes = [obstypes]  # convert to list

        # test if obstype is known
        for obstype in obstypes:
            if obstype not in self.obstypes:
                sys.exit(f"{obstype} is not a known obstype of the Modeldata instance.")

            # test if the obstype is mapped in the era5 hourly dataset
            if "ERA5_hourly" not in self.obstypes[obstype].get_mapped_datasets():
                sys.exit(
                    f"{obstype} has no equivalent mapped band for the ERA5_hourly dataset."
                )

        self.get_gee_dataset_data(
            mapname="ERA5_hourly",
            metadf=metadf,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            obstypes=obstypes,
        )

    # def _interpolate_modeldata(self, to_multiidx):
    #     """Interpolate modeldata in time.

    #     Interpolate the modeldata timeseries, to a given name-datetime
    #     multiindex.

    #     The modeldata will be converted to the timezone of the multiindex.

    #     If no interpolation can be done, Nan values are used.

    #     Parameters
    #     ----------
    #     to_multiidx : pandas.MultiIndex
    #         A name - datetime (tz-aware) multiindex to interpolate the
    #         modeldata timeseries to.

    #     Returns
    #     -------
    #     returndf : pandas.DataFrame
    #         A dataframe with to_multiidx as an index.
    #         The values are the interpolated values.

    #     """
    #     returndf = init_multiindexdf()

    #     recordsdf = init_multiindexdf()
    #     recordsdf.index = to_multiidx
    #     # iterate over stations check to avoid extrapolation is done per stations
    #     for sta in recordsdf.index.get_level_values("name").unique():
    #         sta_recordsdf = xs_save(recordsdf, sta, level="name", drop_level=False)
    #         sta_moddf = xs_save(self.df, sta, level="name", drop_level=False)

    #         if sta_moddf.empty:
    #             logger.warning(f"There are not modeldata records for {sta}!")
    #             # empyt sta_moddg --> not model data --> empyt return
    #             mergedf = sta_recordsdf  # overload the records index
    #             mergedf = mergedf.reindex(
    #                 columns=list(self.df.columns)
    #             )  # add all present obstypes as nan columns

    #         else:
    #             # convert modeldata to timezone of observations
    #             sta_moddf = conv_tz_multiidxdf(
    #                 df=sta_moddf,
    #                 timezone=sta_recordsdf.index.get_level_values("datetime").tz,
    #             )

    #             # check if modeldata is will not be extrapolated !
    #             if min(sta_recordsdf.index.get_level_values("datetime")) < min(
    #                 sta_moddf.index.get_level_values("datetime")
    #             ):
    #                 logger.warning("Modeldata will be extrapolated")
    #             if max(sta_recordsdf.index.get_level_values("datetime")) > max(
    #                 sta_moddf.index.get_level_values("datetime")
    #             ):
    #                 logger.warning("Modeldata will be extrapolated")

    #             # combine model and records
    #             mergedf = sta_recordsdf.merge(
    #                 sta_moddf, how="outer", left_index=True, right_index=True
    #             )

    #             # reset index for time interpolation
    #             mergedf = mergedf.reset_index().set_index("datetime").sort_index()

    #             # interpolate missing modeldata
    #             mergedf = mergedf.drop(columns=["name"])
    #             mergedf.interpolate(method="time", limit_area="inside", inplace=True)
    #             mergedf["name"] = sta

    #             # convert back to multiindex
    #             mergedf = (
    #                 mergedf.reset_index().set_index(["name", "datetime"]).sort_index()
    #             )
    #             # filter only records
    #             mergedf = mergedf.loc[sta_recordsdf.index]

    #         returndf = pd.concat([returndf, mergedf])
    #     return returndf

    # def make_plot(
    #     self,
    #     obstype_model="temp",
    #     dataset=None,
    #     obstype_dataset=None,
    #     stationnames=None,
    #     starttime=None,
    #     endtime=None,
    #     title=None,
    #     show_outliers=True,
    #     show_filled=True,
    #     legend=True,
    #     _ax=None,  # needed for GUI, not recommended use
    # ):
    #     """Plot timeseries of the modeldata.

    #     This function creates a timeseries plot for the Modeldata. When a
    #     metobs_toolkit.Dataset is provided, it is plotted in the same figure.

    #     The line colors represent the timesries for different locations.

    #     Parameters
    #     ----------
    #     obstype_model : string, optional
    #          Fieldname of the Modeldata to visualise. The default is 'temp'.
    #     dataset : metobs_toolkit.Dataset, optional
    #         A Dataset instance with observations plotted in the same figure.
    #         Observations are represented by solid line and modeldata by dashed
    #         lines. The default is None.
    #     obstype_dataset : string, optional
    #         Fieldname of the Dataset to visualise. Only relevent when a dataset
    #         is provided. If None, obsype_dataset = obstype_model. The default
    #         is None.
    #     stationnames : list, optional
    #         A list with stationnames to include in the timeseries. If None is
    #         given, all the stations are used, defaults to None.
    #     starttime : datetime.datetime, optional
    #          Specifiy the start datetime for the plot. If None is given it will
    #          use the start datetime of the dataset, defaults to None.
    #     endtime : datetime.datetime, optional
    #          Specifiy the end datetime for the plot. If None is given it will
    #          use the end datetime of the dataset, defaults to None.
    #     title : string, optional
    #          Title of the figure, if None a default title is generated. The
    #          default is None.
    #     show_outliers : bool, optional
    #          If true the observations labeld as outliers will be included in
    #          the plot. Only relevent when a dataset is provided. The default
    #          is True.
    #     show_filled : bool, optional
    #          If true the filled values for gaps and missing observations will
    #          be included in the plot. Only relevent when a dataset is provided.
    #          The default is True.
    #     legend : bool, optional
    #          If True, a legend is added to the plot. The default is True.

    #     Returns
    #     -------
    #     axis : matplotlib.pyplot.axes
    #          The timeseries axes of the plot is returned.

    #     """
    #     logger.info(f"Make {obstype_model}-timeseries plot of model data")

    #     # Basic test
    #     if obstype_model not in self.df.columns:
    #         logger.warning(
    #             f"{obstype_model} is not foud in the modeldata df (columns = {self.df.columns})."
    #         )
    #         return
    #     if self.df.empty:
    #         logger.warning("The modeldata is empty.")
    #         return
    #     if obstype_dataset is None:
    #         obstype_dataset = obstype_model

    #     if dataset is not None:
    #         if obstype_dataset not in dataset.df.columns:
    #             logger.warning(f"{obstype_dataset} is not foud in the Dataframe df.")
    #             return

    #     model_df = self.df

    #     # ------ filter model ------------

    #     # Filter on obstype
    #     model_df = model_df[[obstype_model]]

    #     # Subset on stationnames
    #     if stationnames is not None:
    #         model_df = model_df[
    #             model_df.index.get_level_values("name").isin(stationnames)
    #         ]

    #     # Subset on start and endtime
    #     model_df = multiindexdf_datetime_subsetting(model_df, starttime, endtime)

    #     #  -------- Filter dataset (if available) -----------
    #     if dataset is not None:
    #         # combine all dataframes
    #         mergedf = dataset.get_full_status_df()

    #         # subset to obstype
    #         mergedf = xs_save(mergedf, obstype_dataset, level="obstype")

    #         # Subset on stationnames
    #         if stationnames is not None:
    #             mergedf = mergedf[
    #                 mergedf.index.get_level_values("name").isin(stationnames)
    #             ]

    #         # Subset on start and endtime
    #         mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

    #     # Generate ylabel
    #     y_label = self.obstypes[obstype_model].get_plot_y_label(mapname=self.modelname)

    #     # Generate title
    #     title = f"{self.modelname}"
    #     if dataset is not None:
    #         title = f"{title} and {self.obstypes[obstype_dataset].name} observations."

    #     # make plot of the observations
    #     if dataset is not None:
    #         # make plot of the observations
    #         _ax, col_map = timeseries_plot(
    #             mergedf=mergedf,
    #             title=title,
    #             ylabel=y_label,
    #             colorby="name",
    #             show_legend=legend,
    #             show_outliers=show_outliers,
    #             show_filled=show_filled,
    #             settings=dataset.settings,
    #             _ax=_ax,
    #         )

    #         # Make plot of the model on the previous axes
    #         ax, col_map = model_timeseries_plot(
    #             df=model_df,
    #             obstype=obstype_model,
    #             title=title,
    #             ylabel=y_label,
    #             settings=self._settings,
    #             show_primary_legend=False,
    #             add_second_legend=True,
    #             _ax=_ax,
    #             colorby_name_colordict=col_map,
    #         )

    #     else:
    #         # Make plot of model on empty axes
    #         ax, _colmap = model_timeseries_plot(
    #             df=model_df,
    #             obstype=obstype_model,
    #             title=title,
    #             ylabel=y_label,
    #             settings=self._settings,
    #             show_primary_legend=legend,
    #             add_second_legend=False,
    #             _ax=_ax,
    #         )

    #     return ax
