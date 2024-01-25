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

from metobs_toolkit.df_helpers import (
    init_multiindexdf,
    conv_tz_multiidxdf,
    xs_save,
    multiindexdf_datetime_subsetting,
)

from metobs_toolkit.landcover_functions import connect_to_gee, gee_extract_timeseries

from metobs_toolkit.plotting_functions import model_timeseries_plot, timeseries_plot

# from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.obstypes import Obstype as Obstype_class
from metobs_toolkit.obstype_modeldata import (
    model_obstypes,
    ModelObstype,
    ModelObstype_Vectorfield,
)
from metobs_toolkit.obstype_modeldata import compute_amplitude, compute_angle
from metobs_toolkit.settings import Settings

logger = logging.getLogger(__name__)

# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================


class Modeldata:
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

    def add_gee_dataset(
        self,
        mapname,
        gee_location,
        obstype,
        bandname,
        units,
        scale,
        band_desc=None,
        time_res="1H",
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
        Celcius. This will be implemented in the futur.

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

    def convert_units_to_tlk(self, obstype):
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

        cur_unit = self.obstypes[obstype].get_modelunit(self.modelname)

        if isinstance(self.obstypes[obstype], ModelObstype):
            converted_data = self.obstypes[obstype].convert_to_standard_units(
                input_data=self.df[obstype], input_unit=cur_unit
            )
            # Update the data and the current unit
            self.df[obstype] = converted_data
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

    def exploid_2d_vector_field(self, obstype):
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
        writen to a file and saved on your google drive. In this case, you need
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
            self.df_tz = "UTC"
            # convert to standard units
            for obstype in obstypes:
                self.convert_units_to_tlk(obstype)
                if isinstance(self.obstypes[obstype], ModelObstype_Vectorfield):
                    self.exploid_2d_vector_field(obstype)
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
        writen to a file and saved on your google drive. In this case, you need
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

    def save_modeldata(
        self,
        outputfolder=None,
        filename="saved_modeldata.pkl",
    ):
        """Save a Modeldata instance to a (pickle) file.

        Parameters
        ----------
        outputfolder : str or None, optional
            The path to the folder to save the file. If None, the outputfolder
            from the Settings is used. The default is None.
        filename : str, optional
            The name of the output file. The default is 'saved_modeldata.pkl'.

        Returns
        -------
        None.

        """
        # check if outputfolder is known and exists
        if outputfolder is None:
            outputfolder = self.settings.IO["output_folder"]
            assert (
                outputfolder is not None
            ), "No outputfolder is given, and no outputfolder is found in the settings."

        assert os.path.isdir(outputfolder), f"{outputfolder} is not a directory!"

        # check file extension in the filename:
        if filename[-4:] != ".pkl":
            filename += ".pkl"

        full_path = os.path.join(outputfolder, filename)

        # check if file exists
        assert not os.path.isfile(full_path), f"{full_path} is already a file!"

        with open(full_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

        print(f"Modeldata saved in {full_path}")
        logger.info(f"Modeldata saved in {full_path}")

    def import_modeldata(self, folder_path=None, filename="saved_modeldata.pkl"):
        """Import a modeldata instance from a (pickle) file.

        Parameters
        ----------
        folder_path : str or None, optional
            The path to the folder to save the file. If None, the outputfolder
            from the Settings is used. The default is None.
        filename : str, optional
            The name of the output file. The default is 'saved_modeldata.pkl'.

        Returns
        -------
        metobs_toolkit.Modeldata
            The modeldata instance.

        """
        # check if folder_path is known and exists
        if folder_path is None:
            folder_path = self.settings.IO["output_folder"]
            assert (
                folder_path is not None
            ), "No folder_path is given, and no outputfolder is found in the settings."

        assert os.path.isdir(folder_path), f"{folder_path} is not a directory!"

        full_path = os.path.join(folder_path, filename)

        # check if file exists
        assert os.path.isfile(full_path), f"{full_path} does not exist."

        with open(full_path, "rb") as inp:
            modeldata = pickle.load(inp)

        return modeldata

    def set_model_from_csv(self, csvpath):
        """Import timeseries data that is stored in a csv file.

        The name of the gee dataset the timeseries are coming from must be the
        same as the .modelname attribute of the Modeldata.


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
        # tests ----
        if self.modelname not in self.mapinfo.keys():
            logger.warning(f"{self.modelname} is not found in the gee datasets.")
            return

        # 1. Read csv and set timezone
        df = pd.read_csv(csvpath, sep=",")
        # format datetime
        df["datetime"] = pd.to_datetime(df["datetime"], format="%Y%m%d%H%M%S")
        # (assume all gee dataset are in UTC)
        df["datetime"] = df["datetime"].dt.tz_localize("UTC")

        # 2. Format dataframe
        # format index
        df = df.set_index(["name", "datetime"])
        df = df.sort_index()

        # make a bandname --> tlk name mapper
        bandname_mapper = {}
        for known_obstype in self.obstypes.values():
            bandname_mapper.update(known_obstype.get_bandname_mapper(self.modelname))

        # rename to values to toolkit space
        df = df.rename(columns=bandname_mapper)

        # 3. update attributes
        self.df = df
        self.df_tz = "UTC"

        # 4. Find which obstypes are present
        data_present_obstypes = []
        for col in self.df.columns:
            if col in self.obstypes.keys():
                # column is a regular obstype
                data_present_obstypes.append(col)
            else:
                # check if column represents a vector component
                for known_obs in self.obstypes.values():
                    if isinstance(known_obs, ModelObstype_Vectorfield):
                        comps = [known_obs.get_u_column(), known_obs.get_v_column()]
                        if col in comps:
                            data_present_obstypes.append(known_obs.name)
        data_present_obstypes = list(set(data_present_obstypes))
        # A. scalar obstypes (same name as column)

        # 5. Convert units
        for obstype in data_present_obstypes:
            self.convert_units_to_tlk(obstype)
            if isinstance(self.obstypes[obstype], ModelObstype_Vectorfield):
                self.exploid_2d_vector_field(obstype)

    def interpolate_modeldata(self, to_multiidx):
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
            sta_moddf = xs_save(self.df, sta, level="name", drop_level=False)

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
            mergedf.interpolate(method="time", limit_area="inside", inplace=True)
            mergedf["name"] = sta
            # convert back to multiindex
            mergedf = mergedf.reset_index().set_index(["name", "datetime"]).sort_index()
            # filter only records
            mergedf = mergedf.loc[sta_recordsdf.index]

            returndf = pd.concat([returndf, mergedf])
        return returndf

    def make_plot(
        self,
        obstype_model="temp",
        dataset=None,
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
        dataset : metobs_toolkit.Dataset, optional
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
        if obstype_model not in self.df.columns:
            logger.warning(
                f"{obstype_model} is not foud in the modeldata df (columns = {self.df.columns})."
            )
            return
        if self.df.empty:
            logger.warning("The modeldata is empty.")
            return
        if obstype_dataset is None:
            obstype_dataset = obstype_model

        if dataset is not None:
            if obstype_dataset not in dataset.df.columns:
                logger.warning(f"{obstype_dataset} is not foud in the Dataframe df.")
                return

        model_df = self.df

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
        if dataset is not None:
            # combine all dataframes
            mergedf = dataset.combine_all_to_obsspace()

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
        y_label = self.obstypes[obstype_model].get_plot_y_label(mapname=self.modelname)

        # Generate title
        title = f"{self.modelname}"
        if dataset is not None:
            title = f"{title} and {self.obstypes[obstype_dataset].name} observations."

        # make plot of the observations
        if dataset is not None:
            # make plot of the observations
            _ax, col_map = timeseries_plot(
                mergedf=mergedf,
                title=title,
                ylabel=y_label,
                colorby="name",
                show_legend=legend,
                show_outliers=show_outliers,
                show_filled=show_filled,
                settings=dataset.settings,
                _ax=_ax,
            )

            # Make plot of the model on the previous axes
            ax, col_map = model_timeseries_plot(
                df=model_df,
                obstype=obstype_model,
                title=title,
                ylabel=y_label,
                settings=self._settings,
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
                settings=self._settings,
                show_primary_legend=legend,
                add_second_legend=False,
                _ax=_ax,
            )

        return ax
