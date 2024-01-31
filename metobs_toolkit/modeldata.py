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

# from metobs_toolkit.landcover_functions import connect_to_gee, gee_extract_timeseries

from metobs_toolkit.plotting_functions import model_timeseries_plot, timeseries_plot


from metobs_toolkit.obstypes import Obstype as Obstype_class
from metobs_toolkit.obstype_modeldata import era5_default_model_obstypes
from metobs_toolkit.obstype_modeldata import (
    ModelObstype,
    ModelObstype_Vectorfield,
)
from metobs_toolkit.gee_extractor import GeeExtractor
from metobs_toolkit.obstype_modeldata import compute_amplitude, compute_angle
from metobs_toolkit.settings import Settings

logger = logging.getLogger(__name__)

# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================


class Modeldata:
    """Class holding timeseries data and methods for modeldata."""

    def __init__(self, metadf, extractor=None):
        """Initialize modeldata."""
        self.df = init_multiindexdf()  # Holds the timeseries coming from a model

        metadf = _format_metadf(metadf)
        self.metadf = metadf  # Holds the metadata for where modeldate is extracted
        self.obstypes = []

        self.extractor = extractor
        if isinstance(extractor, GeeExtractor):
            # Default obstypes Dict name: Obstype-instance
            self.obstypes = era5_default_model_obstypes.copy()
        else:
            self.obstypes = []  # only default era5 modelobstypes are defined

        self._settings = Settings()
        self.df_tz = "UTC"  # the timezone of the datetimes stored in the df

    def __str__(self):
        """Print overview information of the modeldata."""
        if self.df.empty:
            return f"Empty Modeldata instance linked to {self.extractor}."
        n_stations = self.df.index.get_level_values("name").unique().shape[0]
        obstypes = self.df.columns.to_list()
        startdt = self.df.index.get_level_values("datetime").min()
        enddt = self.df.index.get_level_values("datetime").max()
        data_units = [self.obstypes[col].get_standard_unit() for col in self.df.columns]

        return f"Modeldata instance linked to {self.extractor}, containing: \n \
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

    def _get_present_obstypes(self):
        """Return a list of ModelObstypes that are present in the df."""
        if self.df.empty:
            return []
        return [
            self.obstypes[col] for col in self.df.columns if col in self.obstypes.keys()
        ]

    def add_obstype(self, Modelobstype):
        """
        TODO: update docstring
        Add a new Observation type for the current Modeldata.


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
        if isinstance(Modelobstype, ModelObstype_Vectorfield):
            assert (
                not Modelobstype.name in self.obstypes.keys()
            ), f"{Modelobstype.name} already present in the known ModelObstypes: {self.obstypes}"
            # add Modelobstype
            self.obstypes[Modelobstype.name] = Modelobstype
            logger.info(
                f"{Modelobstype.name} added to the known observation types as ModelObstype_Vectorfield."
            )

        elif isinstance(Modelobstype, ModelObstype):
            assert (
                not Modelobstype.name in self.obstypes.keys()
            ), f"{Modelobstype.name} already present in the known Modelobstypes: {self.obstypes}"
            # add Modelobstype
            self.obstypes[Modelobstype.name] = Modelobstype
            logger.info(
                f"{Modelobstype.name} added to the known observation types as ModelObstype."
            )

        else:
            if isinstance(Modelobstype, Obstype_class):
                sys.exit(
                    f"{Modelobstype} is not an instance of metobs_toolkit.Modelobstype_modeldata.ModelModelobstype. Convert the Modelobstype to a ModelModelobstype first."
                )

            sys.exit(
                f"{Modelobstype} is not an instance of metobs_toolkit.Modelobstype_modeldata.ModelModelobstype."
            )

    # =============================================================================
    # Importers
    # =============================================================================

    def import_from_gee(
        self, target_obstypes, start_utc, end_utc, gdrive_filename="era5_data"
    ):
        """Extract timeseries of a gee dataset.

        This method can only be used when a GeeExtractor is set as the extractor.

        This import takes care of unit conversions and exploiting of vector fields
        to scalar fields (components, direction and amplitude).


        Parameters
        ----------
        target_obstypes : ModelObstype, ModelObstype_Vectorfield, name of a
            known obstype, or a list of these.
            The observations to extract from the extractor.
        startdt_utc : datetime.datetime
            Start datetime of the timeseries in UTC.
        enddt_utc : datetime.datetime
            Last datetime of the timeseries in UTC.
        gdrive_filename : str, optional
            Filename to use when the data is writen to a csv file on your
            Google Drive. The default is 'era5_data'.

        Returns
        -------
        None.

        Note
        ------
        When extracting large amounts of data, the timeseries data will be
        writen to a file and saved on your google drive. In this case, you need
        to provide the Modeldata with the data using the .import_gee_data_from_csv()
        method.

        """
        logger.info(f"Importing {target_obstypes} from gee.")
        # tests if gee is set as extractor
        assert isinstance(
            self.extractor, GeeExtractor
        ), f"Importing a gee dataset from csv is only valid if the extractor is a GeeExtractor. The current extractor is {self.extractor}"

        # Cast to a list of obstypes
        if isinstance(target_obstypes, Obstype_class):
            target_obstypes_list = [target_obstypes]

        elif isinstance(target_obstypes, str):
            assert (
                target_obstypes in self.obstypes.keys()
            ), f"{target_obstypes} not found in known ModelObtypes: {self.obstypes}"
            target_obstypes_list = [self.obstypes[target_obstypes]]

        elif isinstance(target_obstypes, list):
            target_obstypes_list = []
            for obs in target_obstypes:
                if isinstance(obs, str):
                    assert (
                        obs in self.obstypes.keys()
                    ), f"{obs} not found in known ModelObtypes: {self.obstypes}"
                    target_obstypes_list.append(self.obstypes[obs])
                elif isinstance(obs, Obstype_class):
                    target_obstypes_list.append(obs)
                else:
                    sys.exit(f"{obs} not a valid input type for target_obstypes.")
        else:
            sys.exit(f"{target_obstypes} not a valid input type for target_obstypes.")
        logger.debug(
            f"Following obstypes will be extracted from gee: {target_obstypes_list}"
        )
        assert bool(target_obstypes_list), "No target_obstypes are found."

        # get bandnmaes for target obstypes
        bandnames_map = {}  # keys are bandnames, values are the obstype names
        for obs in target_obstypes_list:
            bandnames_map.update(obs.get_bandname_mapper())

        # Extract timeseries
        df = self.extractor.extract_timeseries(
            metadf=self.metadf,
            bandnames=list(bandnames_map.keys()),
            start_utc=start_utc,
            end_utc=end_utc,
            gdrive_filename=gdrive_filename,
        )
        self.df = df
        if not self.df.empty:
            self.df = self._update_df_attribute_from_gee_df(
                rawdf=self.df, expected_obstypes=target_obstypes_list
            )

    def import_gee_data_from_csv(self, csvpath):
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
        # tests if gee is set as extractor
        assert isinstance(
            self.extractor, GeeExtractor
        ), f"Importing a gee dataset from csv is only valid if the extractor is a GeeExtractor. The current extractor is {self.extractor}"

        # 1. Read csv and set timezone
        df = pd.read_csv(csvpath, sep=",")
        # renamte name column
        df = df.rename(columns={"feature_idx": "name"})
        # format datetime
        df["datetime"] = pd.to_datetime(df["datetime"], format="%Y%m%d%H%M%S")
        # (assume all gee dataset are in UTC)
        df["datetime"] = df["datetime"].dt.tz_localize("UTC")

        # 2. Format dataframe
        # format index
        df = df.set_index(["name", "datetime"])
        df = df.sort_index()

        # test if all names in the df corresponds with names of the metadf
        assert set(df.index.get_level_values("name")) == set(
            self.metadf.index
        ), f"The station names of the csv file are not equal to those stored in the metadf attribute:  {set(df.index.get_level_values('name'))} != {set(self.metadf.index)}"

        # Find wicht obstypes are found in the dataframe
        found_obstypes = []
        for obs in self.obstypes.values():
            bandnames = list(obs.get_bandname_mapper().keys())
            for bandname in bandnames:
                if bandname in df.columns:
                    found_obstypes.append(obs)
        found_obstypes = list(set(found_obstypes))

        self.df = self._update_df_attribute_from_gee_df(
            rawdf=df, expected_obstypes=found_obstypes
        )

    def _update_df_attribute_from_gee_df(self, rawdf, expected_obstypes):
        """Rename columns, convert to std units and create aggregated fields
        if vector field detected."""

        for obs in expected_obstypes:
            if isinstance(obs, ModelObstype):
                # rename column
                rawdf = rawdf.rename(columns=obs.get_bandname_mapper())
                # convert to standard units
                rawdf[obs.name] = obs.convert_to_standard_units(
                    input_data=rawdf[obs.name], input_unit=obs.get_bandunit()
                )
            elif isinstance(obs, ModelObstype_Vectorfield):
                # rename columns
                rawdf = rawdf.rename(columns=obs.get_bandname_mapper())

                # ------ Amplitude  --------
                # calculate amplitude field
                ampseries, amp_obstype = compute_amplitude(
                    modelobs_vectorfield=obs, df=rawdf
                )
                self.obstypes[amp_obstype.name] = amp_obstype

                # convert amplitudes to standard units
                ampseries = amp_obstype.convert_to_standard_units(
                    input_data=ampseries,  # test this
                    input_unit=amp_obstype.get_bandunit(),
                )

                # add amplitude clumn
                rawdf[amp_obstype.name] = ampseries

                # ------ Direction  --------
                # calculate angle field
                angleseries, angle_obstype = compute_angle(
                    modelobs_vectorfield=obs, df=rawdf
                )
                self.obstypes[angle_obstype.name] = angle_obstype

                # convert angles to standard units
                # - not needed because angles are always created in degrees -

                # add amplitude clumn
                rawdf[angle_obstype.name] = angleseries

                # ----- Keep the components --------
                # we can keep the components, but we need to create two new obstypes
                # since we will interpret them as scalar 2d fields
                u_comp_obstype, v_comp_obstype = obs.create_the_scalar_modelobstypes()
                self.obstypes[u_comp_obstype.name] = u_comp_obstype
                self.obstypes[v_comp_obstype.name] = v_comp_obstype

            else:
                sys.exit("Something went wrong. Report this error.")
        return rawdf

    # =============================================================================
    # IO
    # =============================================================================
    def save_modeldata(
        self, outputfolder=None, filename="saved_modeldata.pkl", overwrite=False
    ):
        """Save a Modeldata instance to a (pickle) file.

        Parameters
        ----------
        outputfolder : str or None, optional
            The path to the folder to save the file. If None, the outputfolder
            from the Settings is used. The default is None.
        filename : str, optional
            The name of the output file. The default is 'saved_modeldata.pkl'.
        overwrite : bool, optional
            If the target file already exists, it will be overwritten if True,
            else an error is thrown. The default is False

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
        if os.path.isfile(full_path):
            if overwrite:
                # remove preexisting file in advance
                os.remove(full_path)
            else:
                sys.exit(f"{full_path} is already a file!")

        with open(full_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

        print(f"Modeldata saved in {full_path}")
        logger.info(f"Modeldata saved in {full_path}")

    # =============================================================================
    # Plotters
    # =============================================================================

    # =============================================================================
    # Helpers
    # =============================================================================
    def sample_data_as(self, target, interp_method="time", **kwargs):
        """Resample the modeldata to a target.

        This methods will convert the modeldata to the same index as present
        in the target_df. Thus this method will transform timezones (if needed)
        and will interpolate to the missing timestamps (required for
        timeresolution mismatch between target_df and the model observations).

        This method is only applicable when target_df and self.df share the
        same levels of the index.

        Parameters
        ----------
        target : pandas.DataFrame, or pandas.index
            The pandas dataframe with the target index to convert the modeldata
            to.
        interp_method : str, optional
            The interpolation method used by pandas.DataFrame.interpolate(). The default is 'time'.
        **kwargs:
            Additional keyword arguments passed to pandas.DataFrame.interpolate().
        Returns
        -------
        target_modeldf : pandas.DataFrame
            Modeldata resampled to target index.

        """
        assert (
            not self.df.empty
        ), f"Sample_data_as not possible because no modeldata found."

        if isinstance(target, pd.DataFrame):
            target_idx = target.index
        elif isinstance(target, pd.Index):
            target_idx = target

        assert set(target_idx.names) == set(
            self.df.index.names
        ), f"sample data to target not possible since index levels are not the same for modeldata and target dataframe."

        # Convert model to tz of the target
        target_tz = target_idx.get_level_values("datetime").tz
        target_modeldf = (
            self.df.reset_index()
            .set_index("datetime")
            .tz_convert(target_tz)
            .reset_index()
            .set_index(["name", "datetime"])
        )

        # get target idices not in modeldf and convert to dataframe
        tg_idx_new = target_idx[~target_idx.isin(target_modeldf.index)]
        tg_idx_new_df = pd.DataFrame(index=tg_idx_new)

        # add the missing idices to the modeldata with Nans as values
        target_modeldf = pd.concat([target_modeldf, tg_idx_new_df]).sort_index()

        # iterpolate the nan values (groupby name first)
        target_modeldf = target_modeldf.reset_index().set_index("datetime")

        target_modeldf = target_modeldf.groupby("name", group_keys=True).apply(
            lambda group: group.interpolate(
                method=interp_method, limit_area="inside", axes="columns", **kwargs
            )
        )

        target_modeldf = target_modeldf.drop(columns=["name"])
        # subset to target idices
        target_modeldf = target_modeldf.loc[target_idx]

        return target_modeldf

    # def _conv_to_timezone(self, tzstr):
    #     """Convert the timezone of the datetime index of the df attribute.

    #     Parameters
    #     ----------
    #     tzstr : str
    #         TImezonstring from the pytz module.

    #     Returns
    #     -------
    #     None.

    #     """
    #     # get tzstr by datetimindex.tz.zone

    #     df = self.df
    #     df["datetime_utc"] = df.index.get_level_values("datetime").tz_convert(tzstr)
    #     df = df.reset_index()
    #     df = df.drop(columns=["datetime"])
    #     df = df.rename(columns={"datetime_utc": "datetime"})
    #     df = df.set_index(["name", "datetime"])
    #     self.df = df
    #     self.df_tz = tzstr

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
        obstype_model : string or ModelObstype, optional
            The observation to plot. The default is 'temp'.
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

        # type casters
        if isinstance(obstype_model, str):
            assert (
                obstype_model in self.obstypes
            ), f"{obstype_model} not in the known obstypes of the Modeldata: {self.obstypes}"
            obstype_model = self.obstypes[obstype_model]
        if isinstance(obstype_dataset, str):
            assert not (
                dataset is None
            ), "You cannot provide a obstype_dataset without providing a dataset"
            assert (
                obstype_dataset in dataset.obstypes
            ), f"{obstype_dataset} not in the known obstypes of the Dataset: {dataset.obstypes}"
            obstype_dataset = dataset.obstypes[obstype_dataset]
            assert (
                obstype_dataset.name in dataset.df.columns
            ), f"{obstype_dataset} not found in the Modeldata data: {dataset.df.columns}"

        assert not self.df.empty, "Make plot not possible for empty Modeldata"

        # Check of obstype is present in the df
        assert (
            obstype_model.name in self.df.columns
        ), f"{obstype_model} not found in the Modeldata data: {self.df.columns}"

        # Basic test
        if obstype_dataset is None:
            obstype_dataset = obstype_model

        if dataset is not None:
            assert (
                obstype_dataset.name in dataset.df.columns
            ), f"{obstype_dataset} is not foud in the Dataframe df."

        model_df = self.df

        # ------ filter model ------------

        # Filter on obstype
        model_df = model_df[[obstype_model.name]]

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
            mergedf = xs_save(mergedf, obstype_dataset.name, level="obstype")

            # Subset on stationnames
            if stationnames is not None:
                mergedf = mergedf[
                    mergedf.index.get_level_values("name").isin(stationnames)
                ]

            # Subset on start and endtime
            mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)

        # Generate ylabel
        y_label = obstype_model.get_plot_y_label()

        # Generate title
        title = f"{self.extractor.usage} data for {obstype_dataset.name} from {self.extractor.location}"
        if dataset is not None:
            title = f" {title} \n and {obstype_dataset.name} observations."

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
                obstypename=obstype_model.name,
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
                obstypename=obstype_model.name,
                title=title,
                ylabel=y_label,
                settings=self._settings,
                show_primary_legend=legend,
                add_second_legend=False,
                _ax=_ax,
            )

        return ax


def import_modeldata(target_pkl_file):
    """Import a modeldata instance from a (pickle) file.

    Parameters
    ----------
    target_pkl_file : str,
        The path to the target pkl file to import.
    Returns
    -------
    metobs_toolkit.Modeldata
        The modeldata instance.

    """
    # # check if folder_path is known and exists
    # if folder_path is None:
    #     folder_path = self.settings.IO["output_folder"]
    #     assert (
    #         folder_path is not None
    #     ), "No folder_path is given, and no outputfolder is found in the settings."

    # assert os.path.isdir(folder_path), f"{folder_path} is not a directory!"

    # full_path = os.path.join(folder_path, filename)

    # check if file exists
    assert os.path.isfile(target_pkl_file), f"{target_pkl_file} does not exist."

    with open(target_pkl_file, "rb") as inp:
        modeldata = pickle.load(inp)

    return modeldata


def _format_metadf(metadf):
    _ = _is_metadf_valid(metadf)

    if metadf.index.name is None:
        if not isinstance(metadf.index, pd.MultiIndex):
            metadf.index.name = "name"
    return metadf


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
