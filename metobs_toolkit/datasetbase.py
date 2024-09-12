#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 14 11:28:20 2024

@author: thoverga
"""

import copy
import logging
import numpy as np
import pandas as pd

import datetime as datetimemodule
from metobs_toolkit.printing import dataset_string_repr
from metobs_toolkit.settings_files.default_formats_settings import label_def
from metobs_toolkit.qc_checks import duplicate_timestamp_check
from metobs_toolkit.df_helpers import (
    # multiindexdf_datetime_subsetting,
    # fmt_datetime_argument,
    # init_multiindex,
    # init_multiindexdf,
    # init_triple_multiindexdf,
    empty_outliers_df,
    metadf_to_gdf,
    # conv_applied_qc_to_df,
    # get_freqency_series,
    get_likely_frequency,
    _simplify_time,
    # value_labeled_doubleidxdf_to_triple_idxdf,
    xs_save,
    # concat_save,
)
from metobs_toolkit.template import Template
from metobs_toolkit.settings import Settings
from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.modeldata import default_datasets
from metobs_toolkit.gap import find_gaps

logger = logging.getLogger(__name__)


# %%
class DatasetBase(object):
    """Abstract class for holding data and metadata."""

    def __init__(self):
        """Construct all the necessary attributes for Dataset object."""
        logger.info("Initialise dataset")

        # Dataset with 'good' observations
        self.df = pd.DataFrame()

        # Dataset with outlier observations (same structure as df, and all
        # these records must be present (with nan values) in the df)
        self.outliersdf = empty_outliers_df()

        # for interpretation of qc effectifness
        self._applied_qc = pd.DataFrame(columns=["obstype", "checkname"])

        # Dataset with metadata (static)
        self.metadf = pd.DataFrame()

        # dictionary storing present observationtypes
        self.obstypes = copy.copy(tlk_obstypes)  # init with all tlk obstypes
        self.settings = copy.deepcopy(Settings())

        # Gaps are stored as a list of Gap()
        self.gaps = None

        # Template
        self.template = Template()

        # GEE datasets defenitions
        self.gee_datasets = copy.deepcopy(default_datasets)

    # =============================================================================
    # Specials
    # =============================================================================
    def __str__(self):
        """Represent as text."""

        return dataset_string_repr(self)

    def __repr__(self):
        """Info representation."""
        class_name = type(self).__name__
        return f"Instance of {class_name} at {hex(id(self))}"

    # =============================================================================
    #   Argument checkers
    # =============================================================================

    def _timedelta_arg_check(self, timedeltaarg, none_is_none=True):
        if (none_is_none) & (timedeltaarg is None):
            return None

        if isinstance(timedeltaarg, datetimemodule.timedelta):
            dt = pd.Timedelta(timedeltaarg)
        elif isinstance(timedeltaarg, str):
            dt = pd.Timedelta(timedeltaarg)
        elif isinstance(timedeltaarg, pd.Timedelta):
            dt = dt
        else:
            MetobsDatasetBaseError(
                f"{timedeltaarg} could not be interpreted as a Timedelta, convert it to a pd.Timedelta()."
            )

        if dt == pd.Timedelta(0):
            logger.warning(f"A {dt} is given as an argument for a timedelta.")

        return dt

    def _datetime_arg_check(self, datetimearg, none_is_none=True):
        """Formats a datetime given by a user in an argument.

        Conversion to a pandas.Timestamp in the tz of the records.

        """

        if (none_is_none) & (datetimearg is None):
            return None
        # check type and cast to pd.Timestamp
        if isinstance(datetimearg, datetimemodule.datetime):
            dt = pd.Timestamp(datetimearg)
        elif isinstance(datetimearg, pd.Timestamp):
            dt = datetimearg
        else:
            raise MetobsDatasetBaseError(
                f"{datetimearg} is not in a datetime format (datetime.datetime or pandas.Timestamp)."
            )

        # check timezone and make tz-awer
        if dt.tz is None:
            # tz naive timestamp
            dt = dt.tz_localize(tz=self._get_tz())
        else:
            # tz aware timestamp --> convert to tz of records
            dt = dt.tz_convert(tz=self._get_tz())

        return dt

    # =============================================================================
    #     Setters
    # =============================================================================
    def _set_df(
        self,
        df,
        apply_structure_checks=True,
        apply_dup_checks=True,
        tz_aware_check=True,
    ):
        # TODO: run simple checks
        if apply_structure_checks:
            # Test index order
            if not list(df.index.names) == ["name", "obstype", "datetime"]:
                raise MetobsDatasetBaseError(
                    f"A dataframe is being set as Dataset.df with wrong df.index.names: {df.index.names}"
                )
            if not list(df.columns) == ["value"]:
                raise MetobsDatasetBaseError(
                    f"A dataframe is being set as Dataset.df with wrong df.columns: {list(df.columns)}"
                )
        if apply_dup_checks:
            if df.index.duplicated().any():
                raise MetobsDatasetBaseError(
                    f"A dataframe is being set as Dataset.df with duplicates in the index: \n {df.loc[df.index.duplicated()]}"
                )
        if tz_aware_check:
            if df.index.get_level_values("datetime").tz is None:
                raise MetobsDatasetBaseError(
                    f'A dataframe is being set as Dataset.df with a timezone-unaware datetimeindex: {df.index.get_level_values("datetime")}'
                )

        # set attr
        self.df = df

    def _set_metadf(self, metadf):
        if not metadf.index.name == "name":
            raise MetobsDatasetBaseError(
                f"A dataframe is being set as Dataset.metadf with wrong df.index.name: {metadf.index.name}"
            )

        if "lat" in metadf.columns:
            if (
                metadf["lat"]
                .dropna()
                .apply(lambda x: ((x < -90.0) | ((x > 90.0))))
                .any()
            ):
                raise MetobsDatasetBaseError(
                    f"The latitude coordinates in the metadata are not all in [-90; 90] range: {metadf['lat']}"
                )
        if "lon" in metadf.columns:
            if (
                metadf["lon"]
                .dropna()
                .apply(lambda x: ((x < 0.0) | ((x > 180.0))))
                .any()
            ):
                raise MetobsDatasetBaseError(
                    f"The longitude coordinates in the metadata are not all in [0; 180] range: {metadf['lon']}"
                )

        self.metadf = metadf

    def _set_outliersdf(
        self, outliersdf, apply_structure_checks=True, apply_dup_checks=True
    ):
        # TODO: run simple checks
        if apply_structure_checks:
            # Test index order
            if not list(outliersdf.index.names) == ["name", "obstype", "datetime"]:
                raise MetobsDatasetBaseError(
                    f"A dataframe is being set as Dataset.outliersdf with wrong df.index.names: {outliersdf.index.names}"
                )
            if not list(outliersdf.columns) == ["value", "label"]:
                raise MetobsDatasetBaseError(
                    f"A dataframe is being set as Dataset.outliersdf with wrong df.columns: {list(outliersdf.columns)}"
                )
        # set attr
        self.outliersdf = outliersdf

    def _set_obstypes(self, obstypes):
        # TODO: run simple checks
        self.obstypes = obstypes.copy()

    def _set_settings(self, settings):
        # TODO: run simple checks
        self.settings = copy.deepcopy(settings)  # deep copy?

    def _set_gaps(self, gapslist):
        # TODO: run simple checks
        self.gaps = gapslist

    def _set_gee_dataset(self, geedatasetlist):
        # clear all metadata and extracted timeseries form the geedatasets
        gee_dataset_dict = {}
        for gee_dataset in geedatasetlist:
            gee_dataset._clear_data()
            gee_dataset_dict[gee_dataset.name] = gee_dataset

        self.gee_datasets = gee_dataset_dict

    def _append_to_applied_qc(self, obstypename, checkname):
        """Add an observationtype to the _applied_qc.

        The applied qc is mainly used to save the order of applied checks,
        to validate the effectiveness of one check.

        Parameters
        ----------
        obstypename : str or list of str
            The names of the obstypes that are checked.
        checkname : str
            The name of the check (must be a key in label_def).

        Returns
        -------
        None.

        """

        if checkname not in label_def.keys():
            raise MetobsDatasetBaseError(f"{checkname} is not a known checkname.")
        # add it to the applied checks
        if isinstance(obstypename, str):
            obstypes = [obstypename]
        else:
            obstypes = list(obstypename)

        self._applied_qc = pd.concat(
            [
                self._applied_qc,
                pd.DataFrame(
                    data={"obstype": obstypes, "checkname": [checkname] * len(obstypes)}
                ),
            ]
        )

    # =============================================================================
    # Getters
    # =============================================================================
    def _get_tz(self):
        # IF no data --> tz is UTC (to work without data for gee fun)
        if self.df.empty:
            return "UTC"
        else:
            return self.df.index.get_level_values("datetime").tz

    def _get_present_obstypes(self):
        """Get all present obstypenames in the df.

        (This does not guarantee that these obstypesnames are known obstypes!)

        Returns
        -------
        list
            All present obstypenames.

        """
        return self.df.index.get_level_values("obstype").unique().to_list()

    def _get_timestamps_info(
        self,
        freq_estimation_method="median",
        freq_simplify_tolerance="2min",
        origin_simplify_tolerance="5min",
    ):
        """Get details of the time resolution for each station.

        It is assumed that ideally, records are perfect periodically. In
        practice this is not always the case. We can specify an ideal set
        of timestamps by a start, and end frequency.

        This method estimates, for each station:

            - Frequency: the frequency of a station is the highest frequency
            detected for all its observationtypes. The frequency estimate
            of an obsertype is done by a 'mean' or 'highest' approach on
            consecutive records. One can specify a tolerance to simplify
            the frequency estimates.

            - start timestamp: The first timestamp registered by a station,
            overall its observationtypes. One can specify a tolerance to
            simplify the origin.

            - end timestamp: The latest timestemp registerd by a station,
            over all its observationtypes. One can specify a tolerance to
            simplify it.

        This information is added to the metadf attribute.

        Note
        -----
        The assumption is made that all obstypes per station have the same start
        and end timestamp.

        """

        df = self.df

        # 1. Find frequencies per station ()
        if pd.Timedelta(freq_simplify_tolerance).seconds < 1:
            freq_simplify = False
        else:
            freq_simplify = True

        freqs_dict = {sta: [] for sta in df.index.get_level_values("name")}

        for groupidx, groupdf in df.reset_index().groupby(["name", "obstype"]):
            # calculate the frequecy of each station, and each of its obstypes
            groupdf = groupdf.set_index("datetime")

            freqs_dict[groupidx[0]].append(
                get_likely_frequency(
                    timestamps=groupdf.index,
                    method=freq_estimation_method,
                    simplify=freq_simplify,
                    max_simplify_error=freq_simplify_tolerance,
                )
            )

        # The target frequency of a station is set as the minimum of its obstypes
        freqs_target = {sta: min(val) for sta, val in freqs_dict.items()}

        # 2. Find origins and last timestamps

        origin_dict = {}
        last_timestamp_dict = {}
        for sta in df.index.get_level_values("name").unique():
            stadatetimes = df.xs(sta, level="name").index.get_level_values("datetime")

            # find origin
            naive_origin = stadatetimes.min()
            sta_origin = _simplify_time(
                time=naive_origin, max_simplyfi_error=origin_simplify_tolerance
            )
            origin_dict[sta] = sta_origin

            # find last timestamp
            naive_last = stadatetimes.max()
            last_timestamp = pd.Timestamp(
                sta_origin
                + (
                    int((naive_last - sta_origin) / freqs_target[sta])
                    * freqs_target[sta]
                )
            )
            last_timestamp_dict[sta] = last_timestamp

        # 3. update metadf
        self.metadf["dataset_resolution"] = pd.Series(freqs_target)
        self.metadf["dt_start"] = pd.Series(origin_dict)
        self.metadf["dt_end"] = pd.Series(last_timestamp_dict)

    # =============================================================================
    # Checks
    # =============================================================================

    def _data_is_required_check(self):
        """Raise an error when there is no observational data found."""
        if self.df.empty:
            if self.metadf.empty:
                raise MetobsDatasetBaseError(
                    f"There is no observational data stored in {self}"
                )
            else:
                raise MetobsDatasetBaseError(
                    f'There is no observational data stored in {self} \n(This method is not a "metadata-only" method.)'
                )
        return

    def _are_all_present_obstypes_knonw(self):
        """Raise an error if there are unknown obstypes in the df."""
        present_obstypes = self._get_present_obstypes()
        unknown_obs_cols = [
            obs_pres
            for obs_pres in present_obstypes
            if obs_pres not in self.obstypes.keys()
        ]
        if len(unknown_obs_cols) > 0:
            raise MetobsDatasetBaseError(
                f"The following observationtypes are found in the data, but are unknown obstypes: {unknown_obs_cols}"
            )

    # =============================================================================
    # Construction methods for creating a Dataset
    # =============================================================================

    def _construct_dataset(
        self,
        df,
        freq_estimation_method,
        freq_estimation_simplify_tolerance,
        origin_simplify_tolerance,
        timestamp_tolerance,
        use_metadata,
        # fixed_freq_series=None,
        # update_full_metadf=True,
    ):
        """Construct the Dataset class from a IO dataframe.

        1. Set the dataframe and metadataframe attributes
        2. Drop stations that have Nan as name.
        2. Find the duplicates (remove them from observations +  add them to outliers)
        3. Convert the values to standard units + update the observationtypes (some template specific attributes)
        5. Find gaps in the records (duplicates are excluded from the gaps)
        6. Get a frequency estimate per station
        7. Initiate the gaps (find missing records)
        8. Add the missing records to the dataframe


        Parameters
        ----------
        df : pandas.dataframe
            The dataframe containing the input observations and metadata.
        freq_estimation_method : 'highest' or 'median'
            Select which method to use for the frequency estimation. If
            'highest', the highest appearing frequency is used. If 'median', the
            median of the appearing frequencies is used.
        freq_estimation_simplify : bool
            If True, the likely frequency is converted to round hours or round minutes.
            The "freq_estimation_simplify_error' is used as a constraint. If the constraint is not met,
            the simplification is not performed.
        freq_estimation_simplify_error : Timedelta or str, optional
            The tolerance string or object represents the maximum translation in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1H', is one hour.
        fixed_freq_series : pandas.Series or None, optional
            If you do not want the frequencies to be recalculated, you can pass the
            frequency series to update the metadf["dataset_resolution"]. If None, the frequencies will be estimated. The default is None.
        update_full_metadf : bool, optional
            If True, the full Dataset.metadf will be updated. If False, only the frequency columns in the Dataset.metadf will be updated. The default is True.


        Returns
        -------
        None.

        """
        # Set the df attribute
        self._construct_df(dataframe=df)

        # Set the metadf attribute
        self._construct_metadf(dataframe=df, use_metadata=use_metadata)

        # Apply QC on Nan and duplicates (needed before unit conversion and gapcreation)
        # Remove nan names
        self._remove_nan_names()

        # Convert to numeric --> "invalid check' will be triggered if not possible
        self._to_num_and_invalid_check()
        self._append_to_applied_qc(
            obstypename=self._get_present_obstypes(), checkname="invalid_input"
        )

        # Remove duplicates (needed in order to convert the units and find gaps)
        df, outliersdf = duplicate_timestamp_check(
            df=self.df,
            checks_settings=self.settings.qc["qc_check_settings"],
        )
        self._set_df(df=df)
        self._update_outliersdf(outliersdf)
        self._append_to_applied_qc(
            obstypename=self._get_present_obstypes(), checkname="duplicated_timestamp"
        )

        # self._covert_timestamps_to_utc()

        # Check observation types and convert units if needed.
        self._setup_of_obstypes_and_units()

        # find the start, end timestamps and frequency for each station + write it to the metadf
        self._get_timestamps_info(
            freq_estimation_method=freq_estimation_method,
            freq_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
        )

        # Convert the records to clean equidistanced records for both the df and outliersdf
        self._construct_equi_spaced_records(
            timestamp_mapping_tolerance=timestamp_tolerance
        )

        # # Find gaps on Import resolution
        gaps = find_gaps(
            df=self.df,
            metadf=self.metadf,
            outliersdf=self.outliersdf,
            obstypes=self.obstypes,
        )
        self._set_gaps(gaps)

    def _setup_of_obstypes_and_units(self):
        """Function to set up all attributes related to observation types and
        convert to standard units."""
        # Check if all present observation types are known.
        present_obstypes = self._get_present_obstypes()

        # check if all present obstypes (in the df), are linked to a knonw Obstype
        self._are_all_present_obstypes_knonw()

        # Found that it is approx 70 times faster convert obstype per obstype,
        # add them to a list and concat them, than using the .loc method to
        # assign the converted values directly to the df attribute

        subdf_list = []
        for present_obs in present_obstypes:
            # Convert the units to the toolkit standards (if unit is known)
            input_unit = self.template._get_input_unit_of_tlk_obstype(present_obs)

            # locate the specific obstype records
            obstype_values = xs_save(self.df, present_obs, "obstype", drop_level=False)[
                "value"
            ]
            # Convert to standard unit and replace them in the df attribute
            subdf_list.append(
                self.obstypes[present_obs].convert_to_standard_units(
                    input_data=obstype_values, input_unit=input_unit
                )
            )
            # Update the description of the obstype
            self.obstypes[present_obs].set_description(
                desc=self.template._get_description_of_tlk_obstype(present_obs)
            )

            # Update the original name of the obstype (used for titles in plots)
            self.obstypes[present_obs].set_original_name(
                columnname=self.template._get_original_obstype_columnname(present_obs)
            )

            # Update the original unit of the obstype (not an application yet)
            self.obstypes[present_obs].set_original_unit(input_unit)

        df = pd.concat(subdf_list).to_frame().sort_index()
        self._set_df(df)

    def _to_num_and_invalid_check(self):
        # 8. map to numeric dtypes
        # When converting to numeric, this overrules the invalid check.
        checkname = "invalid_input"
        df = self.df

        # 1 subset to the records with Nan values --> do not check these, just add them back in the end
        nandf = df[~df["value"].notnull()]
        # 2 Get the other subet with values not nan (can be numerics and strings) --> filter out the strings
        to_checkdf = df[df["value"].notnull()]

        # 3 Convert to numeric
        to_checkdf["value"] = pd.to_numeric(to_checkdf["value"], errors="coerce")

        # 4 All the Nan's in the to_checkdf are outliers triggerd as 'invalid'
        invalid_records = to_checkdf[~to_checkdf["value"].notnull()]
        # add the label of "invalid check' to it
        invalid_records["label"] = label_def[checkname]["label"]
        # special case: duplicates in the invalid records
        invalid_records = invalid_records[~invalid_records.index.duplicated()]

        # 5. Combine the df's back to one
        totaldf = pd.concat([nandf, to_checkdf]).sort_index()

        # Set attributes
        # Note that at this point, the duplicated check is not performed yet.
        self._set_df(totaldf, apply_dup_checks=False)
        self._update_outliersdf(invalid_records)
        return

    def _construct_df(self, dataframe):
        """Fill the df attribute

        The dataframe is wide with data and metadata combined. This method will
        subset and format it to a long structured df to be set as the df attribute.
        """

        # subset to name, datetime and all obstypes
        df_cols = self.template._get_all_mapped_data_cols_in_tlk_space()
        # name and datetime are already index of dataframe, so drop them from df_cols
        df_cols = list(set(df_cols) - set(["name", "datetime"]))
        # subset the column
        df = dataframe.loc[:, df_cols]

        # convert the wide df to a long format
        triple_df = df.stack(future_stack=True)
        # rename the last level of the index to obstype
        triple_df.index.rename(names="obstype", level=-1, inplace=True)
        # rename the series to value
        triple_df.rename("value", inplace=True)
        # fix index order
        triple_df = triple_df.reorder_levels(
            ["name", "obstype", "datetime"]
        ).sort_index()
        # sort by index (by name --> datetime --> obstype)
        triple_df.sort_index(inplace=True)
        # convert to frame
        # TODO is this needed?
        triple_df = triple_df.to_frame()

        # set the attribute
        self._set_df(
            df=triple_df, apply_dup_checks=False
        )  # duplicate have yet to be removed

    def _remove_nan_names(self):
        """if the name is Nan, remove these records from df, and metadf (before)
        they end up in the gaps and missing observations."""

        if np.nan in self.df.index.get_level_values("name"):
            logger.warning(
                f'Following observations are not linked to a station name and will be removed: {xs_save(self.df, np.nan, "name")}'
            )
            self.df = self.df[~self.df.index.get_level_values("name").isna()]

        if np.nan in self.metadf.index:
            logger.warning(
                f"Following station will be removed from the Dataset {self.metadf[self.metadf.index.isna()]}"
            )
            self.metadf = self.metadf[~self.metadf.index.isna()]

    def _construct_metadf(self, dataframe, use_metadata):
        """Fill the metadf attribute.

        The dataframe is wide with data and metadata combined. This method will
        subset and format the data and set the metadf attribute.

        use_metadata is a bool, if True map the metadata normally. If False,
        the user has not specialized a metadata file --³ create a minimal
        metadf. (Because a template can hold mapping of metadata, but
        if the user does not specify the metadata file --> the template should
        not be used)

        """
        if use_metadata:
            meta_cols = list(self.template._get_metadata_column_map().values())

            metadf = (
                dataframe.reset_index()
                .loc[:, meta_cols]
                .drop_duplicates()
                .set_index("name")
            )
        else:
            # Construct a minimal metadf, compatible with df
            metadf = pd.DataFrame(
                index=dataframe.index.get_level_values("name").unique()
            )

        # Construct columns that are required
        if "lat" not in metadf.columns:
            metadf["lat"] = np.nan
        if "lon" not in metadf.columns:
            metadf["lon"] = np.nan

        # Convert to geopandas dataframe
        metadf = metadf_to_gdf(metadf)

        # set the attribute
        self._set_metadf(metadf=metadf)

    def _setup_of_obstypes_and_units(self):
        """Function to set up all attributes related to observation types and
        convert to standard units."""
        # Check if all present observation types are known.
        present_obstypes = self._get_present_obstypes()

        # check if all present obstypes (in the df), are linked to a knonw Obstype
        self._are_all_present_obstypes_knonw()

        # Found that it is approx 70 times faster convert obstype per obstype,
        # add them to a list and concat them, than using the .loc method to
        # assign the converted values directly to the df attribute

        subdf_list = []
        for present_obs in present_obstypes:
            # Convert the units to the toolkit standards (if unit is known)
            input_unit = self.template._get_input_unit_of_tlk_obstype(present_obs)

            # locate the specific obstype records
            obstype_values = xs_save(self.df, present_obs, "obstype", drop_level=False)[
                "value"
            ]
            # Convert to standard unit and replace them in the df attribute
            subdf_list.append(
                self.obstypes[present_obs].convert_to_standard_units(
                    input_data=obstype_values, input_unit=input_unit
                )
            )
            # Update the description of the obstype
            self.obstypes[present_obs].set_description(
                desc=self.template._get_description_of_tlk_obstype(present_obs)
            )

            # Update the original name of the obstype (used for titles in plots)
            self.obstypes[present_obs].set_original_name(
                columnname=self.template._get_original_obstype_columnname(present_obs)
            )

            # Update the original unit of the obstype (not an application yet)
            self.obstypes[present_obs].set_original_unit(input_unit)

        df = pd.concat(subdf_list).to_frame().sort_index()
        self._set_df(df)

    def _construct_equi_spaced_records(
        self, timestamp_mapping_tolerance="4min", direction="nearest"
    ):
        """
        Convert the records to regular records.

        This is done by creating a target list of timestamps (per station, per obstype)
        , and map the records (in self.df) to the target records. The target
        records are constructed by 'dt_start', 'dt_end' and 'dataset_resolution'
        columns present in the metadf.

        The mapping to the target is done by a nearest-merge and respecting a
        time_mapping_tolerance.

        The references to the outliers are mapped to the target as well. This is
        to ensure that the outliersdf and df attributes are compatible (= each
        Nan in the .df attribute must be present in the outliersdf, or covered
        by a Gap).

        The self.df and self.outliersdf are updated.


        Parameters
        ----------
        timestamp_mapping_tolerance : Timedelta or str
            The tolerance string or object representing the maximum translation
            (in time) to map a timestamp to a target timestamp.
            Ex: '5min' is 5 minutes.
        direction : 'backward', 'forward', or 'nearest' (default)
            Whether to search for prior, subsequent, or closest matches for
            mapping to ideal timestamps.


        Returns
        -------
        None

        Warning
        -----------
        This method will corrupt the outliers and gaps, thus they are initialized.
        Typically, gaps are located after this method.

        Note
        -------
        It can happen that the same original timestamp is mapped to multiple
        target timestamps if the nearest and tolerance conditions are met! This is
        not perse problematic, but this could affect the performance of repetitions QC checks.
        By setting the tolerance substantially smaller than the frequency, this
        phenomenon can be avoided.

        """

        # Construct a target index, by making dtranges as defined by a start, end
        # freq and tz as defined in the metadf
        trg_df_list = []
        for groupidx, groupd in self.df.reset_index().groupby(["name", "obstype"]):
            staname = groupidx[0]
            target_dtrange = pd.date_range(
                start=self.metadf.loc[staname, "dt_start"],
                end=self.metadf.loc[staname, "dt_end"],
                freq=self.metadf.loc[staname, "dataset_resolution"],
                tz=self._get_tz(),
            )

            multi_idx = pd.MultiIndex.from_arrays(
                arrays=[
                    [staname] * len(target_dtrange),
                    [groupidx[1]] * len(target_dtrange),
                    target_dtrange,
                ],
                names=["name", "obstype", "datetime"],
            )
            trg = pd.DataFrame(index=multi_idx)
            trg_df_list.append(trg)

        trg_df = pd.concat(trg_df_list).sort_index()

        # merge

        # allert! this merge asof can reduce the data amount, since some timestamps,
        # that are mapped to the same target timestamp are rejected (only the closest
        # one survives) and if there is no mapping candidate within the tolerance
        # the timestamp record is not used (and thus removed)
        trg_df["trg_datetime"] = trg_df.index.get_level_values("datetime")
        df = self.df

        df["label"] = "ok"

        # remove gaps from the df! because:
        # there records are Nan, and so a lot of Nan records are mapped to
        # the ideal timestamps

        # AND, after this method is applied, gaps are always recalculated

        df = df.dropna(subset=["value"])

        # Add a label column and concat the outliers to the df. In that way,
        # we can reconstruct the outliersdf again within the same clean resolutionspace
        if not self.outliersdf.empty:
            df = pd.concat([df, self.outliersdf])

        df = df[
            ~df.index.duplicated(keep="last")
        ]  # remove duplicated records, to keep the outliers
        df = df.sort_index()

        df["obs_datetime"] = df.index.get_level_values("datetime")

        dtmapping = pd.merge_asof(
            left=trg_df.reset_index().sort_values("datetime"),
            right=df.reset_index().sort_values("datetime"),
            by=["name", "obstype"],
            # suffixes=('_obs', '_trg'),
            left_on="datetime",
            right_on="datetime",
            tolerance=pd.Timedelta(timestamp_mapping_tolerance),
            direction=direction,
        )
        # Note: merge_asof is a left-merge under the hood. Because left == target,
        # there are no duplicates in the dtmapping by ['name', 'obstype', 'datetime']
        # and it could safely be set as index

        # Note2: It can happen that the same original timestamp is mapped to multiple
        # target timestamps, if the nearest and tolerance conditions are met! This is
        # not perse problematic, but this could affect the performance of repetitions QC checks.
        # By setting the tolerance substantially smaller than the frequency, this
        # phenomena can be avoided.

        # set index and sort
        dtmapping = dtmapping.set_index(["name", "obstype", "datetime"]).sort_index()

        # split the records in 'ok', outliers and Gaps the ouliers again

        # ------ find gaps -------
        # All timestamps, that could not have been mapped from the original records/outliers
        # are missing --> these are the gap records

        # Note: More pracktical to have a seperate find_gaps() function
        # gaprecords = dtmapping.loc[dtmapping['label'].isna()]

        # ------- find outliers -------
        # collect all the outliers that are sucsesfully mapped to the ideal records
        outliersdf = dtmapping.loc[
            dtmapping["label"].isin(self.outliersdf["label"].unique())
        ]
        outliersdf = outliersdf[["value", "label"]]

        # ------- find records --------
        # because the df, contains all the records (also the gap and outliers),
        # no subsetting is needed. BUT, make sure that the values of outlierrecords
        # are set to Nan!
        recordsdf = dtmapping
        recordsdf.loc[recordsdf["label"] != "ok", "value"] = np.nan
        recordsdf = recordsdf[["value"]]

        # Set attributes
        self._set_df(recordsdf)
        self._set_outliersdf(
            outliersdf
        )  # do not use the append_to_outliers, but overwrite them.
        return


# =============================================================================
# Errors
# =============================================================================


class MetobsDatasetBaseError(Exception):
    """Exception raised for errors in the datasetbase."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
