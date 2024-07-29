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


from metobs_toolkit.df_helpers import (
    # multiindexdf_datetime_subsetting,
    # fmt_datetime_argument,
    # init_multiindex,
    # init_multiindexdf,
    # init_triple_multiindexdf,
    empty_outliers_df,
    # metadf_to_gdf,
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


logger = logging.getLogger(__name__)


# %%
class _DatasetBase(object):
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

    # =============================================================================
    # Specials
    # =============================================================================
    def __str__(self):
        """Represent as text."""
        if self.df.empty:
            if self._istype == "Dataset":
                return "Empty instance of a Dataset."
            elif self._istype == "Station":
                return "Empty instance of a Station."
            else:
                return "Empty instance of a Analysis."

        add_info = ""
        n_stations = self.df.index.get_level_values("name").unique().shape[0]
        n_obs_tot = self.df["value"].count()
        n_outl = self.outliersdf.shape[0]
        startdt = self.df.index.get_level_values("datetime").min()
        enddt = self.df.index.get_level_values("datetime").max()

        if (not self.metadf["lat"].isnull().all()) & (
            not self.metadf["lon"].isnull().all()
        ):
            add_info += "    *Coordinates are available for all stations."

        return (
            f"{self._istype} instance containing: \n \
    *{n_stations} stations \n \
    *{self.df.index.get_level_values('obstype').unique().to_list()} observation types present\n \
    *{n_obs_tot} observation records (not Nan's) \n \
    *{n_outl} records labeled as outliers \n \
    *{len(self.gaps)} gaps \n \
    *records range: {startdt} --> {enddt} (total duration:  {enddt - startdt}) \n \
    *time zone of the records: {str(self._get_tz())} \n "
            + add_info
        )

    def __repr__(self):
        """Info representation."""
        class_name = type(self).__name__
        return f"Instance of {class_name} at {hex(id(self))}"

    def __add__(self, other, gapsize=None):
        """Addition of two Datasets."""
        # important !!!!!

        # the toolkit makes a new dataframe, and assumes the df from self and other
        # to be the input data.
        # This means that missing obs, gaps, invalid and duplicated records are
        # being looked for in the concatenation of both dataset, using their current
        # resolution !

        new = Dataset()
        self_obstypes = self.df.columns.to_list().copy()
        #  ---- df ----

        # check if observation of self are also in other
        assert all([(obs in other.df.columns) for obs in self_obstypes])
        # subset obstype of other to self
        other.df = other.df[self.df.columns.to_list()]

        # remove duplicate rows
        common_indexes = self.df.index.intersection(other.df.index)
        other.df = other.df.drop(common_indexes)

        # set new df
        new.df = concat_save([self.df, other.df])
        new.df = new.df.sort_index()

        #  ----- outliers df ---------

        other_outliers = other.outliersdf.reset_index()
        other_outliers = other_outliers[other_outliers["obstype"].isin(self_obstypes)]
        other_outliers = other_outliers.set_index(["name", "datetime", "obstype"])
        new.outliersdf = concat_save([self.outliersdf, other_outliers])
        new.outliersdf = new.outliersdf.sort_index()

        #  ------- Gaps -------------
        # Gaps have to be recaluculated using a frequency assumtion from the
        # combination of self.df and other.df, thus NOT the native frequency if
        # their is a coarsening allied on either of them.
        new.gaps = []

        # ---------- missing ---------
        # Missing observations have to be recaluculated using a frequency assumtion from the
        # combination of self.df and other.df, thus NOT the native frequency if
        # their is a coarsening allied on either of them.
        new.missing_obs = None

        # ---------- metadf -----------
        # Use the metadf from self and add new rows if they are present in other
        new.metadf = concat_save([self.metadf, other.metadf])
        new.metadf = new.metadf.drop_duplicates(keep="first")
        new.metadf = new.metadf.sort_index()

        # ------- specific attributes ----------

        # Template (units and descritpions) are taken from self
        new.template = self.template

        # Inherit Settings from self
        new.settings = copy.deepcopy(self.settings)

        # Applied qc:
        # TODO:  is this oke to do?
        new._applied_qc = pd.DataFrame(columns=["obstype", "checkname"])

        # set init_dataframe to empty
        # NOTE: this is not necesarry but users will use this method when they
        # have a datafile that is to big. So storing and overloading a copy of
        # the very big datafile is invalid for these cases.
        new.input_df = pd.DataFrame()

        # ----- Apply IO QC ---------
        # Apply only checks that are relevant on records in between self and other
        # OR
        # that are dependand on the frequency (since the freq of the .df is used,
        # which is not the naitive frequency if coarsening is applied on either. )

        # missing and gap check
        if gapsize is None:
            gapsize = new.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"]

        # note gapsize is now defined on the frequency of self
        new.missing_obs, new.gaps = missing_timestamp_and_gap_check(
            df=new.df,
            gapsize_n=self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"],
        )

        # duplicate check
        new.df, dup_outl_df = duplicate_timestamp_check(
            df=new.df,
            checks_info=new.settings.qc["qc_checks_info"],
            checks_settings=new.settings.qc["qc_check_settings"],
        )

        if not dup_outl_df.empty:
            new.update_outliersdf(add_to_outliersdf=dup_outl_df)

        # update the order and which qc is applied on which obstype
        checked_obstypes = list(self.obstypes.keys())

        checknames = ["duplicated_timestamp"]  # KEEP order

        new._applied_qc = concat_save(
            [
                new._applied_qc,
                conv_applied_qc_to_df(
                    obstypes=checked_obstypes, ordered_checknames=checknames
                ),
            ],
            ignore_index=True,
        )

        return new

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
        """Formats a datetime given by a useser in an argument."""

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
    #     attribute setters
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
        # TODO: run simple checks
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

    def _append_to_applied_qc(self, obstypename, checkname):
        """Add record to _applied_qc.
        The applied qc is mainly used to save the order of applied checks,
        to validate the effectivenes of one check.

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

        if checkname not in self.settings.label_def.keys():
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
        return self.df.index.get_level_values("datetime").tz

    def _get_present_obstypes(self):
        """Get all present obstypenames in the df.

        (This does not guarantee that these obstypesnames are knonw obstypes!)

        Returns
        -------
        list
            All present obstypenames.

        """
        return self.df.index.get_level_values("obstype").unique().to_list()

    def _get_timestamps_info(
        self,
        freq_estimation_method="median",
        freq_simplify_tolerance="2T",
        origin_simplify_tolerance="5T",
    ):
        """Get details of the time resolution for each station.

        It is assumed that ideally, records are perfect periodically. In
        pracktice this is not always the case. We can specify an ideal set
        of timestamps by a start, end and frequency.

        This methods estimates, for each station:

            - Frequency: the frequency of a station is the highest frequency
            detected for all it's observationtypes. The frequency estimate
            of an obsertype is done by a 'mean' or 'highest' approach on
            consecutive records. One can specify a tolerance to simplyfy
            the frequency estimates.

            - start timestamp: The first timestemp registerd by a station,
            over all its observationtypes. One can specify a tolerance to
            simplify the origin.

            - end timestamp: The latest timestemp registerd by a station,
            over all its observationtypes. One can specify a tolerance to
            simplify it.

        This information is added to the metadf attribute.

        Note
        -----
        The assumtion is made that all obstypes per station have the same start
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
    # Data attribute updaters
    # =============================================================================

    def _remove_nan_names(self):
        """if the name is Nan, remove these records from df, and metadf (before)
        # they end up in the gaps and missing obs"""

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

    def construct_equi_spaced_records(
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
        target timestamps, if the nearest and tolerance conditions are met! This is
        not perse problematic, but this could affect the performance of repetitions QC checks.
        By setting the tolerance substantially smaller than the frequency, this
        phenomena can be avoided.

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
                tz=self.metadf.loc[staname, "dt_start"].tz,
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

        # Add a label column and concat the outliers to the df. In that way,
        # we can reconstruct the outliersdf again within the same clean resolutionspace
        df["label"] = "ok"
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
    # Base fuctions
    # =============================================================================

    def get_full_status_df(
        self,
        return_as_wide=True,
    ):
        """
        TODO: docstring

        Make one dataframe with all observations and their labels.

        Combine all observations, outliers, missing observations and gaps into
        one Dataframe. All observation types are combined an a label is added
        in a serperate column.

        When gaps and missing records are updated from outliers one has to choice
        to represent these records as outliers or gaps. There can not be duplicates
        in the return dataframe.

        By default the observation values of the outliers are saved, one can
        choice to use these values or NaN's.
        following checks!

        Parameters
        ----------
        repr_outl_as_nan : bool, optional
            If True, Nan's are use for the values of the outliers. The
            default is False.
        overwrite_outliers_by_gaps_and_missing : Bool, optional
            If True, records that are labeld as gap/missing and outlier are
            labeled as gaps/missing. This has only effect when the gaps/missing
            observations are updated from the outliers. The default is True.

         Returns
         ---------
         combdf : pandas.DataFrame()
            A dataframe containing a continious time resolution of records, where each
            record is labeld.

        Examples
        --------
        .. code-block:: python

            >>> import metobs_toolkit
            >>>
            >>> # Import data into a Dataset
            >>> dataset = metobs_toolkit.Dataset()
            >>> dataset.update_settings(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.import_data_from_file()
            >>> dataset.coarsen_time_resolution(freq='1h')
            >>>
            >>> # Apply quality control on the temperature observations
            >>> dataset.apply_quality_control(obstype='temp') #Using the default QC settings
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                 *10080 observation records
                 *1676 records labeled as outliers
                 *0 gaps
                 *3 missing observations
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
            >>>
            >>> # Combine all records to one dataframe in Observation-resolution
            >>> overview_df = dataset.get_full_status_df()
            >>> overview_df.head(12)
                                                                    value  ... toolkit_representation
            name      datetime                  obstype                    ...
            vlinder01 2022-09-01 00:00:00+00:00 humidity        65.000000  ...            observation
                                                temp            18.800000  ...            observation
                                                wind_direction  65.000000  ...            observation
                                                wind_speed       1.555556  ...            observation
                      2022-09-01 01:00:00+00:00 humidity        65.000000  ...            observation
                                                temp            18.400000  ...            observation
                                                wind_direction  55.000000  ...            observation
                                                wind_speed       1.416667  ...            observation
                      2022-09-01 02:00:00+00:00 humidity        68.000000  ...            observation
                                                temp            17.100000  ...            observation
                                                wind_direction  45.000000  ...            observation
                                                wind_speed       1.583333  ...            observation
            <BLANKLINE>
            [12 rows x 3 columns]

        """
        # TODO: label values from settings not hardcoding

        # TODO: use the repr_outl_as_nan argumenten here
        # =============================================================================
        # Stack observations and outliers
        # =============================================================================
        df = self.df
        # note: df is a pointer, and adding these colmns will add them
        # also in the self.df
        df["label"] = "ok"
        df["toolkit_representation"] = "observation"

        # =============================================================================
        # Stack outliers
        # =============================================================================

        outliersdf = self.outliersdf
        outliersdf["toolkit_representation"] = "outlier"

        combdf = pd.concat([df, outliersdf])  # combine the two

        # Since outliers are present records in the df (as NaN's) we introduce
        # duplicats in the index of combdf. We drop the duplicates and keep,
        # the records comming from outliersdf (=last)

        combdf = combdf[~combdf.index.duplicated(keep="last")]

        # =============================================================================
        # Stack gaps
        # =============================================================================

        gapsdf = (
            self._get_gaps_df_for_stacking()
        )  # get a gapdf in the long (similar as outliersdf) structure
        # map labels to known labels (must have a color def in the settings)
        if not gapsdf.empty:
            gapsdf["label"] = gapsdf["fill_method"].replace(
                {"not filled": self.settings.label_def["regular_gap"]["label"]}
            )

        gapsdf = gapsdf[["value", "label"]]

        gapsdf["toolkit_representation"] = "gap"

        combdf = pd.concat([combdf, gapsdf])  # combine

        # Since gaps are present records in the df (as NaN's, because of the
        # ideal freq structure in the df) we introduce
        # duplicats in the index of combdf. We drop the duplicates and keep,
        # the records comming from outliersdf (=last)

        combdf = combdf[~combdf.index.duplicated(keep="last")]
        # =============================================================================
        # Formatting the combineddf
        # =============================================================================

        assert (
            not combdf.index.duplicated().any()
        ), "Duplicates found in the combdf --> report bug."

        # for some reason the dtype of the datetime index-level is 'obstype' and
        # thus not a datetimeindex. This must be fixed
        combdf = combdf.reset_index()
        combdf["datetime"] = pd.to_datetime(combdf["datetime"])
        combdf = combdf.set_index(["name", "obstype", "datetime"]).sort_index()

        if return_as_wide:
            combdf = combdf.unstack(level="obstype").reorder_levels(
                order=[1, 0], axis=1
            )

        # pointer issue
        self.df = self.df[["value"]]
        self.outliersdf = self.outliersdf[["value", "label"]]
        return combdf


# =============================================================================
# Errors
# =============================================================================


class MetobsDatasetBaseError(Exception):
    """Exception raised for errors in the datasetbase."""

    pass
