#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This module contains the Gap class and all its methods.

A Gap contains all information and methods of a data-gap.
"""
import numpy as np
import pandas as pd

import logging
import datetime

import metobs_toolkit.gap_filling as gap_filling

from metobs_toolkit.df_helpers import (
    init_multiindexdf,
    xs_save,
    concat_save,
)

from metobs_toolkit.obstypes import Obstype as Obstype_class


from metobs_toolkit.settings_files.default_formats_settings import (
    label_def,
)  # to assign labels

logger = logging.getLogger(__name__)


# =============================================================================
# Gap class

# a gap is a sequence of repeting missing obs
# =============================================================================


class Gap:
    """Gap class holds all gap information and methods for gaps."""

    def __init__(self, name, startdt, enddt, obstype, records_freq):
        """
        Initiate Gap object.


        Parameters
        ----------
        name : String
            Station name where the gap occurs.
        startdt : datetime.datetime
            Start datetime of the gap (included).
        enddt : datetime.datetime
            End datetime of the gap (included)
        obstype : metobs_toolkit.Obstype
           The corresponding observation type of the gap.
        records_freq : pd.Timedelta()
           The assumed frequency of records for this observationtype.

        Returns
        -------
        None

        Note
        ------
        Gaps are automatically located in your Dataset. In common practices,
        a user does not have to initiate a Gap.

        """
        # Set key attributes
        self.name = str(name)

        assert isinstance(
            startdt, datetime.datetime
        ), f"{startdt} not a datetime instance."
        assert isinstance(enddt, datetime.datetime), f"{enddt} not a datetime instance."
        assert enddt >= startdt, f"Enddt: {enddt} is not after startdt: {startdt}."

        self.startdt = startdt
        self.enddt = enddt
        self.duration = enddt - startdt

        assert isinstance(
            obstype, Obstype_class
        ), f"{obstype} is not a metobs_toolkit.Obstype."
        self.obstype = obstype

        assert isinstance(
            records_freq, pd.Timedelta
        ), f"{records_freq} is not a pd.Timedelta."
        self.records_freq = records_freq

        # gapdf holds all the records in the gap (depending on the freq of the dataset!)
        self.gapdf = (
            init_multiindexdf()
        )  # columns: obstype.name, obstype.name"_fill', "fill_method", "msg"
        # anchordf holds all the records outside the gap, used by fill methods.
        self.anchordf = (
            init_multiindexdf()
        )  # columns: obstype.name, obstype.name"_fill', "fill_method", "msg"

        # Construct derived attributes and initiate them
        self._setup_of_gapdf()

    # =============================================================================
    # Specials
    # =============================================================================

    def __str__(self):
        return f"{self.obstype.name}-gap of {self.name} for {self.startdt} --> {self.enddt}, duration: {self.duration} or {self.gapdf.shape[0]} records."

    def __repr__(self):
        return self.__str__()

    # =============================================================================
    # Setup (methods automatically applied after creation)
    # =============================================================================
    def _setup_of_gapdf(self):
        """Setup of the gapdf attribute from the core attributes."""
        # Compute missing records
        missing_records = pd.date_range(
            start=self.startdt,
            end=self.enddt,
            freq=self.records_freq,
            tz=self.startdt.tz,
        )

        # ----  Construct gapdf ------
        gapdf = init_multiindexdf()
        # Create index
        gapdf.index = pd.MultiIndex.from_arrays(
            arrays=[
                [self.name] * len(missing_records),
                missing_records,
            ],
            names=["name", "datetime"],
        )

        # Set columns with default values
        gapdf[f"{self.obstype.name}"] = np.nan  # fill with nan values
        gapdf[f"{self.obstype.name}_fill"] = np.nan  # fill with nan values
        gapdf["fill_method"] = "not filled"  # on gap creation, will be overwritten
        gapdf["msg"] = None  # Nonedataset

        self._set_gapdf(gapdf)

    # =============================================================================
    # Getters
    # =============================================================================

    def _get_gapfill_status(self):
        """Returns a label on the status of gap filling."""
        # Warning! do not change these labels (or replace all hardcode occurences)
        # if not self.does_gap_holds_missing_records():
        #     return "Gap does not exist in observation space"
        if self.gapdf[f"{self.obstype.name}_fill"].isnull().all():
            return "Unfilled gap"
        if self.gapdf[f"{self.obstype.name}_fill"].isnull().any():
            return "Partially filled gap"
        else:
            return "Filled gap"

    def get_info(self):
        """Print out detailed info about the Gap.

        Returns
        -------
        None.

        """

        print("---- Gap info -----")
        print(
            "(Note: gaps start and end are defined on the frequency estimation of the native dataset.)"
        )
        print(f"  * Gap for station: {self.name}")
        print(f"  * Start gap: {self.startdt}")
        print(f"  * End gap: {self.enddt}")
        print(f"  * Duration gap: {self.duration}")
        print(f"  * For {self.obstype.name}")
        print(f"  * Gapfill status >>>> {self._get_gapfill_status()}")
        # if not self.does_gap_holds_missing_records():
        #     print(
        #         "The gap is unreal (=no missing records at the current frequency resolution)."
        #     )

        if not self.anchordf.empty:
            print("---- Anchor Data Frame ----- ")
            print(self.anchordf)

        print("---- Gap Data Frame -----")
        print(self.gapdf)

    # =============================================================================
    # Setters
    # =============================================================================

    def _set_gapdf(self, gapdf):
        """Setter for the gapdf attribute."""
        if list(gapdf.index.names) != ["name", "datetime"]:
            raise MetobsGapError(f"The gapdf does not have a correct index: {gapdf}")

        try:
            if gapdf.index.get_level_values("datetime").tz is None:
                raise MetobsGapError(
                    f"The gapdf-datetimeindex is tz-naive. {gapdf.index}"
                )
            pass
        except:
            raise MetobsGapError(
                f"The gapdf-datetimeindex is not a datetime dtype. {gapdf.index}"
            )
        if list(gapdf.columns) != [
            f"{self.obstype.name}",
            f"{self.obstype.name}_fill",
            "fill_method",
            "msg",
        ]:
            raise MetobsGapError(f"The gapdf does not have a correct columns: {gapdf}")

        self.gapdf = gapdf

    def _set_anchordf(self, anchordf):
        """Setter for the anchordf."""
        if list(anchordf.index.names) != ["name", "datetime"]:
            raise MetobsGapError(
                f"The anchordf does not have a correct index: {anchordf}"
            )
        try:
            if anchordf.index.get_level_values("datetime").tz is None:
                raise MetobsGapError(
                    f"The gapdf-datetimeindex is tz-naive. {anchordf.index}"
                )
            pass
        except:
            raise MetobsGapError(
                f"The gapdf-datetimeindex is not a datetime dtype. {anchordf.index}"
            )

        if list(anchordf.columns) != [self.obstype.name, "fill_method", "msg"]:
            raise MetobsGapError(
                f"The anchordf does not have a correct columns: {anchordf}"
            )
        self.anchordf = anchordf

    # =============================================================================
    # Checks
    # =============================================================================

    def _can_be_filled(self, overwrite_fill):
        """Returns True if no fill values are present."""

        if self.gapdf[f"{self.obstype.name}_fill"].isnull().all():
            return True
        elif overwrite_fill:
            logger.warning(f"The filling of {self}, will be overwritten!")
            self._setup_of_gapdf()  # reset gapfill
            return True

        else:
            return False

    # =============================================================================
    #         Gapfill
    # =============================================================================

    def interpolate(
        self,
        Dataset,
        method="time",
        max_consec_fill=10,
        n_leading_anchors=1,
        n_trailing_anchors=1,
        max_lead_to_gap_distance=None,
        max_trail_to_gap_distance=None,
        method_kwargs={},
    ):
        """Fill a Gap using interpolation.

        The Gap (data) attributes will be updated.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset or metobs_toolkit.Station
            The dataset that holds the records of the same station as the Gap
            is from.
        method : str, optional
            Interpolation technique to use. See pandas.DataFrame.interpolate
            'method' argument for possible values. Make sure that
            `n_leading_anchors`, `n_trailing_anchors` and `method_kwargs` are
            set accordingly to the method. The default is "time".
        max_consec_fill : int, optional
            The maximum number of consecutive missing records to fill. The default is 10.
        n_leading_anchor : int, optional
            The number of leading anchors to use for the interpolation.
            Higher-order polynomial interpolation techniques require multiple leading
            anchors. The default is 1.
        n_trailing_anchor : int, optional
            The number of trailing anchors to use for the interpolation. Higher
            order polynomial interpolation techniques require multiple trailing
            anchors. The default is 1.
        max_lead_to_gap_distance : str or pandas.Timedelta, optional
            The maximum time difference between the start of the gap and a
            suitable lead (= the good record to start the interpolation from).
            If None, the first occurring good records before the start of the gap
            is used. The default is None.
        max_trail_to_gap_distance : str or pandas.Timedelta, optional
            The maximum time difference between the end of the gap and a
            suitable trail (= the good record to end the interpolation on).
            If None, the first occurring good records after the gap
            is used. The default is None.
        method_kwargs: dict, optional
            A dictionary of kwargs passed to pandas.Dataframe.interpolate(). In
            practice, extra arguments for specific interpolation methods are
            put in method_kwargs. The default is {}.

        Returns
        -------
        None

        See Also
        --------
        Gap: The Gap class.
        Dataset.interpolate_gaps: Interpolate all gaps in a Dataset.
        Gap.raw_model_gapfill: Raw modeldata gapfill method.
        Gap.debias_model_gapfill: Debiased modeldata gapfill method.
        Gap.diurnal_debias_model_gapfill: Diurnal debiased modeldata gapfill method.
        Gap.weighted_diurnal_debias_model_gapfill: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description :

        1. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        2. Find a leading (the last observations before the gap) record and a trailing record (the last observation after the gap).
        3. Check if the leading and trailing records fulfill the criteria of maximum time difference.
        4. Using the leading and trailing record an interpolation is applied to fill the missing records. A maximum consecutive fill threshold is applied, if exceeded the fill values are Nan's.
        5. The gap is updated with the interpolated values

        Note
        -------
        The impact of `max_consec_fill` depends highly on the resolution
        of your records.

        Note
        ------
        If you want to use a higher-order method of interpolation, make sure to
        increase the `n_leading_anchors` and `n_trailing_anchors` accordingly.

        Examples
        ----------
        See ``Dataset.interpolate_gaps`` for examples.

        """
        gapdf = self.gapdf
        obsname = self.obstype.name

        # 1. Get leading and trailing
        # Create anchor periods
        anchordf, lead_msg, trail_msg = (
            gap_filling._create_anchor_df_for_leading_trailing_periods_by_size(
                Gap=self,
                Dataset=Dataset,
                n_lead_records=n_leading_anchors,
                n_trail_records=n_trailing_anchors,
                max_lead_duration=max_lead_to_gap_distance,
                max_trail_duration=max_trail_to_gap_distance,
            )
        )
        # Update attributes
        self._set_anchordf(anchordf)

        gapdf["fill_method"] = label_def["interpolated_gap"]["label"]

        # check if gapfill can proceed
        if (lead_msg != "ok") & (trail_msg != "ok"):
            logger.warning(
                f"Cannot fill {self}, because leading and trailing records are not valid."
            )
            print(
                f"Warning! Cannot fill {self}, because leading and trailing periods are not valid."
            )

            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["fill_method"] = label_def["failed_interpolation_gap"]["label"]
            gapdf["msg"] = f"{lead_msg} and {trail_msg}"

            self._set_gapdf(gapdf)
            return
        elif (lead_msg != "ok") & (trail_msg == "ok"):
            logger.warning(f"Cannot fill {self}, because leading record is not valid.")
            # print(f"Warning! Cannot fill {self}, because leading record is not valid.")

            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["fill_method"] = label_def["failed_interpolation_gap"]["label"]
            gapdf["msg"] = f"{lead_msg}"
            self._set_gapdf(gapdf)
            return
        elif (lead_msg == "ok") & (trail_msg != "ok"):
            logger.warning(f"Cannot fill {self}, because trailing record is not valid.")
            # print(f"Warning! Cannot fill {self}, because trailing record is not valid.")

            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["fill_method"] = label_def["failed_interpolation_gap"]["label"]
            gapdf["msg"] = f"{trail_msg}"
            self._set_gapdf(gapdf)
            return
        else:
            pass

        # 3. Combine the anchors with the observations
        tofilldf = pd.concat([gapdf[[obsname]], self.anchordf[[obsname]]])

        # 4. Replace the NaN's (GAPFILLING)
        # make shure only datetimes are in the index, and sorted by datetimes
        tofilldf = tofilldf.reset_index().set_index("datetime").sort_index()

        # Interpolate series
        tofilldf[obsname] = tofilldf[obsname].interpolate(
            method=method,
            limit=max_consec_fill,  # Maximum number of consecutive NaNs to fill. Must be greater than 0.
            limit_area="inside",
            **method_kwargs,
        )
        # revert to multiindex
        tofilldf = tofilldf.reset_index().set_index(["name", "datetime"]).sort_index()

        # 5. Get messages for all filled records
        # filter to gap records
        tofilldf = tofilldf.loc[gapdf.index]

        # messager
        def _msg_interp(row, obstypename):
            if np.isnan(row[obstypename]):
                return "Permitted_by_max_consecutive_fill"
            else:
                if bool(method_kwargs):
                    return f"ok (method: {method}, with {method_kwargs})"
                else:
                    return f"ok (method: {method})"

        if (lead_msg == "ok") & (
            trail_msg == "ok"
        ):  # labels are result of interpolation
            tofilldf["msg"] = tofilldf.apply(
                _msg_interp, axis="columns", args=([obsname])
            )
        elif (lead_msg == "ok") & (
            trail_msg != "ok"
        ):  # labels are result of invalid trailing
            tofilldf["msg"] = "Invalid trailing anchor"
        elif (lead_msg != "ok") & (
            trail_msg == "ok"
        ):  # labels are result of invalid leading
            tofilldf["msg"] = "Invalid leading anchor"
        elif (lead_msg != "ok") & (
            trail_msg != "ok"
        ):  # labels are result of invalid trailing
            tofilldf["msg"] = "Invalid leading and trailing anchor"
        else:
            raise MetobsGapError(
                "unforseen situation of lead and trail messages in interpolation method."
            )

        # 6. Update the gapdf
        gapdf[f"{obsname}_fill"] = tofilldf[obsname]  # update values
        gapdf["msg"] = tofilldf["msg"]  # give the same msg to all

        # 7. Add FAIL labels if needed
        gapdf.loc[gapdf[f"{obsname}_fill"].isna(), "fill_method"] = label_def[
            "failed_interpolation_gap"
        ]["label"]

        self._set_gapdf(gapdf)

        return

    # =============================================================================
    # Gapfill using modeldata
    # =============================================================================

    def raw_model_gapfill(self, Dataset, Model):
        """Fill the Gap with raw modeldata.

        The Gap (data) attributes will be updated.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset or metobs_toolkit.Station
            The dataset that holds the records of the same station as the Gap
            is from.
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.

        Returns
        -------
        None.

        See Also
        --------
        Gap: The Gap class.
        Dataset.fill_gaps_with_raw_modeldata: Equivalent for all gaps in a Dataset.
        GeeDynamicModelData: The Gee Model data (timeseries).
        Dataset.get_modeldata: Method for creating a modeldata from a dataset.
        Gap.interpolate: Interpolate gap.
        Gap.raw_model_gapfill: Raw modeldata gapfill method.
        Gap.debias_model_gapfill: Debiased modeldata gapfill method.
        Gap.diurnal_debias_model_gapfill: Diurnal debiased modeldata gapfill method.
        Gap.weighted_diurnal_debias_model_gapfill: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description of the raw modeldata gap fill:

        1. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        2. The modeldata is interpolated (in time) to the missing records.
        3. The gap is updated with the interpolated values from the modeldata.

        Examples
        ----------
        See ``Dataset.fill_gaps_with_raw_modeldata`` for examples.

        """
        obsname = self.obstype.name

        # 1. Get leading and trailing info
        # 2. Update the anchordf
        # 3. Combine the anchors with the observations
        label = label_def["raw_modeldata_fill"]["label"]
        fail_label = label_def["failed_raw_modeldata_fill"]["label"]

        # add the gap period
        filldf = Model._interpolate_modeldata(self.gapdf.index)
        filldf = filldf[[obsname]]  # get only relevant obstype
        filldf = filldf.rename(columns={obsname: "modelvalues"})

        # combine the modeldata and gapdata in one df
        filldf = pd.merge(
            left=self.gapdf, right=filldf, how="left", left_index=True, right_index=True
        )

        # 4. Fill the gap
        fill_values, msg = gap_filling._raw_modeldata_fill(filldf)

        # 5. create a gapdf
        gapdf = filldf
        gapdf = gapdf[[f"{obsname}", f"{obsname}_fill", "fill_method", "msg"]]

        gapdf[f"{obsname}_fill"] = fill_values
        gapdf["msg"] = msg

        gapdf["fill_method"] = label_def["raw_modeldata_fill"]["label"]
        gapdf.loc[gapdf[f"{obsname}_fill"].isna(), "fill_method"] = label_def[
            "failed_raw_modeldata_fill"
        ]["label"]

        self._set_gapdf(gapdf)

    def debias_model_gapfill(
        self,
        Dataset,
        Model,
        leading_period_duration="24h",
        min_leading_records_total=60,
        trailing_period_duration="24h",
        min_trailing_records_total=60,
    ):
        """Fill Gap with debiased modeldata.

        The gap will be updated with fill values using the Modeldata. The
        Modeldata is interpolated (in time) to the missing records,
        and corrected with a bias-correction. The bias is estimated by making use
        of a leading (before the gap) and trailing (after the gap) period.

        The Gap (data) attributes will be updated.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset or metobs_toolkit.Station
            The dataset that holds the records of the same station as the Gap
            is from.
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        leading_period_duration : str or pandas.Timedelta, optional
            The duration of the leading period. The default is "24h".
        min_leading_records_total : int, optional
            The minimum number of good records in the leading period. The default
            is 60.
        trailing_period_duration : str or pandas.Timedelta, optional
            The duration of the trailing period. The default is "24h".
        min_trailing_records_total : int, optional
            The minimum number of good records in the trailing period. The
            default is 60.

        Returns
        -------
        None.

        See Also
        --------
        Gap: The Gap class.
        Dataset.fill_gaps_with_debiased_modeldata: Equivalent for all gaps in a Dataset.
        GeeDynamicModelData: The Gee Model data (timeseries).
        Dataset.get_modeldata: Method for creating a modeldata from a dataset.
        Gap.interpolate: Interpolate gap.
        Gap.raw_model_gapfill: Raw modeldata gapfill method.
        Gap.diurnal_debias_model_gapfill: Diurnal debiased modeldata gapfill method.
        Gap.weighted_diurnal_debias_model_gapfill: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description:

        1. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        2. The good observations of the leading and trailing periods are selected and checked if they fulfill the conditions.
        3. The modeldata is interpolated (in time) to the missing records, the leading, and the trailing period.
        4. By combining the leading and trailing periods of both records and modeldata, a bias is calculated.
        5. The gap is updated with the interpolated modeldata, corrected by the calculated bias.

        Examples
        ----------
        See ``Dataset.fill_gaps_with_debiased_modeldata`` for examples.

        """

        # 1. Get leading and trailing info
        obsname = self.obstype.name

        # Create anchor periods
        anchordf = gap_filling._create_anchor_df_for_leading_trailing_periods(
            Gap=self,
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            trailing_period_duration=trailing_period_duration,
        )
        anchordf["msg"] = "_init_"
        # Check validity of anchors
        # 1. Check leading period total sample size and add msg
        if (
            anchordf[anchordf["fill_method"] == "leading period"].shape[0]
            < min_leading_records_total
        ):
            msg_lead = f"minimum records ({min_leading_records_total}) for leading period not met."
        else:
            msg_lead = "ok"

        anchordf.loc[anchordf["fill_method"] == "leading period", "msg"] = msg_lead

        # 1B. Check trailing period total sample size and add msg
        if (
            anchordf[anchordf["fill_method"] == "trailing period"].shape[0]
            < min_trailing_records_total
        ):
            msg_trail = f"minimum records ({min_trailing_records_total}) for trailing period not met."
        else:
            msg_trail = "ok"

        anchordf.loc[anchordf["fill_method"] == "trailing period", "msg"] = msg_lead

        # 2. Update the Anchor attribute
        self._set_anchordf(anchordf)

        # 3.Test if anchor period is valid, if not write default attribute
        if (msg_lead != "ok") & (msg_trail != "ok"):
            logmsg = f"Cannot fill {self}, because leading and trailing periods are not valid."
            gapmsg = f"{msg_lead} and {msg_trail}"
            err = True
        elif (msg_lead != "ok") & (msg_trail == "ok"):
            logmsg = f"Cannot fill {self}, because leading period is not valid."
            gapmsg = f"{msg_lead}"
            err = True
        elif (msg_lead == "ok") & (msg_trail != "ok"):
            logmsg = f"Cannot fill {self}, because trailing period is not valid."
            gapmsg = f"{msg_trail}"
            err = True
        else:
            err = False

        gapdf = self.gapdf
        if err:
            logger.warning(logmsg)
            # print("Warning! ", logmsg)
            gapdf[f"{obsname}_fill"] = np.nan
            gapdf["msg"] = gapmsg
            gapdf["fill_method"] = label_def["failed_debias_modeldata_fill"]["label"]
        else:

            # 4. combine learning and gap period
            filldf = gap_filling._combine_learning_and_gap_to_one_df(
                Modeldata=Model,
                anchordf=self.anchordf,
                gapdf=gapdf,
                obsname=self.obstype.name,
            )

            # 5. Fill the gap
            fill_values, msg = gap_filling._debias_modeldata_fill(filldf)

            # 6. Update gapdf attribute
            gapdf[f"{obsname}_fill"] = fill_values
            gapdf["msg"] = msg
            gapdf["fill_method"] = label_def["debias_modeldata_fill"]["label"]
            gapdf.loc[gapdf[f"{obsname}_fill"].isna(), "fill_method"] = label_def[
                "failed_debias_modeldata_fill"
            ]["label"]

        # 7. Update gapdf attribute
        self._set_gapdf(gapdf)

    def diurnal_debias_model_gapfill(
        self,
        Dataset,
        Model,
        leading_period_duration="24h",
        min_debias_sample_size=6,
        trailing_period_duration="24h",
    ):
        """Fill the Gap with (diurnal debiased) model data.

        The gap will be updated with fill values using the
        Modeldata. The Modeldata is interpolated (in time) to the missing records,
        and corrected with a bias-correction. Multiple biasses are computed, one
        for each timestamp present in the missing records, by using a leading
        (before the gap) and trailing (after the gap) period. Each bias is
        computed at each timestamp, thus computing a diurnal-bias-cycle.

        The Gap (data) attributes will be updated.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset or metobs_toolkit.Station
            The dataset that holds the records of the same station as the Gap
            is from.
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        leading_period_duration : str or pandas.Timedelta, optional
            The duration of the leading period. The default is "24h".
        min_debias_sample_size : int, optional
            The minimum number of good records to
            calculate a diurnal bias. The default is 6.
        trailing_period_duration : str or pandas.Timedelta, optional
            The duration of the trailing period. The default is "24h".

        Returns
        -------
        None.

        See Also
        --------
        Gap: The Gap class.
        Dataset.fill_gaps_with_diurnal_debiased_modeldata: Equivalent for all gaps in a Dataset.
        GeeDynamicModelData: The Gee Model data (timeseries).
        Dataset.get_modeldata: Method for creating a modeldata from a dataset.
        Gap.interpolate: Interpolate gap.
        Gap.raw_model_gapfill: Raw modeldata gapfill method.
        Gap.debias_model_gapfill: Debiased modeldata gapfill method.
        Gap.weighted_diurnal_debias_model_gapfill: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description of the linear gap fill:

        1. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good observations of the leading and trailing periods are selected and grouped per timestamp.
        3. Each group (corresponding to a timestamp) is checked if they fulfill the conditions.
        4. The modeldata is interpolated (in time) to the missing records, the leading, and the trailing period.
        5. A bias for each group is computed by combining the corresponding leading and trailing groups.
        6. The gap is updated with the interpolated modeldata, corrected by the calculated bias corresponding to the specific timestamp.

        Notes
        -------
        This method requires inter-day records. The timestamps for which the
        biases are computed, are the same timestamps as found in the records.

        Examples
        ----------
        See ``Dataset.fill_gaps_with_diurnal_debiased_modeldata`` for examples.

        """
        # 1. Get leading and trailing info

        # Create anchor periods
        anchordf = gap_filling._create_anchor_df_for_leading_trailing_periods(
            Gap=self,
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            trailing_period_duration=trailing_period_duration,
        )

        # Check validity of anchors and label them
        anchordf = gap_filling._label_anchors_for_diurnal_gapfilling(
            anchordf=anchordf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
            min_anchors_per_diurnal_timestamp=min_debias_sample_size,
        )

        self._set_anchordf(anchordf)

        # 4. combine learning and gap period
        filldf = gap_filling._combine_learning_and_gap_to_one_df(
            Modeldata=Model,
            anchordf=self.anchordf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
        )

        # 5. apply blacklisting (remove anchors that do not fulfill criteria)
        blacklist = self.anchordf[self.anchordf["msg"] != "ok"].index
        filldf = filldf[~filldf.index.isin(blacklist)]

        # 6. Fill the gap
        fill_values, msg = gap_filling._diurnal_debias_modeldata_fill(filldf)

        # 7. Update gapdf attribute
        gapdf = self.gapdf
        gapdf[f"{self.obstype.name}_fill"] = fill_values
        gapdf["msg"] = msg
        gapdf["fill_method"] = label_def["diurnal_debias_modeldata_fill"]["label"]
        gapdf.loc[gapdf[f"{self.obstype.name}_fill"].isna(), "fill_method"] = label_def[
            "failed_diurnal_debias_modeldata_fill"
        ]["label"]

        self._set_gapdf(gapdf)

    def weighted_diurnal_debias_model_gapfill(
        self,
        Dataset,
        Model,
        leading_period_duration="48h",
        min_lead_debias_sample_size=3,
        trailing_period_duration="48h",
        min_trail_debias_sample_size=3,
    ):
        """Fill the gap with (weighted-diurnal-debiased) modeldata.

        The gap will be updated with fill values using the modeldata. The
        modeldata is interpolated (in time) to the missing records,
        and corrected with a bias-correction. Multiple biasses are computed, one
        for each timestamp present in the missing records, by using a leading
        (before the gap) and trailing (after the gap) period. Each bias is
        computed at each timestamp, thus computing a diurnal-bias-cycle.

        The modeldata values, used for filling the gaps are corrected by a
        weighted sum of the diurnal biases as they are computed for the leading
        and trailing period. The weights represent the normalized distance (in
        time) to the leading and trailing period.


        The Gap (data) attributes will be updated.

        Parameters
        ----------
        Dataset : metobs_toolkit.Dataset or metobs_toolkit.Station
            The dataset that holds the records of the same station as the Gap
            is from.
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        leading_period_duration : str or pandas.Timedelta, optional
            The duration of the leading period. The default is "48h".
        min_lead_debias_sample_size : int, optional
            The minimum number of good records in the leading period, to
            calculate a diurnal bias. The default is 2.
        trailing_period_duration : str or pandas.Timedelta, optional
            The duration of the trailing period. The default is "48h".
        min_trail_debias_sample_size : int, optional
            The minimum number of good records in the trailing period, to
            calculate a diurnal bias. The default is 2.

        Returns
        -------
        None.

        See Also
        --------
        Gap: The Gap class.
        Dataset.fill_gaps_with_weighted_diurnal_debias_modeldata: Equivalent for all gaps in a Dataset.
        GeeDynamicModelData: The Gee Model data (timeseries).
        Dataset.get_modeldata: Method for creating a modeldata from a dataset.
        Gap.interpolate: Interpolate gap.
        Gap.raw_model_gapfill: Raw modeldata gapfill method.
        Gap.debias_model_gapfill: Debiased modeldata gapfill method.
        Gap.diurnal_debias_model_gapfill: Diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description:

        1. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        2. The good observations of the leading and trailing periods are selected and grouped per timestamp.
        3. Each group (corresponding to a timestamp) is checked if they fulfill the conditions.
        4. The modeldata is interpolated (in time) to the missing records, the leading, and the trailing period.
        5. A bias for each group is computed for the leading and trailing groups separatly.
        6. Two weights are assigned to each missing record, that is the normalized distance to the leading and trailing period respectively.
        7. The gap is updated with the interpolated modeldata, corrected by the weighted sum of calculated bias corresponding to the specific timestamp.

        Notes
        -------
        This method requires inter-day records. The timestamps for which the
        biases are computed, are the same timestamps as found in the records.

        Examples
        ----------
        See ``Dataset.fill_gaps_with_weighted_diurnal_debias_modeldata`` for examples.

        """

        obsname = self.obstype.name
        # Create anchor periods
        # TODO: make shure datetime index is datetime dtype
        anchordf = gap_filling._create_anchor_df_for_leading_trailing_periods(
            Gap=self,
            Dataset=Dataset,
            leading_period_duration=leading_period_duration,
            trailing_period_duration=trailing_period_duration,
        )

        # Check validity of anchors and label them
        lead_anchorsdf = anchordf[anchordf["fill_method"] == "leading period"]
        lead_anchorsdf = gap_filling._label_anchors_for_diurnal_gapfilling(
            anchordf=lead_anchorsdf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
            min_anchors_per_diurnal_timestamp=min_lead_debias_sample_size,
        )

        trail_anchorsdf = anchordf[anchordf["fill_method"] == "trailing period"]
        trail_anchorsdf = gap_filling._label_anchors_for_diurnal_gapfilling(
            anchordf=trail_anchorsdf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
            min_anchors_per_diurnal_timestamp=min_trail_debias_sample_size,
        )

        anchordf = pd.concat([lead_anchorsdf, trail_anchorsdf])
        self._set_anchordf(anchordf)

        # 4. combine learning and gap period
        filldf = gap_filling._combine_learning_and_gap_to_one_df(
            Modeldata=Model,
            anchordf=self.anchordf,
            gapdf=self.gapdf,
            obsname=self.obstype.name,
        )

        # 5. apply blacklisting (remove anchors that do not fulfill criteria)
        blacklist = self.anchordf[self.anchordf["msg"] != "ok"].index
        filldf = filldf[~filldf.index.isin(blacklist)]

        # 6. Fill the gap
        fill_values, msg = gap_filling._weighted_diurnal_debias_modeldata(filldf)

        # 7. Update gapdf attribute
        gapdf = self.gapdf
        gapdf[f"{obsname}_fill"] = fill_values
        gapdf["msg"] = msg
        gapdf["fill_method"] = label_def["weighted_diurnal_debias_modeldata_fill"][
            "label"
        ]
        gapdf.loc[gapdf[f"{self.obstype.name}_fill"].isna(), "fill_method"] = label_def[
            "failed_weighted_diurnal_debias_modeldata_fill"
        ]["label"]

        self._set_gapdf(gapdf)

    # =============================================================================
    #  Get Anchor methods
    # =============================================================================
    def _get_leading_period(
        self,
        observations_series,
        leading_period_duration="24h",
    ):
        """Get the leading period of a Gap

        The leading period are the (good) records in advance of a gap.

        Parameters
        ----------
        observation_series : pandas.Series
            The observation values for the specific station and obstype. This
            is thus a pandas.Series with datetime as index. Outliers can be
            represented by nan's.
        leading_period_duration : str, optional
            A timedelta string, represents the maximum duration of the
            leading period.
        Returns
        -------
        tuple : (name, timestamps, obs-value)
            All information on the leading records is stored in a tuple.

        """
        # sta_obs = xs_save(Dataset.df, self.name, level="name")
        # sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = observations_series.dropna()  # remove the Nans (from qc)
        sta_obs_records = sta_obs.index

        period_start = self.startdt - pd.Timedelta(leading_period_duration)
        leading_period = sta_obs_records[
            ((sta_obs_records < self.startdt) & (sta_obs_records >= period_start))
        ]

        if leading_period.empty:  # No records found
            return tuple(
                [
                    self.name,
                    pd.DatetimeIndex([]),
                    np.nan,
                ]
            )

        else:
            return tuple(
                [
                    self.name,
                    leading_period,
                    sta_obs.loc[leading_period].to_list(),
                ]
            )

    def _get_trailing_period(self, observations_series, trailing_period_duration="24h"):
        """Get the trailing period of a Gap

        The trailing period are the (good) records just after a gap.

        Parameters
        ----------
        observation_series : pandas.Series
            The observation values for the specific station and obstype. This
            is thus a pandas.Series with datetime as index. Outliers can be
            represented by nan's.
        trailing_period_duration : str, optional
            A timedelta string, represents the maximum duration of the
            trailing period.
        Returns
        -------
        tuple : (name, timestamps, obs-value)
            All information on the trailing record is stored in a tuple.
        """
        # sta_obs = xs_save(Dataset.df, self.name, level="name")
        # sta_obs = sta_obs[[self.obstype.name]]  # subset to the gap variable
        sta_obs = observations_series.dropna()  # remove the Nans (from qc)
        sta_obs_records = sta_obs.index

        period_end = self.enddt + pd.Timedelta(trailing_period_duration)
        trailing_period = sta_obs_records[
            ((sta_obs_records > self.enddt) & (sta_obs_records <= period_end))
        ]

        if trailing_period.empty:  # No records found
            return tuple(
                [
                    self.name,
                    pd.DatetimeIndex([]),
                    [],
                ]
            )

        else:
            return tuple(
                [
                    self.name,
                    trailing_period,
                    sta_obs.loc[trailing_period].to_list(),
                ]
            )


# =============================================================================
# Find gaps and missing values
# =============================================================================
def get_station_gaps(gapslist, name):
    """Subset to gaps of a single station by name.

    Parameters
    ----------
    gapslist : list of gaps
        The list of gaps to subset.
    name : String
        Name of the station to extract a Gaps_collection from.

    Returns
    -------
    list
        All gaps of the specified station.

    """
    return [gap for gap in gapslist if gap.name == name]


# =============================================================================
# Gap finders
# =============================================================================
def find_gaps(df, metadf, outliersdf, obstypes):
    """Find missing records and create Gaps of them.

    Gaps are scanned for per station. The records (of a station) are assumed
    to have a perfect frequency, which is defined in the metadf 'dataset_resolution'.
    Each record is tested if it occurs exactly as expected after the previous
    record. If that is not the case, then a gap is located in between these records.

    Outliers are temporarily added to the records, to scan for gaps.


    Parameters
    ----------
    df : pandas.DataFrame()
        The df-attribute of Dataset holding the good records.
    metadf : pandas.DataFrame()
        The metadf-attribute of the Dataset holding the 'dataset_resolution' column.
    outliersdf : pandas.DataFrame()
        The outliersdf-attribute of the Dataset holding the outliers. This is
        needed because an outlier does not count (by default) as a gap.
    obstypes : dict of metobs_toolkit.Obstype's
        The obstypes-attribute of a Dataset.

    Returns
    -------
    gap_list : list
        A list Gap's .

    """

    gap_list = []
    # Replace the nans in the df, which are assigned as outliers, to a default value
    df.loc[outliersdf.index, "value"] = (
        -999
    )  # WARNING! df is a pointer, so this affects the self.df of the dataset !!!
    # make sure to revert this!!

    # Subset to the missing records
    missing_records = df[df["value"].isnull()]

    # Group by station and obstype to locate and define gaps
    for groupidx, groupdf in missing_records.reset_index().groupby(["name", "obstype"]):
        staname = groupidx[0]
        obsname = groupidx[1]
        sta_res = metadf.loc[staname, "dataset_resolution"]

        # calculate differences (at level of obstype) because else
        # gaps are computed wrong !! (DO NOT TRY TO REDUCE COMPUTATIONS ON THIS STEP! )
        groupdf["diff"] = groupdf["datetime"] - groupdf["datetime"].shift()

        # find groups of consecutive missing records
        groupdf["gap_def"] = ((groupdf["diff"] != sta_res)).cumsum()
        gapgroups = groupdf.groupby("gap_def")

        for _idx, gapgroup in gapgroups:
            # Construct gap
            gap = Gap(
                name=staname,
                startdt=gapgroup["datetime"].min(),
                enddt=gapgroup["datetime"].max(),
                obstype=obstypes[obsname],
                records_freq=sta_res,
            )

            # # Compute missing records
            # missing_records = pd.date_range(
            #     start=gap.startdt, end=gap.enddt, freq=sta_res, tz=gap.startdt.tz
            # )
            # gap._construct_gapdf(missing_datetime_records=missing_records)

            # add gap to list
            gap_list.append(gap)

    # Make sure to revert all outlier values to Nan (because they are present as -999)
    df.loc[outliersdf.index, "value"] = np.nan

    return gap_list


# =============================================================================
# Errors
# =============================================================================
class MetobsGapError(Exception):
    """Exception raised for errors in the Gap() class"""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
