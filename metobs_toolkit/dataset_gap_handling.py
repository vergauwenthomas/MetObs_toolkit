#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:33:07 2024

@author: thoverga
"""

import logging

import pandas as pd

logger = logging.getLogger(__name__)


from metobs_toolkit.gap import (
    find_gaps,
)


from metobs_toolkit.df_helpers import (
    # multiindexdf_datetime_subsetting,
    # fmt_datetime_argument,
    # init_multiindex,
    # init_multiindexdf,
    empty_outliers_df,
    # init_triple_multiindexdf,
    # metadf_to_gdf,
    # get_freqency_series,
    # value_labeled_doubleidxdf_to_triple_idxdf,
    xs_save,
    concat_save,
)

logger = logging.getLogger(__name__)


class DatasetGapCore:
    """Extension on the metobs_toolkit.Dataset class with gap related methods"""

    # =============================================================================
    # Get info and details
    # =============================================================================

    def _get_gaps_df_for_stacking(self):
        """
        Construct a long (name-obstype-datetime) df columns from Gap.gapdf.
        This is used when stacking df's together in the update_to_one_df.

        """

        # TODO: docstring
        if not bool(self.gaps):
            return empty_outliers_df()  # empyt df, same structure as outliersdf

        gapdflist = []
        for gap in self.gaps:
            # subset to relevant colunmns
            longgapdf = gap.gapdf

            # #rename
            longgapdf = longgapdf.rename(columns={f"{gap.obstype.name}_fill": "value"})

            # To long format
            longgapdf["obstype"] = gap.obstype.name
            longgapdf = (
                longgapdf.reset_index()
                .set_index(["name", "obstype", "datetime"])
                .sort_index()
            )

            # Drop the obstype column (contains only Nan's) --> not extra information
            # and ennoying when concatting multiple obstypes
            longgapdf = longgapdf.drop(columns=[gap.obstype.name])

            # add it to the list
            gapdflist.append(longgapdf)

        return pd.concat(gapdflist)

    def get_gaps_fill_df(self):
        """Construct a Dataframe with all gap-records and info

        This method constructs a dataframe that contains all the present gaps,
        exploited in records, and info. When the gaps are filled, the filled value


        Returns
        -------
        gapsdf : TYPE
            DESCRIPTION.

        """

        # TODO: docstring
        gapsdf = self._get_gaps_df_for_stacking()
        gapsdf = gapsdf.sort_index()
        return gapsdf

    # =============================================================================
    # Update gaps
    # =============================================================================

    def convert_outliers_to_gaps(self):
        """Convert all outliers to gaps.

        This method will convert all outliers to gaps, so that they can be filled.


        Returns
        -------
        None.

        Warning
        ------
        Information of the value and QC flag of the outliers will be lost.

        """
        if self.outliersdf.empty:
            logger.warning("No outliers are found to convert!")
            return

        if bool(self.gaps):
            logger.warning("The current gaps will be removed and new gaps are formed!")

        outliersdf = empty_outliers_df()

        # Since the outliers are Nan values in the dataset.df, and the find_gaps() method
        # looks for missing records AND invalid values (thus also nans), this method could be applied

        newgaps = find_gaps(
            df=self.df,
            metadf=self.metadf,
            outliersdf=outliersdf,
            obstypes=self.obstypes,
        )

        # update attributes
        self.gaps = newgaps  # overwrite previous gaps !
        self._set_outliersdf(outliersdf)

    # =============================================================================
    # Fill gaps
    # =============================================================================

    def interpolate_gaps(
        self,
        obstype="temp",
        overwrite_fill=False,
        method="time",
        max_consec_fill=10,
        n_leading_anchors=1,
        n_trailing_anchors=1,
        max_lead_to_gap_distance=None,
        max_trail_to_gap_distance=None,
        method_kwargs={},
    ):
        """Fill all the gaps using interpolation.

        The gaps of the Datasetinstance will be updated with fill values.

        Parameters
        ----------
        obstype : str, optional
            The observationtype to fill the gaps of. The default is "temp".
        overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwitten. The default is False.
        method : str, optional
            Intepolation technique to use. See pandas.DataFrame.interpolate
            'method' argument for possible values. Make shure that
            `n_leading_anchors`, `n_trailing_anchors` and `method_kwargs` are
            set accordingly to the method. The default is "time".
        max_consec_fill : int, optional
            The maximum number of consecutive missing records to fill. The default is 10.
        n_leading_anchor : int, optional
            The number of leading anchors to use for the interpolation. Higher
            order polynomial interpolation techniques require multiple leading
            anchors. The default is 1.
        n_trailing_anchor : int, optional
            The number of trailing anchors to use for the interpolation. Higher
            order polynomial interpolation techniques require multiple trailing
            anchors. The default is 1.
        max_lead_to_gap_distance : str or pandas.Timedelta, optional
            The maximum time difference between the start of the gap and a
            suitable lead (= the good record to start the interpolation from).
            If None, the first occuring good records before the start of the gap
            is used. The default is None.
        max_trail_to_gap_distance : str or pandas.Timedelta, optional
            The maximum time difference between the end of the gap and a
            suitable trail (= the good record to end the interpolation on).
            If None, the first occuring good records after the gap
            is used. The default is None.
        method_kwargs: dict, optional
            A dictionary of kwargs passed to pandas.Dataframe.interpolate(). In
            pracktice, extra arguments for specific interpolation methods are
            put in method_kwargs. The default is {}.

        Returns
        -------
        None.

        Notes
        -----
        A schematic description of the linear gap fill:

        1. Iterate over all gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. Find a leading (the last observations before the gap) record and a trailing record (the last observation after the gap).
        4. Check if the leading and trailing records fulfill the critea of maximum timediffernece.
        4. By using the leading and trailing record an interpolation is applied to fill the missing records. A maximum consecutive fill threshold is applied, if exceeded the fill values are Nan's.
        5. The gap is updated with the interpolated values

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        fill_gaps_with_raw_modeldata: Gapfill method using raw modeldata


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
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                  *28 stations
                  *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                  *10080 observation records
                  *0 records labeled as outliers
                  *2 gaps
                  *1473 missing observations
                  *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                  *time zone of the records: UTC
                  *Coordinates are available for all stations.
            >>>
            >>> #Update the gapfill settings (else the defaults are used)
            >>> dataset.update_gap_and_missing_fill_settings(gap_interpolation_max_consec_fill=35)
            >>>
            >>> # Fill the gaps
            >>> dataset.fill_gaps_linear(obstype='temp')
                                                      temp   temp_final_label
            name      datetime
            vlinder05 2022-09-06 21:00:00+00:00  21.378710  gap_interpolation
                      2022-09-06 22:00:00+00:00  21.357419  gap_interpolation
                      2022-09-06 23:00:00+00:00  21.336129  gap_interpolation
                      2022-09-07 00:00:00+00:00  21.314839  gap_interpolation
                      2022-09-07 01:00:00+00:00  21.293548  gap_interpolation
                      2022-09-07 02:00:00+00:00  21.272258  gap_interpolation
                      2022-09-07 03:00:00+00:00  21.250968  gap_interpolation
                      2022-09-07 04:00:00+00:00  21.229677  gap_interpolation
                      2022-09-07 05:00:00+00:00  21.208387  gap_interpolation
                      2022-09-07 06:00:00+00:00  21.187097  gap_interpolation
                      2022-09-07 07:00:00+00:00  21.165806  gap_interpolation
                      2022-09-07 08:00:00+00:00  21.144516  gap_interpolation
                      2022-09-07 09:00:00+00:00  21.123226  gap_interpolation
                      2022-09-07 10:00:00+00:00  21.101935  gap_interpolation
                      2022-09-07 11:00:00+00:00  21.080645  gap_interpolation
                      2022-09-07 12:00:00+00:00  21.059355  gap_interpolation
                      2022-09-07 13:00:00+00:00  21.038065  gap_interpolation
                      2022-09-07 14:00:00+00:00  21.016774  gap_interpolation
                      2022-09-07 15:00:00+00:00  20.995484  gap_interpolation
                      2022-09-07 16:00:00+00:00  20.974194  gap_interpolation
                      2022-09-07 17:00:00+00:00  20.952903  gap_interpolation
                      2022-09-07 18:00:00+00:00  20.931613  gap_interpolation
                      2022-09-07 19:00:00+00:00  20.910323  gap_interpolation
                      2022-09-07 20:00:00+00:00  20.889032  gap_interpolation
                      2022-09-07 21:00:00+00:00  20.867742  gap_interpolation
                      2022-09-07 22:00:00+00:00  20.846452  gap_interpolation
                      2022-09-07 23:00:00+00:00  20.825161  gap_interpolation
                      2022-09-08 00:00:00+00:00  20.803871  gap_interpolation
                      2022-09-08 01:00:00+00:00  20.782581  gap_interpolation
                      2022-09-08 02:00:00+00:00  20.761290  gap_interpolation
                      2022-09-08 03:00:00+00:00  20.740000  gap_interpolation
                      2022-09-08 04:00:00+00:00  20.718710  gap_interpolation
                      2022-09-08 05:00:00+00:00  20.697419  gap_interpolation
                      2022-09-08 06:00:00+00:00  20.676129  gap_interpolation
                      2022-09-08 07:00:00+00:00  20.654839  gap_interpolation
            >>> dataset.get_gaps_info()
            Gap for vlinder05 with:...


        Note
        -------
        The inpact of `max_consec_fill` is highly depending on the resolution
        of your records.

        """
        max_lead_to_gap_distance = self._timedelta_arg_check(max_lead_to_gap_distance)
        max_trail_to_gap_distance = self._timedelta_arg_check(max_trail_to_gap_distance)

        # TODO logging
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                if not gap._can_be_filled(overwrite_fill):
                    logger.warning(
                        f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                    )
                    continue
                else:
                    logger.debug(f"filling {gap} with {method} interpolation.")
                    gap.interpolate_gap(
                        Dataset=self,
                        method=method,
                        max_consec_fill=max_consec_fill,
                        n_leading_anchors=n_leading_anchors,
                        n_trailing_anchors=n_trailing_anchors,
                        max_lead_to_gap_distance=max_lead_to_gap_distance,
                        max_trail_to_gap_distance=max_trail_to_gap_distance,
                        method_kwargs=method_kwargs,
                    )

    def fill_gaps_with_raw_modeldata(
        self, Modeldata, obstype="temp", overwrite_fill=False
    ):
        """Fill all the gaps using raw modeldata.

        The gaps of the Datasetinstance will be updated with fill values by
        directly interpolating Modeldata to the missing records.


        Parameters
        ----------
        Modeldata : metobs_toolkit.Modeldata
            The modeldata that is used to fill the gaps records. The modeldata
            must be compatibel with the Dataset to fill the gaps (names,
            obstypes, period)
        obstype : str, optional
            The observationtype to fill the gaps of. The default is "temp".
        overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwitten. The default is False.

        Returns
        -------
        None.

        Notes
        -----
        A schematic description of the raw modeldata gap fill:

        1. Iterate over all gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The modeldata is interpolated (in time) to the missing records.
        4. The gap is updated with the interpolated values from the modeldata.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gapfill method.
        fill_gaps_with_diurnal_debiased_modeldata: Diurnal debiased modeldata gapfill method.


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
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                  *28 stations
                  *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                  *10080 observation records
                  *0 records labeled as outliers
                  *2 gaps
                  *1473 missing observations
                  *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                  *time zone of the records: UTC
                  *Coordinates are available for all stations.
            >>>
            >>> #Update the gapfill settings (else the defaults are used)
            >>> dataset.update_gap_and_missing_fill_settings(gap_interpolation_max_consec_fill=35)
            >>>
            >>> # Fill the gaps
            >>> dataset.fill_gaps_with_raw_modeldata(obstype='temp')
                                                      temp   temp_final_label
            name      datetime
            vlinder05 2022-09-06 21:00:00+00:00  21.378710  gap_interpolation
                      2022-09-06 22:00:00+00:00  21.357419  gap_interpolation
                      2022-09-06 23:00:00+00:00  21.336129  gap_interpolation
                      2022-09-07 00:00:00+00:00  21.314839  gap_interpolation
                      2022-09-07 01:00:00+00:00  21.293548  gap_interpolation
                      2022-09-07 02:00:00+00:00  21.272258  gap_interpolation
                      2022-09-07 03:00:00+00:00  21.250968  gap_interpolation
                      2022-09-07 04:00:00+00:00  21.229677  gap_interpolation
                      2022-09-07 05:00:00+00:00  21.208387  gap_interpolation
                      2022-09-07 06:00:00+00:00  21.187097  gap_interpolation
                      2022-09-07 07:00:00+00:00  21.165806  gap_interpolation
                      2022-09-07 08:00:00+00:00  21.144516  gap_interpolation
                      2022-09-07 09:00:00+00:00  21.123226  gap_interpolation
                      2022-09-07 10:00:00+00:00  21.101935  gap_interpolation
                      2022-09-07 11:00:00+00:00  21.080645  gap_interpolation
                      2022-09-07 12:00:00+00:00  21.059355  gap_interpolation
                      2022-09-07 13:00:00+00:00  21.038065  gap_interpolation
                      2022-09-07 14:00:00+00:00  21.016774  gap_interpolation
                      2022-09-07 15:00:00+00:00  20.995484  gap_interpolation
                      2022-09-07 16:00:00+00:00  20.974194  gap_interpolation
                      2022-09-07 17:00:00+00:00  20.952903  gap_interpolation
                      2022-09-07 18:00:00+00:00  20.931613  gap_interpolation
                      2022-09-07 19:00:00+00:00  20.910323  gap_interpolation
                      2022-09-07 20:00:00+00:00  20.889032  gap_interpolation
                      2022-09-07 21:00:00+00:00  20.867742  gap_interpolation
                      2022-09-07 22:00:00+00:00  20.846452  gap_interpolation
                      2022-09-07 23:00:00+00:00  20.825161  gap_interpolation
                      2022-09-08 00:00:00+00:00  20.803871  gap_interpolation
                      2022-09-08 01:00:00+00:00  20.782581  gap_interpolation
                      2022-09-08 02:00:00+00:00  20.761290  gap_interpolation
                      2022-09-08 03:00:00+00:00  20.740000  gap_interpolation
                      2022-09-08 04:00:00+00:00  20.718710  gap_interpolation
                      2022-09-08 05:00:00+00:00  20.697419  gap_interpolation
                      2022-09-08 06:00:00+00:00  20.676129  gap_interpolation
                      2022-09-08 07:00:00+00:00  20.654839  gap_interpolation
            >>> dataset.get_gaps_info()
            Gap for vlinder05 with:...

        """

        # check if modeldata has the obstype
        if obstype not in Modeldata.df.columns:
            raise MetobsDatasetGapHandlingError(
                f"{obstype} is not a present observationtype in {Modeldata}."
            )

        # TODO logging
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                if not gap._can_be_filled(overwrite_fill):
                    logger.warning(
                        f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                    )
                    continue
                else:
                    logger.debug(f"filling {gap} with Raw modeldata")
                    gap.raw_model_gapfill(Dataset=self, Modeldata=Modeldata)

    def fill_gaps_with_debiased_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="24h",
        min_leading_records_total=60,
        trailing_period_duration="24h",
        min_trailing_records_total=60,
    ):
        """Fill all the gaps using debiased modeldata.


        The gaps of the Datasetinstance will be updated with fill values using
        Modeldata. The Modeldata is interpolated (in time) to the missing records,
        and corrected with a bias-correction. The bias is estimated by making use
        of a leading (before the gap) and trailing (after the gap) period.

        Parameters
        ----------
        Modeldata : metobs_toolkit.Modeldata
            The modeldata that is used to fill the gaps records. The modeldata
            must be compatibel with the Dataset to fill the gaps (names,
            obstypes, period)
        obstype : str, optional
            The observationtype to fill the gaps of. The default is "temp".
        overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwitten. The default is False.
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

        Notes
        -----
        A schematic description of the debiased modeldata gap fill:

        1. Iterate over all gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good obervations of the leading and trailing period are selected and checked if they fulfill the conditions.
        4. The modeldata is interpolated (in time) to the missing records, the leading, and trailing period.
        5. By combining the leading and trailing period both records and modeldata, a bias is calculated.
        6. The gap is updated with the interpolated modeldata, corrected by the calculated bias.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_diurnal_debiased_modeldata: Diurnal debiased modeldata gapfill method.


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
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                  *28 stations
                  *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                  *10080 observation records
                  *0 records labeled as outliers
                  *2 gaps
                  *1473 missing observations
                  *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                  *time zone of the records: UTC
                  *Coordinates are available for all stations.
            >>>
            >>> #Update the gapfill settings (else the defaults are used)
            >>> dataset.update_gap_and_missing_fill_settings(gap_interpolation_max_consec_fill=35)
            >>>
            >>> # Fill the gaps
            >>> dataset.fill_gaps_with_debiased_modeldata(obstype='temp')
                                                      temp   temp_final_label
            name      datetime
            vlinder05 2022-09-06 21:00:00+00:00  21.378710  gap_interpolation
                      2022-09-06 22:00:00+00:00  21.357419  gap_interpolation
                      2022-09-06 23:00:00+00:00  21.336129  gap_interpolation
                      2022-09-07 00:00:00+00:00  21.314839  gap_interpolation
                      2022-09-07 01:00:00+00:00  21.293548  gap_interpolation
                      2022-09-07 02:00:00+00:00  21.272258  gap_interpolation
                      2022-09-07 03:00:00+00:00  21.250968  gap_interpolation
                      2022-09-07 04:00:00+00:00  21.229677  gap_interpolation
                      2022-09-07 05:00:00+00:00  21.208387  gap_interpolation
                      2022-09-07 06:00:00+00:00  21.187097  gap_interpolation
                      2022-09-07 07:00:00+00:00  21.165806  gap_interpolation
                      2022-09-07 08:00:00+00:00  21.144516  gap_interpolation
                      2022-09-07 09:00:00+00:00  21.123226  gap_interpolation
                      2022-09-07 10:00:00+00:00  21.101935  gap_interpolation
                      2022-09-07 11:00:00+00:00  21.080645  gap_interpolation
                      2022-09-07 12:00:00+00:00  21.059355  gap_interpolation
                      2022-09-07 13:00:00+00:00  21.038065  gap_interpolation
                      2022-09-07 14:00:00+00:00  21.016774  gap_interpolation
                      2022-09-07 15:00:00+00:00  20.995484  gap_interpolation
                      2022-09-07 16:00:00+00:00  20.974194  gap_interpolation
                      2022-09-07 17:00:00+00:00  20.952903  gap_interpolation
                      2022-09-07 18:00:00+00:00  20.931613  gap_interpolation
                      2022-09-07 19:00:00+00:00  20.910323  gap_interpolation
                      2022-09-07 20:00:00+00:00  20.889032  gap_interpolation
                      2022-09-07 21:00:00+00:00  20.867742  gap_interpolation
                      2022-09-07 22:00:00+00:00  20.846452  gap_interpolation
                      2022-09-07 23:00:00+00:00  20.825161  gap_interpolation
                      2022-09-08 00:00:00+00:00  20.803871  gap_interpolation
                      2022-09-08 01:00:00+00:00  20.782581  gap_interpolation
                      2022-09-08 02:00:00+00:00  20.761290  gap_interpolation
                      2022-09-08 03:00:00+00:00  20.740000  gap_interpolation
                      2022-09-08 04:00:00+00:00  20.718710  gap_interpolation
                      2022-09-08 05:00:00+00:00  20.697419  gap_interpolation
                      2022-09-08 06:00:00+00:00  20.676129  gap_interpolation
                      2022-09-08 07:00:00+00:00  20.654839  gap_interpolation
            >>> dataset.get_gaps_info()
            Gap for vlinder05 with:...

        """

        # check if modeldata has the obstype
        if obstype not in Modeldata.df.columns:
            raise MetobsDatasetGapHandlingError(
                f"{obstype} is not a present observationtype in {Modeldata}."
            )

        # TODO logging
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                if not gap._can_be_filled(overwrite_fill):
                    logger.warning(
                        f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                    )
                    continue
                else:
                    logger.debug(f"filling {gap} with Debiased modeldata")
                    gap.debias_model_gapfill(
                        Dataset=self,
                        Modeldata=Modeldata,
                        leading_period_duration=leading_period_duration,
                        min_leading_records_total=min_leading_records_total,
                        trailing_period_duration=trailing_period_duration,
                        min_trailing_records_total=min_trailing_records_total,
                    )

    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="24h",
        min_debias_sample_size=6,
        trailing_period_duration="24h",
    ):
        # TODO docstring

        # check if modeldata has the obstype
        if obstype not in Modeldata.df.columns:
            raise MetobsDatasetGapHandlingError(
                f"{obstype} is not a present observationtype in {Modeldata}."
            )

        # TODO logging
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                if not gap._can_be_filled(overwrite_fill):
                    logger.warning(
                        f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                    )
                    continue
                else:
                    logger.debug(f"filling {gap} with Diurnal debiased modeldata")
                    gap.diurnal_debias_model_gapfill(
                        Dataset=self,
                        Modeldata=Modeldata,
                        leading_period_duration=leading_period_duration,
                        min_debias_sample_size=min_debias_sample_size,
                        trailing_period_duration=trailing_period_duration,
                    )

    def fill_gaps_with_weighted_diurnal_debias_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="48h",
        min_lead_debias_sample_size=2,
        trailing_period_duration="48h",
        min_trail_debias_sample_size=2,
    ):
        """Fill all the gaps using diurnal-debiased modeldata.


         The gaps of the Datasetinstance will be updated with fill values using
         Modeldata. The Modeldata is interpolated (in time) to the missing records,
         and corrected with a bias-correction. Multiple biasses are computed, one
         for each timestamp present in the missing records, by using a leading
         (before the gap) and trailing (after the gap) period. Each bias is
         computed at each timestamp, thus computing a diurnal-bias-cycle.


         Parameters
         ----------
         Modeldata : metobs_toolkit.Modeldata
             The modeldata that is used to fill the gaps records. The modeldata
             must be compatibel with the Dataset to fill the gaps (names,
             obstypes, period)
         obstype : str, optional
             The observationtype to fill the gaps of. The default is "temp".
         overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwitten. The default is False.
         leading_period_duration : str or pandas.Timedelta, optional
             The duration of the leading period. The default is "48h".
         min_lead_debias_sample_size : int, optional
             The minimum number of good records in the leading period, to
             calculate a diurnal bias of. The default is 2.
         trailing_period_duration : str or pandas.Timedelta, optional
             The duration of the trailing period. The default is "48h".
         min_trail_debias_sample_size : int, optional
             The minimum number of good records in the trailing period, to
             calculate a diurnal bias of. The default is 2.


         Returns
         -------
         None.


        Notes
        -----
        A schematic description of the linear gap fill:

        1. Iterate over all gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good obervations of the leading and trailing period are selected and grouped per timestamp.
        4. Each group (coresponding to a timestamp) is checked if they fulfill the conditions.
        5. The modeldata is interpolated (in time) to the missing records, the leading, and trailing period.
        6. A bias for each group is computed by combining the corresponding leading and trailing groups.
        7. The gap is updated with the interpolated modeldata, corrected by the calculated bias corresponding to the specific timestamp.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gapfill method.


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
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                  *28 stations
                  *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                  *10080 observation records
                  *0 records labeled as outliers
                  *2 gaps
                  *1473 missing observations
                  *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                  *time zone of the records: UTC
                  *Coordinates are available for all stations.
            >>>
            >>> #Update the gapfill settings (else the defaults are used)
            >>> dataset.update_gap_and_missing_fill_settings(gap_interpolation_max_consec_fill=35)
            >>>
            >>> # Fill the gaps
            >>> dataset.fill_gaps_with_weighted_diurnal_debias_modeldata(obstype='temp')
                                                      temp   temp_final_label
            name      datetime
            vlinder05 2022-09-06 21:00:00+00:00  21.378710  gap_interpolation
                      2022-09-06 22:00:00+00:00  21.357419  gap_interpolation
                      2022-09-06 23:00:00+00:00  21.336129  gap_interpolation
                      2022-09-07 00:00:00+00:00  21.314839  gap_interpolation
                      2022-09-07 01:00:00+00:00  21.293548  gap_interpolation
                      2022-09-07 02:00:00+00:00  21.272258  gap_interpolation
                      2022-09-07 03:00:00+00:00  21.250968  gap_interpolation
                      2022-09-07 04:00:00+00:00  21.229677  gap_interpolation
                      2022-09-07 05:00:00+00:00  21.208387  gap_interpolation
                      2022-09-07 06:00:00+00:00  21.187097  gap_interpolation
                      2022-09-07 07:00:00+00:00  21.165806  gap_interpolation
                      2022-09-07 08:00:00+00:00  21.144516  gap_interpolation
                      2022-09-07 09:00:00+00:00  21.123226  gap_interpolation
                      2022-09-07 10:00:00+00:00  21.101935  gap_interpolation
                      2022-09-07 11:00:00+00:00  21.080645  gap_interpolation
                      2022-09-07 12:00:00+00:00  21.059355  gap_interpolation
                      2022-09-07 13:00:00+00:00  21.038065  gap_interpolation
                      2022-09-07 14:00:00+00:00  21.016774  gap_interpolation
                      2022-09-07 15:00:00+00:00  20.995484  gap_interpolation
                      2022-09-07 16:00:00+00:00  20.974194  gap_interpolation
                      2022-09-07 17:00:00+00:00  20.952903  gap_interpolation
                      2022-09-07 18:00:00+00:00  20.931613  gap_interpolation
                      2022-09-07 19:00:00+00:00  20.910323  gap_interpolation
                      2022-09-07 20:00:00+00:00  20.889032  gap_interpolation
                      2022-09-07 21:00:00+00:00  20.867742  gap_interpolation
                      2022-09-07 22:00:00+00:00  20.846452  gap_interpolation
                      2022-09-07 23:00:00+00:00  20.825161  gap_interpolation
                      2022-09-08 00:00:00+00:00  20.803871  gap_interpolation
                      2022-09-08 01:00:00+00:00  20.782581  gap_interpolation
                      2022-09-08 02:00:00+00:00  20.761290  gap_interpolation
                      2022-09-08 03:00:00+00:00  20.740000  gap_interpolation
                      2022-09-08 04:00:00+00:00  20.718710  gap_interpolation
                      2022-09-08 05:00:00+00:00  20.697419  gap_interpolation
                      2022-09-08 06:00:00+00:00  20.676129  gap_interpolation
                      2022-09-08 07:00:00+00:00  20.654839  gap_interpolation
            >>> dataset.get_gaps_info()
            Gap for vlinder05 with:...

        """

        # check if modeldata has the obstype
        if obstype not in Modeldata.df.columns:
            raise MetobsDatasetGapHandlingError(
                f"{obstype} is not a present observationtype in {Modeldata}."
            )

        # TODO logging
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                if not gap._can_be_filled(overwrite_fill):
                    logger.warning(
                        f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                    )
                    continue
                else:
                    logger.debug(
                        f"filling {gap} with Weighted diurnal debiased modeldata"
                    )
                    gap.weighted_diurnal_debias_model_gapfill(
                        Dataset=self,
                        Modeldata=Modeldata,
                        leading_period_duration=leading_period_duration,
                        min_lead_debias_sample_size=min_lead_debias_sample_size,
                        trailing_period_duration=trailing_period_duration,
                        min_trail_debias_sample_size=min_trail_debias_sample_size,
                    )

    def fill_gaps_with_weighted_diurnal_debias_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="48h",
        min_lead_debias_sample_size=2,
        trailing_period_duration="48h",
        min_trail_debias_sample_size=2,
    ):
        """Fill all the gaps using weighted-diurnal-debiased modeldata.


         The gaps of the Datasetinstance will be updated with fill values using
         Modeldata. The Modeldata is interpolated (in time) to the missing records,
         and corrected with a bias-correction. Multiple biasses are computed, one
         for each timestamp present in the missing records, by using a leading
         (before the gap) and trailing (after the gap) period. Each bias is
         computed at each timestamp, thus computing a diurnal-bias-cycle.

         The modeldata values, used for filling the gaps are corrected by a
         weighted sum of the diurnal biases as they are computed for the leading
         and trailing period. The weights represent the normalized distance (in
         time) to the leading and trailing period.


         Parameters
         ----------
         Modeldata : metobs_toolkit.Modeldata
             The modeldata that is used to fill the gaps records. The modeldata
             must be compatibel with the Dataset to fill the gaps (names,
             obstypes, period)
         obstype : str, optional
             The observationtype to fill the gaps of. The default is "temp".
         overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwitten. The default is False.
         leading_period_duration : str or pandas.Timedelta, optional
             The duration of the leading period. The default is "48h".
         min_lead_debias_sample_size : int, optional
             The minimum number of good records in the leading period, to
             calculate a diurnal bias of. The default is 2.
         trailing_period_duration : str or pandas.Timedelta, optional
             The duration of the trailing period. The default is "48h".
         min_trail_debias_sample_size : int, optional
             The minimum number of good records in the trailing period, to
             calculate a diurnal bias of. The default is 2.


         Returns
         -------
         None.


        Notes
        -----
        A schematic description of the linear gap fill:

        1. Iterate over all gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good obervations of the leading and trailing period are selected and grouped per timestamp.
        4. Each group (coresponding to a timestamp) is checked if they fulfill the conditions.
        5. The modeldata is interpolated (in time) to the missing records, the leading, and trailing period.
        5. A bias for each group is computed for the leading and trailing groups seperatly.
        6. Two weights are assigned to each missing record, that is the normalized distance to the leading an trailing period respectively.
        7. The gap is updated with the interpolated modeldata, corrected by the weighted sum of calculated bias corresponding to the specific timestamp.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gapfill method.


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
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                  *28 stations
                  *['temp', 'humidity', 'wind_speed', 'wind_direction'] observation types
                  *10080 observation records
                  *0 records labeled as outliers
                  *2 gaps
                  *1473 missing observations
                  *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                  *time zone of the records: UTC
                  *Coordinates are available for all stations.
            >>>
            >>> #Update the gapfill settings (else the defaults are used)
            >>> dataset.update_gap_and_missing_fill_settings(gap_interpolation_max_consec_fill=35)
            >>>
            >>> # Fill the gaps
            >>> dataset.fill_gaps_with_weighted_diurnal_debias_modeldata(obstype='temp')
                                                      temp   temp_final_label
            name      datetime
            vlinder05 2022-09-06 21:00:00+00:00  21.378710  gap_interpolation
                      2022-09-06 22:00:00+00:00  21.357419  gap_interpolation
                      2022-09-06 23:00:00+00:00  21.336129  gap_interpolation
                      2022-09-07 00:00:00+00:00  21.314839  gap_interpolation
                      2022-09-07 01:00:00+00:00  21.293548  gap_interpolation
                      2022-09-07 02:00:00+00:00  21.272258  gap_interpolation
                      2022-09-07 03:00:00+00:00  21.250968  gap_interpolation
                      2022-09-07 04:00:00+00:00  21.229677  gap_interpolation
                      2022-09-07 05:00:00+00:00  21.208387  gap_interpolation
                      2022-09-07 06:00:00+00:00  21.187097  gap_interpolation
                      2022-09-07 07:00:00+00:00  21.165806  gap_interpolation
                      2022-09-07 08:00:00+00:00  21.144516  gap_interpolation
                      2022-09-07 09:00:00+00:00  21.123226  gap_interpolation
                      2022-09-07 10:00:00+00:00  21.101935  gap_interpolation
                      2022-09-07 11:00:00+00:00  21.080645  gap_interpolation
                      2022-09-07 12:00:00+00:00  21.059355  gap_interpolation
                      2022-09-07 13:00:00+00:00  21.038065  gap_interpolation
                      2022-09-07 14:00:00+00:00  21.016774  gap_interpolation
                      2022-09-07 15:00:00+00:00  20.995484  gap_interpolation
                      2022-09-07 16:00:00+00:00  20.974194  gap_interpolation
                      2022-09-07 17:00:00+00:00  20.952903  gap_interpolation
                      2022-09-07 18:00:00+00:00  20.931613  gap_interpolation
                      2022-09-07 19:00:00+00:00  20.910323  gap_interpolation
                      2022-09-07 20:00:00+00:00  20.889032  gap_interpolation
                      2022-09-07 21:00:00+00:00  20.867742  gap_interpolation
                      2022-09-07 22:00:00+00:00  20.846452  gap_interpolation
                      2022-09-07 23:00:00+00:00  20.825161  gap_interpolation
                      2022-09-08 00:00:00+00:00  20.803871  gap_interpolation
                      2022-09-08 01:00:00+00:00  20.782581  gap_interpolation
                      2022-09-08 02:00:00+00:00  20.761290  gap_interpolation
                      2022-09-08 03:00:00+00:00  20.740000  gap_interpolation
                      2022-09-08 04:00:00+00:00  20.718710  gap_interpolation
                      2022-09-08 05:00:00+00:00  20.697419  gap_interpolation
                      2022-09-08 06:00:00+00:00  20.676129  gap_interpolation
                      2022-09-08 07:00:00+00:00  20.654839  gap_interpolation
            >>> dataset.get_gaps_info()
            Gap for vlinder05 with:...

        """

        # TODO docstring

        # check if modeldata has the obstype
        if obstype not in Modeldata.df.columns:
            raise MetobsDatasetGapHandlingError(
                f"{obstype} is not a present observationtype in {Modeldata}."
            )

        # TODO logging
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                if not gap._can_be_filled(overwrite_fill):
                    logger.warning(
                        f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                    )
                    continue
                else:
                    logger.debug(
                        f"filling {gap} with Weighted diurnal debiased modeldata"
                    )
                    gap.weighted_diurnal_debias_model_gapfill(
                        Dataset=self,
                        Modeldata=Modeldata,
                        leading_period_duration=leading_period_duration,
                        min_lead_debias_sample_size=min_lead_debias_sample_size,
                        trailing_period_duration=trailing_period_duration,
                        min_trail_debias_sample_size=min_trail_debias_sample_size,
                    )

    def find_gap(self, stationname, obstype, in_gap_timestamp):
        """Find a specific gap

        This method helps to find a specific gap in your Dataset.

        Parameters
        ----------
        stationname : str
            The name of the station to look for a gap.
        obstype : str
            The obstypename of the gap.
        in_gap_timestamp : datetime.datetime or pandas.Timestamp
            A timestamp that is located in the gap. If a datetime naive timestamp
            is given, it is assumed to be in the same timzone as the records.

        Returns
        -------
        gap : metobs_toolkit.Gap or None
            The specific gap. If no gap was found that matches the arguments,
            None is returned


        """

        in_gap_timestamp = self._datetime_arg_check(in_gap_timestamp)

        for gap in self.gaps:
            if gap.name != str(stationname):
                continue
            if gap.obstype.name != str(obstype):
                continue
            if not (gap.startdt <= in_gap_timestamp):
                continue

            if not (gap.enddt >= in_gap_timestamp):
                continue

            return gap

        logger.warning(
            f"No gap is found with: ({stationname}) - ({obstype}) - ({in_gap_timestamp})"
        )
        return None


class MetobsDatasetGapHandlingError(Exception):
    """Exception raised for errors in the datasetgaphandling."""

    pass
