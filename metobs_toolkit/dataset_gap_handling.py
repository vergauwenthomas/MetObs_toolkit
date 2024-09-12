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
    empty_outliers_df,
    xs_save,
    concat_save,
)

logger = logging.getLogger(__name__)


class DatasetGapCore:
    """Extension on the metobs_toolkit.Dataset class with gap-related methods"""

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
        gapsdf : pandas.DataFrame
            A dataframe with ['name', 'obstype', 'datetime'] as index, and
            [value, fill_method, msg] as columns. If a gap is not filled, the
            default values are used.

        See Also
        -----------
        Dataset.get_full_status_df : Combine all records, outliers and gaps in one Dataframe.

        Examples
        --------

        We start by creating a Dataset, and importing data.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> print(dataset)
        Dataset instance containing:
             *28 stations
             *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
             *483828 observation records (not Nan's)
             *0 records labeled as outliers
             *8 gaps
             *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
             *time zone of the records: UTC
             *Coordinates are available for all stations.
             *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

        We can combine all the gaps in one dataframe:

        >>> gapdf = dataset.get_gaps_fill_df()
        >>> gapdf
                                                            value fill_method   msg
        name      obstype        datetime
        vlinder02 humidity       2022-09-10 17:10:00+00:00    NaN  not filled  None
                                 2022-09-10 17:15:00+00:00    NaN  not filled  None
                                 2022-09-10 17:45:00+00:00    NaN  not filled  None
                  temp           2022-09-10 17:10:00+00:00    NaN  not filled  None
                                 2022-09-10 17:15:00+00:00    NaN  not filled  None
        ...                                                   ...         ...   ...
                  wind_direction 2022-09-10 17:15:00+00:00    NaN  not filled  None
                                 2022-09-10 17:45:00+00:00    NaN  not filled  None
                  wind_speed     2022-09-10 17:10:00+00:00    NaN  not filled  None
                                 2022-09-10 17:15:00+00:00    NaN  not filled  None
                                 2022-09-10 17:45:00+00:00    NaN  not filled  None
        <BLANKLINE>
        [12 rows x 3 columns]

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
        Information on the value and QC flag of the outliers will be lost.

        See Also
        -----------
        Dataset.get_full_status_df : Combine all records, outliers and gaps in one Dataframe.
        Dataset.get_gaps_fill_df: Construct a Dataframe with all gap-records and info

        Examples
        --------

        We start by creating a Dataset, and importing data.

        >>> import metobs_toolkit
        >>>
        >>> #Create your Dataset
        >>> dataset = metobs_toolkit.Dataset() #empty Dataset
        >>> dataset.import_data_from_file(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )

        We can now apply quality control on the records.

        >>> dataset.apply_quality_control('temp')
        >>> print(dataset)
        Dataset instance containing:
             *28 stations
             *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
             *447513 observation records (not Nan's)
             *36315 records labeled as outliers
             *8 gaps
             *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
             *time zone of the records: UTC
             *Coordinates are available for all stations.
             *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

        Now we convert the outliers to gaps, so that they can be filled.

        >>> dataset.convert_outliers_to_gaps()
        >>> print(dataset)
        Dataset instance containing:
             *28 stations
             *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
             *447513 observation records (not Nan's)
             *0 records labeled as outliers
             *1697 gaps
             *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:55:00+00:00 (total duration:  14 days 23:55:00)
             *time zone of the records: UTC
             *Coordinates are available for all stations.
             *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']


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

        The gaps of the Dataset will be updated with fill values.

        Parameters
        ----------
        obstype : str, optional
            The observationtype to fill the gaps. The default is "temp".
        overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gap fill
            values will be overwritten. The default is False.
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
            The number of trailing anchors to use for the interpolation.
            Higher-order polynomial interpolation techniques require multiple trailing
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
        None.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gapfill method.
        fill_gaps_with_diurnal_debiased_modeldata: Diurnal debiased modeldata gapfill method.
        fill_gaps_with_weighted_diurnal_debiased_modeldata: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description:

        1. Iterate over all the gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. Find a leading (the last observations before the gap) record and a trailing record (the last observation after the gap).
        4. Check if the leading and trailing records fulfill the criteria of maximum timedifference.
        5. By using the leading and trailing record an interpolation is applied to fill the missing records. A maximum consecutive fill threshold is applied, if exceeded the fill values are Nan's.
        6. The gap is updated with the interpolated values

        Note
        -------
        The impact of `max_consec_fill` is highly dependent on the resolution
        of your records.

        Note
        ------
        If you want to use a higher-order method of interpolation, make sure to
        increase the `n_leading_anchors` and `n_trailing_anchors` accordingly.


        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            To reduce the data for this example, we coarsen the data to hourly records.

            >>> dataset.coarsen_time_resolution(freq='1h')

            To create some gaps, we apply quality control first and then convert
            the outliers to gaps.

            >>> dataset.apply_quality_control('temp')
            >>> dataset.convert_outliers_to_gaps()
            >>> print(dataset)
            Dataset instance containing:
                 *28 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *38644 observation records (not Nan's)
                 *0 records labeled as outliers
                 *89 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            As we can see, we now have a dataset with gaps (for temperature). It is
            often handy to combine all present gaps into one pandas Dataframe, for
            inspection.

            >>> comb_gap_df = dataset.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value fill_method   msg
            name      obstype datetime
            vlinder01 temp    2022-09-02 16:00:00+00:00    NaN  not filled  None
                              2022-09-02 17:00:00+00:00    NaN  not filled  None
                              2022-09-02 18:00:00+00:00    NaN  not filled  None
                              2022-09-02 19:00:00+00:00    NaN  not filled  None
                              2022-09-02 20:00:00+00:00    NaN  not filled  None
            ...                                            ...         ...   ...
            vlinder28 temp    2022-09-15 04:00:00+00:00    NaN  not filled  None
                              2022-09-15 05:00:00+00:00    NaN  not filled  None
                              2022-09-15 06:00:00+00:00    NaN  not filled  None
                              2022-09-15 07:00:00+00:00    NaN  not filled  None
                              2022-09-15 08:00:00+00:00    NaN  not filled  None
            <BLANKLINE>
            [1676 rows x 3 columns]


            Now we are going to fill the gaps using interpolation.

            >>> dataset.interpolate_gaps(
            ...                 obstype="temp",
            ...                 overwrite_fill=False,
            ...                 method="time", #prefered over linear
            ...                 max_consec_fill=5) # thus interpolate maximum 5 hours (for hourly records)

            We can inspect the filled gaps by plotting or by using `get_gaps_fill_df()` method.

            >>> dataset.make_plot(obstype='temp', colorby='label', title='After interpolation of the gaps')
            <Axes: title={'center': 'After interpolation of the gaps'}, xlabel='datetime', ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            >>> comb_gap_df = dataset.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value           fill_method                       msg
            name      obstype datetime
            vlinder01 temp    2022-09-02 16:00:00+00:00  25.92         interpolation         ok (method: time)
                              2022-09-02 17:00:00+00:00  25.04         interpolation         ok (method: time)
                              2022-09-02 18:00:00+00:00  24.15         interpolation         ok (method: time)
                              2022-09-02 19:00:00+00:00  23.27         interpolation         ok (method: time)
                              2022-09-02 20:00:00+00:00  22.39         interpolation         ok (method: time)
            ...                                            ...                   ...                       ...
            vlinder28 temp    2022-09-15 04:00:00+00:00  12.98         interpolation         ok (method: time)
                              2022-09-15 05:00:00+00:00  13.80         interpolation         ok (method: time)
                              2022-09-15 06:00:00+00:00  14.62         interpolation         ok (method: time)
                              2022-09-15 07:00:00+00:00    NaN  failed interpolation  Permitted_by_max_cons...
                              2022-09-15 08:00:00+00:00    NaN  failed interpolation  Permitted_by_max_cons...
            <BLANKLINE>
            [1676 rows x 3 columns]



            More advanced interpolation methods can be used. As an example we use
            spline interpolation on the gaps of a single station.

            >>> dataset.get_station('vlinder02').interpolate_gaps(
            ...             obstype="temp",
            ...             overwrite_fill=True,
            ...             method="cubicspline",
            ...             max_consec_fill=7,
            ...             n_leading_anchors=3,#depends on method requirements
            ...             n_trailing_anchors=3, #depends on method requirements
            ...             max_lead_to_gap_distance='4h',
            ...             max_trail_to_gap_distance='4h',
            ...             method_kwargs={'bc_type':'not-a-knot'})
            >>> dataset.get_station('vlinder02').get_gaps_fill_df()
                                                         value           fill_method                       msg
            name      obstype datetime
            vlinder02 temp    2022-09-02 16:00:00+00:00  27.76         interpolation  ok (method: cubicspli...
                              2022-09-02 17:00:00+00:00  27.67         interpolation  ok (method: cubicspli...
                              2022-09-02 18:00:00+00:00  27.11         interpolation  ok (method: cubicspli...
                              2022-09-02 19:00:00+00:00  26.18         interpolation  ok (method: cubicspli...
                              2022-09-02 20:00:00+00:00  24.99         interpolation  ok (method: cubicspli...
            ...                                            ...                   ...                       ...
                              2022-09-09 05:00:00+00:00  14.79         interpolation  ok (method: cubicspli...
                              2022-09-09 06:00:00+00:00  15.11         interpolation  ok (method: cubicspli...
                              2022-09-09 07:00:00+00:00  15.47         interpolation  ok (method: cubicspli...
                              2022-09-09 08:00:00+00:00    NaN  failed interpolation  Permitted_by_max_cons...
                              2022-09-09 09:00:00+00:00    NaN  failed interpolation  Permitted_by_max_cons...
            <BLANKLINE>
            [51 rows x 3 columns]

            >>> dataset.get_station('vlinder02').make_plot(obstype='temp', colorby='label')
            <Axes: title={'center': 'Temperatuur of vlinder02'}, xlabel='datetime', ylabel='temp (Celsius)'>



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
                    gap.interpolate(
                        Dataset=self,
                        method=method,
                        max_consec_fill=max_consec_fill,
                        n_leading_anchors=n_leading_anchors,
                        n_trailing_anchors=n_trailing_anchors,
                        max_lead_to_gap_distance=max_lead_to_gap_distance,
                        max_trail_to_gap_distance=max_trail_to_gap_distance,
                        method_kwargs=method_kwargs,
                    )

    def fill_gaps_with_raw_modeldata(self, Model, obstype="temp", overwrite_fill=False):
        """Fill all the gaps using raw modeldata.

        The gaps of the Datasetinstance will be updated with fill values by
        directly interpolating Modeldata to the missing records.


        Parameters
        ----------
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        obstype : str, optional
            The observationtype to fill the gaps. The default is "temp".
        overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwritten. The default is False.

        Returns
        -------
        None.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gap filled records and info.
        metobs_toolkit.GeeDynamicModelData: The Gee Model data (timeseries).
        get_modeldata: Method for creating a modeldata from a dataset.
        interpolate_gaps: Fill gaps by interpolation.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gap fill method.
        fill_gaps_with_diurnal_debiased_modeldata: Diurnal debiased modeldata gap fill method.
        fill_gaps_with_weighted_diurnal_debiased_modeldata: Weighted diurnal debiased modeldata gap fill method.

        Notes
        -----
        A schematic description of the raw modeldata gap fill:

        1. Iterate over all the gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The modeldata is interpolated (in time) to the missing records.
        4. The gap is updated with the interpolated values from the modeldata.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            To reduce the data for this example, we coarsen the data to hourly records
            and focus on the records of a single station. This example is will also
            work on a full Dataset with multiple stations.

            >>> dataset.coarsen_time_resolution(freq='1h')
            >>> sta = dataset.get_station('vlinder05')

            To create some gaps, we apply quality control first and then convert
            the outliers to gaps.

            >>> sta.apply_quality_control('temp')
            >>> sta.convert_outliers_to_gaps()
            >>> print(sta)
            Station instance containing:
                 *1 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *1160 observation records (not Nan's)
                 *0 records labeled as outliers
                 *7 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            As we can see, we now have a Station (or Dataset) with gaps (for temperature). It is
            often handy to combine all present gaps into one pandas Dataframe, for
            inspection.

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value fill_method   msg
            name      obstype datetime
            vlinder05 temp    2022-09-01 00:00:00+00:00    NaN  not filled  None
                              2022-09-01 01:00:00+00:00    NaN  not filled  None
                              2022-09-01 02:00:00+00:00    NaN  not filled  None
                              2022-09-01 03:00:00+00:00    NaN  not filled  None
                              2022-09-01 04:00:00+00:00    NaN  not filled  None
            ...                                            ...         ...   ...
                              2022-09-15 19:00:00+00:00    NaN  not filled  None
                              2022-09-15 20:00:00+00:00    NaN  not filled  None
                              2022-09-15 21:00:00+00:00    NaN  not filled  None
                              2022-09-15 22:00:00+00:00    NaN  not filled  None
                              2022-09-15 23:00:00+00:00    NaN  not filled  None
            <BLANKLINE>
            [280 rows x 3 columns]

            Since we will use Modeldata to fill the gap, we first need to import
            modeldata. We can do this by importing ERA5 (Land) data directly from
            the Google Earth Engine by using the `Dataset.get_modeldata()` method.

            >>> era5_data = sta.get_modeldata(Model=sta.gee_datasets['ERA5-land'],
            ...                               obstypes=['temp'])

            (When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is ERA5_hourly)
            >>> # For large datafiles, the modeldata is writen to a csv file. See Dataset.get_modeldata() for more info.
            >>> print(era5_data)
            GeeDynamicModelData instance of ERA5-land with modeldata

            Now we are going to fill the gaps with this raw modeldata.

            >>> sta.fill_gaps_with_raw_modeldata(
            ...                    Model=era5_data,
            ...                    obstype="temp")

            We can inspect the filled gaps by plotting or by using `get_gaps_fill_df()` method.

            >>> sta.make_plot(obstype='temp', colorby='label', title='After filling the gaps')
            <Axes: title={'center': 'After filling the gaps'}, xlabel='datetime', ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value         fill_method                       msg
            name      obstype datetime
            vlinder05 temp    2022-09-01 00:00:00+00:00  18.22  raw modeldata fill  Modelvalue: 18.22, wi...
                              2022-09-01 01:00:00+00:00  17.68  raw modeldata fill  Modelvalue: 17.68, wi...
                              2022-09-01 02:00:00+00:00  17.31  raw modeldata fill  Modelvalue: 17.31, wi...
                              2022-09-01 03:00:00+00:00  16.74  raw modeldata fill  Modelvalue: 16.74, wi...
                              2022-09-01 04:00:00+00:00  16.42  raw modeldata fill  Modelvalue: 16.42, wi...
            ...                                            ...                 ...                       ...
                              2022-09-15 19:00:00+00:00  15.20  raw modeldata fill  Modelvalue: 15.20, wi...
                              2022-09-15 20:00:00+00:00  14.58  raw modeldata fill  Modelvalue: 14.58, wi...
                              2022-09-15 21:00:00+00:00  14.02  raw modeldata fill  Modelvalue: 14.02, wi...
                              2022-09-15 22:00:00+00:00  13.48  raw modeldata fill  Modelvalue: 13.48, wi...
                              2022-09-15 23:00:00+00:00  13.26  raw modeldata fill  Modelvalue: 13.26, wi...
            <BLANKLINE>
            [280 rows x 3 columns]

        """
        # Check if the Model has the compatible data
        _check_if_model_can_be_used(trg_obstypename=obstype, Model=Model)

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
                    gap.raw_model_gapfill(Dataset=self, Model=Model)

    def fill_gaps_with_debiased_modeldata(
        self,
        Model,
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
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        obstype : str, optional
            The observationtype to fill the gaps. The default is "temp".
        overwrite_fill : bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values will be overwritten. The default is False.
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
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        get_modeldata: Method for creating a modeldata from a dataset.
        interpolate_gaps: Fill gaps by interpolation.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_diurnal_debiased_modeldata: Diurnal debiased modeldata gapfill method.
        fill_gaps_with_weighted_diurnal_debiased_modeldata: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description of the debiased modeldata gap fill:

        1. Iterate over all the gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good observations of the leading and trailing periods are selected and checked if they fulfill the conditions.
        4. The modeldata is interpolated (in time) to the missing records, the leading, and the trailing period.
        5. By combining the leading and trailing periods both records and modeldata, a bias is calculated.
        6. The gap is updated with the interpolated modeldata, corrected by the calculated bias.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            To reduce the data for this example, we coarsen the data to hourly records
            and focus on the records of a single station. This example is will also
            work on a full Dataset with multiple stations.

            >>> dataset.coarsen_time_resolution(freq='1h')
            >>> sta = dataset.get_station('vlinder05')

            To create some gaps, we apply quality control first and then convert
            the outliers to gaps.

            >>> sta.apply_quality_control('temp')
            >>> sta.convert_outliers_to_gaps()
            >>> print(sta)
            Station instance containing:
                 *1 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *1160 observation records (not Nan's)
                 *0 records labeled as outliers
                 *7 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            As we can see, we now have a Station (or Dataset) with gaps (for temperature). It is
            often handy to combine all present gaps into one pandas Dataframe, for
            inspection.

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value fill_method   msg
            name      obstype datetime
            vlinder05 temp    2022-09-01 00:00:00+00:00    NaN  not filled  None
                              2022-09-01 01:00:00+00:00    NaN  not filled  None
                              2022-09-01 02:00:00+00:00    NaN  not filled  None
                              2022-09-01 03:00:00+00:00    NaN  not filled  None
                              2022-09-01 04:00:00+00:00    NaN  not filled  None
            ...                                            ...         ...   ...
                              2022-09-15 19:00:00+00:00    NaN  not filled  None
                              2022-09-15 20:00:00+00:00    NaN  not filled  None
                              2022-09-15 21:00:00+00:00    NaN  not filled  None
                              2022-09-15 22:00:00+00:00    NaN  not filled  None
                              2022-09-15 23:00:00+00:00    NaN  not filled  None
            <BLANKLINE>
            [280 rows x 3 columns]

            Since we will use Modeldata to fill the gap, we first need to import
            modeldata. We can do this by importing ERA5 (Land) data directly from
            the Google Earth Engine by using the `Dataset.get_modeldata()` method.

            >>> era5_data = sta.get_modeldata(Model=sta.gee_datasets['ERA5-land'],
            ...                               obstypes=['temp'])

            (When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is ERA5_hourly)
            >>> # For large datafiles, the modeldata is writen to a csv file. See Dataset.get_modeldata() for more info.
            >>> print(era5_data)
            GeeDynamicModelData instance of ERA5-land with modeldata

            Now we are going to fill the gaps with this debiased modeldata. Do this
            by specifying a leading and trailing period (by duration and the minimum
            number of records), and using the `Dataset.fill_gaps_with_debiased_modeldata()` method.

            >>> sta.fill_gaps_with_debiased_modeldata(
            ...                Model=era5_data,
            ...                obstype="temp",
            ...                overwrite_fill=False,
            ...                leading_period_duration="24h",
            ...                min_leading_records_total=10, #higly depending on resolution
            ...                trailing_period_duration="24h",
            ...                min_trailing_records_total=10) #higly depending on resolution

            We can inspect the filled gaps by plotting or by using `get_gaps_fill_df()` method.

            >>> sta.make_plot(obstype='temp', colorby='label', title='After filling the gaps')
            <Axes: title={'center': 'After filling the gaps'}, xlabel='datetime', ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df[:10]
                                                             value               fill_method                       msg
            name      obstype datetime
            vlinder05 temp    2022-09-01 00:00:00+00:00    NaN  failed debiased model...  minimum records (10) ...
                              2022-09-01 01:00:00+00:00    NaN  failed debiased model...  minimum records (10) ...
                              2022-09-01 02:00:00+00:00    NaN  failed debiased model...  minimum records (10) ...
                              2022-09-01 03:00:00+00:00    NaN  failed debiased model...  minimum records (10) ...
                              2022-09-01 04:00:00+00:00    NaN  failed debiased model...  minimum records (10) ...
                              2022-09-01 05:00:00+00:00    NaN  failed debiased model...  minimum records (10) ...
                              2022-09-01 20:00:00+00:00  19.60   debiased modeldata fill  Modelvalue: 20.20 cor...
                              2022-09-01 21:00:00+00:00  18.87   debiased modeldata fill  Modelvalue: 19.47 cor...
                              2022-09-01 22:00:00+00:00  18.17   debiased modeldata fill  Modelvalue: 18.78 cor...
                              2022-09-01 23:00:00+00:00  17.63   debiased modeldata fill  Modelvalue: 18.23 cor...

        """

        # Check if the Model has the compatible data
        _check_if_model_can_be_used(trg_obstypename=obstype, Model=Model)

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
                        Model=Model,
                        leading_period_duration=leading_period_duration,
                        min_leading_records_total=min_leading_records_total,
                        trailing_period_duration=trailing_period_duration,
                        min_trailing_records_total=min_trailing_records_total,
                    )

    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        Model,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="24h",
        min_debias_sample_size=6,
        trailing_period_duration="24h",
    ):
        """Fill all the gaps using diurnal debiased modeldata.


        The gaps in the Dataset will be updated with fill values using
        Modeldata. The Modeldata is interpolated (in time) to the missing records,
        and corrected with a bias-correction. Multiple biasses are computed, one
        for each timestamp present in the missing records, by using a leading
        (before the gap) and trailing (after the gap) period. Each bias is
        computed at each timestamp, thus computing a diurnal-bias-cycle.


        Parameters
        ----------
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        obstype : str, optional
            The observationtype to fill the gaps. The default is "temp".
        overwrite_fill : bool, optional
           If a gap has already filled values, the interpolation of this gap
           is skipped if overwrite_fill is False. If set to True, the gapfill
           values will be overwritten. The default is False.
        leading_period_duration : str or pandas.Timedelta, optional
            The duration of the leading period. The default is "24h".
        min_debias_sample_size : int, optional
            The minimum number of good records to
            calculate a diurnal bias of. The default is 6.
        trailing_period_duration : str or pandas.Timedelta, optional
            The duration of the trailing period. The default is "24h".


        Returns
        -------
        None.

        See Also
        --------
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        get_modeldata: Method for creating a modeldata from a dataset.
        interpolate_gaps: Fill gaps by interpolation.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gapfill method.
        fill_gaps_with_weighted_diurnal_debiased_modeldata: Weighted diurnal debiased modeldata gapfill method.

        Notes
        -----
        A schematic description of the linear gap fill:

        1. Iterate over all the gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good observations of the leading and trailing periods are selected and grouped per timestamp.
        4. Each group (corresponding to a timestamp) is checked if they fulfill the conditions.
        5. The modeldata is interpolated (in time) to the missing records, the leading, and the trailing period.
        6. A bias for each group is computed by combining the corresponding leading and trailing groups.
        7. The gap is updated with the interpolated modeldata, corrected by the calculated bias corresponding to the specific timestamp.

        Notes
        -------
        This method requires inter-day records. The timestamps for which the
        biases are computed, are the same timestamps as found in the records.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            To reduce the data for this example, we coarsen the data to hourly records
            and focus on the records of a single station. This example is will also
            work on a full Dataset with multiple stations.

            >>> dataset.coarsen_time_resolution(freq='1h')
            >>> sta = dataset.get_station('vlinder02')

            To create some gaps, we apply quality control first and then convert
            the outliers to gaps.

            >>> sta.apply_quality_control('temp')
            >>> sta.convert_outliers_to_gaps()
            >>> print(sta)
            Station instance containing:
                 *1 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *1389 observation records (not Nan's)
                 *0 records labeled as outliers
                 *3 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            As we can see, we now have a Station (or Dataset) with gaps (for temperature). It is
            often handy to combine all present gaps into one pandas Dataframe, for
            inspection.

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value fill_method   msg
            name      obstype datetime
            vlinder02 temp    2022-09-02 16:00:00+00:00    NaN  not filled  None
                              2022-09-02 17:00:00+00:00    NaN  not filled  None
                              2022-09-02 18:00:00+00:00    NaN  not filled  None
                              2022-09-02 19:00:00+00:00    NaN  not filled  None
                              2022-09-02 20:00:00+00:00    NaN  not filled  None
            ...                                            ...         ...   ...
                              2022-09-09 05:00:00+00:00    NaN  not filled  None
                              2022-09-09 06:00:00+00:00    NaN  not filled  None
                              2022-09-09 07:00:00+00:00    NaN  not filled  None
                              2022-09-09 08:00:00+00:00    NaN  not filled  None
                              2022-09-09 09:00:00+00:00    NaN  not filled  None
            <BLANKLINE>
            [51 rows x 3 columns]



            Since we will use Modeldata to fill the gap, we first need to import
            modeldata. We can do this by importing ERA5 (Land) data directly from
            the Google Earth Engine by using the `Dataset.get_modeldata()` method.

            >>> era5_data = sta.get_modeldata(Model=sta.gee_datasets['ERA5-land'],
            ...                               obstypes=['temp'])

            (When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is ERA5_hourly)
            >>> # For large datafiles, the modeldata is writen to a csv file. See Dataset.get_modeldata() for more info.
            >>> print(era5_data)
            GeeDynamicModelData instance of ERA5-land with modeldata

            Now we are going to fill the gaps with this diurnal debiased modeldata.
            Do this by specifying a leading and trailing period (by duration and
            the minimum number of records), and using the
            `Dataset.fill_gaps_with_diurnal_debiased_modeldata()` method.

            >>> sta.fill_gaps_with_diurnal_debiased_modeldata(
            ...                Model=era5_data,
            ...                obstype="temp",
            ...                overwrite_fill=False,
            ...                leading_period_duration="24h",
            ...                min_debias_sample_size=2,
            ...                trailing_period_duration="24h")

            We can inspect the filled gaps by plotting or by using `get_gaps_fill_df()` method.

            >>> sta.make_plot(obstype='temp', colorby='label', title='After filling the gaps')
            <Axes: title={'center': 'After filling the gaps'}, xlabel='datetime', ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value               fill_method                       msg
            name      obstype datetime
            vlinder02 temp    2022-09-02 16:00:00+00:00  25.65  diurnal debiased mode...  Modelvalue: 25.65 cor...
                              2022-09-02 17:00:00+00:00  25.11  diurnal debiased mode...  Modelvalue: 25.17 cor...
                              2022-09-02 18:00:00+00:00  24.03  diurnal debiased mode...  Modelvalue: 23.93 cor...
                              2022-09-02 19:00:00+00:00  21.49  diurnal debiased mode...  Modelvalue: 21.87 cor...
                              2022-09-02 20:00:00+00:00  21.84  diurnal debiased mode...  Modelvalue: 21.23 cor...
            ...                                            ...                       ...                       ...
                              2022-09-09 05:00:00+00:00    NaN  failed diurnal debias...  Modelvalue: 15.52 can...
                              2022-09-09 06:00:00+00:00    NaN  failed diurnal debias...  Modelvalue: 15.73 can...
                              2022-09-09 07:00:00+00:00    NaN  failed diurnal debias...  Modelvalue: 15.74 can...
                              2022-09-09 08:00:00+00:00    NaN  failed diurnal debias...  Modelvalue: 16.12 can...
                              2022-09-09 09:00:00+00:00    NaN  failed diurnal debias...  Modelvalue: 16.84 can...
            <BLANKLINE>
            [51 rows x 3 columns]

        """

        # Check if the Model has the compatible data
        _check_if_model_can_be_used(trg_obstypename=obstype, Model=Model)

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
                        Model=Model,
                        leading_period_duration=leading_period_duration,
                        min_debias_sample_size=min_debias_sample_size,
                        trailing_period_duration=trailing_period_duration,
                    )

    def fill_gaps_with_weighted_diurnal_debias_modeldata(
        self,
        Model,
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
        Model : metobs_toolkit.GeeDynamicModelData
            The model that is used to fill the gaps records. The modeldata
            must be compatible (same metadata and `ModelObstype` equivalent
            of obstype) to fill the gaps.
        obstype : str, optional
            The observationtype to fill the gaps. The default is "temp".
        overwrite_fill : bool, optional
           If a gap has already filled values, the interpolation of this gap
           is skipped if overwrite_fill is False. If set to True, the gapfill
           values will be overwritten. The default is False.
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
        get_gaps_fill_df: Get an overview dataframe of gapfilled records and info.
        metobs_toolkit.Modeldata: Modeldata class.
        get_modeldata: Method for creating a modeldata from a dataset.
        interpolate_gaps: Fill gaps by interpolation.
        fill_gaps_with_raw_modeldata: Raw modeldata gapfill method.
        fill_gaps_with_debiased_modeldata: Debiased modeldata gapfill method.
        fill_gaps_with_diurnal_debiased_modeldata: Diurnal debiased modeldata gapfill method.


        Notes
        -----
        A schematic description of the linear gap fill:

        1. Iterate over all the gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. The good observations of the leading and trailing periods are selected and grouped per timestamp.
        4. Each group (corresponding to a timestamp) is checked if they fulfill the conditions.
        5. The modeldata is interpolated (in time) to the missing records, the leading, and the trailing period.
        6. A bias for each group is computed for the leading and trailing groups seperatly.
        7. Two weights are assigned to each missing record, that is the normalized distance to the leading and trailing period respectively.
        8. The gap is updated with the interpolated modeldata, corrected by the weighted sum of calculated bias corresponding to the specific timestamp.

        Notes
        -------
        This method requires inter-day records. The timestamps for which the
        biases are computed, are the same timestamps as found in the records.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data.

            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )

            To reduce the data for this example, we coarsen the data to hourly records
            and focus on the records of a single station. This example is will also
            work on a full Dataset with multiple stations.

            >>> dataset.coarsen_time_resolution(freq='1h')
            >>> sta = dataset.get_station('vlinder02')

            To create some gaps, we apply quality control first and then convert
            the outliers to gaps.

            >>> sta.apply_quality_control('temp')
            >>> sta.convert_outliers_to_gaps()
            >>> print(sta)
            Station instance containing:
                 *1 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *1389 observation records (not Nan's)
                 *0 records labeled as outliers
                 *3 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']

            As we can see, we now have a Station (or Dataset) with gaps (for temperature). It is
            often handy to combine all present gaps into one pandas Dataframe, for
            inspection.

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                         value fill_method   msg
            name      obstype datetime
            vlinder02 temp    2022-09-02 16:00:00+00:00    NaN  not filled  None
                              2022-09-02 17:00:00+00:00    NaN  not filled  None
                              2022-09-02 18:00:00+00:00    NaN  not filled  None
                              2022-09-02 19:00:00+00:00    NaN  not filled  None
                              2022-09-02 20:00:00+00:00    NaN  not filled  None
            ...                                            ...         ...   ...
                              2022-09-09 05:00:00+00:00    NaN  not filled  None
                              2022-09-09 06:00:00+00:00    NaN  not filled  None
                              2022-09-09 07:00:00+00:00    NaN  not filled  None
                              2022-09-09 08:00:00+00:00    NaN  not filled  None
                              2022-09-09 09:00:00+00:00    NaN  not filled  None
            <BLANKLINE>
            [51 rows x 3 columns]



            Since we will use Modeldata for filling the gap, we first need to import
            modeldata. We can do this by importing ERA5 (Land) data directly from
            the Google Earth Engine by using the `Dataset.get_modeldata()` method.

            >>> era5_data = sta.get_modeldata(Model=sta.gee_datasets['ERA5-land'],
            ...                               obstypes=['temp'])


            (When using the .set_model_from_csv() method, make sure the modelname of your Modeldata is ERA5_hourly)
            >>> # For large datafiles, the modeldata is writen to a csv file. See Dataset.get_modeldata() for more info.
            >>> print(era5_data)
            GeeDynamicModelData instance of ERA5-land with modeldata


            Now we are going to fill the gaps with weighted diurnal debiased modeldata.
            Do this by specifying a leading and trailing period (by duration and
            the minimum number of records), and using the
            `Dataset.fill_gaps_with_diurnal_debiased_modeldata()` method.

            >>> sta.fill_gaps_with_weighted_diurnal_debias_modeldata(
            ...                    Model=era5_data,
            ...                    obstype="temp",
            ...                    overwrite_fill=False,
            ...                    leading_period_duration="48h",
            ...                    min_lead_debias_sample_size=1,
            ...                    trailing_period_duration="48h",
            ...                    min_trail_debias_sample_size=1)

            We can inspect the filled gaps by plotting or by using `get_gaps_fill_df()` method.

            >>> sta.make_plot(obstype='temp', colorby='label', title='After filling the gaps')
            <Axes: title={'center': 'After filling the gaps'}, xlabel='datetime', ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            >>> comb_gap_df = sta.get_gaps_fill_df()
            >>> comb_gap_df
                                                        value               fill_method                       msg
            name      obstype datetime
            vlinder02 temp    2022-09-02 16:00:00+00:00  25.92  Weighted diurnal debi...  Modelvalue: 25.65 cor...
                              2022-09-02 17:00:00+00:00  25.13  Weighted diurnal debi...  Modelvalue: 25.17 cor...
                              2022-09-02 18:00:00+00:00  24.10  Weighted diurnal debi...  Modelvalue: 23.93 cor...
                              2022-09-02 19:00:00+00:00  21.69  Weighted diurnal debi...  Modelvalue: 21.87 cor...
                              2022-09-02 20:00:00+00:00  21.32  Weighted diurnal debi...  Modelvalue: 21.23 cor...
            ...                                            ...                       ...                       ...
                              2022-09-09 05:00:00+00:00  16.04  Weighted diurnal debi...  Modelvalue: 15.52 cor...
                              2022-09-09 06:00:00+00:00  15.85  Weighted diurnal debi...  Modelvalue: 15.73 cor...
                              2022-09-09 07:00:00+00:00    NaN  failed Weighted diurn...  Modelvalue: 15.74 can...
                              2022-09-09 08:00:00+00:00    NaN  failed Weighted diurn...  Modelvalue: 16.12 can...
                              2022-09-09 09:00:00+00:00    NaN  failed Weighted diurn...  Modelvalue: 16.84 can...
            <BLANKLINE>
            [51 rows x 3 columns]

        """

        # Check if the Model has the compatible data
        _check_if_model_can_be_used(trg_obstypename=obstype, Model=Model)

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
                        Model=Model,
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
            is given, it is assumed to be in the same timezone as the records.

        Returns
        -------
        gap : metobs_toolkit.Gap or None
            The specific gap. If no gap was found that matches the arguments,
            None is returned

        See Also
        ---------
        metobs_toolkit.Gap: The gap class.
        Gap.get_info: Print out details of the gap.
        convert_outliers_to_gaps: Convert outliers to gaps.

        Examples
        --------

        .. plot::
            :context: close-figs

            We start by creating a Dataset, and importing data. To create some gaps
            in the Dataset, we apply quality control and convert the outliers to gaps.


            >>> import metobs_toolkit
            >>>
            >>> #Create your Dataset
            >>> dataset = metobs_toolkit.Dataset() #empty Dataset
            >>> dataset.import_data_from_file(
            ...                         input_data_file=metobs_toolkit.demo_datafile,
            ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
            ...                         template_file=metobs_toolkit.demo_template,
            ...                         )
            >>> dataset.coarsen_time_resolution(freq='1h') #to reduce data amount in this example
            >>> dataset.apply_quality_control(obstype='temp')
            >>> dataset.convert_outliers_to_gaps()
            >>> print(dataset)
            Dataset instance containing:
                 *28 stations
                 *['humidity', 'temp', 'wind_direction', 'wind_speed'] observation types present
                 *38644 observation records (not Nan's)
                 *0 records labeled as outliers
                 *89 gaps
                 *records range: 2022-09-01 00:00:00+00:00 --> 2022-09-15 23:00:00+00:00 (total duration:  14 days 23:00:00)
                 *time zone of the records: UTC
                 *Coordinates are available for all stations.
                 *Known GEE datasets for: ['lcz', 'altitude', 'worldcover', 'ERA5-land']


            As we can see, we now have a Dataset with gaps (for temperature). The
            gaps are stored as a list of Gap's

            >>> dataset.gaps[:5] # print only the first 5 gaps
            [temp-gap of vlinder01 for 2022-09-02 16:00:00+00:00 --> 2022-09-03 01:00:00+00:00, duration: 0 days 09:00:00 or 10 records., temp-gap of vlinder01 for 2022-09-07 07:00:00+00:00 --> 2022-09-08 14:00:00+00:00, duration: 1 days 07:00:00 or 32 records., temp-gap of vlinder01 for 2022-09-09 01:00:00+00:00 --> 2022-09-09 09:00:00+00:00, duration: 0 days 08:00:00 or 9 records., temp-gap of vlinder02 for 2022-09-02 16:00:00+00:00 --> 2022-09-03 01:00:00+00:00, duration: 0 days 09:00:00 or 10 records., temp-gap of vlinder02 for 2022-09-07 07:00:00+00:00 --> 2022-09-08 14:00:00+00:00, duration: 1 days 07:00:00 or 32 records.]

            It is often handy to combine all present gaps into one pandas Dataframe, for
            inspection.

            >>> comb_gap_df = dataset.get_gaps_fill_df()
            >>> comb_gap_df
                                                          value fill_method   msg
            name      obstype datetime
            vlinder01 temp    2022-09-02 16:00:00+00:00    NaN  not filled  None
                              2022-09-02 17:00:00+00:00    NaN  not filled  None
                              2022-09-02 18:00:00+00:00    NaN  not filled  None
                              2022-09-02 19:00:00+00:00    NaN  not filled  None
                              2022-09-02 20:00:00+00:00    NaN  not filled  None
            ...                                            ...         ...   ...
            vlinder28 temp    2022-09-15 04:00:00+00:00    NaN  not filled  None
                              2022-09-15 05:00:00+00:00    NaN  not filled  None
                              2022-09-15 06:00:00+00:00    NaN  not filled  None
                              2022-09-15 07:00:00+00:00    NaN  not filled  None
                              2022-09-15 08:00:00+00:00    NaN  not filled  None
            <BLANKLINE>
            [1676 rows x 3 columns]


            If you want all details or need to interact directly with a specific
            Gap, you can find the Gap (from the .gaps attribute) by using the
            `Dataset.find_gap()` method. In practice, one can make a timeseries
            plot, identify where the gap-of-interest is, and then use this method.

            >>> dataset.make_plot(obstype='temp', colorby='label')
            <Axes: title={'center': 'Temperatuur for all stations. '}, xlabel='datetime', ylabel='temp (Celsius)'>

        .. plot::
            :context: close-figs

            ... Localize the gap-of-interest ...

            >>> import datetime
            >>> your_favorite_gap = dataset.find_gap(
            ...                         stationname='vlinder02',
            ...                         obstype='temp',
            ...                         in_gap_timestamp=datetime.datetime(2022, 9, 9, 4)) #2022-09-9 4:00 UTC
            >>> print(your_favorite_gap)
            temp-gap of vlinder02 for 2022-09-09 01:00:00+00:00 --> 2022-09-09 09:00:00+00:00, duration: 0 days 08:00:00 or 9 records.

            If you want more info on the gap:

            >>> your_favorite_gap.get_info()
            ---- Gap info -----
            (Note: gaps start and end are defined on the frequency estimation of the native dataset.)
              * Gap for station: vlinder02
              * Start gap: 2022-09-09 01:00:00+00:00
              * End gap: 2022-09-09 09:00:00+00:00
              * Duration gap: 0 days 08:00:00
              * For temp
              * Gapfill status >>>> Unfilled gap
            ---- Gap Data Frame -----
                                                 temp  temp_fill fill_method   msg
            name      datetime
            vlinder02 2022-09-09 01:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 02:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 03:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 04:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 05:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 06:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 07:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 08:00:00+00:00   NaN        NaN  not filled  None
                      2022-09-09 09:00:00+00:00   NaN        NaN  not filled  None


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


def _check_if_model_can_be_used(trg_obstypename, Model):
    """Raise error if Model cannot be used to fill for target obstype."""
    # check if modeldata has the obstype
    if trg_obstypename not in Model.modelobstypes.keys():
        raise MetobsDatasetGapHandlingError(
            f"{Model} does not have a known ModelObstype equivalent of {trg_obstypename}."
        )
    # Check if modeldata has timeseries
    if Model.modeldf.empty:
        raise MetobsDatasetGapHandlingError(f"{Model} does not have any modeldata.")
    # Check if obstype is present in the modeldata
    if trg_obstypename not in Model.modeldf.columns:
        raise MetobsDatasetGapHandlingError(
            f"{Model} does not have modeldata for {trg_obstypename}"
        )
    return


# =============================================================================
# Errors
# =============================================================================


class MetobsDatasetGapHandlingError(Exception):
    """Exception raised for errors in the datasetgaphandling."""

    pass


# =============================================================================
# Docstring test
# =============================================================================
if __name__ == "__main__":
    from metobs_toolkit.doctest_fmt import setup_and_run_doctest

    setup_and_run_doctest()
