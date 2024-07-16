#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:33:07 2024

@author: thoverga
"""

import logging
import sys
import pandas as pd

logger = logging.getLogger(__name__)


from metobs_toolkit import Dataset

from metobs_toolkit.missingobs import Missingob_collection

from metobs_toolkit.gap import (
    Gap,
    remove_gaps_from_obs,
    remove_gaps_from_outliers,
    missing_timestamp_and_gap_check,
    get_gaps_indx_in_obs_space,
    get_station_gaps,
    apply_interpolate_gaps,
    make_gapfill_df,
    apply_debias_era5_gapfill,
    gaps_to_df,
)


from metobs_toolkit.df_helpers import (
    multiindexdf_datetime_subsetting,
    fmt_datetime_argument,
    init_multiindex,
    init_multiindexdf,
    init_triple_multiindexdf,
    metadf_to_gdf,
    conv_applied_qc_to_df,
    get_freqency_series,
    value_labeled_doubleidxdf_to_triple_idxdf,
    xs_save,
    concat_save,
)

logger = logging.getLogger(__name__)


class Dataset(Dataset):
    """Extension on the metobs_toolkit.Dataset class with gap related methods"""

    # =============================================================================
    # Get info and details
    # =============================================================================

    # =============================================================================
    # Update gaps
    # =============================================================================

    def update_gaps_and_missing_from_outliers(self, obstype="temp", n_gapsize=None):
        """Interpret the outliers as missing observations.

        If there is a sequence
        of these outliers for a station, larger than n_gapsize than this will
        be interpreted as a gap.

        The outliers are not removed.

        Parameters
        ----------
        obstype : str, optional
            Use the outliers on this observation type to update the gaps and
            missing timestamps. The default is 'temp'.
        n_gapsize : int, optional
            The minimum number of consecutive missing observations to define
            as a gap. If None, n_gapsize is taken from the settings defenition
            of gaps. The default is None.

        Returns
        -------
        None.

        Note
        -------
        Gaps and missing observations resulting from an outlier on a specific
        obstype, are assumed to be gaps/missing observation for all obstypes.

        Note
        ------
        Be aware that n_gapsize is used for the current resolution of the Dataset,
        this is different from the gap check applied on the inported data, if
        the dataset is coarsend.

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

        """
        if n_gapsize is None:
            n_gapsize = self.settings.gap["gaps_settings"]["gaps_finder"]["gapsize_n"]
            if (
                not self.metadf["assumed_import_frequency"]
                .eq(self.metadf["dataset_resolution"])
                .all()
            ):
                logger.info(
                    f"The defenition of the gapsize (n_gapsize = {n_gapsize}) \
                                will have another effect on the update of the gaps and missing \
                                    timestamps because coarsening is applied and the defenition \
                                    of the gapsize is not changed."
                )

        # combine to one dataframe
        mergedf = self.combine_all_to_obsspace()
        mergedf = xs_save(mergedf, obstype, level="obstype")

        # ignore labels
        possible_outlier_labels = [
            vals["outlier_flag"] for vals in self.settings.qc["qc_checks_info"].values()
        ]

        # create groups when the final label changes
        persistance_filter = ((mergedf["label"].shift() != mergedf["label"])).cumsum()
        grouped = mergedf.groupby(["name", persistance_filter])

        # locate new gaps by size of consecutive the same final label per station
        group_sizes = grouped.size()
        large_groups = group_sizes[group_sizes > n_gapsize]

        # find only groups with final label as an outlier
        gaps = []
        # new_gapsdf = pd.DataFrame()
        new_gaps_idx = init_multiindex()
        for group_idx in large_groups.index:
            groupdf = grouped.get_group(group_idx)
            group_final_label = groupdf["label"].iloc[0]
            if group_final_label not in possible_outlier_labels:
                # no gap candidates
                continue
            else:
                gap = Gap(
                    name=groupdf.index.get_level_values("name")[0],
                    startdt=groupdf.index.get_level_values("datetime").min(),
                    enddt=groupdf.index.get_level_values("datetime").max(),
                )

                gaps.append(gap)
                new_gaps_idx = new_gaps_idx.union(groupdf.index, sort=False)

        # add all the outliers, that are not in the new gaps to the new missing obs
        new_missing_obs = mergedf[mergedf["label"].isin(possible_outlier_labels)].index
        new_missing_obs = new_missing_obs.drop(new_gaps_idx.to_numpy(), errors="ignore")

        # to series
        missing_obs_series = (
            new_missing_obs.to_frame()
            .reset_index(drop=True)
            .set_index("name")["datetime"]
        )
        # Create missing obs
        new_missing_collection = Missingob_collection(missing_obs_series)

        # update self
        self.gaps.extend(gaps)
        self.missing_obs = self.missing_obs + new_missing_collection

        # remove outliers that are converted to gaps
        self.outliersdf = remove_gaps_from_outliers(
            gaplist=gaps, outldf=self.outliersdf
        )

        # remove outliers that are converted to missing obs
        self.outliersdf = self.missing_obs.remove_missing_from_outliers(self.outliersdf)

    # =============================================================================
    # Fill gaps
    # =============================================================================

    def fill_gaps_era5(
        self, modeldata, method="debias", obstype="temp", overwrite_fill=False
    ):
        """Fill the gaps using a diurnal debiased modeldata approach.


        Parameters
        ----------
        modeldata : metobs_toolkit.Modeldata
            The modeldata to use for the gapfill. This model data should the required
            timeseries to fill all gaps present in the dataset.
        method : 'debias', optional
            Specify which method to use. The default is 'debias'.
        obstype : String, optional
            Name of the observationtype you want to apply gap filling on. The
            modeldata must contain this observation type as well. The
            default is 'temp'.
        overwrite_fill: bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values and info will be overwitten. The default is False.

        Returns
        -------
        Gapfilldf : pandas.DataFrame
            A dataframe containing all gap filled values and the use method.

        Notes
        -----
        A schematic description of the fill_gaps_era5 method:

        1. Modeldata is converted to the timezone of the observations.
        2. Iterate over all gaps.
            * The gap is converted into a set of missing records (depending on the time resolution of the observations).
            * Find a leading and trailing period. These periods are a subset
              of observations respectively before and after the gap. The size
              of these subsets is set by a target size (in records) and a minimum
              size (in records). If the subset of observations is smaller than
              the corresponding minimum size, the gap cannot be filled.
            * Modeldata, for the corresponding station and observation type, is extracted for the leading and trailing period.
            * By comparing the model data with the observations of the
              leading and trailing period, and grouping all records to their
              timestamp (i.g. diurnal categories), biasses are computed.
            * Modeldata for the missing records is extracted.
            * Weights ([0;1]) are computed for each gap record, representing
              the normalized distance (in time), to the beginning and end of
              the gap.
            * The modeldata at the missing records is then corrected by
              a weighted sum of the leading and trailing biases at the
              corresponding timestamp. In general, this means that the diurnal
              trend of the observations is restored as well as possible.
        3. The gap is updated with the interpolated values (metobs_toolkit.Gap.gapfill_df)

        Note
        -------
        A scientific publication on the performance of this technique is expected.

        Examples
        --------
        .. code-block:: python

            import metobs_toolkit

            your_dataset = metobs_toolkit.Dataset()
            your_dataset.update_settings(
                input_data_file=metobs_toolkit.demo_datafile, # path to the data file
                input_metadata_file=metobs_toolkit.demo_metadatafile,
                template_file=metobs_toolkit.demo_template,
            )
            # Specify the gap defenition
            your_dataset.update_qc_settings(gapsize_in_records = 20)

            #Update the gapsize BEFORE importing the data
            your_dataset.import_data_from_file()

            #Update the settings (definition of the period to calculate biases for)
            your_dataset.update_gap_and_missing_fill_settings(
                                                              gap_debias_prefered_leading_period_hours=24,
                                                              gap_debias_prefered_trailing_period_hours=24,
                                                              gap_debias_minimum_leading_period_hours=6,
                                                              gap_debias_minimum_trailing_period_hours=6,
                                                              )
            #(As a demonstration, we will fill the gaps of a single station. The following functions can also be
            # directly applied to the dataset.)
            your_station = your_dataset.get_station('vlinder05')


            #Get ERA5 modeldata at the location of your stations and period.
            ERA5_modeldata = your_station.get_modeldata(modelname='ERA5_hourly',
                                                        obstype='temp')

            #Use the debias method to fill the gaps
            gapfill_df = your_station.fill_gaps_era5(modeldata=ERA5_modeldata,
                                                      obstype='temp')

        """
        # check if modeldata is available
        if modeldata is None:
            logger.warning(
                "The dataset has no modeldate. Use the set_modeldata() function to add modeldata."
            )
            return None
        # check if obstype is present in eramodel
        assert (
            obstype in modeldata.df.columns
        ), f"{obstype} is not present in the modeldate: {modeldata}"
        # check if all station are present in eramodeldata
        # stations = self.gaps.to_df().index.unique().to_list()
        stations = list(set([gap.name for gap in self.gaps]))
        assert all(
            [sta in modeldata.df.index.get_level_values("name") for sta in stations]
        ), "Not all stations with gaps are in the modeldata!"

        if method == "debias":
            fill_settings_debias = self.settings.gap["gaps_fill_settings"][
                "model_debias"
            ]

            apply_debias_era5_gapfill(
                gapslist=self.gaps,
                dataset=self,
                eraModelData=modeldata,
                obstype=obstype,
                debias_settings=fill_settings_debias,
                overwrite_fill=overwrite_fill,
            )

            # get fill df
            filldf = make_gapfill_df(self.gaps)
        else:
            sys.exit(f"{method} not implemented yet")

        # update attribute
        self.gapfilldf = filldf

        return filldf

    def fill_gaps_automatic(
        self,
        modeldata,
        obstype="temp",
        max_interpolate_duration_str=None,
        overwrite_fill=False,
    ):
        """Fill the gaps by using linear interpolation or debiased modeldata.

        This method serves as a triage to select the gaps to be filled with
        linear interpolation and those to be filled using a diurnal debias
        gapfill. When the duration of a gap is smaller or equal than
        max_interpolation_duration, the linear interpolation method is applied
        else the debiased modeldata method.

        For a detailed description of these methods, we refer to the
        corresponding metobs_toolkit.Dataset.fill_gaps_linear() and
        metobs_toolkit.Dataset.fill_gaps_era5().

        Parameters
        ----------
        modeldata : metobs_toolkit.Modeldata
            The modeldata to use for the gapfill. This model data should the required
            timeseries to fill all gaps present in the dataset.
        obstype : String, optional
            Name of the observationtype you want to apply gap filling on. The
            modeldata must contain this observation type as well. The
            default is 'temp'.
        max_interpolate_duration_str : Timedelta or str, optional
            Maximum duration to apply interpolation for gapfill when using the
            automatic gapfill method. Gaps with longer durations will be filled
            using debiased modeldata. The default is None.
        overwrite_fill: bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values and info will be overwitten. The default is False.

        Returns
        -------
        comb_df : pandas.DataFrame
            A dataframe containing all the filled records.

        Examples
        --------
        .. code-block:: python

            import metobs_toolkit

            your_dataset = metobs_toolkit.Dataset()
            your_dataset.update_settings(
                input_data_file=metobs_toolkit.demo_datafile, # path to the data file
                input_metadata_file=metobs_toolkit.demo_metadatafile,
                template_file=metobs_toolkit.demo_template,
            )
            # Specify the gap defenition
            your_dataset.update_qc_settings(gapsize_in_records = 20)

            #Update the gapsize BEFORE importing the data
            your_dataset.import_data_from_file()

            #Update the settings (definition of the period to calculate biases for)
            your_dataset.update_gap_and_missing_fill_settings(
                                                              gap_debias_prefered_leading_period_hours=24,
                                                              gap_debias_prefered_trailing_period_hours=24,
                                                              gap_debias_minimum_leading_period_hours=6,
                                                              gap_debias_minimum_trailing_period_hours=6,
                                                              )
            #(As a demonstration, we will fill the gaps of a single station. The following functions can also be
            # directly applied to the dataset.)
            your_station = your_dataset.get_station('vlinder05')


            #Get ERA5 modeldata at the location of your stations and period.
            ERA5_modeldata = your_station.get_modeldata(modelname='ERA5_hourly',
                                                        obstype='temp')

            #Use the debias method to fill the gaps
            gapfill_df = your_station.fill_gaps_automatic(modeldata=ERA5_modeldata,
                                                          max_interpolate_duration_str='6h', # <6 hours will be interpolated
                                                          obstype='temp')

        """
        #  ----------- Validate ----------------------------------------

        # check if modeldata is available
        if modeldata is None:
            logger.warning(
                "The dataset has no modeldate. Use the set_modeldata() function to add modeldata."
            )
            return None

        # check if obstype is present in eramodel
        assert (
            obstype in modeldata.df.columns
        ), f"{obstype} is not present in the modeldate: {modeldata}"

        # check if all station are present in eramodeldata
        # stations = self.gaps.to_df().index.unique().to_list()
        stations = list(set([gap.name for gap in self.gaps]))
        assert all(
            [sta in modeldata.df.index.get_level_values("name") for sta in stations]
        ), "Not all stations with gaps are in the modeldata!"

        if max_interpolate_duration_str is None:
            max_interpolate_duration_str = self.settings.gap["gaps_fill_settings"][
                "automatic"
            ]["max_interpolation_duration_str"]

        #  ------------select the method to apply gapfill per gap ----------
        interpolate_gaps = []
        debias_gaps = []

        for gap in self.gaps:
            if gap.duration <= pd.to_timedelta(max_interpolate_duration_str):
                interpolate_gaps.append(gap)
            else:
                debias_gaps.append(gap)

        # 1   ---------------Fill by interpolation ---------------------

        fill_settings_interp = self.settings.gap["gaps_fill_settings"]["linear"]

        apply_interpolate_gaps(
            gapslist=interpolate_gaps,
            obsdf=self.df,
            outliersdf=self.outliersdf,
            dataset_res=self.metadf["dataset_resolution"],
            gapfill_settings=self.settings.gap["gaps_fill_info"],
            obstype=obstype,
            method=fill_settings_interp["method"],
            max_consec_fill=fill_settings_interp["max_consec_fill"],
            overwrite_fill=overwrite_fill,
        )

        filldf_interp = make_gapfill_df(interpolate_gaps)

        # 2  --------------  Fill by debias -----------------------------

        fill_settings_debias = self.settings.gap["gaps_fill_settings"]["model_debias"]

        apply_debias_era5_gapfill(
            gapslist=debias_gaps,
            dataset=self,
            eraModelData=modeldata,
            obstype=obstype,
            debias_settings=fill_settings_debias,
            overwrite_fill=overwrite_fill,
        )

        # add label column
        filldf_debias = make_gapfill_df(debias_gaps)

        # combine both fill df's
        comb_df = concat_save([filldf_interp, filldf_debias])

        # update attr
        self.gapfilldf = comb_df

        return comb_df

    def fill_gaps_linear(self, obstype="temp", overwrite_fill=False):
        """Fill the gaps using linear interpolation.

        The gapsfilldf attribute of the Datasetinstance will be updated if
        the gaps are not filled yet or if overwrite_fill is set to True.

        Parameters
        ----------
        obstype : string, optional
            Fieldname to visualise. This can be an observation or station
            attribute. The default is 'temp'.
        overwrite_fill: bool, optional
            If a gap has already filled values, the interpolation of this gap
            is skipped if overwrite_fill is False. If set to True, the gapfill
            values and info will be overwitten. The default is False.

        Returns
        -------
        gapfilldf : pandas.DataFrame
            A dataframe containing all the filled records.

        Notes
        -----
        A schematic description of the linear gap fill:

        1. Iterate over all gaps.
        2. The gap is converted into a set of missing records (depending on the time resolution of the observations).
        3. Find a leading (the last observations before the gap) record and a trailing record (the last observation after the gap).
        4. By using the leading and trailing record an interpolation is applied to fill the missing records. A maximum consecutive fill threshold is applied, if exceeded the fill values are Nan's.
        5. The gap is updated with the interpolated values (metobs_toolkit.Gap.gapfill_df)

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

        """
        # TODO logging
        fill_settings = self.settings.gap["gaps_fill_settings"]["linear"]

        # fill gaps
        apply_interpolate_gaps(
            gapslist=self.gaps,
            obsdf=self.df,
            outliersdf=self.outliersdf,
            dataset_res=self.metadf["dataset_resolution"],
            gapfill_settings=self.settings.gap["gaps_fill_info"],
            obstype=obstype,
            method=fill_settings["method"],
            max_consec_fill=fill_settings["max_consec_fill"],
            overwrite_fill=overwrite_fill,
        )

        # get gapfilldf
        gapfilldf = make_gapfill_df(self.gaps)

        # update attr
        self.gapfilldf = gapfilldf

        return gapfilldf

    def fill_missing_obs_linear(self, obstype="temp"):
        """Interpolate missing observations.

        Fill in the missing observation rectords using interpolation. The
        missing_fill_df attribute of the Dataset will be updated.

        Parameters
        ----------
        obstype : string, optional
            Fieldname to visualise. This can be an observation or station
            attribute. The default is 'temp'.

        Returns
        -------
        None.

        Notes
        -----
        A schematic description of the linear fill of missing observations:

        1. Iterate over all missing observations.
        2. The missing observations are converted into a set of missing records (depending on the time resolution of the observations).
        3. Find a leading (the last observations before the missing observation) record and a trailing record (the last observation after the missing observation).
        4. By using the leading and trailing records, interpolation is applied to fill the missing records.
        5. The missing record is updated with the interpolated values (metobs_toolkit.Gap.gapfill_df).

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
            >>> # Fill the missing observations
            >>> dataset.fill_missing_obs_linear(obstype='temp')
            >>> dataset.missing_obs.get_info()
            -------- Missing observations info --------
            (Note: missing observations are defined on the frequency estimation of the native dataset.)
              * 1473 missing observations
              * For 28 stations
              * Missing observations are filled with interpolate for:
                temp:
                                                        temp
            name      datetime
            vlinder01 2022-09-08 08:00:00+00:00  18.630303
                      2022-09-07 23:00:00+00:00  17.512121
                      2022-09-08 00:00:00+00:00  17.636364
                      2022-09-08 02:00:00+00:00  17.884848
                      2022-09-08 03:00:00+00:00  18.009091
            ...

        """
        # TODO logging
        fill_settings = self.settings.missing_obs["missing_obs_fill_settings"]["linear"]
        fill_info = self.settings.missing_obs["missing_obs_fill_info"]

        # fill missing obs
        self.missing_obs.interpolate_missing(
            obsdf=self.df,
            resolutionseries=self.metadf["dataset_resolution"],
            obstype=obstype,
            method=fill_settings["method"],
        )

        missing_fill_df = self.missing_obs.fill_df

        missing_fill_df[obstype + "_" + fill_info["label_columnname"]] = fill_info[
            "label"
        ]["linear"]

        # Update attribute

        self.missing_fill_df = missing_fill_df

    def get_gaps_df(self):
        """
        List all gaps into an overview dataframe.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with stationnames as index, and the start, end and duretion
            of the gaps as columns.

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
            >>> dataset.get_gaps_df()
                                      start_gap                   end_gap        duration
            name
            vlinder05 2022-09-06 21:00:00+00:00 2022-09-13 06:00:00+00:00 6 days 09:00:00
            vlinder05 2022-09-13 20:00:00+00:00 2022-09-15 23:00:00+00:00 2 days 03:00:00


        """
        return gaps_to_df(self.gaps)

    def get_gaps_info(self):
        """Print out detailed information of the gaps.

        Returns
        -------
        None.

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
            >>> dataset.get_gaps_info()
            Gap for vlinder05 with:
            ---- Gap info -----
            (Note: gaps are defined on the frequency estimation of the native dataset.)
              * Start gap: 2022-09-06 21:00:00+00:00
              * End gap: 2022-09-13 06:00:00+00:00
              * Duration gap: 6 days 09:00:00
            ---- Gap fill info -----
            (No gapfill applied)
            Gap for vlinder05 with:
            ---- Gap info -----
            (Note: gaps are defined on the frequency estimation of the native dataset.)
              * Start gap: 2022-09-13 20:00:00+00:00
              * End gap: 2022-09-15 23:00:00+00:00
              * Duration gap: 2 days 03:00:00
            ---- Gap fill info -----
            (No gapfill applied)

        """
        if bool(self.gaps):
            # there are gaps
            for gap in self.gaps:
                gap.get_info()
        else:
            # no gaps
            print("There are no gaps.")

    def get_missing_obs_info(self):
        """Print out detailed information of the missing observations.

        Returns
        -------
        None.

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
            >>> dataset.get_missing_obs_info()
            -------- Missing observations info --------
            (Note: missing observations are defined on the frequency estimation of the native dataset.)
              * 1473 missing observations
              * For 28 stations
              * The missing observations are not filled.
            (More details on the missing observation can be found in the .series and .fill_df attributes.)

        """
        # empty obs protector in the .get_info method.
        self.missing_obs.get_info()
