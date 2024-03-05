#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 11:46:29 2024

@author: thoverga
"""

import logging
import pandas as pd


from metobs_toolkit import Dataset
from metobs_toolkit.df_helpers import xs_save

logger = logging.getLogger(__name__)


class Dataset(Dataset):
    """Extension on the metobs_toolkit.Dataset class with gapfill techniques"""

    def get_gaps_fill_df(self):
        """
        Create a dataframe with information on how all gaps are filled.

        Returns
        -------
        pandas.DataFrame
            A DataFrame with ['name', 'datetime', 'obstype'] as index, ['fill',
            'fill_method', 'msg'] as columns. Nan and None values are used
            for records that have not been gapfilled.

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
            >>> dataset.coarsen_time_resolution(freq='1H')
            >>>
            >>> # Apply quality control on the temperature observations
            >>> dataset.apply_quality_control(obstype='temp') #Using the default QC settings
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'radiation_temp', 'pressure', 'pressure_at_sea_level', 'precip', 'precip_sum', 'wind_speed', 'wind_gust', 'wind_direction'] observation types
                 *10080 observation records
                 *235 records labeled as outliers
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
        gaps_infodf_list = []
        for gap in self.gaps:
            gapfilldf = pd.concat([gap.gapdf, gap.anchordf])
            # Create a triple index
            gapfilldf["obstype"] = gap.obstype.name
            gapfilldf = gapfilldf.rename(columns={f"{gap.obstype.name}_fill": "fill"})
            gapfilldf = gapfilldf.drop(columns=[gap.obstype.name])
            gapfilldf = gapfilldf.reset_index().set_index(
                ["name", "datetime", "obstype"]
            )
            # add a gapid column
            gapfilldf["gap_ID"] = f"{gap.name};{gap.startdt};{gap.enddt}"

            gaps_infodf_list.append(gapfilldf)

        return pd.concat(gaps_infodf_list).sort_index()

    # def fill_gaps_automatic(
    #     self,
    #     modeldata,
    #     obstype="temp",
    #     max_interpolate_duration_str=None,
    #     overwrite_fill=False,
    # ):
    #     """Fill the gaps by using linear interpolation or debiased modeldata.

    #     The method that is applied to perform the gapfill will be determined by
    #     the duration of the gap.

    #     When the duration of a gap is smaller or equal than
    #     max_interpolation_duration, the linear interpolation method is applied
    #     else the debiased modeldata method.

    #     Parameters
    #     ----------
    #     modeldata : metobs_toolkit.Modeldata
    #         The modeldata to use for the gapfill. This model data should the required
    #         timeseries to fill all gaps present in the dataset.
    #     obstype : String, optional
    #         Name of the observationtype you want to apply gap filling on. The
    #         modeldata must contain this observation type as well. The
    #         default is 'temp'.
    #     max_interpolate_duration_str : Timedelta or str, optional
    #         Maximum duration to apply interpolation for gapfill when using the
    #         automatic gapfill method. Gaps with longer durations will be filled
    #         using debiased modeldata. The default is None.
    #     overwrite_fill: bool, optional
    #         If a gap has already filled values, the interpolation of this gap
    #         is skipped if overwrite_fill is False. If set to True, the gapfill
    #         values and info will be overwitten. The default is False.

    #     Returns
    #     -------
    #     comb_df : TYPE
    #         gapfilldf : pandas.DataFrame
    #             A dataframe containing all the filled records.

    #     Examples
    #     --------
    #     .. code-block:: python

    #         import metobs_toolkit

    #         your_dataset = metobs_toolkit.Dataset()
    #         your_dataset.update_settings(
    #             input_data_file=metobs_toolkit.demo_datafile, # path to the data file
    #             input_metadata_file=metobs_toolkit.demo_metadatafile,
    #             template_file=metobs_toolkit.demo_template,
    #         )
    #         # Specify the gap defenition
    #         your_dataset.update_qc_settings(gapsize_in_records = 20)

    #         #Update the gapsize BEFORE importing the data
    #         your_dataset.import_data_from_file()

    #         #Update the settings (definition of the period to calculate biases for)
    #         your_dataset.update_gap_and_missing_fill_settings(
    #                                                           gap_debias_prefered_leading_period_hours=24,
    #                                                           gap_debias_prefered_trailing_period_hours=24,
    #                                                           gap_debias_minimum_leading_period_hours=6,
    #                                                           gap_debias_minimum_trailing_period_hours=6,
    #                                                           )
    #         #(As a demonstration, we will fill the gaps of a single station. The following functions can also be
    #         # directly applied to the dataset.)
    #         your_station = your_dataset.get_station('vlinder05')

    #         #Get ERA5 modeldata at the location of your stations and period.
    #         ERA5_modeldata = your_station.get_modeldata(modelname='ERA5_hourly',
    #                                                     obstype='temp')

    #         #Use the debias method to fill the gaps
    #         gapfill_df = your_station.fill_gaps_automatic(modeldata=ERA5_modeldata,
    #                                                       max_interpolate_duration_str='6H', # <6 hours will be interpolated
    #                                                       obstype='temp')

    #     """
    #     #  ----------- Validate ----------------------------------------

    #     # check if modeldata is available
    #     if modeldata is None:
    #         logger.warning(
    #             "The dataset has no modeldate. Use the set_modeldata() function to add modeldata."
    #         )
    #         return None

    #     # check if obstype is present in eramodel
    #     assert (
    #         obstype in modeldata.df.columns
    #     ), f"{obstype} is not present in the modeldate: {modeldata}"

    #     # check if all station are present in eramodeldata
    #     # stations = self.gaps.to_df().index.unique().to_list()
    #     stations = list(set([gap.name for gap in self.gaps]))
    #     assert all(
    #         [sta in modeldata.df.index.get_level_values("name") for sta in stations]
    #     ), "Not all stations with gaps are in the modeldata!"

    #     if max_interpolate_duration_str is None:
    #         max_interpolate_duration_str = self.settings.gap["gaps_fill_settings"][
    #             "automatic"
    #         ]["max_interpolation_duration_str"]

    #     #  ------------select the method to apply gapfill per gap ----------
    #     interpolate_gaps = []
    #     debias_gaps = []

    #     for gap in self.gaps:
    #         if gap.duration <= pd.to_timedelta(max_interpolate_duration_str):
    #             interpolate_gaps.append(gap)
    #         else:
    #             debias_gaps.append(gap)

    #     # 1   ---------------Fill by interpolation ---------------------

    #     fill_settings_interp = self.settings.gap["gaps_fill_settings"]["linear"]

    #     apply_interpolate_gaps(
    #         gapslist=interpolate_gaps,
    #         obsdf=self.df,
    #         outliersdf=self.outliersdf,
    #         dataset_res=self.metadf["dataset_resolution"],
    #         gapfill_settings=self.settings.gap["gaps_fill_info"],
    #         obstype=obstype,
    #         method=fill_settings_interp["method"],
    #         max_consec_fill=fill_settings_interp["max_consec_fill"],
    #         overwrite_fill=overwrite_fill,
    #     )

    #     filldf_interp = make_gapfill_df(interpolate_gaps)

    #     # 2  --------------  Fill by debias -----------------------------

    #     fill_settings_debias = self.settings.gap["gaps_fill_settings"]["model_debias"]

    #     apply_debias_era5_gapfill(
    #         gapslist=debias_gaps,
    #         dataset=self,
    #         eraModelData=modeldata,
    #         obstype=obstype,
    #         debias_settings=fill_settings_debias,
    #         overwrite_fill=overwrite_fill,
    #     )

    #     # add label column
    #     filldf_debias = make_gapfill_df(debias_gaps)

    #     # combine both fill df's
    #     comb_df = concat_save([filldf_interp, filldf_debias])

    #     # update attr
    #     # self.gapfilldf = comb_df

    #     return comb_df

    def interpolate_gaps(
        self,
        obstype="temp",
        overwrite_fill=False,
        method="time",
        max_consec_fill=10,
        max_lead_to_gap_distance=None,
        max_trail_to_gap_distance=None,
    ):
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
            >>> dataset.coarsen_time_resolution(freq='1H')
            >>>
            >>> # Apply quality control on the temperature observations
            >>> dataset.apply_quality_control(obstype='temp') #Using the default QC settings
            >>>
            >>> # Interpret the outliers as missing/gaps
            >>> dataset.update_gaps_and_missing_from_outliers(obstype='temp')
            >>> dataset
            Dataset instance containing:
                 *28 stations
                 *['temp', 'humidity', 'radiation_temp', 'pressure', 'pressure_at_sea_level', 'precip', 'precip_sum', 'wind_speed', 'wind_gust', 'wind_direction'] observation types
                 *10080 observation records
                 *235 records labeled as outliers
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
        for gap in self.gaps:
            # filter the gaps to those of target obstype
            if gap.obstype.name == str(obstype):
                # check if gap is filled
                status = gap._get_gapfill_status()
                if status == "Gap does not exist in observation space":
                    logger.info(f"{gap} cannot be filled because it is not a real one")
                    continue
                elif status != "Unfilled gap":
                    # already filled values
                    if not overwrite_fill:
                        logger.warning(
                            f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                        )
                        print(
                            f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite_fill}."
                        )
                        continue
                else:
                    logger.debug(f"filling {gap} with {method} interpolation.")
                    gap.interpolate_gap(
                        Dataset=self,
                        method=method,
                        max_consec_fill=max_consec_fill,
                        max_lead_to_gap_distance=max_lead_to_gap_distance,
                        max_trail_to_gap_distance=max_trail_to_gap_distance,
                    )

        infodf = self.get_gaps_fill_df()
        return xs_save(infodf, obstype, "obstype").sort_index()

    def fill_gaps_with_raw_modeldata(
        self, Modeldata, obstype="temp", overwrite_fill=False
    ):

        # Get all gaps that are candidates for the filling
        to_fill_list = [gap for gap in self.gaps if gap.obstype.name == obstype]

        # Return None if there are no candidates
        if len(to_fill_list) == 0:
            logger.warning("There are no gaps to fill.")
            print("There are no gaps to fill.")
            return

        for gap in to_fill_list:
            fillbool = _can_gap_be_filled(gap, obstype, overwrite_fill, Modeldata)
            if fillbool:
                # fill the gap
                logger.debug(f"filling {gap} with raw modeldata: {Modeldata} ")
                gap.raw_model_gapfill(Dataset=self, Modeldata=Modeldata)
            else:
                # gap could not be filled
                continue

        infodf = self.get_gaps_fill_df()
        if obstype is None:
            return infodf.sort_index()
        else:
            return xs_save(infodf, obstype, "obstype").sort_index()

    def fill_gaps_with_debias_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="24H",
        min_leading_records_total=5,
        trailing_period_duration="24H",
        min_trailing_records_total=5,
    ):

        # Get all gaps that are candidates for the filling
        to_fill_list = [gap for gap in self.gaps if gap.obstype.name == obstype]

        # Return None if there are no candidates
        if len(to_fill_list) == 0:
            logger.warning("There are no gaps to fill.")
            print("There are no gaps to fill.")
            return

        for gap in to_fill_list:
            fillbool = _can_gap_be_filled(gap, obstype, overwrite_fill, Modeldata)
            if fillbool:
                # fill the gap
                logger.debug(f"filling {gap} with debiased modeldata: {Modeldata} ")
                gap.debias_model_gapfill(
                    Dataset=self,
                    Modeldata=Modeldata,
                    leading_period_duration=leading_period_duration,
                    min_leading_records_total=min_leading_records_total,
                    trailing_period_duration=trailing_period_duration,
                    min_trailing_records_total=min_trailing_records_total,
                )
            else:
                # gap could not be filled
                continue

        infodf = self.get_gaps_fill_df()
        if obstype is None:
            return infodf.sort_index()
        else:
            return xs_save(infodf, obstype, "obstype").sort_index()

    def fill_gaps_with_diurnal_debias_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="48H",
        min_debias_sample_size=3,
        trailing_period_duration="48H",
    ):

        # Get all gaps that are candidates for the filling
        to_fill_list = [gap for gap in self.gaps if gap.obstype.name == obstype]

        # Return None if there are no candidates
        if len(to_fill_list) == 0:
            logger.warning("There are no gaps to fill.")
            print("There are no gaps to fill.")
            return

        for gap in to_fill_list:
            fillbool = _can_gap_be_filled(gap, obstype, overwrite_fill, Modeldata)
            if fillbool:
                # fill the gap
                logger.debug(
                    f"filling {gap} with diurnal debiased modeldata: {Modeldata} "
                )
                gap.diurnal_debias_model_gapfill(
                    Dataset=self,
                    Modeldata=Modeldata,
                    leading_period_duration=leading_period_duration,
                    min_debias_sample_size=min_debias_sample_size,
                    trailing_period_duration=trailing_period_duration,
                )
            else:
                # gap could not be filled
                continue

        infodf = self.get_gaps_fill_df()
        if obstype is None:
            return infodf.sort_index()
        else:
            return xs_save(infodf, obstype, "obstype").sort_index()

    def fill_gaps_with_weighted_diurnal_debias_modeldata(
        self,
        Modeldata,
        obstype="temp",
        overwrite_fill=False,
        leading_period_duration="48H",
        min_lead_debias_sample_size=2,
        trailing_period_duration="48H",
        min_trail_debias_sample_size=2,
    ):

        # Get all gaps that are candidates for the filling
        to_fill_list = [gap for gap in self.gaps if gap.obstype.name == obstype]

        # Return None if there are no candidates
        if len(to_fill_list) == 0:
            logger.warning("There are no gaps to fill.")
            print("There are no gaps to fill.")
            return

        for gap in to_fill_list:
            fillbool = _can_gap_be_filled(gap, obstype, overwrite_fill, Modeldata)
            if fillbool:
                # fill the gap
                logger.debug(
                    f"filling {gap} with weighted diurnal debiased modeldata: {Modeldata} "
                )
                gap.weighted_diurnal_debias_model_gapfill(
                    Dataset=self,
                    Modeldata=Modeldata,
                    leading_period_duration=leading_period_duration,
                    min_lead_debias_sample_size=min_lead_debias_sample_size,
                    trailing_period_duration=trailing_period_duration,
                    min_trail_debias_sample_size=min_trail_debias_sample_size,
                )
            else:
                # gap could not be filled
                continue

        infodf = self.get_gaps_fill_df()
        if obstype is None:
            return infodf.sort_index()
        else:
            return xs_save(infodf, obstype, "obstype").sort_index()


def _can_gap_be_filled(gap, obstypename, overwrite, Modeldata):
    status = gap._get_gapfill_status()
    if status == "Gap does not exist in observation space":
        logger.info(f"{gap} cannot be filled because it is not a real one")
        return False
    if status != "Unfilled gap":
        # already filled values
        if not overwrite:
            logger.warning(
                f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite}."
            )
            print(
                f"{gap} cannot be filled because it already contains filled values, and overwrite fill is {overwrite}."
            )
            return False
    if not (gap.obstype.name in Modeldata.df.columns):
        logger.warning(
            f"{gap} cannot be filled because {gap.obstype.name} has no values in {Modeldata}."
        )
        return False

    return True
