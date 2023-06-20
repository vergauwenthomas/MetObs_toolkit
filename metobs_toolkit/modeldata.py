#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:50:17 2023

@author: thoverga
"""

import pandas as pd
import logging

from metobs_toolkit.df_helpers import init_multiindexdf, conv_tz_multiidxdf, xs_save, multiindexdf_datetime_subsetting

from metobs_toolkit.landcover_functions import connect_to_gee, gee_extract_timeseries

from metobs_toolkit.plotting_functions import model_timeseries_plot, timeseries_plot

from metobs_toolkit.convertors import convert_to_toolkit_units

from metobs_toolkit.settings import Settings

logger = logging.getLogger(__name__)

# =============================================================================
# Class Model data (collection of external model data)
# =============================================================================


class Modeldata:
    def __init__(self, modelname):
        self.df = init_multiindexdf
        self.modelname = modelname

        self._settings = Settings()
        self.mapinfo = self._settings.gee["gee_dataset_info"]

    def __repr__(self):
        return f"ModelData instance: {self.modelname} model data of {list(self.df.columns)}"

    def __str__(self):
        return self.__repr__()


    def _conv_to_timezone(self, tzstr):
        # get tzstr by datetimindex.tz.zone


        df=self.df
        df['datetime_utc'] = df.index.get_level_values('datetime').tz_convert(tzstr)
        df = df.reset_index()
        df = df.drop(columns=['datetime'])
        df = df.rename(columns={'datetime_utc': 'datetime'})
        df = df.set_index(['name', 'datetime'])
        self.df = df

    def get_ERA5_data(self, metadf, startdt, enddt, obstype="temp"):
        # startdt and enddt IN UTC FORMAT!!!!!

        # Subset metadf to stations with coordinates
        no_coord_meta = metadf[metadf[['lat','lon']].isna().any(axis=1)]
        if not no_coord_meta.empty:
            print(f'WARNING. Following stations do not have coordinates, and thus no modeldata extraction is possible: {no_coord_meta.index.to_list()}')
            metadf = metadf[~metadf[['lat','lon']].isna().any(axis=1)]


        era_mapinfo = self.mapinfo["ERA5_hourly"]
        # Connect to Gee
        connect_to_gee()

        # Get data using GEE
        df = gee_extract_timeseries(
            metadf=metadf,
            mapinfo=era_mapinfo,
            startdt=startdt,
            enddt=enddt,
            obstype=obstype,
            latcolname="lat",
            loncolname="lon",
        )
        if not df.empty:
            # Convert to toolkit units
            df[obstype], _tlk_unit = convert_to_toolkit_units(
                data=df[obstype], data_unit=era_mapinfo["band_of_use"][obstype]["units"]
            )

        self.df = df
        self.modelname = "ERA5_hourly"

    def set_model_from_csv(
        self,
        csvpath,
        modelname="ERA5_hourly",
        convert_units=True,
        obstype="temp",
        datatimezone="UTC",
    ):
        df = pd.read_csv(csvpath, sep=",")
        # format datetime
        df["datetime"] = pd.to_datetime(df["datetime"], format="%Y%m%d%H%M%S")
        df["datetime"] = df["datetime"].dt.tz_localize(datatimezone)

        # format index
        df = df.set_index(["name", "datetime"])
        df = df.sort_index()

        # rename to values to toolkit space
        bandname = self.mapinfo[modelname]["band_of_use"][obstype]["name"]
        df = df.rename(columns={bandname: obstype})

        # convert units
        if convert_units:
            df[obstype], _tlk_unit = convert_to_toolkit_units(
                data=df[obstype],
                data_unit=self.mapinfo[modelname]["band_of_use"][obstype]["units"],
            )
        df = df[obstype].to_frame()
        self.df = df
        self.modelname = modelname

    def interpolate_modeldata(self, to_multiidx, obstype="temp"):
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
                print("Extrapolation")
            if max(sta_recordsdf.index.get_level_values("datetime")) > max(
                sta_moddf.index.get_level_values("datetime")
            ):
                print("Extrapolation")

            # combine model and records
            mergedf = sta_recordsdf.merge(
                sta_moddf, how="outer", left_index=True, right_index=True
            )

            # reset index for time interpolation
            mergedf = mergedf.reset_index().set_index("datetime").sort_index()

            # interpolate missing modeldata
            mergedf[obstype].interpolate(
                method="time", limit_area="inside", inplace=True
            )
            # convert back to multiindex
            mergedf = mergedf.reset_index().set_index(["name", "datetime"]).sort_index()
            # filter only records
            mergedf = mergedf.loc[sta_recordsdf.index]

            returndf = pd.concat([returndf, mergedf])
        return returndf

    def make_plot(self, obstype="temp", dataset = None, stationnames=None,
        starttime=None, endtime=None, title=None, show_outliers=True,
        show_filled=True, legend=True,
        _ax=None, #needed for GUI, not recommended use
        ):

        """
        This function creates a timeseries plot for the Modeldata. When a
        metobs_toolkit.Dataset is provided, it is plotted in the same figure.

        The line colors represent the timesries for different locations.



        Parameters
        ----------
        obstype : string, optional
             Fieldname to visualise. This can be an observation or station
             attribute. The default is 'temp'.
        dataset : metobs_toolkit.Dataset, optional
            A Dataset instance with observations plotted in the same figure.
            Observations are represented by solid line and modeldata by dashed
            lines. The default is None.
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


        logger.info(f"Make {obstype}-timeseries plot of model data")

        # Basic test
        if obstype not in self.df.columns:
            print(f'ERROR: {obstype} is not foud in the modeldata df.')
            return
        if self.df.empty:
            print('ERROR: The modeldata is empty.')
            return
        if (not dataset is None):
            if (obstype not in dataset.df.columns):
                print(f'ERROR: {obstype} is not foud in the Dataframe df.')
                return


        model_df = self.df

        # ------ filter model ------------

        # Filter on obstype
        model_df = model_df[[obstype]]

        # Subset on stationnames
        if not stationnames is None:
            model_df = model_df[model_df.index.get_level_values('name').isin(stationnames)]

        # Subset on start and endtime
        model_df = multiindexdf_datetime_subsetting(model_df, starttime, endtime)


        #  -------- Filter dataset (if available) -----------
        if not dataset is None:
            # combine all dataframes
            mergedf = dataset.combine_all_to_obsspace()

            # subset to obstype
            mergedf = xs_save(mergedf, obstype, level='obstype')

            # Subset on stationnames
            if not stationnames is None:
                mergedf = mergedf[mergedf.index.get_level_values('name').isin(stationnames)]

            # Subset on start and endtime
            mergedf = multiindexdf_datetime_subsetting(mergedf, starttime, endtime)


        # Generate ylabel
        try:
            model_true_field_name = self.mapinfo[self.modelname]['band_of_use'][obstype]['name']
        except KeyError:
            print (f'No model field name found for {obstype} in {self}.')
            model_true_field_name = 'Unknown fieldname'
        y_label = f'{model_true_field_name}'

        if not dataset is None:
            dataset_obs_orig_name = dataset.data_template[obstype]['orig_name']
            units = dataset.data_template[obstype]['units']

            y_label = f'{y_label} \n {dataset_obs_orig_name} ({units})'


        #Generate title
        title = f'{self.modelname} : {model_true_field_name}'
        if not dataset is None:
            title = f'{title} and {dataset_obs_orig_name} observations.'


        # make plot of the observations
        if not dataset is None:
            # make plot of the observations
            ax, col_map = timeseries_plot(mergedf=mergedf,
                                 title=title,
                                 ylabel=y_label,
                                 colorby='name',
                                 show_legend=legend,
                                 show_outliers=show_outliers,
                                 show_filled=show_filled ,
                                 settings = dataset.settings)


            # Make plot of the model on the previous axes
            ax, col_map = model_timeseries_plot(
                                    df=model_df,
                                    obstype=obstype,
                                    title=title,
                                    ylabel=y_label,
                                    settings = self._settings,
                                    show_primary_legend=False,
                                    add_second_legend=True,
                                    _ax = ax,
                                    colorby_name_colordict=col_map)


        else:

            # Make plot of model on empty axes
            ax, _colmap = model_timeseries_plot(
                    df=model_df,
                    obstype=obstype,
                    title=title,
                    ylabel=y_label,
                    settings = self._settings,
                    show_primary_legend=legend,
                    add_second_legend=False,
                    _ax = None
                    )

        return ax