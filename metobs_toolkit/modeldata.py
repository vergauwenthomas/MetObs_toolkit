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
    
    def make_plot(
        self,
        dataset,
        stationnames=None,
        obstype="temp",
        starttime=None,
        endtime=None,
        title=None,
        colorby='label',
        legend=True,
        _ax=None, #needed for GUI, not recommended use
        ):
        """
        This function creates a timeseries plot for the model data. The variable observation type
        is plotted for all stationnames from a starttime to an endtime.

        All styling attributes are extracted from the Settings.

        Parameters
        ----------

        stationnames : list, optional
            A list with stationnames to include in the timeseries. If None is given, all the stations are used, defaults to None.
        obstype : string, optional
             Fieldname to visualise. This can be an observation or station
             attribute. The default is 'temp'.
        starttime : datetime.datetime, optional
             Specifiy the start datetime for the plot. If None is given it will use the start datetime of the dataset, defaults to None.
        endtime : datetime.datetime, optional
             Specifiy the end datetime for the plot. If None is given it will use the end datetime of the dataset, defaults to None.
        title : string, optional
             Title of the figure, if None a default title is generated. The default is None.
        y_label : string, optional
             y-axes label of the figure, if None a default label is generated. The default is None.
        legend : bool, optional
             I True, a legend is added to the plot. The default is True.


        Returns
        -------
        axis : matplotlib.pyplot.axes
             The timeseries axes of the plot is returned.

        """
        if stationnames is None:
            logger.info(f"Make {obstype}-timeseries plot of model data for all stations")
        else:
            logger.info(f"Make {obstype}-timeseries plot of model data for {stationnames}")
        
        model_df = self.df
        
        # subset model to stationnames in dataset/station
        if not dataset is None:
            model_df = model_df.reset_index()
            
            
            if dataset._istype == 'Dataset':
                # dataset object
                model_df = model_df.loc[model_df['name'].isin(
                    dataset.df.index.get_level_values('name').unique())]
            else:
                # station object
                model_df = model_df.loc[model_df['name'] == dataset.name]
                
            model_df = model_df.set_index(['name', 'datetime'])
        
        # Subset on stationnames
        if not stationnames is None:
            model_df = model_df.reset_index()
            model_df = model_df.loc[model_df['name'].isin(stationnames)]
            model_df = model_df.set_index(['name', 'datetime'])
            dataset.df = dataset.df.reset_index()
            dataset.df = dataset.df.loc[dataset.df['name'].isin(stationnames)]
            dataset.df = dataset.df.set_index(['name', 'datetime'])


        # ylabel tips
        try:
            model_true_field = self.mapinfo[self.modelname]['band_of_use'][obstype]
        except KeyError:
            print (f'No model data available for {obstype}.')
            
        try:
            obs_true_field = dataset.data_template[obstype]
        except KeyError:
            print (f'No observation data available for {obstype}.')
        
        y_label = f'{model_true_field["name"]} ({model_true_field["units"]}) \n {obs_true_field["orig_name"]} ({obs_true_field["units"]})'
       
        # Subset on start and endtime
        model_df = multiindexdf_datetime_subsetting(model_df, starttime, endtime)

        # fig, ax = plt.subplots()
        stationlist = sorted(list(model_df.index.get_level_values("name").unique()))
        for i in range(5, len(stationlist), 5):
            stationlist.insert(i, "\n")
        
        title = f'{model_true_field["name"]}/{obs_true_field["orig_name"]} for [{", ".join(stationlist)}]'
  
        # check if df is empty
        if model_df.empty:
            print(f'No model data available')
       
        ax = dataset.make_plot(stationnames=stationnames, starttime=starttime, endtime=endtime, title=title, y_label=y_label, legend=False, _ax=None)
       
        # Make plot
        ax = model_timeseries_plot(
            df=model_df,
            obstype=obstype,
            title=title,
            ylabel=y_label,
            show_legend=legend,
            settings = dataset.settings,
            _ax = ax
        )

        return ax