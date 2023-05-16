#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 13:50:17 2023

@author: thoverga
"""

import pandas as pd

from metobs_toolkit.df_helpers import init_multiindexdf, conv_tz_multiidxdf

from metobs_toolkit.landcover_functions import connect_to_gee, gee_extract_timeseries

from metobs_toolkit.convertors import convert_to_toolkit_units

from metobs_toolkit.settings import Settings

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
            sta_recordsdf = recordsdf.xs(sta, level="name", drop_level=False)
            sta_moddf = self.df.xs(sta, level="name", drop_level=False)

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
