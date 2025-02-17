import logging
import numpy as np
import pandas as pd
from typing import Literal
from matplotlib.pyplot import Axes

from metobs_toolkit.site import Site
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_datetime_arg,
    fmt_timedelta_arg,
)
import metobs_toolkit.backend_collection.timeseries_plotting as plotting

from metobs_toolkit.backend_collection.errorclasses import *


from metobs_toolkit.modeldata import GeeStaticDataset, GeeDynamicDataset
from metobs_toolkit.modeldata import default_datasets as default_gee_datasets
from metobs_toolkit.modeltimeseries import ModelTimeSeries


logger = logging.getLogger(__file__)


class Station:
    def __init__(self, stationname: str, site: Site, all_sensor_data: list):

        # dimension attributes
        self._name = str(stationname)
        self._site = site
        self.obsdata = self._set_stationdata(
            all_sensor_data
        )  # obstypename : SensorData

        # Extra extracted data
        self._modeldata = {}  # dict of ModelTimeSeries

    @property
    def name(self):
        return str(self._name)

    @property
    def site(self):
        return self._site

    @property
    def modeldata(self):
        return dict(self._modeldata)

    @property
    def df(self):
        # return dataframe with ['datetime', 'obstype'] as index and 'value' as single column.
        concatdf = pd.concat([sensor.df for sensor in self.obsdata.values()])

        # sort by datetime
        concatdf.sort_index(inplace=True)
        return concatdf

    @property
    def outliersdf(self):
        concatlist = []
        for sensordata in self.obsdata.values():
            stadf = sensordata.outliersdf[["value", "label"]].reset_index()
            stadf["obstype"] = sensordata.obstype.name
            concatlist.append(stadf.set_index(["datetime", "obstype"]))

        combdf = pd.concat(concatlist)
        combdf.sort_index(inplace=True)
        return combdf

    @property
    def gapsdf(self) -> pd.DataFrame:

        concatlist = []
        for sensordata in self.obsdata.values():
            stadf = sensordata.gapsdf[["value", "label"]].reset_index()
            stadf["obstype"] = sensordata.obstype.name
            concatlist.append(stadf.set_index(["datetime", "obstype"]))

        combdf = pd.concat(concatlist)
        combdf.sort_index(inplace=True)
        return combdf

    @property
    def metadf(self) -> pd.DataFrame:
        return self.site.metadf

    @property
    def start_datetime(self):
        return min([sensdata.start_datetime for sensdata in self.obsdata.values()])

    @property
    def end_datetime(self):
        return max([sensdata.end_datetime for sensdata in self.obsdata.values()])

    def add_to_modeldata(self, new_modeltimeseries, force_update=False):
        if not isinstance(new_modeltimeseries, ModelTimeSeries):
            raise MetObsWrongType(
                f"{new_modeltimeseries} is not an instance of ModelTimeSeries."
            )

        # Test if there is already modeldata for the same obstype available
        if (new_modeltimeseries.obstype.name in self.modeldata) & (not force_update):
            raise MetObsDataAlreadyPresent(
                f"There is already an modeltimeseries instance represinting {new_modeltimeseries.obstype.name}, and force_update is False."
            )

        self._modeldata.update({new_modeltimeseries.obstype.name: new_modeltimeseries})

    # ------------------------------------------
    #    Specials
    # ------------------------------------------
    def get_info(self, printout: bool = True) -> str | None:

        infostr = f"Station {self.name} with: \n"
        infostr += "\n --- Sensor Data ---\n"
        for sensordata in self.obsdata.values():
            infostr += sensordata.get_info(printout=False)

        infostr += "\n--- Site ---\n"
        infostr += self.site.get_info(printout=False)

        infostr += "\n--- Modeldata ---\n"
        if bool(self.modeldata):
            for modeldataseries in self.modeldata.values():
                infostr += modeldataseries.get_info(printout=False)
        else:
            infostr += "  (no modeldata timeseries present)  \n"

        if printout:
            print(infostr)
        else:
            return infostr

    def resample(
        self,
        target_freq,
        shift_tolerance=pd.Timedelta("4min"),
        origin=None,
        direction="nearest",
    ):

        target_freq = fmt_timedelta_arg(target_freq)
        shift_tolerance = fmt_timedelta_arg(shift_tolerance)

        for sensor in self.obsdata.values():
            sensor.resample(
                target_freq=target_freq,
                shift_tolerance=shift_tolerance,
                origin=origin,
                direction=direction,
            )

    # ------------------------------------------
    #    Modeldata extraction
    # ------------------------------------------

    def get_gee_point_data(
        self, geestaticdataset, overwrite: bool = True, initialize_gee: bool = True
    ):
        if not isinstance(geestaticdataset, GeeStaticDataset):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        value = self.site.get_gee_point_metadata(
            geestaticdataset=geestaticdataset, initialize_gee=initialize_gee
        )
        if overwrite:
            self.site.set_geedata(geestaticdataset.name, value)

        return value

    def get_lcz_from_gee(self, overwrite=True, initialize_gee=True):
        return self.get_gee_point_data(
            geestaticdataset=default_gee_datasets["lcz"],
            overwrite=overwrite,
            initialize_gee=initialize_gee,
        )

    def get_altitude_from_gee(self, overwrite=True, initialize_gee=True):
        return self.get_gee_point_data(
            geestaticdataset=default_gee_datasets["altitude"],
            overwrite=overwrite,
            initialize_gee=initialize_gee,
        )

    def get_gee_pointbuffer_data(
        self,
        geestaticdataset,
        buffers=[100],
        aggregate=False,
        overwrite: bool = True,
        initialize_gee: bool = True,
    ):

        if not isinstance(geestaticdataset, GeeStaticDataset):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        nesteddict = self.site.get_gee_point_buffer_fractions(
            geestaticdataset=geestaticdataset,
            buffers=buffers,
            aggregate=aggregate,
            initialize_gee=initialize_gee,
        )
        if overwrite:
            for bufferrad, fractions in nesteddict.items():
                self.site.set_gee_buffered_frac_data(buffer=bufferrad, data=fractions)
        return nesteddict

    def get_landcover_fractions(self, buffers=[100], aggregate=False, overwrite=True):
        self.get_gee_pointbuffer_data(
            geestaticdataset=default_gee_datasets["worldcover"],
            buffers=buffers,
            aggregate=aggregate,
            overwrite=overwrite,
        )

    def get_gee_timeseries_data(
        self,
        geedynamicdataset,
        startdt_utc=None,
        enddt_utc=None,
        target_obstypes=["temp"],
        get_all_bands=False,
        drive_filename=None,
        drive_folder="gee_timeseries_data",
        force_direct_transfer=False,
        force_to_drive=False,
    ):

        # Check geedynamic dataset
        if not isinstance(geedynamicdataset, GeeDynamicDataset):
            raise ValueError(
                f"geedynamicdataset should be an isntance of GeeDynamicDataset, not {type(geedynamicdataset)}"
            )

        # Format datetime arguments
        if startdt_utc is None:
            startdt_utc = self.start_datetime.tz_convert("UTC")
        else:
            startdt_utc = fmt_datetime_arg(startdt_utc, input_tz="UTC")

        if enddt_utc is None:
            enddt_utc = self.end_datetime.tz_convert("UTC")
        else:
            enddt_utc = fmt_datetime_arg(enddt_utc, input_tz="UTC")

        # check if target_obstypes are mapped to bands
        for obst in target_obstypes:
            if obst not in geedynamicdataset.modelobstypes.keys():
                raise MetObsMetadataNotFound(
                    f"{obst} is not a known modelobstype of {geedynamicdataset}."
                )

        # create specific name for the file that might be written to Drive
        if drive_filename is None:
            drive_filename = (
                f"{geedynamicdataset.name}_timeseries_data_of_{self.name}.csv"
            )

        df = geedynamicdataset.extract_timeseries_data(
            metadf=self.metadf,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            obstypes=target_obstypes,
            get_all_bands=get_all_bands,
            drive_filename=drive_filename,
            drive_folder=drive_folder,
            force_direct_transfer=force_direct_transfer,
            force_to_drive=force_to_drive,
        )
        if df is None:
            print("No data is returned by the GEE api request.")
            logger.warning("No data is returned by the GEE api request.")
            return

        # Create ModelTimeSeries instances
        for modelobscol in df.columns:
            if modelobscol in geedynamicdataset.modelobstypes.keys():
                modeltimeseries = ModelTimeSeries(
                    site=self.site,
                    datarecords=df[modelobscol].to_numpy(),
                    timestamps=df.index.get_level_values("datetime").to_numpy(),
                    obstype=geedynamicdataset.modelobstypes[modelobscol],
                    datadtype=np.float32,
                    timezone="UTC",
                    modelname=geedynamicdataset.name,
                    modelvariable=geedynamicdataset.modelobstypes[
                        modelobscol
                    ].get_modelband(),
                )
                # todo: duplicacy check
                self._modeldata.append(modeltimeseries)
            else:
                logger.info(
                    f"Skip {modelobscol} for creating a ModelTimeeries because of unknown obstype."
                )

        return df

    # ------------------------------------------
    #    QC checks (value based)
    # ------------------------------------------
    def gross_value_check(
        self,
        target_obstype: str = "temp",
        lower_threshold: float = -15.0,
        upper_threshold: float = 39.0,
    ):

        # argument validity checks
        self._obstype_is_known_check(target_obstype)

        self.obsdata[target_obstype].gross_value_check(
            lower_threshold=lower_threshold, upper_threshold=upper_threshold
        )

    def persistence_check(
        self,
        target_obstype: str = "temp",
        timewindow: str | pd.Timedelta = pd.Timedelta("60min"),
        min_records_per_window: int = 5,
    ):

        # argument checks
        self._obstype_is_known_check(target_obstype)
        timewindow = fmt_timedelta_arg(timewindow)

        # apply check on the sensordata
        self.obsdata[target_obstype].persistence_check(
            timewindow=timewindow, min_records_per_window=min_records_per_window
        )

    def repetitions_check(
        self, target_obstype: str = "temp", max_N_repetitions: int = 5
    ):

        # argument checks
        self._obstype_is_known_check(target_obstype)

        # apply check on the sensordata
        self.obsdata[target_obstype].repetitions_check(
            max_N_repetitions=max_N_repetitions
        )

    def step_check(
        self,
        target_obstype: str = "temp",
        max_increase_per_second: int | float = 8.0 / 3600.0,
        max_decrease_per_second: int | float = -10.0 / 3600.0,
    ):

        # argument checks
        self._obstype_is_known_check(target_obstype)

        # apply check on the sensordata
        self.obsdata[target_obstype].step_check(
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

    def window_variation_check(
        self,
        target_obstype: str = "temp",
        timewindow: pd.Timedelta = pd.Timedelta("1h"),
        min_records_per_window: int = 3,
        max_increase_per_second: int | float = 8.0 / 3600,
        max_decrease_per_second: int | float = -10.0 / 3600,
    ):

        # argument checks
        self._obstype_is_known_check(target_obstype)
        timewindow = fmt_timedelta_arg(timewindow)

        # apply check on the sensordata
        self.obsdata[target_obstype].window_variation_check(
            timewindow=timewindow,
            min_records_per_window=min_records_per_window,
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

    # ------------------------------------------
    #    Getters
    # ------------------------------------------

    # ------------------------------------------
    #    Setters
    # ------------------------------------------

    def _set_stationdata(self, all_sensor_data: list) -> dict:
        return {
            sensor_data.obstype.name: sensor_data for sensor_data in all_sensor_data
        }

    def _set_tz(self, current_tz):
        pass

    # ------------------------------------------
    #    Plotting
    # ------------------------------------------

    def make_plot_of_modeldata(
        self,
        obstype: str = "temp",
        linecolor=None,
        ax=None,
        figkwargs: dict = {},
        title: str | None = None,
    ) -> Axes:

        # test if the obstype has modeldata
        self._obstype_has_modeldata_check(obstype)

        ax = self.modeldata[obstype].make_plot(
            linecolor=linecolor, ax=ax, figkwargs=figkwargs, title=title
        )
        return ax

    def make_plot(
        self,
        obstype: str = "temp",
        colorby: Literal["station", "label"] = "label",
        show_modeldata=False,
        linecolor=None,
        show_outliers=True,
        show_gaps=True,
        ax=None,
        figkwargs: dict = {},
        title: str | None = None,
    ) -> Axes:

        # test if obstype have sensordata
        self._obstype_is_known_check(obstype)

        # Define linecolor (needed here if modeldata is added )
        if linecolor is None:
            linecolor = plotting.create_station_color_map(["dummy"])["dummy"]

        # first the metadata (so that styling attributes are overwritten by the sensordata)
        if show_modeldata:
            ax = self.make_plot_of_modeldata(
                ax=ax,
                obstype=obstype,
                linecolor=linecolor,
                figkwargs=figkwargs,
                title=None,
            )

        # add the records to the axes
        ax = self.obsdata[obstype].make_plot(
            colorby=colorby,
            linecolor=linecolor,
            show_outliers=show_outliers,
            show_gaps=show_gaps,
            ax=ax,
            figkwargs=figkwargs,
            title=title,
        )
        return ax

    # ------------------------------------------
    #    Commons
    # ------------------------------------------

    def _obstype_is_known_check(self, obstype: str) -> None:
        if obstype not in self.obsdata.keys():
            raise MetObsSensorDataNotFound(
                f"{self} does not hold {obstype} sensordata. The present sensordata is: {list(self.obsdata.keys())}"
            )

    def _obstype_has_modeldata_check(self, obstype: str) -> None:
        if obstype not in self.modeldata.keys():
            raise MetObsObstypeNotFound(
                f"There is no {obstype} - modeldata present for {self}"
            )
