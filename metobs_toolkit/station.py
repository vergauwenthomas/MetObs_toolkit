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
import metobs_toolkit.plot_collection as plotting

from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.settings_collection import label_def
from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
)
from metobs_toolkit.geedatasetmanagers import default_datasets as default_gee_datasets
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

    def __eq__(self, other):
        if not isinstance(other, Station):
            return False
        return (
            self.name == other.name
            and self.site == other.site
            and self.obsdata == other.obsdata
            and self._modeldata == other._modeldata
        )

    def __repr__(self):
        return f"Station instance of {self.name}"

    @property
    def name(self):
        return str(self._name)

    @property
    def site(self):
        return self._site

    @property
    def sensordata(self):
        return dict(self.obsdata)

    @property
    def df(self):
        # return dataframe with ['datetime', 'obstype'] as index and 'value' as single column.
        concatdf = save_concat(([sensor.df for sensor in self.obsdata.values()]))

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

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:

            combdf = pd.DataFrame(
                columns=["value", "label"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )
        return combdf

    @property
    def gapsdf(self) -> pd.DataFrame:

        concatlist = []
        for sensordata in self.obsdata.values():
            stadf = sensordata.gapsdf.reset_index()
            stadf["obstype"] = sensordata.obstype.name
            concatlist.append(stadf.set_index(["datetime", "obstype"]))

        combdf = save_concat(concatlist)
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )

        return combdf

    @property
    def metadf(self) -> pd.DataFrame:
        return self.site.metadf

    @property
    def modeldatadf(self) -> pd.DataFrame:
        concatlist = []
        for modeldata in self.modeldata.values():
            df = (
                modeldata.df.assign(obstype=modeldata.obstype.name)
                .assign(
                    details=f"{modeldata.modelname}:{modeldata.modelvariable} converted from {modeldata.obstype.model_unit} -> {modeldata.obstype.std_unit}"
                )
                .reset_index()
                .set_index(["datetime", "obstype"])
            )
            concatlist.append(df)
        combdf = save_concat(concatlist)
        combdf = combdf[["value", "details"]]
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "details"],
                index=pd.MultiIndex(
                    levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
                ),
            )
        return combdf

    @property
    def start_datetime(self):
        return min([sensdata.start_datetime for sensdata in self.obsdata.values()])

    @property
    def end_datetime(self):
        return max([sensdata.end_datetime for sensdata in self.obsdata.values()])

    @property
    def modeldata(self):
        return self._modeldata

    def add_to_modeldata(self, new_modeltimeseries, force_update=False):
        if not isinstance(new_modeltimeseries, ModelTimeSeries):
            raise MetObsWrongType(
                f"{new_modeltimeseries} is not an instance of ModelTimeSeries."
            )

        # Test if there is already modeldata for the same obstype available
        if (new_modeltimeseries.obstype.name in self.sensordata) & (not force_update):
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
        if bool(self.sensordata):
            for modeldataseries in self.sensordata.values():
                infostr += modeldataseries.get_info(printout=False)
        else:
            infostr += "  (no modeldata timeseries present)  \n"

        if printout:
            print(infostr)
        else:
            return infostr

    def sync_records(self, **kwargs):
        raise MetObsStationClassError(
            "sync_records() can only be applied on metobs.Dataset instances, not on metobs.Station instances."
        )

    def resample(
        self,
        target_freq,
        target_obstype: str | None = None,
        shift_tolerance=pd.Timedelta("4min"),
        origin=None,
        origin_simplify_tolerance=pd.Timedelta("4min"),
    ):
        # NOTE: if target_obstype is none, all sensors are resampled to the same freq!

        # format arguments
        target_freq = fmt_timedelta_arg(target_freq)
        shift_tolerance = fmt_timedelta_arg(shift_tolerance)

        if target_obstype is None:
            for sensor in self.obsdata.values():
                sensor.resample(
                    target_freq=target_freq,
                    shift_tolerance=shift_tolerance,
                    origin=origin,
                    origin_simplify_tolerance=origin_simplify_tolerance,
                )
        else:
            # check if target obstype is knonw
            self._obstype_is_known_check(target_obstype)

            # resample
            self.obsdata[target_obstype].resample(
                target_freq=target_freq,
                shift_tolerance=shift_tolerance,
                origin=origin,
                origin_simplify_tolerance=origin_simplify_tolerance,
            )

    def convert_outliers_to_gaps(self, all_observations=True, obstype="temp"):
        if all_observations:
            for sensor in self.obsdata.values():
                sensor.convert_outliers_to_gaps()

        else:
            self._obstype_is_known_check(obstype)
            self.obsdata[obstype].convert_outliers_to_gaps()

    def _rename(self, targetname):
        # Note: Not for users, one could accidentaly rename to another station in the dataset.
        # So --> only accecible as method in the dataset, that checks this possible error.

        # rename all
        self._name = str(targetname)
        self._site._stationname = str(targetname)
        for sensordat in self.obsdata.values():
            sensordat._rename(targetname)

    # ------------------------------------------
    #    Modeldata extraction
    # ------------------------------------------

    def get_gee_point_data(
        self, geestaticdataset, overwrite: bool = True, initialize_gee: bool = True
    ):
        if not isinstance(geestaticdataset, GEEStaticDatasetManager):
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

        if not isinstance(geestaticdataset, GEEStaticDatasetManager):
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
        if not isinstance(geedynamicdataset, GEEDynamicDatasetManager):
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
                    ].model_band,
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

    def get_qc_stats(self, target_obstype="temp", make_plot=True):
        # argument checks
        self._obstype_is_known_check(target_obstype)
        # get freq statistics
        qc_df = self.obsdata[target_obstype].get_qc_freq_statistics()

        if make_plot:
            plotdf = qc_df.reset_index().drop(columns=["name"]).set_index("qc_check")

            fig = plotting.qc_overview_pies(df=plotdf)
            fig.suptitle(
                f"QC frequency statistics of {target_obstype} on Station level: {self.stationname}."
            )
            return fig
        else:
            return qc_df

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
        linestyle="--",
    ) -> Axes:

        # test if the obstype has modeldata
        self._obstype_has_modeldata_check(obstype)

        # Create new axes if needed
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        plotdf = (
            self.modeldatadf.xs(obstype, level="obstype", drop_level=False)
            .assign(name=self.name)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        plotdf = plotdf[["value"]]
        plotdf["label"] = label_def["goodrecord"][
            "label"
        ]  # Just so that they are plotted as lines

        # Define linecolor (needed here if modeldata is added )
        if linecolor is None:
            colormap = plotting.create_categorical_color_map([self.name])
        else:
            colormap = {self.name: linecolor}
        ax = plotting.plot_timeseries_color_by_station(
            plotdf=plotdf,
            colormap=colormap,
            show_outliers=False,  # will not be used,
            show_gaps=False,  # will not be used
            ax=ax,
            linestyle=linestyle,
            legend_prefix=f"{self.modeldata[obstype].modelname}:{self.modeldata[obstype].modelvariable}@",
        )
        # Styling
        obstypeinstance = self.modeldata[obstype].obstype

        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{obstypeinstance.name} data for station {self.name}"
            )
        else:
            plotting.set_title(ax, title)

        # Set ylabel
        plotting.set_ylabel(ax, obstypeinstance._get_plot_y_label())

        # Set xlabel
        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax)

        # Add legend
        plotting.set_legend(ax)

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

        # Create new axes if needed
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        if show_modeldata:
            if linecolor is None:
                colormap = plotting.create_categorical_color_map([self.name])
            else:
                colormap = {self.name: linecolor}

            ax = self.make_plot_of_modeldata(
                obstype=obstype,
                linecolor=linecolor,
                ax=ax,
                figkwargs=figkwargs,
                title=title,
            )

        # Create plotdf
        plotdf = (
            self.df.xs(obstype, level="obstype", drop_level=False)
            .assign(name=self.name)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        if colorby == "station":
            # Define linecolor (needed here if modeldata is added )
            if linecolor is None:
                colormap = plotting.create_categorical_color_map([self.name])
            else:
                colormap = {self.name: linecolor}
            ax = plotting.plot_timeseries_color_by_station(
                plotdf=plotdf,
                colormap=colormap,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
                linestyle="-",
            )

        elif colorby == "label":
            ax = plotting.plot_timeseries_color_by_label(
                plotdf=plotdf,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
            )

        else:
            raise ValueError(
                f'colorby is either "station" or "label" but not {colorby}'
            )

        # Styling
        obstypeinstance = self.obsdata[obstype].obstype

        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{obstypeinstance.name} data for station {self.name}"
            )
        else:
            plotting.set_title(ax, title)

        # Set ylabel
        plotting.set_ylabel(ax, obstypeinstance._get_plot_y_label())

        # Set xlabel
        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax)

        # Add legend
        plotting.set_legend(ax)

        return ax

    # ------------------------------------------
    #    Gapfilling
    # ------------------------------------------
    def fill_gaps_with_raw_modeldata(self, target_obstype: str, overwrite_fill=False):

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.obsdata[target_obstype].fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries, method="raw", overwrite_fill=overwrite_fill
        )

    def fill_gaps_with_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration=pd.Timedelta("24h"),
        min_leading_records_total: int = 60,
        trailing_period_duration=pd.Timedelta("24h"),
        min_trailing_records_total: int = 60,
        overwrite_fill=False,
    ):
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.obsdata[target_obstype].fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries,
            method="debiased",
            method_kwargs={
                "leading_period_duration": leading_period_duration,
                "min_leading_records_total": min_leading_records_total,
                "trailing_period_duration": trailing_period_duration,
                "min_trailing_records_total": min_trailing_records_total,
            },
            overwrite_fill=overwrite_fill,
        )

    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration=pd.Timedelta("24h"),
        trailing_period_duration=pd.Timedelta("24h"),
        min_debias_sample_size: int = 6,
        overwrite_fill=False,
    ):
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.obsdata[target_obstype].fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries,
            method="diurnal_debiased",
            method_kwargs={
                "leading_period_duration": leading_period_duration,
                "trailing_period_duration": trailing_period_duration,
                "min_debias_sample_size": min_debias_sample_size,
            },
            overwrite_fill=overwrite_fill,
        )

    def fill_gaps_with_weighted_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration=pd.Timedelta("24h"),
        trailing_period_duration=pd.Timedelta("24h"),
        min_lead_debias_sample_size: int = 2,
        min_trail_debias_sample_size: int = 2,
        overwrite_fill=False,
    ):
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # Get modeltimeseries
        if target_obstype not in self.modeldata:
            raise MetObsModelDataError(
                f"No Modeldata found for {target_obstype} in {self}"
            )

        modeltimeseries = self.modeldata[target_obstype]

        # fill the gaps
        self.obsdata[target_obstype].fill_gap_with_modeldata(
            modeltimeseries=modeltimeseries,
            method="weighted_diurnal_debiased",
            method_kwargs={
                "leading_period_duration": leading_period_duration,
                "trailing_period_duration": trailing_period_duration,
                "min_lead_debias_sample_size": min_lead_debias_sample_size,
                "min_trail_debias_sample_size": min_trail_debias_sample_size,
            },
            overwrite_fill=overwrite_fill,
        )

    def interpolate_gaps(
        self,
        target_obstype: str,
        method: str = "time",
        max_consec_fill: int = 10,
        n_leading_anchors: int = 1,
        n_trailing_anchors: int = 1,
        max_lead_to_gap_distance: pd.Timedelta | None = None,
        max_trail_to_gap_distance: pd.Timedelta | None = None,
        overwrite_fill=False,
        method_kwargs={},
    ):
        # special formatters
        max_lead_to_gap_distance = fmt_timedelta_arg(max_lead_to_gap_distance)
        max_trail_to_gap_distance = fmt_timedelta_arg(max_trail_to_gap_distance)
        # obstype check
        self._obstype_is_known_check(obstype=target_obstype)

        # interpolate all the gaps
        self.obsdata[target_obstype].interpolate_gaps(
            overwrite_fill=overwrite_fill,
            method=method,
            max_consec_fill=max_consec_fill,
            n_leading_anchors=n_leading_anchors,
            n_trailing_anchors=n_trailing_anchors,
            max_lead_to_gap_distance=max_lead_to_gap_distance,
            max_trail_to_gap_distance=max_trail_to_gap_distance,
            method_kwargs=method_kwargs,
        )

    # ------------------------------------------
    #    Commons
    # ------------------------------------------

    def _obstype_is_known_check(self, obstype: str) -> None:
        if obstype not in self.obsdata.keys():
            raise MetObsSensorDataNotFound(
                f"{self} does not hold {obstype} sensordata. The present sensordata is: {list(self.obsdata.keys())}"
            )

    def _obstype_has_modeldata_check(self, obstype: str) -> None:
        if obstype not in self.sensordata.keys():
            raise MetObsObstypeNotFound(
                f"There is no {obstype} - modeldata present for {self}"
            )
