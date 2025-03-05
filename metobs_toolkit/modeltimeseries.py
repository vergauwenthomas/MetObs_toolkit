import logging
import pandas as pd
import numpy as np


from matplotlib.pyplot import Axes

from metobs_toolkit.obstypes import Obstype
import metobs_toolkit.plot_collection.timeseries_plotting as plotting


class ModelTimeSeries:
    """Class equivalent of SensorData for model-based timeseries at one location."""

    def __init__(
        self,
        site,
        datarecords: np.ndarray,
        timestamps: np.ndarray,
        obstype: Obstype,
        datadtype=np.float32,
        timezone="UTC",
        modelname=None,
        modelvariable=None,
    ):

        self.site = site
        self.obstype = obstype

        # Data
        data = pd.Series(
            data=pd.to_numeric(datarecords, errors="coerce").astype(datadtype),
            index=pd.DatetimeIndex(data=timestamps, tz=timezone),
            name=obstype.name,
        )
        self.series = data

        # model metadata
        self.modelname = modelname
        self.modelvariable = modelvariable

    # ------------------------------------------
    #    Specials
    # ------------------------------------------
    @property
    def df(self):
        # get all records
        df = (
            self.series.to_frame()
            .rename(columns={self.obstype.name: "value", self.stationname: "value"})
            .assign(model=self.modelname)
        )
        return df

    @property
    def stationname(self) -> str:
        """
        Return the name of the station this SensorData belongs to.

        Returns
        -------
        str
            station name
        """
        return self.site.stationname

    @property
    def tz(self):
        """
        Return the (up-to-date) timezone of the stored timestamps.

        Returns
        -------
        str
            Timezone of the stored timestamps.
        """

        return self.series.index.tz

    @property
    def start_datetime(self) -> pd.Timestamp:
        """
        Return the start datetime of the series.

        Returns
        -------
        pd.Timestamp
            Start datetime of the series.
        """
        return self.series.index.min()

    @property
    def end_datetime(self) -> pd.Timestamp:
        """
        Return the end datetime of the series.

        Returns
        -------
        pd.Timestamp
            End datetime of the series.
        """
        return self.series.index.max()

    @property
    def freq(self) -> pd.Timedelta:
        """
        Return the frequency of the series.

        Returns
        -------
        pd.Timedelta
            Frequency of the series.
        """
        freq = pd.infer_freq(self.series.index)
        if freq is None:
            raise ValueError("Frequency could not be computed.")
        # note: sometimes 'h' is returned, and this gives issues, so add a 1 in front
        if not freq[0].isdigit():
            freq = "1" + freq
        return pd.Timedelta(freq)

    def get_info(self, printout: bool = True):
        infostr = ""
        infostr += f"{self.obstype.name} modeldata at location of {self.stationname}:\n"
        infostr += (
            f"  * origin {self.modelname} -> variable/band: {self.modelvariable}\n"
        )
        infostr += f"  * from {self.start_datetime} --> {self.end_datetime}\n"
        infostr += f"  * assumed frequency: {self.freq}\n"
        infostr += f"  * Number of records: {self.series.shape[0]}\n"
        if printout:
            print(infostr)
        else:
            return infostr

    def make_plot(
        self, linecolor=None, ax=None, figkwargs: dict = {}, title: str | None = None
    ) -> Axes:

        # define figure
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        # Define a color
        if linecolor is None:
            # create a new color
            color = plotting.create_station_color_map(["dummy"])["dummy"]
        else:
            color = linecolor

        legendname = f"{self.modelname}:{self.modelvariable}@{self.stationname}"

        ax = plotting.add_lines_to_axes(
            ax=ax,
            series=self.series,
            legend_label=legendname,
            linestyle="--",
            color=color,
        )

        # Add Styling attributes
        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{self.obstype.name} data for station {self.stationname}"
            )
        else:
            plotting.set_title(ax, title)
        # Set ylabel
        plotting.set_ylabel(ax, self.obstype._get_plot_y_label())

        # Set xlabel
        plotting.set_xlabel(ax, f"Timestamps (in {self.tz})")

        # Add legend
        plotting.set_legend(ax)

        return ax
