import logging
from typing import Union
import pandas as pd
import numpy as np

from matplotlib.pyplot import Axes

from metobs_toolkit.backend_collection.df_helpers import to_timedelta
import metobs_toolkit.backend_collection.printing_collection as printing
from metobs_toolkit.obstypes import Obstype
import metobs_toolkit.plot_collection.timeseries_plotting as plotting


logger = logging.getLogger("<metobs_toolkit>")

class ModelTimeSeries:
    """Class for model-based timeseries at one location.

    Parameters
    ----------
    site : object
        The site object representing the location.
    datarecords : np.ndarray
        Array of data records.
    timestamps : np.ndarray
        Array of timestamps corresponding to the data records.
    obstype : Obstype
        The observation type.
    datadtype : type, optional
        Data type for the records, by default np.float32.
    timezone : str, optional
        Timezone for the timestamps, by default "UTC".
    modelname : str, optional
        Name of the model, by default None.
    modelvariable : str, optional
        Name of the model variable, by default None.
    """

    def __init__(
        self,
        site,
        datarecords: np.ndarray,
        timestamps: np.ndarray,
        obstype: Obstype,
        datadtype: type = np.float32,
        timezone: str = "UTC",
        modelname: str = None,
        modelvariable: str = None,
    ):
        self.site = site
        self.obstype = obstype

        # Data
        data = pd.Series(
            data=pd.to_numeric(datarecords, errors="coerce").astype(datadtype),
            index=pd.DatetimeIndex(data=timestamps, tz=timezone, name="datetime"),
            name=obstype.name,
        )
        self.series = data

        # model metadata
        self.modelname = modelname
        self.modelvariable = modelvariable

    # ------------------------------------------
    #    Specials
    # ------------------------------------------
    def __eq__(self, other) -> bool:
        """Check equality with another ModelTimeSeries object."""
        if not isinstance(other, ModelTimeSeries):
            return False
        return (
            self.site == other.site
            and self.obstype == other.obstype
            and self.series.equals(other.series)
            and self.modelname == other.modelname
            and self.modelvariable == other.modelvariable
        )

    @property
    def df(self) -> pd.DataFrame:
        """Return all records as a DataFrame."""
        # get all records
        df = (
            self.series.to_frame()
            .rename(columns={self.obstype.name: "value", self.stationname: "value"})
            .assign(model=self.modelname)
        )
        return df

    @property
    def stationname(self) -> str:
        """Return the name of the station this SensorData belongs to."""
        return self.site.stationname

    @property
    def tz(self) -> str:
        """Return the timezone of the stored timestamps."""
        return self.series.index.tz

    @property
    def start_datetime(self) -> pd.Timestamp:
        """Return the start datetime of the series."""
        return self.series.index.min()

    @property
    def end_datetime(self) -> pd.Timestamp:
        """Return the end datetime of the series."""
        return self.series.index.max()

    @property
    def freq(self) -> pd.Timedelta:
        """Return the frequency of the series."""
        freq = pd.infer_freq(self.series.index)
        if freq is None:
            raise ValueError("Frequency could not be computed.")
        # note: sometimes 'h' is returned, and this gives issues, so add a 1 in front
        return to_timedelta(freq)

    def _get_info_core(self, nident_root=1) -> dict:
        infostr = ""
        infostr += printing.print_fmt_line(f"Origin {self.modelname} -> variable/band: {self.modelvariable}", nident_root)
        infostr += printing.print_fmt_line(f"From {self.start_datetime} --> {self.end_datetime}", nident_root)
        infostr += printing.print_fmt_line(f"Assumed frequency: {self.freq}", nident_root)
        infostr += printing.print_fmt_line(f"Number of records: {self.series.shape[0]}", nident_root)
        infostr += printing.print_fmt_line(f"Units are converted from {self.obstype.model_unit} --> {self.obstype.std_unit}", nident_root)

        return infostr
    def get_info(self, printout: bool = True) -> Union[None, str]:
        """
        Print or return information about the ModelTimeSeries.

        Parameters
        ----------
        printout : bool, optional
            If True, print the information. If False, return as string. Default is True.

        Returns
        -------
        None or str
            None if printout is True, otherwise the information string.
        """
        logger.debug(f"{self.__class__.__name__}.get_info called for {self}")
        infostr = ""
        infostr += printing.print_fmt_title('General info of ModelTimeSeries')
        infostr += printing.print_fmt_line(f"{self.obstype.name} model data at location of {self.stationname}")
        infostr += self._get_info_core(nident_root=1)
        if printout:
            print(infostr)
        else:
            return infostr

    def make_plot(
        self,
        linecolor: str = None,
        ax: Union[Axes, None] = None,
        figkwargs: dict = {},
        title: Union[str, None] = None,
    ) -> Axes:
        """
        Create a plot of the model time series.

        Parameters
        ----------
        linecolor : str, optional
            Color of the line, by default None.
        ax : Axes, optional
            Matplotlib Axes to plot on, by default None.
        figkwargs : dict, optional
            Additional keyword arguments for figure creation, by default {}.
        title : str or None, optional
            Title for the plot, by default None.

        Returns
        -------
        Axes
            The matplotlib Axes with the plot.
        """
        logger.debug(f"{self.__class__.__name__}.make_plot called for {self}")
        # define figure
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        # Define a color
        if linecolor is None:
            # create a new color
            color = plotting.create_categorical_color_map(["dummy"])["dummy"]
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
