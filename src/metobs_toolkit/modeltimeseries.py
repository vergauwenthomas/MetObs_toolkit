import logging
from typing import Union
import pandas as pd
import numpy as np

from matplotlib.pyplot import Axes

from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.df_helpers import (
    to_timedelta,
    convert_to_numeric_series,
)
from metobs_toolkit.xrconversions import modeltimeseries_to_xr
import metobs_toolkit.backend_collection.printing_collection as printing
from metobs_toolkit.obstypes import ModelObstype
from metobs_toolkit.plot_collection.general_functions import (
    create_axes,
    set_legend,
    set_xlabel,
    set_ylabel,
    set_title,
    create_categorical_color_map,
)
from metobs_toolkit.plot_collection.timeseries_plotting import add_lines_to_axes

from metobs_toolkit.backend_collection.errorclasses import MetObsUnitUnknown
from metobs_toolkit.backend_collection.loggingmodule import log_entry
from metobs_toolkit.backend_collection.dataframe_constructors import modeltimeseries_df

logger = logging.getLogger("<metobs_toolkit>")


class ModelTimeSeries:
    """Class for model-based timeseries at one location.

    This class stores model data at a specific location and automatically converts
    the data from the model's native units to standard units as defined by the
    observation type. The unit conversion is performed during initialization.

    Parameters
    ----------
    site : object
        The site object representing the location.
    datarecords : np.ndarray
        Array of data records in the model's native units. These will be
        automatically converted to standard units during initialization.
    timestamps : np.ndarray
        Array of timestamps corresponding to the data records.
    modelobstype : ModelObstype
        The observation type containing unit information. Must have a model_unit
        attribute set for unit conversion to work properly.
    datadtype : type, optional
        Data type for the records, by default np.float32.
    timezone : str, optional
        Timezone for the timestamps, by default "UTC".
    modelname : str, optional
        Name of the model, by default None.
    modelvariable : str, optional
        Name of the model variable, by default None.

    Notes
    -----
    The stored data in the `series` attribute will be in standard units as defined
    by `modelobstype.std_unit`, not in the original model units. The conversion is
    performed automatically.

    Raises
    ------
    MetObsUnitUnknown
        If the modelobstype does not have a model_unit set, which is required for
        unit conversion.
    """

    def __init__(
        self,
        site,
        datarecords: np.ndarray,
        timestamps: np.ndarray,
        modelobstype: ModelObstype,
        datadtype: type = np.float32,
        timezone: str = "UTC",
        modelname: str = None,
        modelvariable: str = None,
        _convert_to_standard_units: bool = True,
    ):
        self.site = site

        # Testing the ModelObstype
        self.modelobstype = modelobstype
        if not isinstance(self.modelobstype, ModelObstype):
            raise TypeError(f"Expected ModelObstype, got {type(self.modelobstype)}")
        if self.modelobstype.model_unit is None:
            raise MetObsUnitUnknown(
                f"The model_unit of {self.modelobstype} is not set. Set this using the ModelObstype.model_unit = value syntax."
            )

        # Data
        data = pd.Series(
            data=convert_to_numeric_series(datarecords, datadtype=datadtype).values,
            index=pd.DatetimeIndex(data=timestamps, tz=timezone, name="datetime"),
            name=modelobstype.name,
        )
        if _convert_to_standard_units:
            data = self.modelobstype.convert_to_standard_units(
                input_data=data, input_unit=self.modelobstype.model_unit
            )
        else:
            pass

        self.series = data

        # model metadata
        self.modelname = modelname
        # TODO: is modelvariable not band_unit from the obstype?
        self.modelvariable = modelvariable

    # ------------------------------------------
    #    Specials
    # ------------------------------------------
    def __repr__(self):
        return f"<ModelTimeSeries(id={self._id()},obstype={self.modelobstype.name})>"

    def _id(self) -> str:
        """A physical unique id.

        In the __add__ methods, if the id of two instances differs, adding is
        a regular concatenation.
        """
        return f"{self.site._id()}{self.modelname}_{self.modelvariable}"

    def __eq__(self, other) -> bool:
        """Check equality with another ModelTimeSeries object."""
        if not isinstance(other, ModelTimeSeries):
            return False
        return (
            self.site == other.site
            and self.modelobstype == other.modelobstype
            and self.series.equals(other.series)
            and self.modelname == other.modelname
            and self.modelvariable == other.modelvariable
        )

    def __add__(self, other: "ModelTimeSeries") -> "ModelTimeSeries":
        """
        Combine two ModelTimeSeries objects for the same site and obstype.
        The result contains all unique records, with preference to non-NaN values from 'other'.
        """
        if not isinstance(other, ModelTimeSeries):
            raise TypeError("Can only add ModelTimeSeries to ModelTimeSeries.")
        if self._id() != other._id():
            raise ValueError(
                f"Cannot add ModelTimeSeries for different ID's ({self._id()} != {other._id()})."
            )

        # Align timezones if different
        if self.tz != other.tz:
            other_series = other.series.tz_convert(str(self.tz))
        else:
            other_series = other.series

        # Combine the series, preferring non-NaN from 'other'
        combined_series = self.series.combine_first(other_series)
        combined = ModelTimeSeries(
            site=self.site,
            datarecords=combined_series.values,
            timestamps=combined_series.index.values,
            modelobstype=self.modelobstype + other.modelobstype,
            datadtype=combined_series.dtype,
            timezone=self.tz,
            modelname=self.modelname,
            modelvariable=self.modelvariable,
            _convert_to_standard_units=False,  # !! units are already converted !!
        )
        return combined

    @copy_doc(modeltimeseries_df)
    @property
    def df(self) -> pd.DataFrame:
        return modeltimeseries_df(self)

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

    @copy_doc(modeltimeseries_to_xr)
    @log_entry
    def to_xr(self) -> "xarray.Dataset":
        return modeltimeseries_to_xr(self, fmt_datetime_coordinate=True)

    def _get_info_core(self, nident_root=1) -> dict:
        infostr = ""
        infostr += printing.print_fmt_line(
            f"Modelname {self.modelname} -> variable/band: {self.modelvariable}",
            nident_root,
        )
        infostr += printing.print_fmt_line(
            f"From {self.start_datetime} --> {self.end_datetime}", nident_root
        )
        infostr += printing.print_fmt_line(
            f"Assumed frequency: {self.freq}", nident_root
        )
        infostr += printing.print_fmt_line(
            f"Number of records: {self.series.shape[0]}", nident_root
        )
        infostr += printing.print_fmt_line(
            f"Units are converted from {self.modelobstype.model_unit} --> {self.modelobstype.std_unit}",
            nident_root,
        )

        return infostr

    @log_entry
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
        infostr += printing.print_fmt_title("General info of ModelTimeSeries")
        infostr += printing.print_fmt_line(
            f"{self.modelobstype.name} model data at location of {self.stationname}"
        )
        infostr += self._get_info_core(nident_root=1)
        if printout:
            print(infostr)
        else:
            return infostr

    @log_entry
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
            ax = create_axes(**figkwargs)

        # Define a color
        if linecolor is None:
            # create a new color
            color = create_categorical_color_map(["dummy"])["dummy"]
        else:
            color = linecolor

        legendname = f"{self.modelname}:{self.modelvariable}@{self.stationname}"

        ax = add_lines_to_axes(
            ax=ax,
            series=self.series,
            legend_label=legendname,
            linestyle="--",
            color=color,
        )

        # Add Styling attributes
        # Set title:
        if title is None:
            set_title(
                ax, f"{self.modelobstype.name} data for station {self.stationname}"
            )
        else:
            set_title(ax, title)
        # Set ylabel
        set_ylabel(ax, self.modelobstype._get_plot_y_label())

        # Set xlabel
        set_xlabel(ax, f"Timestamps (in {self.tz})")

        # Add legend
        set_legend(ax)

        return ax
