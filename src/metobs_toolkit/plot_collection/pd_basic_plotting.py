import numpy as np
from matplotlib.pyplot import Axes
from metobs_toolkit.backend_collection.decorators import log_entry


@log_entry
def modeldata_simple_pd_plot(modeltimeseries, **pdplotkwargs) -> Axes:
    """
    A wrapper on the pandas.series.plot function for ModelTimeSeries.

    Parameters
    ----------
    **pdplotkwargs
        Additional keyword arguments passed to pandas.Series.plot().

    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.

    """
    # set name (so that labels match the id)
    plot_series = modeltimeseries.series
    plot_series.name = modeltimeseries._id()
    return plot_series.plot(**pdplotkwargs)


@log_entry
def sensordata_simple_pd_plot(
    sensordata, show_labels: list = ["ok"], **pdplotkwargs
) -> Axes:
    """
    A wrapper on the pandas.series.plot function for SensorData.

    Parameters
    ----------
    show_labels : list of str, optional
        List of quality control labels to include in the plot. Records with
        other labels are converted to NaN to avoid interpolation artifacts.
        Default is ['ok'].
    **pdplotkwargs
        Additional keyword arguments passed to pandas.Series.plot().

    Returns
    -------
    matplotlib.axes.Axes
        The axes object containing the plot.

    Notes
    -----
    * Tip: You can inspect all the present labels (and their counts) by using: `sensordata.df['label'].value_counts()`.

    * The plot excludes records that don't match the specified labels by
      converting them to NaN rather than subsetting, which prevents unwanted
      interpolation between valid data points.

    """

    df = sensordata.df
    # filter by label (do not subset, but convert others to Nan to avoid interpolation)
    df.loc[~df["label"].isin(show_labels), "value"] = np.nan

    # convert to pandas series
    df.index = df.index.droplevel("obstype")  # drop obstype level
    series = df["value"]
    series.name = sensordata._id()
    # make plot
    axs = series.plot(**pdplotkwargs)

    return axs
