import numpy as np
from matplotlib.pyplot import Axes
from metobs_toolkit.backend_collection.loggingmodule import log_entry

@log_entry
def sensordata_simple_pd_plot(sensordata, show_labels: list=['ok'], **pdplotkwargs) -> Axes:
    """
    Create a pandas plot of the sensor data with optional label filtering.

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
    The plot excludes records that don't match the specified labels by 
    converting them to NaN rather than subsetting, which prevents unwanted 
    interpolation between valid data points.
    
    Notes
    -------
    Tip: You can inspect all the present labels (and their counts) by using: `sensordata.df['label'].value_counts()`.
    
    """
    
    df = sensordata.df
    #filter by label (do not subset, but convert others to Nan to avoid interpolation)
    df.loc[~df['label'].isin(show_labels), 'value'] = np.nan
    
    #convert to pandas series
    df.index = df.index.droplevel('obstype') #drop obstype level
    series = df['value']
    series.name = f'{sensordata.stationname}:{sensordata.obstype.name}'
    #make plot
    axs = series.plot(**pdplotkwargs)
    
    return axs