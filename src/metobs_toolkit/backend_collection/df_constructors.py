""" This module is a collection of functions to convert and join
    Metobs classes to pandas DataFrames."""

import pandas as pd

# Rearranged imports: standard libraries, dependencies, local modules
from metobs_toolkit.settings_collection import label_def
from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.loggingmodule import log_entry

# ------------------------------------------
#    Sensor
# ------------------------------------------

@log_entry
def sensor_construct_df(sensordata):
    """Return a DataFrame of the sensor records."""
   
    # get all records
    df = (
        sensordata.series.to_frame()
        .rename(columns={sensordata.obstype.name: "value", sensordata.stationname: "value"})
        .assign(label=label_def["goodrecord"]["label"])
    )

    outliersdf = sensordata.outliersdf[["value", "label"]]

    gapsdf = sensordata.gapsdf[["value", "label"]]

    # concat all together (do not change order)
    to_concat = [df]
    if not outliersdf.empty:
        to_concat.append(outliersdf)
    if not gapsdf.empty:
        to_concat.append(gapsdf)
    df = save_concat((to_concat))
    # remove duplicates
    df = df[~df.index.duplicated(keep="last")].sort_index()

    # add 'obstype' as index
    df = (
        df.assign(obstype=sensordata.obstype.name)
        .reset_index()
        .set_index(["datetime", "obstype"])
    )
    return df

@log_entry
def sensor_construct_outliersdf(sensordata):
    """Return a DataFrame of the outlier records."""
        
    to_concat = []
    for outlierinfo in sensordata.outliers:
        checkname = outlierinfo["checkname"]
        checkdf = outlierinfo["df"]
        checkdf["label"] = label_def[checkname]["label"]
        to_concat.append(checkdf)

    totaldf = save_concat(to_concat)

    if totaldf.empty:
        # return empty dataframe
        totaldf = pd.DataFrame(
            columns=["value", "label"], index=pd.DatetimeIndex([], name="datetime")
        )
    else:
        totaldf.sort_index(inplace=True)

    return totaldf

@log_entry  
def sensor_construct_gapsdf(sensordata):
    """Return a DataFrame of the gap records."""
    to_concat = []
    if bool(sensordata.gaps):
        for gap in sensordata.gaps:
            to_concat.append(gap.df)
        return save_concat((to_concat)).sort_index()
    else:
        return pd.DataFrame(
            columns=["value", "label", "details"],
            index=pd.DatetimeIndex([], name="datetime"),
        )


# ------------------------------------------
#    Station
# ------------------------------------------
@log_entry
def station_construct_df(station):
    """
    Construct a DataFrame representation of the observations.

    Returns
    -------
    pd.DataFrame
        A pandas DataFrame with a single column 'value'.
    """

    # return dataframe with ['datetime', 'obstype'] as index and 'value' as single column.
    concatdf = save_concat(([sensor.df for sensor in station.sensordata.values()]))

    # sort by datetime
    concatdf.sort_index(inplace=True)
    return concatdf

@log_entry
def station_construct_outliersdf(station):
    """
    Construct a DataFrame representation of all the outliers.

    Outliers are the observations that are flagged by the performed quality control.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with two columns ['value', 'label'], representing
        the value and details of the flagged observation.
    """

    concatlist = []
    for sensordata in station.sensordata.values():
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

@log_entry
def station_construct_gapsdf(station):
    """
    Construct a DataFrame representation of all the gaps.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns ['value', 'label', 'details'], representing
        the value, the gap label, and details of the gap record.
    """
    concatlist = []
    for sensordata in station.sensordata.values():
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

@log_entry
def station_construct_modeldatadf(station):
    """
    Construct a DataFrame representation of all the present model data.

    Model data is stored as `ModelTimeSeries` instances, and is set as an attribute of
    a `Station`.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with columns ['value', 'details'], representing
        the value, and details of the corresponding model data.
    """
    concatlist = []
    for modeldata in station.modeldata.values():
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
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["value", "details"],
            index=pd.MultiIndex(
                levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
            ),
        )
    # formatting
    combdf = combdf[["value", "details"]]
    combdf.sort_index(inplace=True)
    return combdf

@log_entry
def station_construct_metadf(station):
    """
    Construct a DataFrame representation of metadata.

    Metadata is the information related to the sensors, that does not change over time.
    The metadata is extracted from the site instance.

    Returns
    -------
    pandas.DataFrame
        A DataFrame with the station names as index, and the metadata as columns.
    """

    return station.site.metadf



# ------------------------------------------
#    Dataset
# ------------------------------------------
@copy_doc(station_construct_df)
@log_entry
def dataset_construct_df(dataset):
    concatlist = []
    for sta in dataset.stations:
        stadf = sta.df.reset_index()
        if stadf.empty:
            continue
        stadf["name"] = sta.name
        concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

    combdf = save_concat((concatlist))
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["value", "label"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["datetime", "obstype", "name"],
            ),
        )
    return combdf

@copy_doc(station_construct_outliersdf)
@log_entry
def dataset_construct_outliersdf(dataset):
    concatlist = []
    for sta in dataset.stations:
        stadf = sta.outliersdf.reset_index()
        stadf["name"] = sta.name
        concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

    combdf = save_concat((concatlist))
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["value", "label"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["datetime", "obstype", "name"],
            ),
        )
    return combdf

@copy_doc(station_construct_gapsdf)
@log_entry
def dataset_construct_gapsdf(dataset):
    concatlist = []
    for sta in dataset.stations:
        stadf = sta.gapsdf.reset_index()
        if stadf.empty:
            continue
        stadf["name"] = sta.name
        concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

    combdf = save_concat((concatlist))
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["value", "label", "details"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["datetime", "obstype", "name"],
            ),
        )
    return combdf

@copy_doc(station_construct_modeldatadf)
@log_entry
def dataset_construct_modeldatadf(dataset):
    concatlist = []
    for sta in dataset.stations:
        stadf = sta.modeldatadf.reset_index()
        if stadf.empty:
            continue
        stadf["name"] = sta.name
        concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

    combdf = save_concat((concatlist))
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["value", "details"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["datetime", "obstype", "name"],
            ),
        )
    return combdf

@copy_doc(station_construct_metadf)
@log_entry
def dataset_construct_metadf(dataset):

    concatlist = []
    for sta in dataset.stations:
        concatlist.append(sta.metadf)
    return save_concat((concatlist)).sort_index()

