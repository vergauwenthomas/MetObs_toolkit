"""Collection of DF constructing functions on various levels
(sensordata, station, dataset) for overviews and summaries of QC checks."""

import pandas as pd
from typing import Union
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.df_helpers import save_concat


def sensordata_qc_overview_df(sensor) -> pd.DataFrame:
    """Build a QC overview DataFrame for a single SensorData object.

    Parameters
    ----------
    sensor : SensorData
        The sensor whose QC results are to be summarised.

    Returns
    -------
    pandas.DataFrame
        Wide DataFrame with a DatetimeIndex.  Columns follow the pattern
        ``('label', checkname)`` and ``('details', checkname)`` for each
        applied check, plus a ``'value'`` column with the raw observation
        values.  Missing checks are filled with ``'Not applied'``.
    """
    # TODO rearange the order of qc columns to reflect the executeion order
    possible_timestamps = sensor.series.index
    qc_before_timecoarsening = ["duplicated_timestamp"]

    to_concat = []
    for qcresult in sensor.outliers:
        checkdf = qcresult.create_outliersdf(
            map_to_basic_labels=False,  # get all flags (ok, outl, notchecked, unmet, saved)
            subset_to_outliers=False,
        )  # Get all flags
        # add checkname to the index
        checkdf["checkname"] = qcresult.checkname
        if qcresult.checkname in qc_before_timecoarsening:
            # Subset to coarsende timestmaps only
            checkdf = checkdf.reindex(possible_timestamps)

        checkdf.set_index("checkname", append=True, inplace=True)
        to_concat.append(checkdf)

    totaldf = save_concat(to_concat)

    if totaldf.empty:
        return pd.DataFrame(
            columns=["value", "label", "details"],
            index=pd.DatetimeIndex([], name="datetime"),
        )

    # Unstack
    totaldf = totaldf.unstack(level="checkname")
    totaldf.fillna("Not applied", inplace=True)

    # add values
    allvals = save_concat(
        [sensor.series, sensor.outliers_values_bin]
    )  # do not sort before removing the duplicates !
    allvals = allvals[~allvals.index.duplicated(keep="last")].sort_index()
    totaldf["value"] = allvals.loc[totaldf.index]

    return totaldf[["value", "label", "details"]]


def station_qc_overview_df(
    station, subset_obstypes: Union[list[str], None] = None
) -> pd.DataFrame:
    """Build a QC overview DataFrame for all sensors of a Station.

    Parameters
    ----------
    station : Station
        The station whose QC results are to be summarised.
    subset_obstypes : list of str or None, optional
        If given, only these observation types are included.  Unknown
        types are silently ignored.  Default is None (all sensors).

    Returns
    -------
    pandas.DataFrame
        Wide DataFrame indexed by ``(datetime, obstype)``.  Column
        structure mirrors :func:`sensordata_qc_overview_df`.
    """
    if subset_obstypes is None:
        sensortargets = station.sensordata.values()
    else:
        sensortargets = []
        for obstype in subset_obstypes:
            if obstype in station.sensordata:
                sensortargets.append(station.get_sensor(obstype))
            else:
                # Log a warning?
                pass

    to_concat = []
    for sensordata in sensortargets:
        stadf = sensordata_qc_overview_df(sensordata).reset_index()
        # add obstype to the index
        if not stadf.empty:
            stadf["obstype"] = sensordata.obstype.name
            stadf = stadf.reset_index().set_index(["datetime", "obstype"])
            to_concat.append(stadf)

    totaldf = save_concat(to_concat)
    totaldf.sort_index(inplace=True)

    if totaldf.empty:
        return pd.DataFrame(
            columns=["value", "label", "details"],
            index=pd.MultiIndex(
                levels=[[], []], codes=[[], []], names=["datetime", "obstype"]
            ),
        )

    return totaldf[["value", "label", "details"]]


def dataset_qc_overview_df(
    dataset,
    subset_stations: Union[list[str], None] = None,
    subset_obstypes: Union[list[str], None] = None,
) -> pd.DataFrame:
    """Build a QC overview DataFrame for all stations in a Dataset.

    Parameters
    ----------
    dataset : Dataset
        The dataset whose QC results are to be summarised.
    subset_stations : list of str or None, optional
        If given, only these station names are included.  Default is None
        (all stations).
    subset_obstypes : list of str or None, optional
        If given, only these observation types are included.  Default is
        None (all observation types).

    Returns
    -------
    pandas.DataFrame
        Wide DataFrame indexed by ``(datetime, obstype, name)``.  Column
        structure mirrors :func:`sensordata_qc_overview_df`.
    """
    if subset_stations is None:
        stationtargets = dataset.stations
    else:
        stationtargets = [
            dataset.get_station(station_name) for station_name in subset_stations
        ]

    to_concat = []
    for station in stationtargets:
        stadf = station_qc_overview_df(
            station, subset_obstypes=subset_obstypes
        ).reset_index()
        # add obstype to the index
        if not stadf.empty:
            stadf["name"] = station.name
            stadf = stadf.reset_index().set_index(["datetime", "obstype", "name"])
            to_concat.append(stadf)

    totaldf = save_concat(to_concat)
    totaldf.sort_index(inplace=True)

    if totaldf.empty:
        return pd.DataFrame(
            columns=["value", "label", "details"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["datetime", "obstype", "name"],
            ),
        )

    return totaldf[["value", "label", "details"]]
