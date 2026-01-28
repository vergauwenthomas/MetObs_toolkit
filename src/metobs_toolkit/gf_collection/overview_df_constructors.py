"""Collection of DF constructing functions on various levels
(sensordata, station, dataset) for overviews and summaries of Gaps."""

import pandas as pd
from typing import Union
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.df_helpers import save_concat

#===============================
# Gap overiview
#===============================

def sensordata_gap_status_overview_df(sensordata) -> pd.DataFrame:
    """
    Create gap status overview DataFrame with one row per gap period.

    Parameters
    ----------
    sensordata : SensorData
        SensorData instance containing gap information.

    Returns
    -------
    pandas.DataFrame
        DataFrame with gap periods indexed by gap start time. Contains columns:

        * gapend : pandas.Timestamp - End time of the gap
        * gapsize : pandas.Timedelta - Duration of the gap
        * label : str - Gap fill status (e.g., 'not filled', 'interpolated')
        * details : str - Gap creation details and methods used

    Notes
    -----
    Unlike gapsdf which lists all missing records, this provides one summary
    row per continuous gap period.
    """

    gap_info_list = []

    if bool(sensordata.gaps):
        for gap in sensordata.gaps:
            gap_df = gap.df

            # Basic gap information using Gap object properties
            gap_start = gap.start_datetime
            gap_end = gap.end_datetime
            gap_size = gap.end_datetime - gap.start_datetime
            gap_label = gap.fillstatus

            # Handle details
            unique_details = gap_df["details"].unique()
            if len(unique_details) == 1:
                gap_details = f"unidetail gap: {unique_details[0]}"
            else:
                gap_details = (
                    f'multi_details gap: {" -- ".join(sorted(unique_details))}'
                )

            # Create gap info dictionary
            gap_info = {
                "gapstart": gap_start,
                "gapend": gap_end,
                "gapsize": gap_size,
                "label": gap_label,
                "details": gap_details,
            }

            gap_info_list.append(gap_info)

        # Create DataFrame from gap info list
        if gap_info_list:
            result_df = (
                pd.DataFrame(gap_info_list)
                .reset_index(drop=True)
                .set_index("gapstart")
                .sort_index()
            )

            return result_df
    else:
        # No gaps present
        return pd.DataFrame(
            columns=["gapend", "gapsize", "label", "details"],
            index=pd.Index([], name="gapstart"),
        )

@copy_doc(sensordata_gap_status_overview_df)
def station_gap_status_overview_df(station) -> pd.DataFrame:
    concatlist = []
    for sensordata in station.sensordata.values():
        stadf = sensordata_gap_status_overview_df(sensordata).reset_index()
        if not stadf.empty:
            stadf["obstype"] = sensordata.obstype.name
            stadf = stadf.set_index(["gapstart", "obstype"])
            concatlist.append(stadf)

    combdf = save_concat(concatlist)
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["gapend", "gapsize", "label", "details"],
            index=pd.MultiIndex(
                levels=[[], []], codes=[[], []], names=["gapstart", "obstype"]
            ),
        )

    return combdf


@copy_doc(sensordata_gap_status_overview_df)
def dataset_gap_status_overview_df(dataset) -> pd.DataFrame:

    concatlist = []
    for sta in dataset.stations:
        stadf = station_gap_status_overview_df(sta).reset_index()
        if stadf.empty:
            continue
        stadf["name"] = sta.name
        concatlist.append(stadf.set_index(["gapstart", "obstype", "name"]))

    combdf = save_concat((concatlist))
    combdf.sort_index(inplace=True)
    if combdf.empty:
        combdf = pd.DataFrame(
            columns=["gapend", "gapsize", "label", "details"],
            index=pd.MultiIndex(
                levels=[[], [], []],
                codes=[[], [], []],
                names=["gapstart", "obstype", "name"],
            ),
        )
    return combdf



#===============================
# QC overiew
#===============================

def sensordata_qc_overview_df(sensor) -> pd.DataFrame:
    #TODO: docstring
    
    #TODO rearange the order of qc columns to reflect the executeion order
    possible_timestamps = sensor.series.index
    qc_before_timecoarsening = ['duplicated_timestamp']
    
    
    to_concat = []
    for qcresult in sensor.outliers:
        checkdf = qcresult.create_outliersdf(subset_to_outliers=False) #Get all flags
        #add checkname to the index
        checkdf["checkname"] = qcresult.checkname
        if qcresult.checkname in qc_before_timecoarsening:
            #Subset to coarsende timestmaps only
            checkdf = checkdf.reindex(possible_timestamps)

        checkdf.set_index("checkname", append=True, inplace=True)
        to_concat.append(checkdf)
        
    totaldf =  save_concat(to_concat)

    if totaldf.empty:
        return pd.DataFrame(columns=['value', 'label', 'details'],
                                index=pd.DatetimeIndex([], name='datetime'))

    #Unstack
    totaldf = totaldf.unstack(level='checkname')
    totaldf.fillna('Not applied', inplace=True)

    #add values
    allvals = pd.concat([sensor.series, sensor.outliers_values_bin]) #do not sort before removing the duplicates !
    allvals = allvals[~allvals.index.duplicated(keep='last')].sort_index()
    totaldf['value'] = allvals.loc[totaldf.index]
    
    return totaldf[['value', 'label', 'details']]


def station_qc_overview_df(station, subset_obstypes:Union[list[str], None] = None) -> pd.DataFrame:
    #TODO: docstring
    
    if subset_obstypes is None:
        sensortargets = station.sensordata.values() 
    else:
        sensortargets = []
        for obstype in subset_obstypes:
            if obstype in station.sensordata:
                sensortargets.append(station.get_sensor(obstype))
            else:
                #Log a warning?
                pass
    
    to_concat = []
    for sensordata in sensortargets:
        stadf = sensordata_qc_overview_df(sensordata).reset_index()
        #add obstype to the index
        if not stadf.empty:
            stadf["obstype"] = sensordata.obstype.name
            stadf = stadf.reset_index().set_index(['datetime', "obstype"])
            to_concat.append(stadf)
    
    
        
    totaldf =  save_concat(to_concat)
    totaldf.sort_index(inplace=True)

    if totaldf.empty:
        return pd.DataFrame(columns=['value', 'label', 'details'],
                            index=pd.MultiIndex(
                levels=[[], []], codes=[[], []], names=["datetime", "obstype"]))

    return totaldf[['value', 'label', 'details']]

def dataset_qc_overview_df(dataset, subset_stations:Union[list[str], None] = None,
                           subset_obstypes:Union[list[str], None] = None) -> pd.DataFrame:
    #TODO: docstring
    if subset_stations is None:
        stationtargets = dataset.stations
    else:
        stationtargets = [dataset.get_station(station_name) for station_name in subset_stations]
    
    to_concat = []
    for station in stationtargets:
        stadf = station_qc_overview_df(station, subset_obstypes=subset_obstypes).reset_index()
        #add obstype to the index
        if not stadf.empty:
            stadf["name"] = station.name
            stadf = stadf.reset_index().set_index(['datetime', "obstype", "name"])
            to_concat.append(stadf)
    
    totaldf =  save_concat(to_concat)
    totaldf.sort_index(inplace=True)

    if totaldf.empty:
        return pd.DataFrame(columns=['value', 'label', 'details'],
                            index=pd.MultiIndex(
                levels=[[], [], []], codes=[[], [], []], names=["datetime", "obstype", "name"]))

    return totaldf[['value', 'label', 'details']]