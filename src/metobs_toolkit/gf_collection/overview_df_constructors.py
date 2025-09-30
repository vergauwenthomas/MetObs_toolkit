"""Collection of DF constructing functions on various levels
(sensordata, station, dataset) for overviews and summaries of Gaps."""

import pandas as pd
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.df_helpers import save_concat


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
