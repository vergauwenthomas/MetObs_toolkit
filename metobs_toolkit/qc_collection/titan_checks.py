import logging


logger = logging.getLogger(__name__)


def create_titanlib_points_dict(obsdf, metadf, obstype):
    """Create a dictionary of titanlib-points.

    Titanlib uses point as dataformats. This method converts the dataframes to
    a dictionary of points.

    Parameters
    ----------
    obsdf : pandas.DataFrame
        Dataset.df
    metadf : pandas.DataFrame
        Dataset.metadf.
    obstype : str
        The observation type to pass to the points.

    Returns
    -------
    points_dict : dict
        The collection of datapoints.

    """
    obs = xs_save(obsdf, obstype, "obstype")
    obs = obs.reset_index()

    # merge metadata
    obs = obs.merge(
        right=metadf[["lat", "lon", "altitude"]],
        how="left",
        left_on="name",
        right_index=True,
    )

    dt_grouper = obs.groupby("datetime")

    points_dict = {}
    for dt, group in dt_grouper:

        check_group = group[~group["value"].isnull()]

        points_dict[dt] = {
            "values": check_group["value"].to_numpy(),
            "names": check_group["name"].to_numpy(),
            "lats": check_group["lat"].to_numpy(),
            "lons": check_group["lon"].to_numpy(),
            "elev": check_group["altitude"].to_numpy(),
            "ignore_names": group[group["value"].isnull()]["name"].to_numpy(),
        }

    return points_dict


def titan_buddy_check(obsdf, metadf, obstype, checks_settings, titan_specific_labeler):
    """Apply the Titanlib buddy check.

    The buddy check compares an observation against its neighbors (i.e. buddies). The check looks for
    buddies in a neighborhood specified by a certain radius. The buddy check flags observations if the
    (absolute value of the) difference between the observations and the average of the neighbors
    normalized by the standard deviation in the circle is greater than a predefined threshold.


    Parameters
    ------------
    obsdf: Pandas.DataFrame
        The dataframe containing the observations
    metadf: Pandas.DataFrame
        The dataframe containing the metadata (e.g. latitude, longitude...)
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'
    checks_settings: Dictionary
        Dictionary with the settings for each check
    titan_specific_labeler: Dictionary
        Dictionary that maps numeric flags to 'ok' or 'outlier' flags for each titan check

    Returns
    ----------
    obsdf: Pandas.DataFrame
        The dataframe containing the unflagged-observations
    outlier_df : Pandas.DataFrame
        The dataframe containing the flagged observations

    """
    try:
        _ = metadf["altitude"]
    except:
        logger.warning("Cannot find altitude of weather stations. Check is skipped!")

    # Create points_dict
    pointsdict = create_titanlib_points_dict(obsdf, metadf, obstype)

    df_list = []
    for dt, point in pointsdict.items():
        obs = list(point["values"])
        titan_points = titanlib.Points(
            np.asarray(point["lats"]),
            np.asarray(point["lons"]),
            np.asarray(point["elev"]),
        )

        num_labels = titanlib.buddy_check(
            titan_points,
            np.asarray(obs),
            np.asarray(
                [checks_settings["radius"]] * len(obs)
            ),  # same radius for all stations
            np.asarray(
                [checks_settings["num_min"]] * len(obs)
            ),  # same min neighbors for all stations
            checks_settings["threshold"],
            checks_settings["max_elev_diff"],
            checks_settings["elev_gradient"],
            checks_settings["min_std"],
            checks_settings["num_iterations"],
            np.full(len(obs), 1),
        )  # check all

        labels = pd.Series(num_labels, name="num_label").to_frame()
        labels["name"] = point["names"]
        labels["datetime"] = dt
        df_list.append(labels)

    checkeddf = pd.concat(df_list)

    # Convert to toolkit format
    outliersdf = checkeddf[checkeddf["num_label"].isin(titan_specific_labeler["outl"])]

    outliersdf = outliersdf.set_index(["name", "datetime"])

    obsdf, outliersdf = make_outlier_df_for_check(
        station_dt_list=outliersdf.index,
        obsdf=obsdf,
        obstype=obstype,
        flag=label_def["titan_buddy_check"]["label"],
    )

    return obsdf, outliersdf


def titan_sct_resistant_check(
    obsdf, metadf, obstype, checks_settings, titan_specific_labeler
):
    """Apply the Titanlib (robust) Spatial-Consistency-Test (SCT).

    The SCT resistant check is a spatial consistency check which compares each observation to what is expected given the other observations in the
    nearby area. If the deviation is large, the observation is removed. The SCT uses optimal interpolation
    (OI) to compute an expected value for each observation. The background for the OI is computed from
    a general vertical profile of observations in the area.

    Parameters
    -------------
    obsdf: Pandas.DataFrame
        The dataframe containing the observations
    metadf: Pandas.DataFrame
        The dataframe containing the metadata (e.g. latitude, longitude...)
    obstype: String, optional
        The observation type that has to be checked. The default is 'temp'.
    checks_settings: Dictionary
        Dictionary with the settings for each check
    titan_specific_labeler: Dictionary
        Dictionary that maps numeric flags to 'ok' or 'outlier' flags for each titan check

    Returns
    ----------
    obsdf: Pandas.DataFrame
        The dataframe containing the unflagged-observations
    outlier_df : Pandas.DataFrame
        The dataframe containing the flagged observations
    """
    import time

    try:
        _ = metadf["altitude"]
    except:
        logger.warning("Cannot find altitude of weather stations. Check is skipped!")

    # Create points_dict
    pointsdict = create_titanlib_points_dict(obsdf, metadf, obstype)

    df_list = []
    for dt, point in pointsdict.items():
        logger.debug(f"sct on observations at {dt}")
        obs = list(point["values"])
        titan_points = titanlib.Points(
            np.asarray(point["lats"]),
            np.asarray(point["lons"]),
            np.asarray(point["elev"]),
        )

        flags, scores = titanlib.sct_resistant(
            points=titan_points,  # points
            values=np.asarray(obs),  # vlues
            obs_to_check=np.full(len(obs), 1),  # obs to check (check all)
            background_values=np.full(len(obs), 0),  # background values
            background_elab_type=titanlib.MedianOuterCircle,  # background elab type
            num_min_outer=checks_settings["num_min_outer"],  # num min outer
            num_max_outer=checks_settings["num_max_outer"],  # num mac outer
            inner_radius=checks_settings["inner_radius"],  # inner radius
            outer_radius=checks_settings["outer_radius"],  # outer radius
            num_iterations=checks_settings["num_iterations"],  # num iterations
            num_min_prof=checks_settings["num_min_prof"],  # num min prof
            min_elev_diff=checks_settings["min_elev_diff"],  # min elev diff
            min_horizontal_scale=checks_settings[
                "min_horizontal_scale"
            ],  # min horizontal scale
            max_horizontal_scale=checks_settings[
                "max_horizontal_scale"
            ],  # max horizontal scale
            kth_closest_obs_horizontal_scale=checks_settings[
                "kth_closest_obs_horizontal_scale"
            ],  # kth closest obs horizontal scale
            vertical_scale=checks_settings["vertical_scale"],  # vertical scale
            value_mina=[
                x - checks_settings["mina_deviation"] for x in obs
            ],  # values mina
            value_maxa=[
                x + checks_settings["maxa_deviation"] for x in obs
            ],  # values maxa
            value_minv=[
                x - checks_settings["minv_deviation"] for x in obs
            ],  # values minv
            value_maxv=[
                x + checks_settings["maxv_deviation"] for x in obs
            ],  # values maxv
            eps2=np.full(len(obs), checks_settings["eps2"]),  # eps2
            tpos=np.full(len(obs), checks_settings["tpos"]),  # tpos
            tneg=np.full(len(obs), checks_settings["tneg"]),  # tneg
            debug=checks_settings["debug"],  # debug
            basic=checks_settings["basic"],
        )  # basic

        logger.debug("Sleeping ... (to avoid segmentaton errors)")
        time.sleep(1)

        labels = pd.Series(flags, name="num_label").to_frame()
        labels["name"] = point["names"]
        labels["datetime"] = dt
        df_list.append(labels)

    checkeddf = pd.concat(df_list)

    # Convert to toolkit format
    outliersdf = checkeddf[checkeddf["num_label"].isin(titan_specific_labeler["outl"])]

    outliersdf = outliersdf.set_index(["name", "datetime"])

    obsdf, outliersdf = make_outlier_df_for_check(
        station_dt_list=outliersdf.index,
        obsdf=obsdf,
        obstype=obstype,
        flag=label_def["titan_sct_resistant_check"]["label"],
    )

    return obsdf, outliersdf
