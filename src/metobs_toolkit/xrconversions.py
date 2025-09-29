from typing import Optional, Dict, Any, Union, Tuple
import logging
import xarray as xr
import numpy as np
import pandas as pd

from metobs_toolkit.settings_collection import __version__
from metobs_toolkit.settings_collection.label_defenitions import label_to_numeric_map


def modeltimeseries_to_xr(
    modeltimeseries: "Modeltimeseries", fmt_datetime_coordinate=True
) -> xr.Dataset:
    """
    Convert a model time series object to an xarray Dataset.

    The returned Dataset contains a single variable named after the
    observation type (e.g. 'temperature'). Its DataArray has three
    dimensions:

    * kind: distinguishes the nature of the stacked data. For model
      time series this contains a single value: 'model'.
    * models: the model name (length 1 here, prepared for concatenation).
    * datetime: timestamps of the model output.

    Attributes on the variable describe the observation type and model
    metadata.

    Parameters
    ----------
    modeltimeseries : ModelTimeseries
        Object holding a time-indexed pandas Series and related metadata.
    fmt_datetime_coordinate : bool, optional
        If True, format the datetime coordinate for CF compliance by removing
        timezone information and adding CF-compliant attributes. Default is True.

    Returns
    -------
    xarray.Dataset
        Dataset with dimensions ('kind', 'models', 'datetime') where
        kind = ['model'].
    """
    attrsdict = {
        "modelobstype_name": modeltimeseries.modelobstype.name,
        "modelobstype_desc": modeltimeseries.modelobstype.description,
        "modelobstype_unit": modeltimeseries.modelobstype.std_unit,
        "modelname": modeltimeseries.modelname,
        "modelvariable": modeltimeseries.modelvariable,
    }
    attrsdict.update(construct_metobs_attr())

    ar = xr.DataArray(
        data=[[modeltimeseries.series.values]],
        coords={
            "datetime": modeltimeseries.series.index.get_level_values("datetime"),
            "models": [modeltimeseries.modelname],
            "kind": ["model"],
        },
        dims=["kind", "models", "datetime"],
        attrs=attrsdict,
    )
    ds = xr.Dataset({modeltimeseries.modelobstype.name: ar})
    if fmt_datetime_coordinate:
        ds = format_datetime_coord(ds)
    return ds


def sensordata_to_xr(
    sensordata: "Sensordata", fmt_datetime_coordinate=True
) -> xr.Dataset:
    """
    Convert sensor observations (including labels) to an xarray Dataset.

    The returned Dataset contains one variable named after the observation
    type (e.g. 'temperature'). Its DataArray has:

    * kind dimension with two entries:

      * 'obs'   -> the measured (and possibly processed) numerical values
      * 'label' -> the associated integer / categorical QC or gap labels

    * datetime dimension with the observation timestamps.

    The 'obs' slice holds the physical observation values. The 'label'
    slice holds the label codes; QC and gap-fill method metadata are stored
    as attributes on that slice (accessible via the DataArray attributes).

    Parameters
    ----------
    sensordata : Sensordata
        Sensor data object containing a DataFrame with 'value' and 'label'
        plus metadata (QC, gap-fill info).
    fmt_datetime_coordinate : bool, optional
        If True, format the datetime coordinate for CF compliance by removing
        timezone information and adding CF-compliant attributes. Default is True.

    Returns
    -------
    xarray.Dataset
        Dataset with variable <obstype> and dimensions ('kind', 'datetime'),
        where kind = ['obs', 'label'].
    """
    df = sensordata.df  # contains obs, outliers and gaps
    varname = sensordata.obstype.name

    # Values
    xr_value = xr.DataArray(
        data=[df["value"].values],
        coords={"datetime": df.index.get_level_values("datetime"), "kind": ["obs"]},
        dims=["kind", "datetime"],
        attrs={},
    )
    # Labels
    # Note: to make the xr.Dataset serializable, we need to use integer labels!!
    present_labels = df["label"].unique()
    df["label_numeric"] = df["label"].map(label_to_numeric_map)

    applied_map = {
        f"Label:{key}": val
        for key, val in label_to_numeric_map.items()
        if key in present_labels
    }

    xr_labels = xr.DataArray(
        data=[df["label_numeric"].values],
        coords={"datetime": df.index.get_level_values("datetime"), "kind": ["label"]},
        dims=["kind", "datetime"],
        attrs={},
    )

    # Attributes
    sensor_attrs = {
        "obstype_name": sensordata.obstype.name,
        "obstype_desc": sensordata.obstype.description,
        "obstype_unit": sensordata.obstype.std_unit,
    }
    # Labels
    sensor_attrs.update(construct_metobs_attr())
    sensor_attrs.update(construct_QC_attr(sensordata))
    sensor_attrs.update(construct_GF_attr(sensordata))
    sensor_attrs.update(applied_map)  # only the applied labels

    # Combine along the type dimension
    xr_comb = xr.concat([xr_value, xr_labels], dim="kind")
    xr_comb.attrs = sensor_attrs

    # Construct dataset
    ds = xr.Dataset({varname: xr_comb})
    if fmt_datetime_coordinate:
        ds = format_datetime_coord(ds)

    return ds


def station_to_xr(station: "Station", fmt_datetime_coordinate=True) -> xr.Dataset:
    """
    Merge all sensor and model data of a station into a single Dataset.

    Each variable (per observation type) preserves its internal 'kind'
    dimension, which may include:

    * 'obs'   : sensor values
    * 'label' : sensor labels
    * 'model' : model time series (if present for that type)

    Datetimes from all contributing sources (sensors and model series) are
    unioned to build a common 'datetime' coordinate; individual variables
    are reindexed onto this union (introducing NaNs where data are absent).

    Station metadata (lat, lon, altitude, LCZ, and any extra data) are added
    as coordinates.

    Parameters
    ----------
    station : Station
        Station object containing sensor and model data.
    fmt_datetime_coordinate : bool, optional
        If True, format the datetime coordinate for CF compliance by removing
        timezone information and adding CF-compliant attributes. Default is True.

    Returns
    -------
    xarray.Dataset
        Dataset with (potential) dimensions ('kind', 'datetime') per
        variable; station metadata as scalar coordinates.
    """

    # --- Create variables per sensor----

    # Construct target sensors
    target_sensors = list(station.sensordata.values())

    # Create xrdataarrays
    station_vars = [
        sensordata_to_xr(sens, fmt_datetime_coordinate=False) for sens in target_sensors
    ]

    # Construct modeldata xr
    modelobs_vars = [
        modeltimeseries_to_xr(modeltimeseries, fmt_datetime_coordinate=False)
        for modeltimeseries in station.modeldata
    ]

    # The 'datetime' coordinate of the observations is not (persee) the
    # same as in the modelobs datasets. This leads to un-mergable dataset. To resolve
    # this we must construct a new datetime (union of all datetimes) and broadcast
    # all datasets to that coordiante

    all_to_join = [*station_vars, *modelobs_vars]
    # Get the union of all datetimes
    all_datetimes = np.unique(
        np.concatenate([pd.to_datetime(ds["datetime"].values) for ds in all_to_join])
    )
    # to datetimes again
    all_datetimes = pd.to_datetime(all_datetimes).tz_localize("UTC")

    # Broadcast
    all_reindexed = []
    for subds in all_to_join:
        all_reindexed.append(subds.reindex(datetime=all_datetimes))
        del subds

    # Create a xr Dataset of all variables
    ds = xr.merge(all_reindexed, combine_attrs="no_conflicts")

    # Station related coordinates
    sta_coords = {
        "lat": station.site.lat,
        "lon": station.site.lon,
        "altitude": station.site.altitude,
        "LCZ": station.site.LCZ,
    }
    extra_data_coords = {key: val for key, val in station.site.extradata.items()}
    sta_coords.update(extra_data_coords)

    ds = ds.assign_coords(sta_coords)
    if fmt_datetime_coordinate:
        ds = format_datetime_coord(ds)
    return ds


def dataset_to_xr(dataset: "Dataset", fmt_datetime_coordinate=True) -> xr.Dataset:
    """
    Concatenate multiple station Datasets into one along a new 'name'
    dimension.

    All per-station variables retain their internal 'kind' dimension
    (e.g. combinations of 'obs', 'label', 'model'). Only variables common
    across stations will align cleanly (xarray merge semantics apply).

    Parameters
    ----------
    dataset : Dataset
        Collection of Station objects.
    fmt_datetime_coordinate : bool, optional
        If True, format the datetime coordinate for CF compliance by removing
        timezone information and adding CF-compliant attributes. Default is True.

    Returns
    -------
    xarray.Dataset
        Multi-station Dataset with a 'name' dimension plus any variable
        dimensions such as ('kind', 'datetime').
    """
    sta_xrdict = {
        sta.name: station_to_xr(sta, fmt_datetime_coordinate=False)
        for sta in dataset.stations
    }
    ds = xr.concat(objs=sta_xrdict.values(), dim="name")
    ds = ds.assign_coords({"name": list(sta_xrdict.keys())})

    if fmt_datetime_coordinate:
        ds = format_datetime_coord(ds)

    return ds


# ------------------------------------------
#    Helper
# ------------------------------------------


# ------------------------------------------
#    Coordinates formatters
# ------------------------------------------
def format_datetime_coord(ds: xr.Dataset) -> xr.Dataset:
    """
    Format a datetime coordinate in an xarray Dataset for CF compliance.

    Removes timezone information from datetime coordinates and adds
    CF-compliant time attributes including 'standard_name', 'long_name',
    and 'timezone' attributes.

    Parameters
    ----------
    ds : xr.Dataset
        The Dataset containing a 'datetime' coordinate to format.

    Returns
    -------
    xr.Dataset
        The Dataset with formatted datetime coordinate.
    """
    if "UTC" in str(ds["datetime"].dtype) or "tz" in str(ds["datetime"].dtype):
        # Convert to pandas datetime index for timezone handling
        dt_index = pd.DatetimeIndex(ds["datetime"].values)

        if dt_index.tz is not None:
            # Store timezone info and convert to naive
            timezone_name = str(dt_index.tz)
            dt_naive = dt_index.tz_localize(None)
        else:
            # Handle datetime64[ns, UTC] format - manually remove timezone
            timezone_name = "UTC"
            # Convert to naive by accessing the underlying datetime64[ns] values
            dt_naive = pd.DatetimeIndex(ds["datetime"].values.astype("datetime64[ns]"))

        # Update coordinate with naive datetime
        ds = ds.assign_coords({"datetime": dt_naive})

        # Add CF-compliant time attributes (but avoid calendar to prevent conflicts)
        ds["datetime"].attrs["timezone"] = timezone_name
        ds["datetime"].attrs["standard_name"] = "time"
        ds["datetime"].attrs["long_name"] = "time"
    else:
        pass

    return ds


# ------------------------------------------
#    Attribute formatters and helpers
# ------------------------------------------
def fmt_attr_value(val: Any) -> Union[str, list, int, float]:
    """
    Format an attribute value for output.

    Parameters
    ----------
    val : Any
        The attribute value to format.

    Returns
    -------
    Any
        The formatted attribute value.
    """
    if isinstance(val, dict):
        raise ValueError(
            f"Nested dictionaries are not allowed in netCDF attributes. Found: {val}"
        )
    elif isinstance(val, type(None)):
        return str(val)
    elif isinstance(val, list):
        # return [fmt_attr_value(v) for v in val]
        return val
    elif isinstance(val, (str, int, float)):
        return val
    else:
        raise ValueError(f"Unsupported attribute type found: {type(val)}")


def construct_metobs_attr() -> Dict:
    """
    Construct metadata attributes for a MetObs dataset.

    Returns
    -------
    dict
        Metadata attributes for the dataset.
    """
    return {
        "MetObs toolkit version": __version__,
    }


def construct_QC_attr(sensordata: "Sensordata") -> Dict:
    """
    Extract quality control check metadata from sensor data.

    Creates attributes documenting all QC checks applied to the sensor data,
    including a list of check names and detailed settings for each check
    as flattened dictionary items.

    Parameters
    ----------
    sensordata : Sensordata
        Sensor data object containing outliers list with QC check metadata.

    Returns
    -------
    dict
        Dictionary with 'QC checks' list and flattened QC settings as
        'QC:{checkname}.{setting}' key-value pairs.
    """
    returndict = {"QC checks": []}
    for qcdict in sensordata.outliers:
        # add name to the list
        returndict["QC checks"].append(qcdict["checkname"])
        # add details as flat dict items
        for key, value in qcdict["settings"].items():
            returndict[f"QC:{qcdict['checkname']}.{key}"] = fmt_attr_value(
                value
            )  # NO NESTED DICT!

    return returndict


def construct_GF_attr(sensordata: "Sensordata") -> Dict[str, Dict[str, Any]]:
    """
    Extract applied gap-fill methods and their settings.

    Parameters
    ----------
    sensordata : Sensordata
        Sensor data containing gap objects with fill settings.

    Returns
    -------
    dict
        Mapping of gap-fill method names to the settings used.
    """
    returndict = {"GF methods": []}
    # NOTE:iteration is done over all the gaps, this is a bit overkill?

    for gap in sensordata.gaps:
        if "applied_gapfill_method" in gap.fillsettings:
            method = gap.fillsettings["applied_gapfill_method"]
            if method not in returndict["GF methods"]:
                returndict["GF methods"].append(method)
                # add details as flat dict items
                for settingname, settingval in gap.fillsettings.items():
                    if settingname == "applied_gapfill_method":
                        continue

                    returndict[f"GF:{method}.{settingname}"] = fmt_attr_value(
                        settingval
                    )  # NO NESTED DICT!

    return returndict
