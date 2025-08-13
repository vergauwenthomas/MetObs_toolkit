from typing import Optional, Dict, Any
import logging
import xarray as xr
import numpy as np
import pandas as pd


def modeltimeseries_to_xr(modeltimeseries: "Modeltimeseries") -> xr.Dataset:
    """
    Convert a model time series object to an xarray Dataset.

    The returned Dataset contains a single variable named after the
    observation type (e.g. 'temperature'). Its DataArray has three
    dimensions:
      - kind: distinguishes the nature of the stacked data. For model
              time series this contains a single value: 'model'.
      - models: the model name (length 1 here, prepared for concatenation).
      - datetime: timestamps of the model output.

    Attributes on the variable describe the observation type and model
    metadata.

    Parameters
    ----------
    modeltimeseries : ModelTimeseries
        Object holding a time-indexed pandas Series and related metadata.

    Returns
    -------
    xarray.Dataset
        Dataset with dimensions ('kind', 'models', 'datetime') where
        kind = ['model'].
    """
    ar = xr.DataArray(
        data=[[modeltimeseries.series.values]],
        coords={
            "datetime": modeltimeseries.series.index.get_level_values("datetime"),
            "models": [modeltimeseries.modelID],
            "kind": ["model"],
        },
        dims=["kind", "models", "datetime"],
        attrs={
            modeltimeseries.modelID: {
                "obstype_name": modeltimeseries.modelobstype.name,
                "obstype_desc": modeltimeseries.modelobstype.description,
                "obstype_unit": modeltimeseries.modelobstype.std_unit,
                "modelnameID": modeltimeseries.modelID,
                "modelvariable": modeltimeseries.modelobstype.model_band,
            }
        },
    )

    return xr.Dataset({modeltimeseries.modelobstype.name: ar})


def sensordata_to_xr(sensordata: "Sensordata") -> xr.Dataset:
    """
    Convert sensor observations (including labels) to an xarray Dataset.

    The returned Dataset contains one variable named after the observation
    type (e.g. 'temperature'). Its DataArray has:
      - kind dimension with two entries:
          'obs'   -> the measured (and possibly processed) numerical values
          'label' -> the associated integer / categorical QC or gap labels
      - datetime dimension with the observation timestamps.

    The 'obs' slice holds the physical observation values. The 'label'
    slice holds the label codes; QC and gap-fill method metadata are stored
    as attributes on that slice (accessible via the DataArray attributes).

    Parameters
    ----------
    sensordata : Sensordata
        Sensor data object containing a DataFrame with 'value' and 'label'
        plus metadata (QC, gap-fill info).

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

    xr_labels = xr.DataArray(
        data=[df["label"].values],
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
    sensor_attrs["QC"] = get_QC_info_in_dict(sensordata)
    sensor_attrs["GF"] = get_GF_info_in_dict(sensordata)

    # Combine along the type dimension
    xr_comb = xr.concat([xr_value, xr_labels], dim="kind")
    xr_comb.attrs = sensor_attrs
    return xr.Dataset({varname: xr_comb})


def station_to_xr(station: "Station") -> xr.Dataset:
    """
    Merge all sensor and model data of a station into a single Dataset.

    Each variable (per observation type) preserves its internal 'kind'
    dimension, which may include:
      - 'obs'   : sensor values
      - 'label' : sensor labels
      - 'model' : model time series (if present for that type)

    Datetimes from all contributing sources (sensors and model series) are
    unioned to build a common 'datetime' coordinate; individual variables
    are reindexed onto this union (introducing NaNs where data are absent).

    Station metadata (lat, lon, altitude, LCZ, and any extra data) are added
    as coordinates.

    Parameters
    ----------
    station : Station
        Station object containing sensor and model data.

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
    station_vars = [sens.to_xr() for sens in target_sensors]

    # Construct modeldata xr
    modelobs_vars = [
        modeltimeseries.to_xr() for modeltimeseries in station.modeldata]
    

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
    return ds


def dataset_to_xr(dataset: "Dataset") -> xr.Dataset:
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

    Returns
    -------
    xarray.Dataset
        Multi-station Dataset with a 'name' dimension plus any variable
        dimensions such as ('kind', 'datetime').
    """
    sta_xrdict = {sta.name: sta.to_xr() for sta in dataset.stations}
    ds = xr.concat(objs=sta_xrdict.values(), dim="name")
    ds = ds.assign_coords({"name": list(sta_xrdict.keys())})
    return ds


# ------------------------------------------
#    Helper
# ------------------------------------------


# ------------------------------------------
#    Attribute formatters and helpers
# ------------------------------------------


def get_QC_info_in_dict(sensordata: "Sensordata") -> Dict[str, Dict[str, Any]]:
    """
    Collect QC check settings applied to a sensor dataset.

    Parameters
    ----------
    sensordata : Sensordata
        Sensor data with recorded QC outlier checks.

    Returns
    -------
    dict
        Mapping of check name to its settings.
    """
    returndict = {}
    for qcdict in sensordata.outliers:
        returndict[qcdict["checkname"]] = {"settings": qcdict["settings"]}
    return returndict


def get_GF_info_in_dict(sensordata: "Sensordata") -> Dict[str, Dict[str, Any]]:
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
    returndict = {}
    # NOTE:iteration is done over all the gaps, this is a bit overkill?

    for gap in sensordata.gaps:
        if "applied_gapfill_method" in gap.fillsettings:
            method = gap.fillsettings["applied_gapfill_method"]
            gapsettings = gap.fillsettings
            del gapsettings["applied_gapfill_method"]  # delete key, so update works
            # create infodict
            gapinfo = {method: gapsettings}
            # UPdate
            returndict.update(gapinfo)
    return returndict
