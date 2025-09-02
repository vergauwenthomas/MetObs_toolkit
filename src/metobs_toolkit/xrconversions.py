from typing import Optional, Dict, Any, Union, Tuple
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
            "models": [modeltimeseries.modelname],
            "kind": ["model"],
        },
        dims=["kind", "models", "datetime"],
        attrs={
            modeltimeseries.modelname: {
                "obstype_name": modeltimeseries.obstype.name,
                "obstype_desc": modeltimeseries.obstype.description,
                "obstype_unit": modeltimeseries.obstype.std_unit,
                "modelname": modeltimeseries.modelname,
                "modelvariable": modeltimeseries.modelvariable,
            }
        },
    )

    return xr.Dataset({modeltimeseries.obstype.name: ar})


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
    modelobs_vars = [modeltimeseries.to_xr() for modeltimeseries in station.modeldata]

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
#    NetCDF serialization helpers
# ------------------------------------------


def flatten_nested_dict(nested_dict: Dict[str, Any], prefix: str = "") -> Dict[str, str]:
    """
    Flatten nested dictionaries into dot-separated keys with string values.
    
    This function converts nested dictionaries to a flat structure that is 
    compatible with netCDF attribute serialization.
    
    Parameters
    ----------
    nested_dict : Dict[str, Any]
        The nested dictionary to flatten.
    prefix : str, optional
        Prefix to add to keys, by default "".
        
    Returns
    -------
    Dict[str, str]
        Flattened dictionary with string values.
        
    Examples
    --------
    >>> nested = {"QC": {"gross_value": {"settings": {"threshold": 10}}}}
    >>> flatten_nested_dict(nested)
    {"QC.gross_value.settings.threshold": "10"}
    """
    flat_dict = {}
    
    for key, value in nested_dict.items():
        new_key = f"{prefix}.{key}" if prefix else key
        
        if isinstance(value, dict):
            if len(value) == 0:
                # Handle empty dictionaries by storing them as a placeholder
                flat_dict[new_key] = "{empty_dict}"
            else:
                # Recursively flatten nested dictionaries
                flat_dict.update(flatten_nested_dict(value, new_key))
        else:
            # Convert all values to strings for netCDF compatibility
            flat_dict[new_key] = str(value)
    
    return flat_dict


def make_dataset_serializable(ds: xr.Dataset) -> xr.Dataset:
    """
    Convert an xarray Dataset to be netCDF serializable.
    
    This function handles:
    - Nested dictionaries in variable attributes (converted to flattened structure)
    - Timezone-aware datetime coordinates (converted to timezone-naive with tz attribute)
    - Mixed-type variables with 'kind' dimension (split into separate obs/label variables)
    - CF convention compliance for metadata
    
    Parameters
    ----------
    ds : xarray.Dataset
        The input Dataset with potentially non-serializable attributes.
        
    Returns
    -------
    xarray.Dataset
        A copy of the Dataset with all attributes made netCDF serializable.
    """
    # Create a copy to avoid modifying the original
    ds_copy = ds.copy(deep=True)
    
    # Handle datetime timezone conversion
    ds_copy = _handle_datetime_timezones(ds_copy)
    
    # Handle mixed-type variables with 'kind' dimension (this also flattens attributes)
    ds_copy = _handle_mixed_type_variables(ds_copy)
    
    # Process global dataset attributes
    new_global_attrs = {}
    for attr_name, attr_val in ds_copy.attrs.items():
        if isinstance(attr_val, dict):
            flattened = flatten_nested_dict(attr_val, attr_name)
            new_global_attrs.update(flattened)
        else:
            new_global_attrs[attr_name] = _make_value_serializable(attr_val)
    
    ds_copy.attrs = new_global_attrs
    
    return ds_copy


def _handle_mixed_type_variables(ds: xr.Dataset) -> xr.Dataset:
    """
    Handle variables with mixed data types along the 'kind' dimension.
    
    Variables that contain both observation data (floats) and label data (strings)
    along a 'kind' dimension are split into separate variables to ensure netCDF
    compatibility.
    
    Parameters
    ----------
    ds : xarray.Dataset
        Dataset with potentially mixed-type variables.
        
    Returns
    -------
    xarray.Dataset
        Dataset with mixed-type variables split into separate obs/label variables.
    """
    new_vars = {}
    vars_to_remove = []
    
    for var_name in ds.data_vars:
        var = ds[var_name]
        
        # Check if this variable has 'kind' dimension with mixed types
        if 'kind' in var.dims and 'kind' in ds.coords:
            kind_values = ds.coords['kind'].values
            
            # Check if we have both 'obs' and 'label' kinds
            if 'obs' in kind_values and 'label' in kind_values:
                # Split into separate variables
                obs_data = var.sel(kind='obs')
                label_data = var.sel(kind='label')
                
                # Prepare flattened attributes from original variable
                flattened_attrs = {}
                for attr_name, attr_val in var.attrs.items():
                    if isinstance(attr_val, dict):
                        # Flatten nested dictionaries
                        flattened = flatten_nested_dict(attr_val, attr_name)
                        flattened_attrs.update(flattened)
                    else:
                        flattened_attrs[attr_name] = _make_value_serializable(attr_val)
                
                # Create obs variable (ensure numeric type)
                obs_var_name = f"{var_name}_obs"
                try:
                    # Try to convert to float64 for consistency
                    obs_values = obs_data.values.astype(float)
                    new_vars[obs_var_name] = xr.DataArray(
                        obs_values,
                        coords={coord: obs_data.coords[coord] for coord in obs_data.coords if coord != 'kind'},
                        dims=[dim for dim in obs_data.dims if dim != 'kind'],
                        attrs={**flattened_attrs, 'data_type': 'observations'}
                    )
                except (ValueError, TypeError):
                    # If conversion fails, keep as object but mark it
                    new_vars[obs_var_name] = obs_data.drop_vars('kind', errors='ignore').assign_attrs(
                        {**flattened_attrs, 'data_type': 'observations'}
                    )
                
                # Create label variable (ensure string type)
                label_var_name = f"{var_name}_label"
                try:
                    # Convert to strings for consistency
                    label_values = label_data.values.astype('U')  # Unicode strings
                    new_vars[label_var_name] = xr.DataArray(
                        label_values,
                        coords={coord: label_data.coords[coord] for coord in label_data.coords if coord != 'kind'},
                        dims=[dim for dim in label_data.dims if dim != 'kind'],
                        attrs={**flattened_attrs, 'data_type': 'labels', 'description': 'Quality control and gap-fill labels'}
                    )
                except (ValueError, TypeError):
                    # If conversion fails, keep as object but mark it
                    new_vars[label_var_name] = label_data.drop_vars('kind', errors='ignore').assign_attrs(
                        {**flattened_attrs, 'data_type': 'labels', 'description': 'Quality control and gap-fill labels'}
                    )
                
                # Mark original variable for removal
                vars_to_remove.append(var_name)
            else:
                # Variable has 'kind' dimension but not mixed types, keep as-is but flatten attributes
                flattened_attrs = {}
                for attr_name, attr_val in var.attrs.items():
                    if isinstance(attr_val, dict):
                        flattened = flatten_nested_dict(attr_val, attr_name)
                        flattened_attrs.update(flattened)
                    else:
                        flattened_attrs[attr_name] = _make_value_serializable(attr_val)
                
                new_vars[var_name] = var.assign_attrs(flattened_attrs)
        else:
            # Variable doesn't have 'kind' dimension, keep as-is but flatten attributes
            flattened_attrs = {}
            for attr_name, attr_val in var.attrs.items():
                if isinstance(attr_val, dict):
                    flattened = flatten_nested_dict(attr_val, attr_name)
                    flattened_attrs.update(flattened)
                else:
                    flattened_attrs[attr_name] = _make_value_serializable(attr_val)
            
            new_vars[var_name] = var.assign_attrs(flattened_attrs)
    
    # Create new dataset
    ds_new = xr.Dataset(new_vars, coords=ds.coords, attrs=ds.attrs)
    
    # Remove 'kind' coordinate if no variables use it anymore
    remaining_vars_with_kind = [name for name in ds_new.data_vars 
                               if 'kind' in ds_new[name].dims]
    if not remaining_vars_with_kind and 'kind' in ds_new.coords:
        ds_new = ds_new.drop_vars('kind')
    
    return ds_new


def _handle_datetime_timezones(ds: xr.Dataset) -> xr.Dataset:
    """
    Handle timezone-aware datetime coordinates for CF compliance.
    
    Convert timezone-aware datetime coordinates to timezone-naive and store
    timezone information as a coordinate attribute following CF conventions.
    
    Parameters
    ----------
    ds : xarray.Dataset
        Input dataset with potentially timezone-aware coordinates.
        
    Returns
    -------
    xarray.Dataset
        Dataset with timezone-naive datetime coordinates and timezone attrs.
    """
    # Process datetime coordinates
    for coord_name, coord in ds.coords.items():
        if coord_name == 'datetime':
            # Check if the coordinate contains timezone info in dtype
            if 'UTC' in str(coord.dtype) or 'tz' in str(coord.dtype):
                try:
                    # Convert to pandas datetime index for timezone handling
                    dt_index = pd.DatetimeIndex(coord.values)
                    
                    if dt_index.tz is not None:
                        # Store timezone info and convert to naive
                        timezone_name = str(dt_index.tz)
                        dt_naive = dt_index.tz_localize(None)
                    else:
                        # Handle datetime64[ns, UTC] format - manually remove timezone
                        timezone_name = 'UTC'
                        # Convert to naive by accessing the underlying datetime64[ns] values
                        dt_naive = pd.DatetimeIndex(coord.values.astype('datetime64[ns]'))
                    
                    # Update coordinate with naive datetime
                    ds = ds.assign_coords({coord_name: dt_naive})
                    
                    # Add CF-compliant time attributes (but avoid calendar to prevent conflicts)
                    ds[coord_name].attrs['timezone'] = timezone_name
                    ds[coord_name].attrs['standard_name'] = 'time'
                    ds[coord_name].attrs['long_name'] = 'time'
                    # Note: calendar attribute is handled by xarray automatically
                        
                except Exception as e:
                    logging.warning(f"Could not process timezone for {coord_name}: {e}")
    
    # Also handle datetime data variables (not just coordinates)
    for var_name in ds.data_vars:
        var = ds[var_name]
        if 'datetime64' in str(var.dtype) and ('UTC' in str(var.dtype) or 'tz' in str(var.dtype)):
            try:
                # Convert timezone-aware datetime values to naive
                dt_values = pd.DatetimeIndex(var.values)
                if dt_values.tz is not None:
                    dt_naive = dt_values.tz_localize(None)
                else:
                    # Handle datetime64[ns, UTC] format manually
                    dt_naive = pd.DatetimeIndex(var.values.astype('datetime64[ns]'))
                
                ds[var_name] = (var.dims, dt_naive, var.attrs)
            except Exception as e:
                logging.warning(f"Could not convert variable {var_name}: {e}")
    
    return ds


def _make_value_serializable(value: Any) -> Union[str, int, float, list, tuple]:
    """
    Convert a value to a netCDF-serializable format.
    
    Parameters
    ----------
    value : Any
        The value to make serializable.
        
    Returns
    -------
    Union[str, int, float, list, tuple]
        Serializable version of the value.
    """
    if isinstance(value, (str, int, float, bool)):
        return value
    elif isinstance(value, (list, tuple)):
        return [_make_value_serializable(v) for v in value]
    elif isinstance(value, np.ndarray):
        return value.tolist()
    elif hasattr(value, '__iter__') and not isinstance(value, (str, bytes)):
        return [_make_value_serializable(v) for v in value]
    else:
        return str(value)


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
