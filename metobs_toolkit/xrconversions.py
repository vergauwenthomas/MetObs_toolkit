from typing import Optional, Dict, Any
import logging
import xarray as xr
import numpy as np
import pandas as pd


def modeltimeseries_to_xr(modeltimeseries: "Modeltimeseries") -> xr.Dataset:
    """
    Convert a model time series object to an xarray Dataset.

    Parameters
    ----------
    modeltimeseries : ModelTimeseries
        Object holding a time-indexed pandas Series and related metadata.

    Returns
    -------
    xarray.Dataset
        Dataset with one variable '<obstype>_modeltimeseries' indexed by
        'models' and 'datetime', including observation and model metadata.
    """
    ar = xr.DataArray(
                data=[modeltimeseries.series.values],
                coords={'datetime': modeltimeseries.series.index.get_level_values('datetime'),
                        'models': [modeltimeseries.modelname]},
                dims=['models', 'datetime'], 
                attrs={
                    'obstype_name': modeltimeseries.obstype.name,
                    'obstype_desc': modeltimeseries.obstype.description,
                    'obstype_unit': modeltimeseries.obstype.std_unit,
                    'modelname': modeltimeseries.modelname,
                    'modelvariable':modeltimeseries.modelvariable,
                    }
                )
    xr_comb = xr.Dataset(data_vars={f'{modeltimeseries.obstype.name}_modeltimeseries': ar})
    return xr_comb



def sensordata_to_xr(sensordata: "Sensordata") -> xr.Dataset:
    """
    Convert sensor observations (including labels) to an xarray Dataset.

    Parameters
    ----------
    sensordata : Sensordata
        Sensor data object containing a DataFrame with 'value' and 'label'
        plus metadata (QC, gap-fill info).

    Returns
    -------
    xarray.Dataset
        Dataset with two variables: '<obstype>' (values) and
        '<obstype>_labels' (label codes) along the 'datetime' dimension.
    """
    df = sensordata.df #contains obs, outliers and gaps

    dict_container={}
    #Values
    xr_value = xr.DataArray(
            data=df['value'].values,
            coords={'datetime': df.index.get_level_values('datetime')},
            dims=['datetime'], 
            attrs={
                'obstype_name': sensordata.obstype.name,
                'obstype_desc': sensordata.obstype.description,
                'obstype_unit': sensordata.obstype.std_unit,
                }
            )
    varname = sensordata.obstype.name
    dict_container[varname] = xr_value
    # labels
    label_attrs = {
        'QC': get_QC_info_in_dict(sensordata),
        'GF': get_GF_info_in_dict(sensordata) }
   

    xr_labels = xr.DataArray(
            data=df['label'].values,
            coords={'datetime': df.index.get_level_values('datetime')},
            dims=['datetime'], 
            attrs=label_attrs,
            )
    varname_labels = f"{varname}_labels"
    dict_container[varname_labels] = xr_labels
    
    xr_comb = xr.Dataset(data_vars=dict_container)
    
    return xr_comb



def station_to_xr(station: "Station", obstype: Optional[str] = None) -> xr.Dataset:
    """
    Merge all sensor and model data of a station into a single Dataset.

    Parameters
    ----------
    station : Station
        Station object containing sensor and model data.
    obstype : str, optional
        If provided, only include the specified observation type.

    Returns
    -------
    xarray.Dataset
        Dataset with unified 'datetime' axis, a 'name' dimension (station),
        and station metadata as coordinates.
    """
    #NOTE: LIMITATION: only usefull for synchronized data

    # --- Create variables per sensor---- 


    #Construct target sensors        
    target_sensors = []
    if obstype is None:
        target_sensors = list(station.sensordata.values())
    else:
        target_sensors.append(station.get_sensor(obstype))
    
    #Create xrdataarrays 
    station_vars = [sens.to_xr() for sens in target_sensors]
   
    #Construct modeldata xr
    modelobs_vars = [modeltimeseries.to_xr() for modeltimeseries in station.modeldata.values()]
    
    #The 'datetime' coordinate of the observations is not (persee) the
    #same as in the modelobs datasets. This leads to un-mergable dataset. To resolve
    #this we must construct a new datetime (union of all datetimes) and broadcast
    #all datasets to that coordiante
    
    all_to_join = [*station_vars, *modelobs_vars]
    # Get the union of all datetimes
    all_datetimes = np.unique(
        np.concatenate([pd.to_datetime(ds['datetime'].values) for ds in all_to_join])
    )
    #to datetimes again
    all_datetimes = pd.to_datetime(all_datetimes).tz_localize('UTC')
    
    #Broadcast
    all_reindexed = []
    for subds in all_to_join:
        all_reindexed.append(subds.reindex(datetime=all_datetimes))
        del subds

    #Create a xr Dataset of all variables
    ds = xr.merge(all_reindexed)
 
    #add the name dimension
    ds = ds.expand_dims('name').assign_coords({'name': [station.name]})

    #Add metadata as coordinates
    #Station related coordinates
    sta_coords = {"lat": ("name", [station.site.lat]),
                "lon": ("name", [station.site.lon]),
                "altitude": ("name", [station.site.altitude]),
                "LCZ": ("name", [station.site.LCZ])}
    extra_data_coords = {key: ('name', [val]) for key, val in station.site.extradata.items()}
    sta_coords.update(extra_data_coords)
    ds = ds.assign_coords(sta_coords)
    return ds


def dataset_to_xr(dataset: "Dataset", obstype: Optional[str] = None) -> xr.Dataset:
    """
    Concatenate multiple station Datasets into one along the 'name' dimension.

    Parameters
    ----------
    dataset : Dataset
        Collection of Station objects.
    obstype : str, optional
        Not currently used (reserved for future filtering).

    Returns
    -------
    xarray.Dataset
        Multi-station Dataset with shared variables and coordinates.
    """
    sta_xrlist = [sta.to_xr() for sta in dataset.stations]
    ds = xr.concat(sta_xrlist, dim='name')
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
        returndict[qcdict['checkname']] = {'settings': qcdict['settings']}
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
    #NOTE:iteration is done over all the gaps, this is a bit overkill?

    for gap in sensordata.gaps:
        if 'applied_gapfill_method' in gap.fillsettings:
            method = gap.fillsettings['applied_gapfill_method']
            gapsettings = gap.fillsettings
            del gapsettings['applied_gapfill_method'] #delete key, so update works
            #create infodict
            gapinfo = {method: gapsettings}
            #UPdate
            returndict.update(gapinfo)
    return returndict