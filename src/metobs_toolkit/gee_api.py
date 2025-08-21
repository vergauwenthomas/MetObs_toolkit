#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions for Google Earth Engine (GEE) interactions.

This module provides utilities for authenticating and interacting with the GEE Python API,
including helpers for converting objects and validating metadata.

@author: thoverga
"""

import os
import json
import logging
from pathlib import Path

import pandas as pd
import ee

from metobs_toolkit.backend_collection.errorclasses import MetObsGEEDatasetError

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger(__file__)


@log_entry
def connect_to_gee(**kwargs) -> None:
    """
    Set up authentication for the use of the GEE Python API.

    For a fresh kernel, without stored credentials, a prompt/browser window
    will appear with further instructions for the authentication.

    Parameters
    ----------
    **kwargs : dict
        Keyword arguments passed to ee.Authenticate(). Used for resetting the GEE connection.

    Returns
    -------
    None

    Notes
    -----

    * This function assumes you have a Google developers account and a project with the
      Google Earth Engine API enabled.
    * During authentication, you may be asked if you want a read-only scope. A read-only
      scope is sufficient for small data transfers, but not for extracting large amounts
      of data (e.g., model data written to Google Drive).
    * If an EEException is thrown, it is likely due to an invalid credential file. You can
      update your credential file and specify a specific authentication method.
    * Example usage:

      .. code-block:: python

          import metobs_toolkit

          metobs_toolkit.connect_to_gee(force=True, # create new credentials
                                        auth_mode='notebook', # 'notebook', 'localhost', 'gcloud', or 'colab'
                                        )

    Raises
    ------
    TypeError
        If kwargs is not a dict.
    """
    logger.info("Entering connect_to_gee function.")

    if "/runner/" in os.getcwd():  # Triggered on GitHub action runner
        auth_on_runner()
        return

    if os.getenv("READTHEDOCS_VIRTUALENV_PATH") is not None:  # Triggered on RTD builds
        auth_on_rtd()
        return

    if bool(kwargs):  # kwargs are always passed by user, so reinitialize
        ee.Authenticate(**kwargs)
        ee.Initialize()
        return

    if not ee.data._credentials:  # check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return


@log_entry
def auth_on_rtd(secret: str = "GEE_SERVICE_ACCOUNT") -> None:
    """
    Authenticate and initialize the GEE API using a service account for ReadTheDocs builds.

    Parameters
    ----------
    secret : str, optional
        The environment variable name containing the service account key JSON,
        by default "GEE_SERVICE_ACCOUNT".

    Raises
    ------
    EnvironmentError
        If the specified environment variable is not set or the key JSON is missing.
    """

    if os.getenv(secret) is None:
        raise EnvironmentError(f"{secret} variable is not set or not present in scope.")
    key_str = r"{}".format(os.getenv(secret))

    # Write the data to a JSON file
    json_data = json.loads(key_str.replace("\n", "\\n"))
    output_path = Path.cwd() / "gee_service_account.json"
    with open(output_path, "w") as json_file:
        json.dump(json_data, json_file, indent=4)

    # Get the credentials
    email = "metobs-service-account@metobs-public.iam.gserviceaccount.com"
    credentials = ee.ServiceAccountCredentials(email=email, key_file=str(output_path))

    # Initiate Google API
    ee.Initialize(credentials)


@log_entry
def auth_on_runner(secret: str = "GEE_SERVICE_ACCOUNT") -> None:
    """
    Authenticate and initialize the GEE API using a service account for CI runners.

    Parameters
    ----------
    secret : str, optional
        The environment variable name containing the service account key JSON,
        by default "GEE_SERVICE_ACCOUNT".

    Raises
    ------
    EnvironmentError
        If the specified environment variable is not set or the key JSON is missing.

    Warning
    -------
    This function is only relevant when there is no stdin available (typical for GitHub runners).
    """
    service_account = "metobs-service-account@metobs-public.iam.gserviceaccount.com"
    key_json = os.getenv(secret)
    if not key_json:
        raise EnvironmentError(f"{secret} secret is not set, are present in scope.")
    credentials = ee.ServiceAccountCredentials(service_account, key_data=key_json)
    ee.Initialize(credentials)


@log_entry
def datetime_to_gee_datetime(datetime_obj) -> ee.Date:
    """
    Convert a Python datetime object to a GEE ee.Date object.

    Parameters
    ----------
    datetime_obj : datetime.datetime
        The datetime object to convert. Assumed to be in UTC.

    Returns
    -------
    ee.Date
        The corresponding GEE ee.Date object.

    Raises
    ------
    TypeError
        If datetime_obj is not a datetime.datetime instance.
    """
    return ee.Date(datetime_obj.replace(tzinfo=None))


@log_entry
def get_ee_obj(model, target_bands: list = [], force_mosaic: bool = None):
    """
    Get an ee.Image or ee.ImageCollection from a GEE object.

    Parameters
    ----------
    model : object
        The model object with attributes for GEE dataset location and bands.
    target_bands : list, optional
        List of bands to select. If empty, uses model.band_of_use.
    force_mosaic : bool, optional
        If True, mosaics the image collection. If None, uses model._is_mosaic.

    Returns
    -------
    ee.Image or ee.ImageCollection
        The resulting GEE object.

    Raises
    ------
    TypeError
        If target_bands is not a list or force_mosaic is not a bool or None.
    MetObsGEEDatasetError
        If no image could be constructed.
    """

    # get the dataset
    if model.is_image:
        obj = ee.Image(model.location).select(model.band_of_use)
        if obj is None:
            raise MetObsGEEDatasetError(
                f"No image returned by: ee.Image({model.location}).select({model.band_of_use}) "
            )
    else:
        obj = ee.ImageCollection(model.location)

        # filter the bands
        if bool(target_bands):
            obj = obj.select(ee.List(target_bands))
        else:
            obj = obj.select(model.band_of_use)

        # mosaic over the collection (== convert imagecollection to image)
        if force_mosaic is None:
            # use attribute of the model
            if model._is_mosaic:
                obj = obj.mosaic()
        elif force_mosaic:
            obj = obj.mosaic()
        else:
            obj = obj

        if obj is None:
            raise MetObsGEEDatasetError(
                "No ee.Image could be constructed. Check the location, band names, and type of the GeeDataset."
            )
    return obj


def _is_eeobj_empty(geeobj) -> bool:
    """
    Check if a GEE object (Image or ImageCollection) is empty.

    Parameters
    ----------
    geeobj : ee.Image or ee.ImageCollection
        The GEE object to check.

    Returns
    -------
    bool
        True if the object is empty, False otherwise.

    Raises
    ------
    TypeError
        If geeobj is not an ee.Image or ee.ImageCollection.
    """
    if isinstance(geeobj, ee.ImageCollection):
        return geeobj.size().getInfo() == 0
    elif isinstance(geeobj, ee.Image):
        return geeobj.bandNames().size().getInfo() == 0


def _is_image(geeobj) -> bool:
    """
    Check if a GEE object is an ee.Image.

    Parameters
    ----------
    geeobj : object
        The object to check.

    Returns
    -------
    bool
        True if geeobj is an ee.Image, False otherwise.
    """
    return isinstance(geeobj, ee.image.Image)


def _is_imagecollection(geeobj) -> bool:
    """
    Check if a GEE object is an ee.ImageCollection.

    Parameters
    ----------
    geeobj : object
        The object to check.

    Returns
    -------
    bool
        True if geeobj is an ee.ImageCollection, False otherwise.
    """
    return isinstance(geeobj, ee.imagecollection.ImageCollection)


def _validate_metadf(metadf: pd.DataFrame) -> bool:
    """
    Test if metadf is valid for GEE extraction.

    Returns True if metadata is suitable for GEE extraction.

    Parameters
    ----------
    metadf : pd.DataFrame
        Metadata dataframe.

    Returns
    -------
    bool
        True if valid, False otherwise.

    Raises
    ------
    TypeError
        If metadf is not a pandas DataFrame.
    """
    if metadf.empty:
        return False
    if metadf["geometry"].x.isnull().values.all():
        return False
    if metadf["geometry"].y.isnull().values.all():
        return False
    try:
        # Just testing if it can be converted
        metadf = metadf.to_crs("epsg:4326")
    except Exception:
        return False

    return True


def _addDate(image) -> ee.Image:
    """
    Add the image datetime as a band.

    Parameters
    ----------
    image : ee.Image
        The image to which the datetime band will be added.

    Returns
    -------
    ee.Image
        The image with an added 'datetime' band.

    Raises
    ------
    TypeError
        If image is not an ee.Image.
    """
    img_date = ee.Date(image.date())
    img_date = ee.Number.parse(img_date.format("YYYYMMddHHmmss"))
    return image.addBands(ee.Image(img_date).rename("datetime"))


def _df_to_features_point_collection(df: pd.DataFrame) -> ee.FeatureCollection:
    """
    Convert a DataFrame to a FeatureCollection row-wise (points).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns 'lon', 'lat', and 'name'.

    Returns
    -------
    ee.FeatureCollection
        FeatureCollection of points with properties.

    Raises
    ------
    TypeError
        If df is not a pandas DataFrame.
    """
    features = []
    for index, row in df.reset_index().iterrows():
        # Construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]])
        # Construct the attributes (properties) for each point
        poi_properties = {"name": row["name"]}
        # Construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def _df_to_features_buffer_collection(
    df: pd.DataFrame, bufferradius: float
) -> ee.FeatureCollection:
    """
    Convert a DataFrame to a FeatureCollection row-wise (buffered points).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with columns 'lon', 'lat', and 'name'.
    bufferradius : float
        Buffer radius in meters.

    Returns
    -------
    ee.FeatureCollection
        FeatureCollection of buffered points with properties.

    Raises
    ------
    TypeError
        If df is not a pandas DataFrame or bufferradius is not a float or int.
    """
    features = []
    for index, row in df.reset_index().iterrows():
        # Construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]]).buffer(
            distance=bufferradius
        )
        # Construct the attributes (properties) for each point
        poi_properties = {"name": row["name"]}
        # Construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


@log_entry
def coordinates_available(
    metadf: pd.DataFrame, latcol: str = "lat", loncol: str = "lon"
) -> bool:
    """
    Test if all coordinates are available in the metadata DataFrame.

    Parameters
    ----------
    metadf : pd.DataFrame
        Metadata DataFrame.
    latcol : str, optional
        Name of the latitude column, by default "lat".
    loncol : str, optional
        Name of the longitude column, by default "lon".

    Returns
    -------
    bool
        True if all coordinates are available, False otherwise.

    Raises
    ------
    TypeError
        If metadf is not a pandas DataFrame or latcol/loncol are not strings.
    """
    if metadf[latcol].isnull().all():
        logger.warning("No coordinates are found!")
        return False
    if metadf[loncol].isnull().all():
        logger.warning("No coordinates are found!")
        return False
    return True


def _estimate_data_size(
    metadf: pd.DataFrame, startdt, enddt, time_res: str, n_bands: int = 1
) -> int:
    """
    Estimate the size of the data to be extracted.

    Parameters
    ----------
    metadf : pd.DataFrame
        Metadata DataFrame.
    startdt : datetime.datetime
        Start datetime.
    enddt : datetime.datetime
        End datetime.
    time_res : str
        Pandas-compatible frequency string (e.g., '1H').
    n_bands : int, optional
        Number of bands, by default 1.

    Returns
    -------
    int
        Estimated number of data points.

    Raises
    ------
    TypeError
        If arguments are not of the correct type.
    """

    datatimerange = pd.date_range(start=startdt, end=enddt, freq=time_res)
    return metadf.shape[0] * len(datatimerange) * n_bands
