#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Functions that are used for GEE interactions.

@author: thoverga
"""

import os
import logging
from time import sleep
import pandas as pd
import ee


logger = logging.getLogger(__name__)


def connect_to_gee(**kwargs):
    """
    Setup authentication for the use of the GEE Python API.

    For a fresh kernel, without stored credentials, a prompt/browser window
    will appear with further instructions for the authentication.


    Parameters
    ----------
    **kwargs : Kwargs passed to ee.Authenticate()
        Kwargs are only used by the user, for resetting the gee connection. See
        the Note below.

    Returns
    -------
    None.

    Note
    ------
    Upon calling, this function assumes you have a Google developers account,
    and a project with the Google Earth Engine API enabled.
    See the * Using Google Earth Engine * page for more info.

    Note
    ------
    During the Authentication, you will be asked if you want a read-only scope.
    A read-only scope is sufficient when the data is transferred directly to your
    machine (small data transfers), but will not be sufficient when extracting
    large amounts of data (typical for extracting Modeldata). This is because
    modeldata is written directly to your Google Drive, and therefore
    the read-only scope is insufficient.

    Note
    ------
    Due to several reasons, an EEExeption may be thrown. This is
    likely because of an invalid credential file. To fix this, you
    can update your credential file, and specify a specific authentication method.
    We found that the "notebook" authentication method works best for most users.

    Here is an example on how to update the credentials:

    .. code-block:: python

        import metobs_toolkit

        metobs_toolkit.connect_to_gee(force=True, #create new credentials
                                      auth_mode='notebook', # 'notebook', 'localhost', 'gcloud' (requires gcloud installed) or 'colab' (works only in colab)
                                      )

    """
    if "/runner/" in os.getcwd():  # Triggered on github action runner
        _auth_on_runner()
        return

    if os.getenv("READTHEDOCS_VIRTUALENV_PATH") is not None:  # Triggered on RTD builds
        _auth_on_rtd()
        return

    if bool(kwargs):  # kwargs are always passed by user, so reinitialize
        ee.Authenticate(**kwargs)
        ee.Initialize()
        return

    if not ee.data._credentials:  # check if ee connection is initialized
        ee.Authenticate()
        ee.Initialize()
    return


def _auth_on_rtd(secret="GEE_SERVICE_ACCOUNT"):
    logger.debug("Entering authentication on RTD funtion")
    service_account = "metobs-service-account@metobs-public.iam.gserviceaccount.com"
    # get variable (string reperesntation)
    key_str = os.getenv(secret)
    if key_str is None:
        raise EnvironmentError(
            f"{secret} variable is not set on RTD, or present in scope."
        )
    print("secret begin found: ", key_str)

    # # convert to json
    # trg_json_file = ".private_key_gee.json"
    # with open(trg_json_file, "w") as key_file:
    #     key_file.write(key_str)

    # authenticae with service account
    credentials = ee.ServiceAccountCredentials(service_account, key_data=key_str)
    ee.Initialize(credentials)


# import json
# #testing

# keystr = """{
#   "type": "service_account",
#   "project_id": "metobs-public",
#   "private_key_id": "3f911bc6ea5a7676b1ddda63ca8f6f4ef7649ec0",
#   "private_key": "-----BEGIN PRIVATE KEY-----\nMIIEvgIBADANBgkqhkiG9w0BAQEFAASCBKgwggSkAgEAAoIBAQC1FsRYveM1sARD\nrQ6UrLPF6HK3bOor7xXPpBXTCGcKmSLCI61wPIDJ6jJhyqGf8RsflZyoHTgsuUqo\nRRv7ipV2FORGrTO58hQ/HJOUPtCxS5e6ChDM2OKNUa1WgVONvQHuPQlTjUQy/rbW\n7Uy+E7NdByveCCBGx37gZCCMsuONRnVTackVoW5PrlksiuTH8zitmokBwxI5cj0X\nan+AjVUiMkwNtZeBLYRn16Zr+cUdilKshP5JcCh/My7G4KMTww8BR2DA+ohyi070\nbGdqbB6FDZ2336tCSUY1D3IqFODdWFwxq7QcDrZmJnfWb/Z9P9eSgPb1VB0MNJOg\nWScgpY3JAgMBAAECggEATT8s+n3l0h0HdKb5tUoGVcHWTZBUQ/F06GIiPSc0bTzt\nqsr1TQ9CEN+qJjT9xPBglZSIgt4T/F/+DNGOIjr3jqtSxSNVEVjGcjWKbo5tD3Qj\ngOSSTg+mdIoG2wPH1IpvrGS0+cMk+GvXKs+HEP3uYRySBeCJhCfNY4LSr7IPh09y\nYTeGVgl3yv7CPdg2txCDSG9J4CE663g/wppssTFnviFw0BxCVhrtaXuM0UMYsxpc\nN3RE1pO+kMey66nVNwYYV3hkUTLb5Daqa9GkQCwqsa4mAFPnthj16VYr5PYRANd6\n6bznJH/n8TggeGsIXlqsl+RIWEsn0lBo5hYVwm1m0QKBgQDtgMWkDuluFiW3EYZU\n1VktlEL44sEZYN4iYI+uMo6ege2KLYM2dWI3So123iLLW4juJAVsdTJvv841T1Z8\no3vj0EVbFB3Bnv7DzSHnrscNpX7DxC5V4S5nj37Y6RZlzUBYYn48nLAaFAaviHCq\n6MvVcGd/WAFkSgFvmpkkXKgLAwKBgQDDMTyGe0+Eq8GUeVoDOcdnuZtg8sSgkFEw\ns8IGtpaRj8+zZKMsiPr/IKadw2Rh1lfX5dx1RN6yE4wPcz9Fv59p7CzhScB9Nek6\nN6CA0LPf2MYgn4nY2+pyNDKGLZlwojPRiwyjn2Rbeii0Vj5KNkCN4svr6AOjgBUT\nVVqo6PTkQwKBgQCg0BuC7VHjZ32cCnKxiFA8y3HZgfgLzgo6rrU61yK4cvM7J3v7\n6Nla5NEKlnhqx4zc6mj9uhEvl2jxscm21R9y7re3ZtSLILQSMhht/mrrc4500aYq\nIjHAj2ntR04SGjsiXXZqZhHbZonsWu8m2kACQnzhvd4bYPy966kb9N0XrwKBgGrx\nmArKlvfGrwotLs8joghHnkUQ3gBm8cTwSDcrZPWdyJuuYOSurG6KMh+wBVjBemGU\n1CQANs90fEGe0CTj53C+cJyCrw8rMXyvK6ZIMrVLbMpE/t6tJxepV6FCbJZUmHWP\nbL5dPzwLgy/DLN+2N8pFX5sCBLRZIbL1CfYPpNVHAoGBAJLuFxnRmM+iS8//VOST\nAy10rVr/U8Ks1LoHqUVoVBrwnnXgmPQsrqqTC1Ju3X0LtalGv3EUp7S8zZ0SVoRd\nJlLVfJWZqRooTjpCM1SNM+Jo1ivsywi/rp7xVc+5868huybnyfJW1ciqnZdhIg+C\nZmjoUmwt1eq6pPIa9sZLbA9L\n-----END PRIVATE KEY-----\n",
#   "client_email": "metobs-service-account@metobs-public.iam.gserviceaccount.com",
#   "client_id": "108556761234228647120",
#   "auth_uri": "https://accounts.google.com/o/oauth2/auth",
#   "token_uri": "https://oauth2.googleapis.com/token",
#   "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
#   "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/metobs-service-account%40metobs-public.iam.gserviceaccount.com",
#   "universe_domain": "googleapis.com"
# }"""

# formatting
# keystr = keystr.replace('\n', '')
# # Convert the keystr to a JSON object
# key_json = json.loads(keystr, strict=False)
# keystr = keystr.replace("'", '"')
# keydict = eval(keystr)
# keydict2 = {
#   "type": "service_account",
#   "project_id": "metobs-public",
#   "private_key_id": "3f911bc6ea5a7676b1ddda63ca8f6f4ef7649ec0",
#   "private_key": "-----BEGIN PRIVATE KEY-----\nMIIEvgIBADANBgkqhkiG9w0BAQEFAASCBKgwggSkAgEAAoIBAQC1FsRYveM1sARD\nrQ6UrLPF6HK3bOor7xXPpBXTCGcKmSLCI61wPIDJ6jJhyqGf8RsflZyoHTgsuUqo\nRRv7ipV2FORGrTO58hQ/HJOUPtCxS5e6ChDM2OKNUa1WgVONvQHuPQlTjUQy/rbW\n7Uy+E7NdByveCCBGx37gZCCMsuONRnVTackVoW5PrlksiuTH8zitmokBwxI5cj0X\nan+AjVUiMkwNtZeBLYRn16Zr+cUdilKshP5JcCh/My7G4KMTww8BR2DA+ohyi070\nbGdqbB6FDZ2336tCSUY1D3IqFODdWFwxq7QcDrZmJnfWb/Z9P9eSgPb1VB0MNJOg\nWScgpY3JAgMBAAECggEATT8s+n3l0h0HdKb5tUoGVcHWTZBUQ/F06GIiPSc0bTzt\nqsr1TQ9CEN+qJjT9xPBglZSIgt4T/F/+DNGOIjr3jqtSxSNVEVjGcjWKbo5tD3Qj\ngOSSTg+mdIoG2wPH1IpvrGS0+cMk+GvXKs+HEP3uYRySBeCJhCfNY4LSr7IPh09y\nYTeGVgl3yv7CPdg2txCDSG9J4CE663g/wppssTFnviFw0BxCVhrtaXuM0UMYsxpc\nN3RE1pO+kMey66nVNwYYV3hkUTLb5Daqa9GkQCwqsa4mAFPnthj16VYr5PYRANd6\n6bznJH/n8TggeGsIXlqsl+RIWEsn0lBo5hYVwm1m0QKBgQDtgMWkDuluFiW3EYZU\n1VktlEL44sEZYN4iYI+uMo6ege2KLYM2dWI3So123iLLW4juJAVsdTJvv841T1Z8\no3vj0EVbFB3Bnv7DzSHnrscNpX7DxC5V4S5nj37Y6RZlzUBYYn48nLAaFAaviHCq\n6MvVcGd/WAFkSgFvmpkkXKgLAwKBgQDDMTyGe0+Eq8GUeVoDOcdnuZtg8sSgkFEw\ns8IGtpaRj8+zZKMsiPr/IKadw2Rh1lfX5dx1RN6yE4wPcz9Fv59p7CzhScB9Nek6\nN6CA0LPf2MYgn4nY2+pyNDKGLZlwojPRiwyjn2Rbeii0Vj5KNkCN4svr6AOjgBUT\nVVqo6PTkQwKBgQCg0BuC7VHjZ32cCnKxiFA8y3HZgfgLzgo6rrU61yK4cvM7J3v7\n6Nla5NEKlnhqx4zc6mj9uhEvl2jxscm21R9y7re3ZtSLILQSMhht/mrrc4500aYq\nIjHAj2ntR04SGjsiXXZqZhHbZonsWu8m2kACQnzhvd4bYPy966kb9N0XrwKBgGrx\nmArKlvfGrwotLs8joghHnkUQ3gBm8cTwSDcrZPWdyJuuYOSurG6KMh+wBVjBemGU\n1CQANs90fEGe0CTj53C+cJyCrw8rMXyvK6ZIMrVLbMpE/t6tJxepV6FCbJZUmHWP\nbL5dPzwLgy/DLN+2N8pFX5sCBLRZIbL1CfYPpNVHAoGBAJLuFxnRmM+iS8//VOST\nAy10rVr/U8Ks1LoHqUVoVBrwnnXgmPQsrqqTC1Ju3X0LtalGv3EUp7S8zZ0SVoRd\nJlLVfJWZqRooTjpCM1SNM+Jo1ivsywi/rp7xVc+5868huybnyfJW1ciqnZdhIg+C\nZmjoUmwt1eq6pPIa9sZLbA9L\n-----END PRIVATE KEY-----\n",
#   "client_email": "metobs-service-account@metobs-public.iam.gserviceaccount.com",
#   "client_id": "108556761234228647120",
#   "auth_uri": "https://accounts.google.com/o/oauth2/auth",
#   "token_uri": "https://oauth2.googleapis.com/token",
#   "auth_provider_x509_cert_url": "https://www.googleapis.com/oauth2/v1/certs",
#   "client_x509_cert_url": "https://www.googleapis.com/robot/v1/metadata/x509/metobs-service-account%40metobs-public.iam.gserviceaccount.com",
#   "universe_domain": "googleapis.com"
# }

# serialize the info into a string
# test_json = json.dumps(keydict, indent=4).replace("'", '"')
# test_json2 = json.dumps(keydict2, indent=4).replace("'", '"')
# Save the JSON string to a file
# trg_json_file = ".private_key_gee.json"
# with open(trg_json_file, "w") as key_file:
#     key_file.write(test_json)

# print(f"JSON saved to {trg_json_file}")
# format:

# test = test.replace('\n','')
# test = test.replace('-----BEGIN PRIVATE KEY-----', '')

# email = keydict['client_email']
# ee.ServiceAccountCredentials(email=email,
#                             key_data=test_json)


# _auth_on_rtd()

print("done")


# ----
def _auth_on_runner(secret="GEE_SERVICE_ACCOUNT"):
    """
    Authenticate and initialize the Google Earth Engine (GEE) API using a service account.



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
    --------
    This function is only relevent when there is no stdinput available (typical github runners)
    """
    logger.debug("Entering authentication on runnen funtion")
    service_account = "metobs-service-account@metobs-public.iam.gserviceaccount.com"
    key_json = os.getenv(secret)
    if not key_json:
        raise EnvironmentError(f"{secret} secret is not set, are present in scope.")
    credentials = ee.ServiceAccountCredentials(service_account, key_data=key_json)
    ee.Initialize(credentials)


# =============================================================================
# Object convertors
# =============================================================================


def _datetime_to_gee_datetime(datetime):
    # ASSUME DATETIME IN UTC!!!
    # utcdt = datetime.astimezone(pytz.utc)#This will c
    logger.debug(datetime.replace(tzinfo=None))
    return ee.Date(datetime.replace(tzinfo=None))


def get_ee_obj(model, target_bands=[], force_mosaic=None):
    """Get an image from a GEE object."""

    # get the dataset
    if model.is_image:
        obj = ee.Image(model.location).select(model.band_of_use)
        if obj is None:
            raise MetObsGEEDatasetError(
                f"No image returend by: ee.Image({model.location}).select({model.band_of_use}) "
            )
    else:
        obj = ee.ImageCollection(model.location)

        # filter the bands
        if bool(target_bands):
            obj = obj.select(ee.List(target_bands))
        else:
            obj = obj.select(model.band_of_use)

        # mosaic over the collectiontion (== convert imagecollection to image)
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
                f"No ee.Image could be constructed. Check the location, bandnames and type of the GeeDataset."
            )
    return obj


# =============================================================================
# Helpers
# =============================================================================


def _is_eeobj_empty(geeobj):
    """Check if a GEE object (Image or ImageCollection) is empty."""
    if isinstance(geeobj, ee.ImageCollection):
        return geeobj.size().getInfo() == 0
    elif isinstance(geeobj, ee.Image):
        return geeobj.bandNames().size().getInfo() == 0
    else:
        raise TypeError("geeobj must be an ee.Image or ee.ImageCollection")


def _is_image(geeobj):
    return isinstance(geeobj, ee.image.Image)


def _is_imagecollection(geeobj):
    return isinstance(geeobj, ee.imagecollection.ImageCollection)


def _validate_metadf(metadf):
    """Test if metadf is valid for GEE extraction.

    Returns True if metadata is suitable for gee extraction.

    :param metadf: metadata dataframe
    :type metadf: pd.DataFrame
    :return: True if oke, else False
    :rtype: Bool

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
    except:
        return False

    return True


def _addDate(image):
    """Add the image datetime as a band."""
    img_date = ee.Date(image.date())
    img_date = ee.Number.parse(img_date.format("YYYYMMddHHmmss"))
    return image.addBands(ee.Image(img_date).rename("datetime"))


def _df_to_features_point_collection(df):
    """Convert a dataframe to a featurecollections row-wise."""
    features = []
    for index, row in df.reset_index().iterrows():
        #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]])
        #     construct the attributes (properties) for each point
        poi_properties = poi_properties = {"name": row["name"]}
        #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def _df_to_features_buffer_collection(df, bufferradius):
    """Convert a dataframe to a featurecollections row-wise."""
    features = []
    for index, row in df.reset_index().iterrows():
        #     construct the geometry from dataframe
        poi_geometry = ee.Geometry.Point([row["lon"], row["lat"]]).buffer(
            distance=bufferradius
        )
        #     construct the attributes (properties) for each point
        poi_properties = poi_properties = {"name": row["name"]}
        #     construct feature combining geometry and properties
        poi_feature = ee.Feature(poi_geometry, poi_properties)
        features.append(poi_feature)

    return ee.FeatureCollection(features)


def coordinates_available(metadf, latcol="lat", loncol="lon"):
    """Test if all coordinates are available."""
    if metadf[latcol].isnull().all():
        logger.warning("No coordinates are found!")
        return False
    if metadf[loncol].isnull().all():
        logger.warning("No coordinates are found!")
        return False
    return True


def _estimate_data_size(metadf, startdt, enddt, time_res, n_bands=1):
    datatimerange = pd.date_range(start=startdt, end=enddt, freq=time_res)

    return metadf.shape[0] * len(datatimerange) * n_bands


class MetObsGEEDatasetError(Exception):
    """Raise when there is an issue with a GEE api call"""

    pass
