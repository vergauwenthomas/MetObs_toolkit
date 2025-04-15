import os
from pathlib import Path
import ee
import json

# set path to repo root
libfolder = Path(str(Path(__file__).resolve())).parent.parent
# hardcode the account name of the service account
service_account = "metobs-service-account@metobs-public.iam.gserviceaccount.com"


# Get the Key (on github action use a secret, locally use the json in the .keys folder)
if "/runner/" in os.getcwd():
    # when executing in docs folder
    key_json = os.getenv("GEE_SERVICE_ACCOUNT")
    print("DEBUG: GEE_SERVICE_ACCOUNT:", key_json)
    if not key_json:
        raise EnvironmentError("GEE_SERVICE_ACCOUNT secret is not set.")
    # key_data = json.loads(key_json)
    credentials = ee.ServiceAccountCredentials(service_account, key_data=key_json)
    ee.Initialize(credentials)

else:
    # key is a json file in de .key direcotry
    if libfolder.joinpath(".keys").exists():
        # when executing in basefolder
        key_file = list(
            libfolder.joinpath(".keys").glob("*.json")
        )  # only 1 json file expected !!
        if len(key_file) == 1:
            key_file = key_file[0]
        elif len(key_file) > 1:
            raise FileNotFoundError(
                f"More than 1 *.json file found {libfolder.joinpath('.keys').absolute()}"
            )
        else:
            raise FileNotFoundError(
                f"No *.json file found {libfolder.joinpath('.keys').absolute()}"
            )
    else:
        raise FileNotFoundError(
            f"No .keys folder found {libfolder.joinpath('.keys').absolute()}"
        )

    credentials = ee.ServiceAccountCredentials(service_account, str(key_file))
    ee.Initialize(credentials)
