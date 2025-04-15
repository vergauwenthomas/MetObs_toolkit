import os
from pathlib import Path
import ee


# set path to repo root
libfolder = Path(str(Path(__file__).resolve())).parent.parent
# hardcode the account name of the service account
service_account = "metobs-service-account@metobs-public.iam.gserviceaccount.com"


class GEE_Authenticator:
    def __init__(self):
        if "/runner/" in os.getcwd():
            self.auth_on_runner()
        else:
            self.auth_locally()

    def auth_on_runner(self, secret="GEE_SERVICE_ACCOUNT"):
        key_json = os.getenv(secret)

        if not key_json:
            raise EnvironmentError(f"{secret} secret is not set, are present in scope.")
        credentials = ee.ServiceAccountCredentials(service_account, key_data=key_json)
        ee.Initialize(credentials)

    def auth_locally(self):
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
