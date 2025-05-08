import logging
from pathlib import Path
import numpy as np
import pandas as pd
import requests
import json
import pickle


# from metobs_toolkit.data_import import read_csv
from abc import ABC, abstractmethod


logger = logging.getLogger(__name__)


class FileReader(ABC):
    def __init__(self, file_path: str, is_url: bool = False):
        # filepath
        if is_url:
            self.file_path = str(file_path)
        else:
            self.file_path = Path(file_path)

        # helping methods
        self.is_url = is_url

    @abstractmethod
    def read(self):
        """Must return a pd.DataFrame representation of the data."""
        pass

    @abstractmethod
    def read_as_local_file(self):
        pass

    @abstractmethod
    def read_as_remote_file(self):
        pass

    def file_exist(self):
        return self.file_path.exists()


class CsvFileReader(FileReader):
    def __init__(self, file_path: str, is_url: bool = False):
        super().__init__(file_path, is_url=is_url)
        self.logger = logging.getLogger(__name__)

    def read(self, **readkwargs):
        if self.is_url:
            data = self.read_as_remote_file(**readkwargs)
        else:
            data = self.read_as_local_file(**readkwargs)

        return data

    def read_as_local_file(self, **readkwargs):
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")
        return read_csv(filepath=self.file_path, **readkwargs)

    def read_as_remote_file(self, **readkwargs):
        # Pandas can read online files with the same syntax
        return self.read_as_local_file(**readkwargs)


class JsonFileReader(FileReader):
    def __init__(self, file_path: str, is_url: bool = False):
        super().__init__(file_path, is_url=is_url)
        self.logger = logging.getLogger(__name__)

    def read(self, **readkwargs):
        if self.is_url:
            self.data = self.read_as_remote_file(**readkwargs)
        else:
            self.data = self.read_as_local_file(**readkwargs)

        self.is_already_read = True

    def read_as_local_file(self, **readkwargs):
        if not str(self.file_path).endswith(".json"):
            raise ImportError(f"The file {self.file_path} is not a JSON file.")
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")

        with open(self.file_path, "r") as file:
            data = json.load(file, **readkwargs)
        return data

    def read_as_remote_file(self, **readkwargs):
        r = requests.get(self.file_path, **readkwargs)
        data = r.json()
        return data


class PickleFileReader(FileReader):
    def __init__(self, file_path: str, is_url: bool = False):
        super().__init__(file_path, is_url=is_url)
        self.logger = logging.getLogger(__name__)

    def read(self, **readkwargs):
        if self.is_url:
            self.data = self.read_as_remote_file(**readkwargs)
        else:
            self.data = self.read_as_local_file(**readkwargs)

        self.is_already_read = True

    def read_as_local_file(self, **readkwargs):
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")

        with open(self.file_path, "rb") as inp:
            obj = pickle.load(inp, **readkwargs)

        return obj

    def read_as_remote_file(self, **readkwargs):
        raise ReferenceError(
            "No remote pickle reader implemented. Please make an issue if you want this functionallity."
        )


# ------------------------------------------
#    functions
# ------------------------------------------
def read_csv(filepath: str, **kwargs):
    logger.debug(f"Reading {filepath} to Dataframe.")

    if "sep" not in kwargs:
        df = read_csv_with_flexible_seperator(filepath, **kwargs)
    else:
        df = pd.read_csv(filepath, **kwargs)

    # Sanity checks
    if df.empty:
        raise EncodingWarning(f"The file {filepath} is read as an empty Dataframe!")
    elif len(df.columns) == 1:
        logger.warning(f"The file {filepath} is read as an single-column-Dataframe!")
    return df


def read_csv_with_flexible_seperator(filepath: str, **kwargs):
    """
    Reads a CSV file into a DataFrame, trying to infer the separator.

    Parameters:
    - filepath: str, path to the CSV file.
    - kwargs: additional keyword arguments to pass to pandas.read_csv.

    Returns:
    - DataFrame containing the data from the CSV file.
    """
    separators = [",", ";", "\t", "|", " "]

    for sep in separators:
        try:
            df = pd.read_csv(filepath, sep=sep, **kwargs)
            if df.shape[1] > 1:  # Assuming a valid CSV has more than one column
                return df
        except Exception as e:
            logger.debug(f"Failed to read {filepath} with separator '{sep}': {e}")

    raise ValueError(f"Could not determine the separator for {filepath}")
