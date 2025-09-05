import logging
from pathlib import Path
import json
import pickle
import requests
from abc import ABC, abstractmethod
from typing import Union

import pandas as pd

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


class FileReader(ABC):
    """
    Abstract base class for file readers.

    Parameters
    ----------
    file_path : str
        Path to the file or URL.
    is_url : bool, optional
        Whether the file_path is a URL (default is False).
    """

    def __init__(self, file_path: str, is_url: bool = False):
        """Initialize FileReader."""
        if is_url:
            self.file_path = str(file_path)
        else:
            self.file_path = Path(file_path)
        self.is_url = is_url

    @abstractmethod
    @log_entry
    def read(self, **kwargs):
        """
        Read the file and return its contents.

        Returns
        -------
        object
            The contents of the file, typically as a pandas DataFrame or other appropriate object.
        """
        pass

    @abstractmethod
    @log_entry
    def read_as_local_file(self, **kwargs):
        """
        Read the file as a local file.

        Returns
        -------
        object
            The contents of the file.
        """
        pass

    @abstractmethod
    @log_entry
    def read_as_remote_file(self, **kwargs):
        """
        Read the file as a remote file (e.g., from a URL).

        Returns
        -------
        object
            The contents of the file.
        """
        pass

    @log_entry
    def file_exist(self) -> bool:
        """
        Check if the file exists.

        Returns
        -------
        bool
            True if the file exists, False otherwise.
        """
        return self.file_path.exists()


class CsvFileReader(FileReader):
    """
    File reader for CSV files.

    Parameters
    ----------
    file_path : str
        Path to the CSV file or URL.
    is_url : bool, optional
        Whether the file_path is a URL (default is False).
    """

    def __init__(self, file_path: str, is_url: bool = False):
        """Initialize CsvFileReader."""
        super().__init__(file_path, is_url=is_url)

    @log_entry
    def read(self, **readkwargs) -> pd.DataFrame:
        """
        Read the CSV file and return its contents as a DataFrame.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the CSV reader.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the CSV data.
        """
        if self.is_url:
            data = self.read_as_remote_file(**readkwargs)
        else:
            data = self.read_as_local_file(**readkwargs)
        return data

    @log_entry
    def read_as_local_file(self, **readkwargs) -> pd.DataFrame:
        """
        Read the CSV file as a local file.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the CSV reader.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the CSV data.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")
        return read_csv(filepath=self.file_path, **readkwargs)

    @log_entry
    def read_as_remote_file(self, **readkwargs) -> pd.DataFrame:
        """
        Read the CSV file as a remote file (URL).

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the CSV reader.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the CSV data.
        """
        # Pandas can read online files with the same syntax
        return self.read_as_local_file(**readkwargs)


class JsonFileReader(FileReader):
    """
    File reader for JSON files.

    Parameters
    ----------
    file_path : str
        Path to the JSON file or URL.
    is_url : bool, optional
        Whether the file_path is a URL (default is False).
    """

    def __init__(self, file_path: str, is_url: bool = False):
        """Initialize JsonFileReader."""
        super().__init__(file_path, is_url=is_url)

    @log_entry
    def read(self, **readkwargs) -> dict:
        """
        Read the JSON file and return its contents as a dictionary.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the JSON loader.

        Returns
        -------
        dict
            Dictionary containing the JSON data.
        """
        if self.is_url:
            self.data = self.read_as_remote_file(**readkwargs)
        else:
            self.data = self.read_as_local_file(**readkwargs)
        self.is_already_read = True

    @log_entry
    def read_as_local_file(self, **readkwargs) -> dict:
        """
        Read the JSON file as a local file.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the JSON loader.

        Returns
        -------
        dict
            Dictionary containing the JSON data.

        Raises
        ------
        ImportError
            If the file is not a JSON file.
        FileNotFoundError
            If the file does not exist.
        """
        if not str(self.file_path).endswith(".json"):
            raise ImportError(f"The file {self.file_path} is not a JSON file.")
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")

        with open(self.file_path, "r") as file:
            data = json.load(file, **readkwargs)
        return data

    @log_entry
    def read_as_remote_file(self, **readkwargs) -> dict:
        """
        Read the JSON file as a remote file (URL).

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the requests.get call.

        Returns
        -------
        dict
            Dictionary containing the JSON data.
        """
        r = requests.get(self.file_path, **readkwargs)
        data = r.json()
        return data


class PickleFileReader(FileReader):
    """
    File reader for pickle files.

    Parameters
    ----------
    file_path : str
        Path to the pickle file.
    is_url : bool, optional
        Whether the file_path is a URL (default is False).
    """

    def __init__(self, file_path: str, is_url: bool = False):
        """Initialize PickleFileReader."""
        super().__init__(file_path, is_url=is_url)

    @log_entry
    def read(self, **readkwargs):
        """
        Read the pickle file and return its contents.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the pickle loader.

        Returns
        -------
        object
            The object loaded from the pickle file.
        """
        if self.is_url:
            self.data = self.read_as_remote_file(**readkwargs)
        else:
            self.data = self.read_as_local_file(**readkwargs)
        self.is_already_read = True

    @log_entry
    def read_as_local_file(self, **readkwargs):
        """
        Read the pickle file as a local file.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the pickle loader.

        Returns
        -------
        object
            The object loaded from the pickle file.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")

        with open(self.file_path, "rb") as inp:
            obj = pickle.load(inp, **readkwargs)

        return obj

    @log_entry
    def read_as_remote_file(self, **readkwargs):
        """
        Not implemented. Raises ReferenceError.

        Raises
        ------
        ReferenceError
            Always raised, as remote pickle reading is not implemented.
        """
        raise ReferenceError(
            "No remote pickle reader implemented. Please make an issue if you want this functionality."
        )


class ParquetFileReader(FileReader):
    """
    File reader for Parquet files.

    Parameters
    ----------
    file_path : str
        Path to the Parquet file or URL.
    is_url : bool, optional
        Whether the file_path is a URL (default is False).
    """

    def __init__(self, file_path: str, is_url: bool = False):
        """Initialize ParquetFileReader."""
        super().__init__(file_path, is_url=is_url)

    @log_entry
    def read(self, **readkwargs) -> pd.DataFrame:
        """
        Read the Parquet file and return its contents as a DataFrame.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the Parquet reader.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the Parquet data.
        """
        if self.is_url:
            data = self.read_as_remote_file(**readkwargs)
        else:
            data = self.read_as_local_file(**readkwargs)
        return data

    @log_entry
    def read_as_local_file(self, **readkwargs) -> pd.DataFrame:
        """
        Read the Parquet file as a local file.

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the Parquet reader.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the Parquet data.

        Raises
        ------
        FileNotFoundError
            If the file does not exist.
        """
        if not self.file_exist():
            raise FileNotFoundError(f"{self.file_path} is not a file.")
        return pd.read_parquet(self.file_path, **readkwargs)

    @log_entry
    def read_as_remote_file(self, **readkwargs) -> pd.DataFrame:
        """Read the Parquet file as a remote file (URL).

        Parameters
        ----------
        **readkwargs
            Additional keyword arguments to pass to the Parquet reader.

        Returns
        -------
        pandas.DataFrame
            DataFrame containing the Parquet data.
        """
        return pd.read_parquet(self.file_path, **readkwargs)


@log_entry
def read_csv(filepath: str, **kwargs) -> pd.DataFrame:
    """
    Read a CSV file into a DataFrame.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.
    **kwargs
        Additional keyword arguments to pass to pandas.read_csv.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the data from the CSV file.

    Raises
    ------
    EncodingWarning
        If the file is read as an empty DataFrame.
    """
    if "sep" not in kwargs:
        df = read_csv_with_flexible_seperator(filepath, **kwargs)
    else:
        df = pd.read_csv(filepath, **kwargs)

    # Sanity checks
    if df.empty:
        raise EncodingWarning(f"The file {filepath} is read as an empty DataFrame!")
    elif len(df.columns) == 1:
        logger.warning(f"The file {filepath} is read as a single-column DataFrame!")
    return df


@log_entry
def read_csv_with_flexible_seperator(filepath: str, **kwargs) -> pd.DataFrame:
    """
    Read a CSV file into a DataFrame, trying to infer the separator.

    Parameters
    ----------
    filepath : str
        Path to the CSV file.
    **kwargs
        Additional keyword arguments to pass to pandas.read_csv.

    Returns
    -------
    pandas.DataFrame
        DataFrame containing the data from the CSV file.

    Raises
    ------
    ValueError
        If the separator could not be determined.
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


@log_entry
def find_suitable_reader(filepath: Union[str, Path], is_url=False) -> FileReader:
    """
    Find a suitable file reader based on file extension.

    Determines the appropriate FileReader subclass (CSV, JSON, or Parquet)
    by examining the file extension and returns an instance of that reader.

    Parameters
    ----------
    filepath : str or Path
        Path to the file or URL.
    is_url : bool, optional
        Whether the filepath is a URL (default is False).

    Returns
    -------
    FileReader
        An instance of the appropriate FileReader subclass.

    Raises
    ------
    FileNotFoundError
        If the local file does not exist.
    NotImplementedError
        If no reader is available for the file extension.
    """

    # to Path
    if isinstance(filepath, str):
        filepath = Path(filepath)

    # 1. Check if the file exists
    if not filepath.exists():
        raise FileNotFoundError(f"{filepath} does not exist.")

    # 2. Get suffix
    suffix = filepath.suffix

    if suffix == ".csv":
        return CsvFileReader(filepath, is_url=is_url)
    elif suffix == ".json":
        return JsonFileReader(filepath, is_url=is_url)
    elif suffix == ".parquet":
        return ParquetFileReader(filepath, is_url=is_url)
    else:
        raise NotImplementedError(
            f"Could not read {filepath}, no reader implemented for file extension {suffix}"
        )
