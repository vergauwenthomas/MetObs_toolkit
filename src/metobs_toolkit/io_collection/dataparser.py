import logging
import pandas as pd

from metobs_toolkit.io_collection.filereaders import FileReader

from metobs_toolkit.template import Template, MetObsTemplateError

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")


class DataParser:
    """
    Parser class for meteorological observation data.

    Parameters
    ----------
    datafilereader : FileReader
        An instance of FileReader to read the data file.
    template : Template
        An instance of Template to validate and map the data.
    """

    def __init__(self, datafilereader: FileReader, template: Template):
        self.filereader = datafilereader
        self.template = template
        self.datadf = pd.DataFrame()  # Metadata in formatted DataFrame style

    @log_entry
    def parse(self, **readkwargs) -> "DataParser":
        """
        Parse the data file and format it according to the template.

        Parameters
        ----------
        **readkwargs : dict
            Additional arguments to pass to the file reader.

        Returns
        -------
        DataParser
            The instance of DataParser with the parsed data.
        """

        if not isinstance(readkwargs, dict):
            raise TypeError("readkwargs must be a dictionary.")

        # Read the raw data file
        rawdf = self.filereader.read(**readkwargs)
        logger.debug("Raw data read successfully.")

        # Check template compatibility
        self.template._data_template_compatibility_test(rawdf.columns)

        # ---- Format data ---------
        # 1. Black label handling
        rawdf = self._rename_blacklabels(rawdf)

        # 2. Set datetime column
        rawdf = _create_datetime_column(rawdf, self.template)

        # 3. Add a dummy name column for single-station cases
        if self.template._is_data_single_station():
            rawdf = self._add_single_station_name(rawdf)

        # 4. Convert wide into long structures
        if not self.template._is_data_long():
            rawdf = wide_to_long(
                df=rawdf, obstypename=self.template._get_wide_obstype()
            )

        # 5. Rename name column
        rawdf = self._set_name(rawdf)

        # 6. Rename columns to toolkit observation types
        rawdf = self._columns_to_tlk_obstypes(rawdf)

        # 7. Subset to mapped columns
        rawdf = self._subset_to_mapped_columns(rawdf)

        # Set data attribute
        self.datadf = rawdf
        logger.info("Data parsing completed.")
        return self

    def _rename_blacklabels(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Rename blacklisted labels in the DataFrame.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The raw DataFrame to process.

        Returns
        -------
        pd.DataFrame
            The DataFrame with renamed columns.
        """

        if not isinstance(rawdf, pd.DataFrame):
            raise TypeError("rawdf must be a pandas DataFrame.")

        blacklist_remapper = self.template._apply_blacklist(
            columns=rawdf.columns, on_data=True
        )
        rawdf.rename(columns=blacklist_remapper, inplace=True)
        return rawdf

    def _add_single_station_name(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Add a dummy name column for single-station data.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The raw DataFrame to process.

        Returns
        -------
        pd.DataFrame
            The DataFrame with the added dummy name column.
        """

        if not isinstance(rawdf, pd.DataFrame):
            raise TypeError("rawdf must be a pandas DataFrame.")

        # Check if there is a column indicating the name of the station that is mapped
        assumed_name_col = list(self.template._get_data_name_map().keys())[0]
        if assumed_name_col is None:
            rawdf["_dummy_name_column"] = (
                self.template._get_single_station_default_name()
            )
            # Add it to the template
            self.template._set_dataname("_dummy_name_column")
        return rawdf

    def _set_name(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Rename the name column and drop rows with NaN names.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The raw DataFrame to process.

        Returns
        -------
        pd.DataFrame
            The DataFrame with renamed name column.
        """

        if not isinstance(rawdf, pd.DataFrame):
            raise TypeError("rawdf must be a pandas DataFrame.")

        rawdf.rename(columns=self.template._get_data_name_map(), inplace=True)
        # Drop rows for which the name is NaN
        todrop = rawdf[rawdf["name"].isnull()]
        if not todrop.empty:
            logger.warning(
                "Records with NaN as station name are found in the data file, these are skipped!"
            )
        rawdf = rawdf[~rawdf["name"].isnull()]

        # Make sure the names are strings
        rawdf["name"] = rawdf["name"].astype(str)
        return rawdf

    def _columns_to_tlk_obstypes(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Rename columns to toolkit observation types.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The raw DataFrame to process.

        Returns
        -------
        pd.DataFrame
            The DataFrame with renamed columns.
        """

        if not isinstance(rawdf, pd.DataFrame):
            raise TypeError("rawdf must be a pandas DataFrame.")

        rawdf.rename(columns=self.template._get_obs_column_map(), inplace=True)
        return rawdf

    def _subset_to_mapped_columns(self, rawdf: pd.DataFrame) -> pd.DataFrame:
        """
        Subset the DataFrame to mapped columns.

        Parameters
        ----------
        rawdf : pd.DataFrame
            The raw DataFrame to process.

        Returns
        -------
        pd.DataFrame
            The DataFrame with only mapped columns.
        """

        if not isinstance(rawdf, pd.DataFrame):
            raise TypeError("rawdf must be a pandas DataFrame.")

        mapped_columns = self.template._get_all_mapped_data_cols_in_tlk_space()

        unmapped_but_present = list(set(rawdf.columns) - set(mapped_columns))
        mapped_but_unpresent = list(set(mapped_columns) - set(rawdf.columns))
        if bool(unmapped_but_present):
            logger.warning(
                f"The following columns are present in the data file, but not in the template! They are skipped!\n {unmapped_but_present}"
            )
        if bool(mapped_but_unpresent):
            logger.warning(
                f"The following variables are mapped in the template, but are not found in the data!\n {mapped_but_unpresent}"
            )
        # mapped and present
        mapped_and_present = list(set(mapped_columns).intersection(set(rawdf.columns)))
        return rawdf[mapped_and_present]

    @log_entry
    def get_df(self) -> pd.DataFrame:
        """Get the parsed DataFrame."""
        logger.info(f"Entering get_df method of {self}.")
        return self.datadf


def _create_datetime_column(df: pd.DataFrame, template: Template) -> pd.DataFrame:
    """
    Use the template to construct a timezone-naive "datetime" column.

    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame to process.
    template : Template
        The template to use for datetime construction.

    Returns
    -------
    pd.DataFrame
        The DataFrame with a "datetime" column.

    Raises
    ------
    TypeError
        If the input types are incorrect.
    MetobsTemplateError
        If required columns are missing or datetime conversion fails.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")
    if not isinstance(template, Template):
        raise TypeError("template must be an instance of Template.")

    template._check_if_datetime_is_mapped()

    if template.timestampinfo["datetimecolumn"] is not None:
        if not (template.timestampinfo["datetimecolumn"] in df.columns):
            raise MetObsTemplateError(
                f'The {template.timestampinfo["datetimecolumn"]} is not found in the columns of the data file: {df.columns}'
            )
        df = df.rename(columns={template.timestampinfo["datetimecolumn"]: "datetime"})

        # NOTE: in parquet files, the datetimecolumn can be a pandas datetime column with tz info.
        # It becomes problematic when the tz in the data, is not the same as the tz set in the template

        # Check if datetime column is already in datetime format
        if pd.api.types.is_datetime64_any_dtype(df["datetime"]):
            # If already datetime, convert to timezone-naive
            if hasattr(df["datetime"].dtype, "tz") and df["datetime"].dt.tz is not None:
                data_tz = df["datetime"].dt.tz
                template_tz = template._get_tz()
                if str(data_tz) != template_tz:
                    raise MetObsTemplateError(
                        f"The timezone of the data ({data_tz}) does not match the template ({template_tz}). Please update the template."
                    )
            return df

        try:
            df["datetime"] = pd.to_datetime(
                df["datetime"], format=template.timestampinfo["fmt"]
            )
        except Exception:
            raise MetObsTemplateError(
                "The timestamps could not be converted to datetimes, check the timestamp format(s) in your template."
            )

    else:
        # By date and time column
        if not (template.timestampinfo["time_column"] in df.columns):
            raise MetObsTemplateError(
                f'The {template.timestampinfo["time_column"]} is not found in the columns of the data file: {df.columns}'
            )
        if not (template.timestampinfo["date_column"] in df.columns):
            raise MetObsTemplateError(
                f'The {template.timestampinfo["date_column"]} is not found in the columns of the data file: {df.columns}'
            )

        df = df.rename(
            columns={
                template.timestampinfo["time_column"]: "_time",
                template.timestampinfo["date_column"]: "_date",
            }
        )
        try:
            df["datetime"] = pd.to_datetime(
                df["_date"] + " " + df["_time"], format=template.timestampinfo["fmt"]
            )

        except Exception:
            raise MetObsTemplateError(
                "The timestamps could not be converted to datetimes, check the timestamp format(s) in your template."
            )

        df = df.drop(columns=["_date", "_time"])

    return df


@log_entry
def wide_to_long(df: pd.DataFrame, obstypename: str) -> pd.DataFrame:
    """
    Convert a wide DataFrame to a long format.

    Parameters
    ----------
    df : pd.DataFrame
        Wide DataFrame with a "datetime" column.
    obstypename : str
        A MetObs observation type name.

    Returns
    -------
    pd.DataFrame
        Long DataFrame in the standard toolkit structure.

    Raises
    ------
    TypeError
        If the input types are incorrect.
    """

    if not isinstance(df, pd.DataFrame):
        raise TypeError("df must be a pandas DataFrame.")
    if not isinstance(obstypename, str):
        raise TypeError("obstypename must be a string.")

    # The DataFrame is assumed to have one datetime column, and the others represent
    # stations with their observation type values
    stationnames = df.columns.to_list()
    stationnames.remove("datetime")

    longdf = pd.melt(
        df,
        id_vars=["datetime"],
        value_vars=stationnames,
        var_name="name",
        value_name=obstypename,
    )
    return longdf
