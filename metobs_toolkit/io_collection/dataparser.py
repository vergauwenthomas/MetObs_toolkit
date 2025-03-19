import logging
import numpy as np
import pandas as pd
from metobs_toolkit.io_collection.filereaders import CsvFileReader
from metobs_toolkit.template import Template, MetobsTemplateError

logger = logging.getLogger(__file__)


class DataParser:
    """Parser class for data"""

    def __init__(self, datafilereader: CsvFileReader, template: Template):
        self.filereader = datafilereader
        self.template = template

        self.datadf = pd.DataFrame()  # metadata in formatted dataframe style

    def parse(self, **readkwargs):
        # Read the raw datafile
        rawdf = self.filereader.read(**readkwargs)
        # check template compatibility
        self.template._data_template_compatibility_test(rawdf.columns)

        # ---- format data ---------
        # 1. black label handling
        rawdf = self._rename_blacklabels(rawdf)

        # 2. Set datetime column
        rawdf = _create_datetime_column(rawdf, self.template)

        # 3. Add a dumy-name column for single-station cases
        if self.template._is_data_single_station():
            rawdf = self._add_single_station_name(rawdf)

        # 4. Convert wide into long structures
        if not self.template._is_data_long():
            rawdf = wide_to_long(
                df=rawdf, obstypename=self.template._get_wide_obstype()
            )

        # 5. rename name column
        rawdf = self._set_name(rawdf)

        # 6. rename columns to toolkit obstypes
        rawdf = self._columns_to_tlk_obstypes(rawdf)

        # 7. Subset to mapped columns
        rawdf = self._subset_to_mapped_columns(rawdf)

        # Set data attribute
        self.datadf = rawdf
        return self

    def _rename_blacklabels(self, rawdf):
        blacklist_remapper = self.template._apply_blacklist(
            columns=rawdf.columns, on_data=True
        )
        rawdf.rename(columns=blacklist_remapper, inplace=True)
        return rawdf

    def _add_single_station_name(self, rawdf):

        # check if there is a column indicating the name of the station that is mapped
        assumed_name_col = list(self.template._get_data_name_map().keys())[0]
        if assumed_name_col is None:
            rawdf["_dummy_name_column"] = (
                self.template._get_single_station_default_name()
            )
            # add it to the template
            self.template._set_dataname("_dummy_name_column")
        return rawdf

    def _set_name(self, rawdf):
        rawdf.rename(columns=self.template._get_data_name_map(), inplace=True)
        # Drop rows for which the name is nan
        todrop = rawdf[rawdf["name"].isnull()]
        if not todrop.empty:
            logger.warning(
                "Records with NaN as stationname are found in the datafile, these are skipped!"
            )
        rawdf = rawdf[~rawdf["name"].isnull()]

        # make sure the names are strings
        rawdf["name"] = rawdf["name"].astype(str)
        return rawdf

    def _columns_to_tlk_obstypes(self, rawdf):
        rawdf.rename(columns=self.template._get_obs_column_map(), inplace=True)
        return rawdf

    def _subset_to_mapped_columns(self, rawdf):
        mapped_columns = self.template._get_all_mapped_data_cols_in_tlk_space()

        unmapped_but_present = list(set(rawdf.columns) - set(mapped_columns))
        mapped_but_unpresent = list(set(mapped_columns) - set(rawdf.columns))
        if bool(unmapped_but_present):
            logger.warning(
                f"The following columns are present in the datafile, but not in the template! They are skipped!\n {unmapped_but_present}"
            )
        if bool(mapped_but_unpresent):
            logger.warning(
                f"The following variables are mapped in the template, but are not found in the data!\n {mapped_but_unpresent}"
            )

        return rawdf[mapped_columns]

    # ------------------------------------------
    #    Getters
    # ------------------------------------------
    def get_df(self) -> pd.DataFrame:
        return self.datadf


# ------------------------------------------
#    Helpers
# ------------------------------------------


def _create_datetime_column(df, template):
    """Use the template to construct a tz-naive "datetime" column."""

    template._check_if_datetime_is_mapped()

    if template.timestampinfo["datetimecolumn"] is not None:
        if not (template.timestampinfo["datetimecolumn"] in df.columns):
            raise MetobsTemplateError(
                f'The {template.timestampinfo["datetimecolumn"]} is not found in the columns of the data file: {df.columns}'
            )
        df = df.rename(columns={template.timestampinfo["datetimecolumn"]: "datetime"})
        try:
            df["datetime"] = pd.to_datetime(
                df["datetime"], format=template.timestampinfo["fmt"]
            )
        except Exception as e:
            raise MetobsTemplateError(
                "The timestamps could not be converted to datetimes, check the timestamp format(s) in your template."
            )

    else:
        # by date and time column
        if not (template.timestampinfo["time_column"] in df.columns):
            raise MetobsTemplateError(
                f'The {template.timestampinfo["time_column"]} is not found in the columns of the data file: {df.columns}'
            )
        if not (template.timestampinfo["date_column"] in df.columns):
            raise MetobsTemplateError(
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

        except Exception as e:
            raise MetobsTemplateError(
                "The timestamps could not be converted to datetimes, check the timestamp format(s) in your template."
            )
            # raise Exception('The timestamps could not be converted to datetimes, check the timestamp format(s) in your template. \n').with_traceback(e.__traceback__)

        df = df.drop(columns=["_date", "_time"])

    return df


def wide_to_long(df, obstypename):
    """Convert a wide dataframe to a long format.

    Convert a wide dataframe that represents obstype-observations to a long
    dataframe (=standard toolkit structure).

    Parameters
    ----------
    df : pandas.DataFrame()
        Wide dataframe with a "datetime" column.
    obstypename : str
        A MetObs obstype name.

    Returns
    -------
    longdf : pandas.DataFrame
        Long dataframe.
    template : dict
        Updated template dictionary.

    """
    # the df is assumed to have one datetime column, and the others represent
    # stations with their obstype values

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
