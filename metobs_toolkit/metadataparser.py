import logging
import numpy as np
import pandas as pd
from metobs_toolkit.filereaders import CsvFileReader
from metobs_toolkit.template import Template

logger = logging.getLogger(__file__)


class MetaDataParser:
    """Parser class for metadata"""

    def __init__(self, metadatafilereader: CsvFileReader, template: Template):
        self.filereader = metadatafilereader
        self.template = template

        # datadf with 'name' as index
        self.datadf = pd.DataFrame()  # metadata in formatted dataframe style

    def parse(self, **readkwargs):
        # Read in the raw metadata
        raw_metadata = self.filereader.read(**readkwargs)
        self.template._metadata_template_compatibility_test(raw_metadata.columns)

        # 1. black label handling
        metadf = self._rename_blacklabels(raw_metadata)

        # 2. name column handling (ID)
        metadf = self._set_name(metadf)

        # 3. Rename all columns to standards
        metadf = self._rename_raw_columns(metadf)

        # 4. Subset to mapped columns
        metadf = self._subset_to_mapped_columns(metadf)

        # set the metadf as data attribute
        self.datadf = metadf

    def _rename_blacklabels(self, rawdf):
        blacklist_mapper = self.template._apply_blacklist(
            columns=rawdf.columns, on_data=False
        )
        rawdf.rename(columns=blacklist_mapper, inplace=True)
        return rawdf

    def _set_name(self, rawdf):
        """Add a name column and set it as index"""
        if (self.template._is_data_single_station()) & (
            self.template.metadata_namemap["name"] is None
        ):
            rawdf["name"] = self.template.single_station_name
        else:
            rawdf.rename(columns=self.template._get_metadata_name_map(), inplace=True)

        # make sure the name column are strings
        rawdf["name"] = rawdf["name"].astype(str)

        # test uniqueness
        if not rawdf["name"].is_unique:
            duplicates = rawdf[rawdf["name"].duplicated(keep=False)]
            logger.warning(
                f'Duplicate names found in metadata: {duplicates["name"].tolist()}'
            )
            rawdf = rawdf.drop_duplicates(subset="name", keep="first")

        # set index
        rawdf.set_index("name", inplace=True)

        return rawdf

    def _rename_raw_columns(self, rawdf):
        metacolmap = self.template._get_metadata_column_map()
        rawdf.rename(columns=metacolmap, inplace=True)
        return rawdf

    def _subset_to_mapped_columns(self, rawdf):
        metacolmap = self.template._get_metadata_column_map()
        relev_columns = list(metacolmap.values())

        # drop name --> that is the index
        relev_columns = [col for col in relev_columns if col != "name"]
        return rawdf[relev_columns]

    def get_station_lon(self, stationname: str) -> float:
        if self._check_stationname_is_known(stationname=stationname):
            if "lon" in self.datadf.columns:
                return float(self.datadf.loc[stationname, "lon"])
            else:
                logger.warning(
                    f"No logitude is found for {stationname} in the metadata!"
                )
                return np.nan
        else:
            # station not in metadata
            return np.nan

    def get_station_lat(self, stationname: str) -> float:
        if self._check_stationname_is_known(stationname=stationname):
            if "lat" in self.datadf.columns:
                return float(self.datadf.loc[stationname, "lat"])
            else:
                logger.warning(
                    f"No latitude is found for {stationname} in the metadata!"
                )
                return np.nan
        else:
            # station not in metadata
            return np.nan

    def get_station_extra_metadata(self, stationname: str) -> dict:
        not_extra_columns = ["lat", "lon"]
        if self._check_stationname_is_known(stationname=stationname):
            extra_info = (
                self.datadf.loc[stationname]
                .drop(labels=not_extra_columns, errors="ignore")
                .to_dict()
            )
            return extra_info
        else:
            # station not in metadata
            return {}

    def _check_stationname_is_known(self, stationname: str) -> bool:
        if stationname in self.datadf.index:
            return True
        else:
            logger.warning(f"{stationname} is not found in the metadata!")
            return False

    # ------------------------------------------
    #    Getters
    # ------------------------------------------

    def get_df(self) -> pd.DataFrame:
        return self.datadf
