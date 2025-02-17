from typing import Literal, Tuple
import os
from pathlib import Path
import logging
import copy
import pickle
import pandas as pd
import numpy as np
from matplotlib.pyplot import Axes
import concurrent.futures


from metobs_toolkit.template import Template, update_known_obstype_with_original_data
from metobs_toolkit.station import Station
from metobs_toolkit.metadataparser import MetaDataParser
from metobs_toolkit.dataparser import DataParser
from metobs_toolkit.filereaders import CsvFileReader, PickleFileReader
from metobs_toolkit.site import Site
from metobs_toolkit.sensordata import SensorData
from metobs_toolkit.backend_collection.printing import dataset_string_repr
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_timedelta_arg,
    fmt_datetime_arg,
)
from metobs_toolkit.obstypes import tlk_obstypes
from metobs_toolkit.obstype_modeldata import ModelObstype
from metobs_toolkit.backend_collection.timeseries_plotting import (
    create_axes,
    create_station_color_map,
)
from metobs_toolkit.qc_collection import toolkit_buddy_check
from metobs_toolkit.backend_collection.docstring_wrapper import copy_doc
from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.modeltimeseries import ModelTimeSeries


from metobs_toolkit.modeldata import (
    GeeStaticDataset,
    GeeDynamicDataset,
    default_datasets,
)
from metobs_toolkit.gee_api import connect_to_gee

logger = logging.getLogger(__file__)


class Dataset:

    def __init__(self):

        self._stations = []  # stationname: Station
        # dictionary storing present observationtypes
        self._obstypes = copy.copy(tlk_obstypes)  # init with all tlk obstypes
        # self._settings = copy.deepcopy(Settings())

        # Gaps are stored as a list of Gap()
        # self.gaps = None

        # Template
        self._template = Template()

        # GEE datasets defenitions
        # self._gee_datasets = copy.deepcopy(default_datasets)

    # ------------------------------------------
    #    specials
    # ------------------------------------------

    def __str__(self):
        """Represent as text."""
        return dataset_string_repr(self)

    def __repr__(self):
        """Info representation."""
        class_name = type(self).__name__
        return f"Instance of {class_name} at {hex(id(self))}"

    def copy(self, deep=True):
        if deep:
            return copy.deepcopy(self)
        return copy.copy(self)

    @property
    def stations(self):
        return self._stations

    @property
    def obstypes(self):
        return self._obstypes

    @property
    def settings(self):
        return self._settings

    @property
    def template(self):
        return self._template

    @property
    def gee_datasets(self):
        return self._gee_datasets

    @stations.setter
    def stations(self, stationlist: list) -> None:
        self._stations = stationlist

    @obstypes.setter
    def obstypes(self, obstypesdict: dict) -> None:
        self._obstypes = obstypesdict

    @property
    def df(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.df.reset_index()
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = pd.concat(concatlist)
        combdf.sort_index(inplace=True)
        return combdf

    @property
    def outliersdf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.outliersdf.reset_index()
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = pd.concat(concatlist)
        combdf.sort_index(inplace=True)
        return combdf

    @property
    def gapsdf(self) -> pd.DataFrame:
        # TODO
        pass

    @property
    def metadf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            concatlist.append(sta.metadf)
        return pd.concat(concatlist).sort_index()

    @property
    def start_datetime(self):
        return min([sta.start_datetime for sta in self.stations])

    @property
    def end_datetime(self):
        return max([sta.end_datetime for sta in self.stations])

    # ------------------------------------------
    #   Extracting data
    # ------------------------------------------

    def get_station(self, stationname: str) -> Station:
        # create lookup by name
        stationlookup = {sta.name: sta for sta in self.stations}
        try:
            station = stationlookup[stationname]
        except KeyError:
            raise MetObsStationNotFound(f"{stationname} is not found in {self}.")

        return station

    def get_info(self, printout=True):

        df = self.df

        present_obstypes = list(df.index.get_level_values("obstype").unique())

        infostr = " --- General Info --- \n\n"
        infostr += "Dataset instance with:\n"
        infostr += f"  *{len(self.stations)} number of stations\n"
        infostr += f"  *{len(present_obstypes)} types of sensor data are present.\n"
        infostr += f'  *Observations from {df.index.get_level_values("datetime").min()} -> {df.index.get_level_values("datetime").max()}\n'

        infostr += "\n --- Observational info ---\n\n"
        # TODO
        pass

        infostr += "\n --- Metadata info ---\n\n"
        # TODO
        pass

        if printout:
            print(infostr)
        else:
            return infostr

    def resample(
        self,
        target_freq,
        shift_tolerance=pd.Timedelta("4min"),
        origin=None,
        direction="nearest",
    ):
        target_freq = fmt_timedelta_arg(target_freq)
        shift_tolerance = fmt_timedelta_arg(shift_tolerance)
        for sta in self.stations:
            sta.resample(
                target_freq=target_freq,
                shift_tolerance=shift_tolerance,
                origin=origin,
                direction=direction,
            )

    # ------------------------------------------
    #   Reading/writing data
    # ------------------------------------------

    def import_gee_data_from_file(
        self,
        filepath: str,
        geedynamicdataset: GeeDynamicDataset,
        force_update=True,
        _force_from_dataframe=None,
    ) -> pd.DataFrame:
        if _force_from_dataframe is None:
            # Reading the data
            reader = CsvFileReader(file_path=filepath)
            data = reader.read_as_local_file()
            force_update = True

            # 1. format the data
            totaldf = geedynamicdataset._format_gee_df_structure(data)

            # 2. convert units
            totaldf = geedynamicdataset._convert_units(totaldf)
            # this totaldf will be returned
        else:
            totaldf = _force_from_dataframe

        # 3. Subset to known obstypes
        known_obstypes = list(geedynamicdataset.modelobstypes.keys())
        cols_to_skip = list(set(totaldf.columns) - set(known_obstypes))
        if bool(cols_to_skip):
            logger.warning(
                f"The folowing columns in the GEE datafile are not present in the known modelobstypes of {geedynamicdataset}: {cols_to_skip}"
            )

        known_and_present = set(totaldf.columns) & set(known_obstypes)
        df = totaldf[list(known_and_present)]

        # 4. add it to the stations
        for sta in self.stations:
            # subset to station data
            if sta.name not in df.index.get_level_values("name"):
                continue  # no data present for the station
            stadf = df.xs(key=sta.name, level="name", drop_level=True)

            # create an object for each variable
            for col in stadf.columns:
                modeltimeseries = ModelTimeSeries(
                    site=sta.site,
                    datarecords=stadf[col].to_numpy(),
                    timestamps=stadf.index.to_numpy(),
                    obstype=geedynamicdataset.modelobstypes[col],
                    timezone="UTC",
                    modelname=geedynamicdataset.name,
                    modelvariable=geedynamicdataset.modelobstypes[col].get_modelband(),
                )
                # add it to the station
                sta.add_to_modeldata(modeltimeseries, force_update=force_update)

        return totaldf

    def save_dataset_to_pkl(
        self, target_folder, filename="saved_dataset.pkl", overwrite=False
    ):

        # Check if the targetfolder exist
        if not Path(target_folder).is_dir():
            raise FileNotFoundError("{target_folder} does not exist")

        # Check filetype extension
        if not filename.endswith(".pkl"):
            filename += ".pkl"

        # Check if target file exist
        target_path = Path(target_folder).joinpath(filename)
        if target_path.exists():
            if overwrite:
                target_path.unlink()
            else:
                raise FileExistsError(
                    f"The file {target_path} already exist. Remove it or set overwrite=True"
                )

        # dump to pickle
        with open(target_path, "wb") as outp:
            pickle.dump(self, outp, pickle.HIGHEST_PROTOCOL)

    def import_data_from_file(
        self,
        input_data_file: str | Path = None,
        input_metadata_file: str | Path = None,
        template_file: str | Path = None,
        freq_estimation_method: Literal["highest", "median"] = "median",
        freq_estimation_simplify_tolerance: str | pd.Timedelta = "2min",
        origin_simplify_tolerance: str | pd.Timedelta = "5min",
        timestamp_tolerance: str | pd.Timedelta = "4min",
        kwargs_data_read: dict = {},
        kwargs_metadata_read: dict = {},
        templatefile_is_url: bool = False,
    ):
        """Read observations from a csv file and fill the Dataset.


        The input data (and metadata) are interpreted by using a template
        (JSON file).

        In order to locate gaps, an ideal set of timestamps is expected. This
        set of timestamps is computed for each station separatly by:

         * Assuming a constant frequency. This frequency is estimated by using
           a freq_estimation_method. If multiple observationtypes are present,
           the assumed frequency is the highest estimated frequency among
           the differnt observationtypes. To simplify the estimated frequency a
           freq_estimation_simplify_error can be specified.
         * A start timestamp (origin) is found for each station. If multiple
           observationtypes are present,
           the start timestamp is the first timestamp among
           the different observationtypes. The start
           timestamp can be simplified by specifying an origin_simplify_tolerance.
         * The last timestamp is found for each station by taking the timestamp
           which is closest and smaller than the latest timestamp found of a station,
           and is an element of the ideal set of timestamps.

        Each present observation record is linked to a timestamp of this ideal set,
        by using a 'nearest' merge. If the timedifference is smaller than the
        timestamp_tolerance, the ideal timestamp is used. Otherwise, the timestamp
        will be interpreted as a (part of a) gap.


        The Dataset attributes are set and the following checks are executed:
                * Duplicate check
                * Invalid input check
                * Find gaps


        Parameters
        ----------
        input_data_file : string, optional
            Path to the input data file with observations. If None, the input
            data path in the settings is used. If provided, the
            template_file must be provided as well. The default is None.
        input_metadata_file : string, optional
            Path to the input metadata file. If None, the input metadata path
            in the settings is used. The default is None
        template_file : string, optional
            Path to the template (JSON) file to be used on the observations
            and metadata. If None, the template path in the settings is used.
            If provided, the input_data_file must be provided as well. The
            default is None.
        freq_estimation_method : 'highest' or 'median', optional
            Select which method to use for the frequency estimation. If
            'highest', the highest appearing frequency is used. If 'median', the
            median of the appearing frequencies is used. The default is 'highest'.
        freq_estimation_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation
            in time to form a simplified frequency estimation.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '2min' (2 minutes).
        origin_simplify_tolerance : Timedelta or str, optional
            The tolerance string or object representing the maximum translation
            in time to apply on the start timestamp to create a simplified timestamp.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '5min' (5 minutes).
        timestamp_tolerance : Timedelta or str, optional
            The tolerance string or object represents the maximum translation
            in time to apply on a timestamp for conversion to an ideal set of timestamps.
            Ex: '5min' is 5 minutes, '1H', is one hour. A zero-tolerance (thus no
            simplification) can be set by '0min'. The default is '4min' (4 minutes).
        kwargs_data_read : dict, optional
            Keyword arguments are collected in a dictionary to pass to the
            pandas.read_csv() function on the data file. The default is {}.
        kwargs_metadata_read : dict, optional
            Keyword arguments are collected in a dictionary to pass to the
            pandas.read_csv() function on the metadata file. The default is {}.
        templatefile_is_url : bool, optional
            If the path to the template file, is a url to an online template file,
            set templatefile_is_url to True. If False, the template_file is
            interpreted as a path.

        Returns
        -------
        None.

        Note
        --------
        In practice, the default arguments will be sufficient for most applications.

        Note
        --------
        If options are present in the template, these will have priority over the arguments of this function.

        Warning
        ---------
        All CSV data files must be in *UTF-8 encoding*. For most CSV files,
        this condition has already been met. To make sure, in Microsoft Excel (or
        similar), you can specify to export as **`CSV UTF-8`**. If you
        encounter an error, mentioning a `"/ueff..."` tag in a CSV file, it is
        often solved by converting the CSV to UTF-8.

        See Also
        --------
        update_settings: Update the (file paths) settings of a Dataset.
        import_only_metadata_from_file: Import metadata without observational data.
        import_dataset: Import a dataset from a pkl file.

        Examples
        --------

        >>> import metobs_toolkit
        >>>
        >>> # Import data into a Dataset
        >>> dataset = metobs_toolkit.Dataset()
        >>> dataset.update_file_paths(
        ...                         input_data_file=metobs_toolkit.demo_datafile,
        ...                         input_metadata_file=metobs_toolkit.demo_metadatafile,
        ...                         template_file=metobs_toolkit.demo_template,
        ...                         )
        >>> dataset.import_data_from_file()

        """
        # Special argschecks
        freq_estimation_simplify_tolerance = fmt_timedelta_arg(
            freq_estimation_simplify_tolerance
        )
        origin_simplify_tolerance = fmt_timedelta_arg(origin_simplify_tolerance)
        timestamp_tolerance = fmt_timedelta_arg(timestamp_tolerance)

        # check input file args
        if (input_data_file is not None) & (template_file is None):
            raise MetObsMissingFile(
                f"A input_data_file is provided, but the template_file is missing."
            )
        if (input_data_file is None) & (template_file is not None):
            raise MetObsMissingFile(
                f"A template_file is provided, but the input_data_file is missing."
            )

        assert template_file is not None, "No templatefile is specified."

        # Read template
        logger.info(f"Reading the templatefile")
        self.template.read_template_from_file(
            jsonpath=template_file, templatefile_is_url=templatefile_is_url
        )

        # Read Datafile
        dataparser = DataParser(
            datafilereader=CsvFileReader(file_path=input_data_file),
            template=self.template,
        )
        dataparser.parse(**kwargs_data_read)  # read and parse to a dataframe

        # Read Metadata
        if input_metadata_file is not None:
            use_metadata = True
            metadataparser = MetaDataParser(
                metadatafilereader=CsvFileReader(file_path=input_metadata_file),
                template=self.template,
            )
            metadataparser.parse(**kwargs_metadata_read)
        else:
            logger.info("No metadatafile is provided.")
            use_metadata = False
            metadataparser = None  # will not be used

        # Add original columnname and units to the known obstypes
        self.obstypes = update_known_obstype_with_original_data(
            known_obstypes=self.obstypes, template=self.template
        )

        # Construct Stations
        stations = createstations(
            data_parser=dataparser,
            metadata_parser=metadataparser,
            use_metadata=use_metadata,
            known_obstypes=self.obstypes,
            timezone=self.template.tz,
            freq_estimation_method=freq_estimation_method,
            freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
            origin_simplify_tolerance=origin_simplify_tolerance,
            timestamp_tolerance=timestamp_tolerance,
        )

        # Set stations attribute
        self.stations = stations

    # ------------------------------------------
    #    Plotting
    # ------------------------------------------

    def make_plot(
        self,
        obstype: str = "temp",
        colorby: Literal["station", "label"] = "label",
        show_modeldata=False,
        show_outliers=True,
        show_gaps=True,
        ax=None,
        figkwargs: dict = {},
        title: str | None = None,
    ) -> Axes:

        if ax is None:
            ax = create_axes(**figkwargs)

        if title is None:
            title = f"{obstype} data of the dataset"

        # Set colors scheme
        if colorby == "station":
            # create color defenitions
            colormap = create_station_color_map(
                catlist=[sta.name for sta in self.stations]
            )
        elif colorby == "label":
            colormap = {sta.name: None for sta in self.stations}
        else:
            raise ValueError(
                f'colorby is either "station" or "label" but not {colorby}'
            )

        for sta in self.stations:
            ax = sta.make_plot(
                obstype=obstype,
                colorby=colorby,
                linecolor=colormap[sta.name],
                show_modeldata=show_modeldata,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
                figkwargs=figkwargs,
                title=title,
            )

        return ax

    def make_gee_plot(
        self,
        geedataset: GeeStaticDataset | GeeDynamicDataset = default_datasets["lcz"],
        timeinstance=None,
        modelobstype: str = None,
        outputfolder=os.getcwd(),
        filename="gee_plot.html",
        save=False,
        vmin=None,
        vmax=None,
        overwrite=False,
    ):
        # check model type
        if isinstance(geedataset, GeeStaticDataset):
            kwargs = dict(
                outputfolder=outputfolder,
                filename=filename,
                save=save,
                vmin=vmin,
                vmax=vmax,
                overwrite=overwrite,
            )

        elif isinstance(geedataset, GeeDynamicDataset):
            timeinstance = fmt_datetime_arg(timeinstance, input_tz="UTC")
            if modelobstype not in geedataset.modelmobstypes:
                raise MetObsObstypeNotFound(
                    f"{modelobstype} is not a known modelobstype of {geedataset}. These are knonw: {geedataset.modelobstypes} "
                )

            kwargs = dict(
                timeinstance=timeinstance,
                modelobstype=modelobstype,
                outputfolder=outputfolder,
                filename=filename,
                save=save,
                vmin=vmin,
                vmax=vmax,
                overwrite=overwrite,
            )

        return geedataset.make_gee_plot(metadf=self.metadf, **kwargs)

    # ------------------------------------------
    #    Gee extracting
    # ------------------------------------------

    def get_static_gee_point_data(
        self,
        geestaticdataset: GeeStaticDataset,
        update_stations: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """simple but slow option is to loop over all stations and get the lcz,
        but this requires N-API calls (N number of stations).

        Faster: construct the metadf with all stations, and get the lcs from one api call
        """
        if not isinstance(geestaticdataset, GeeStaticDataset):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        # initialize gee api
        if initialize_gee:
            connect_to_gee()

        geedf = geestaticdataset.extract_static_point_data(self.metadf)

        # update the station attributes
        varname = geestaticdataset.name
        if update_stations:
            for staname, geedict in geedf.to_dict(orient="index").items():
                self.get_station(staname).site.set_geedata(varname, geedict[varname])
        return geedf

    def get_static_gee_buffer_fraction_data(
        self,
        geestaticdataset: GeeStaticDataset,
        buffers=[100],
        aggregate=False,
        update_stations: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """simple but slow option is to loop over all stations and get the lcz,
        but this requires N-API calls (N number of stations).

        Faster: construct the metadf with all stations, and get the lcs from one api call
        """
        if not isinstance(geestaticdataset, GeeStaticDataset):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdataset)}"
            )

        # initialize gee api
        if initialize_gee:
            connect_to_gee()

        dflist = []
        for bufferradius in buffers:
            geedf = geestaticdataset.extract_static_buffer_frac_data(
                metadf=self.metadf, bufferradius=bufferradius, agg_bool=aggregate
            )
            dflist.append(geedf)

        geedf = pd.concat(dflist)

        # update the station attributes
        if update_stations:
            for staname in geedf.index.get_level_values("name").unique():
                asdict = geedf.loc[staname].to_dict(orient="index")
                for radius, fractions in asdict.items():
                    self.get_station(staname).site.set_gee_buffered_frac_data(
                        buffer=radius, data=fractions
                    )

        return geedf

    def get_lcz(
        self, update_stations: bool = True, initialize_gee: bool = True
    ) -> pd.DataFrame:
        return self.get_static_gee_point_data(
            default_datasets["lcz"],
            update_stations=update_stations,
            initialize_gee=initialize_gee,
        )

    def get_altitude(
        self, update_stations: bool = True, initialize_gee: bool = True
    ) -> pd.DataFrame:
        return self.get_static_gee_point_data(
            default_datasets["altitude"],
            update_stations=update_stations,
            initialize_gee=initialize_gee,
        )

    def get_landcover_fractions(
        self,
        buffers=[100],
        aggregate=False,
        update_stations: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:

        return self.get_static_gee_buffer_fraction_data(
            geestaticdataset=default_datasets["worldcover"],
            buffers=buffers,
            aggregate=aggregate,
            update_stations=update_stations,
            initialize_gee=initialize_gee,
        )

    def get_gee_timeseries_data(
        self,
        geedynamicdataset: GeeDynamicDataset,
        startdt_utc=None,
        enddt_utc=None,
        target_obstypes=["temp"],
        get_all_bands=False,
        drive_filename=None,
        drive_folder="gee_timeseries_data",
        force_direct_transfer=False,
        force_to_drive=False,
    ):

        # Check geedynamic dataset
        if not isinstance(geedynamicdataset, GeeDynamicDataset):
            raise ValueError(
                f"geedynamicdataset should be an isntance of GeeDynamicDataset, not {type(geedynamicdataset)}"
            )

        # Format datetime arguments
        if startdt_utc is None:
            startdt_utc = self.start_datetime.tz_convert("UTC")
        else:
            startdt_utc = fmt_datetime_arg(startdt_utc, input_tz="UTC")

        if enddt_utc is None:
            enddt_utc = self.end_datetime.tz_convert("UTC")
        else:
            enddt_utc = fmt_datetime_arg(enddt_utc, input_tz="UTC")

        # check if target_obstypes are mapped to bands
        for obst in target_obstypes:
            if obst not in geedynamicdataset.modelobstypes.keys():
                raise MetObsMetadataNotFound(
                    f"{obst} is not a known modelobstype of {geedynamicdataset}."
                )

        # create specific name for the file that might be written to Drive
        if drive_filename is None:
            drive_filename = f"{geedynamicdataset.name}_timeseries_data_of_full_dataset_{len(self.stations)}_stations.csv"

        df = geedynamicdataset.extract_timeseries_data(
            metadf=self.metadf,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            obstypes=target_obstypes,
            get_all_bands=get_all_bands,
            drive_filename=drive_filename,
            drive_folder=drive_folder,
            force_direct_transfer=force_direct_transfer,
            force_to_drive=force_to_drive,
        )
        if df is None:
            print("No data is returned by the GEE api request.")
            logger.warning("No data is returned by the GEE api request.")
            return

        # Import the data and add it to the Stations

        # Note: all functionallity of this part is already implemented in
        # import_gee_data_from_file, except for the actual reading/formatting and unit converting.
        # Thus we call that function, by providing it to _force_from_dataframe, these steps are skipped.
        _ = self.import_gee_data_from_file(
            filepath=None,
            geedynamicdataset=geedynamicdataset,
            force_update=True,
            _force_from_dataframe=df,
        )
        return df

    # ------------------------------------------
    #    QC
    # ------------------------------------------
    @copy_doc(Station.gross_value_check)
    def gross_value_check(
        self,
        target_obstype: str = "temp",
        lower_threshold: float = -15.0,
        upper_threshold: float = 39.0,
        use_mp: bool = False,
    ):

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            lower_threshold=lower_threshold,
            upper_threshold=upper_threshold,
        )
        if use_mp:
            # Use multiprocessing generatore (parralelization)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_grossvalue_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            # Use regular generator
            self.stations = list(map(_qc_grossvalue_generatorfunc, func_feed_list))

    @copy_doc(Station.persistence_check)
    def persistence_check(
        self,
        target_obstype: str = "temp",
        timewindow: str | pd.Timedelta = pd.Timedelta("60min"),
        min_records_per_window: int = 5,
        use_mp: bool = True,
    ):
        timewindow = fmt_timedelta_arg(timewindow)
        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            timewindow=timewindow,
            min_records_per_window=min_records_per_window,
        )
        if use_mp:
            # Use multiprocessing generatore (parralelization)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_persistence_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            # Use regular generator
            self.stations = list(map(_qc_persistence_generatorfunc, func_feed_list))

    @copy_doc(Station.repetitions_check)
    def repetitions_check(
        self,
        target_obstype: str = "temp",
        max_N_repetitions: int = 5,
        use_mp: bool = False,
    ):
        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            max_N_repetitions=max_N_repetitions,
        )

        if use_mp:
            # Use multiprocessing generatore (parralelization)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_repetitions_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            # Use regular generator
            self.stations = list(map(_qc_repetitions_generatorfunc, func_feed_list))

    @copy_doc(Station.step_check)
    def step_check(
        self,
        target_obstype: str = "temp",
        max_increase_per_second: int | float = 8.0 / 3600.0,
        max_decrease_per_second: int | float = -10.0 / 3600.0,
        use_mp: bool = True,
    ):

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

        if use_mp:
            # Use multiprocessing generatore (parralelization)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(_qc_step_generatorfunc, func_feed_list)
            self.stations = list(stationgenerator)
        else:
            # Use regular generator
            self.stations = list(map(_qc_step_generatorfunc, func_feed_list))

    @copy_doc(Station.window_variation_check)
    def window_variation_check(
        self,
        target_obstype: str = "temp",
        timewindow: pd.Timedelta = pd.Timedelta("1h"),
        min_records_per_window: int = 3,
        max_increase_per_second: int | float = 0.0022,  # 8./3600
        max_decrease_per_second: int | float = -0.0027,  # -10/3600
        use_mp=True,
    ):

        timewindow = fmt_timedelta_arg(timewindow)

        func_feed_list = _create_qc_arg_set(
            dataset=self,
            target_obstype=target_obstype,
            timewindow=timewindow,
            min_records_per_window=min_records_per_window,
            max_increase_per_second=max_increase_per_second,
            max_decrease_per_second=max_decrease_per_second,
        )

        if use_mp:
            # Use multiprocessing generatore (parralelization)
            with concurrent.futures.ProcessPoolExecutor() as executor:
                stationgenerator = executor.map(
                    _qc_window_var_generatorfunc, func_feed_list
                )
            self.stations = list(stationgenerator)
        else:
            # Use regular generator
            self.stations = list(map(_qc_window_var_generatorfunc, func_feed_list))

    @copy_doc(toolkit_buddy_check)
    def buddy_check(
        self,
        target_obstype: str = "temp",
        buddy_radius: int | float = 10000,
        min_sample_size: int = 4,
        max_alt_diff: int | float = None,
        min_std: int | float = 1.0,
        std_threshold: int | float = 3.1,
        N_iter: int = 2,
        instantanious_tolerance: pd.Timedelta = pd.Timedelta("4min"),
        lapserate: float | None = None,  # -0.0065
        use_mp: bool = True,
    ):

        instantanious_tolerance = fmt_timedelta_arg(instantanious_tolerance)
        if lapserate is not None:
            # test if altitude data is available
            if not all([sta.site.flag_has_altitude() for sta in self.stations]):
                raise MetObsMetadataNotFound(
                    "Not all stations have altitude data, lapserate correction could not be applied."
                )

        qc_kwargs = dict(
            obstype=target_obstype,
            buddy_radius=buddy_radius,
            min_sample_size=min_sample_size,
            max_alt_diff=max_alt_diff,
            min_std=min_std,
            std_threshold=std_threshold,
            N_iter=N_iter,
            instantanious_tolerance=instantanious_tolerance,
            lapserate=lapserate,
            use_mp=use_mp,
        )

        outlierslist, timestamp_map = toolkit_buddy_check(dataset=self, **qc_kwargs)
        # outlierslist is a list of tuples (stationname, datetime) that are outliers
        # timestamp_map is a dict with keys the stationname and values a series to map the syncronized
        # timestamps to the original timestamps

        # convert to a dataframe
        df = pd.DataFrame(data=outlierslist, columns=["name", "datetime"])

        # update all the sensordata
        for station in self.stations:
            if target_obstype in station.obsdata.keys():
                # Get the sensordata object
                sensorddata = station.obsdata[target_obstype]

                # get outlier datetimeindex
                outldt = pd.DatetimeIndex(df[df["name"] == station.name]["datetime"])

                # convert to original timestamps
                dtmap = timestamp_map[station.name]
                outldt = outldt.map(dtmap)

                # update the sensordata
                sensorddata._update_outliers(
                    qccheckname="buddy_check",
                    outliertimestamps=outldt,
                    check_kwargs=qc_kwargs,
                )


def _qc_grossvalue_generatorfunc(input):
    station, kwargs = input  # first element is station, the reset ar kwargs
    station.gross_value_check(**kwargs)
    return station


def _qc_persistence_generatorfunc(input):
    station, kwargs = input  # first element is station, the reset ar kwargs
    station.persistence_check(**kwargs)
    return station


def _qc_repetitions_generatorfunc(input):
    station, kwargs = input  # first element is station, the reset ar kwargs
    station.repetitions_check(**kwargs)
    return station


def _qc_step_generatorfunc(input):
    station, kwargs = input  # first element is station, the reset ar kwargs
    station.step_check(**kwargs)
    return station


def _qc_window_var_generatorfunc(input):
    station, kwargs = input  # first element is station, the reset ar kwargs
    station.window_variation_check(**kwargs)
    return station


# ------------------------------------------
#    Helping functions
# ------------------------------------------


def _create_qc_arg_set(dataset, **qckwargs):
    return [([sta, dict(**qckwargs)]) for sta in dataset.stations]


def createstations(
    data_parser: DataParser,
    metadata_parser: MetaDataParser,
    use_metadata: bool,
    known_obstypes: dict,
    timezone: str,
    freq_estimation_method: str,  #'highest' | 'median',
    freq_estimation_simplify_tolerance: pd.Timedelta | str,
    origin_simplify_tolerance: pd.Timedelta | str,
    timestamp_tolerance: pd.Timedelta | str,
) -> list:
    """
    Create a list of Station objects from parsed data and metadata.
    Args:
        data_parser (DataParser): An object that provides access to the observational data.
        metadata_parser (MetadataParser): An object that provides access to the metadata for the stations.
        known_obstypes (dict): A dictionary mapping observation type names to their corresponding ObsType objects.
        timezone (str): a pytz equivalent string indicating the timezone of the timestamps.
    Returns:
        list: A list of Station objects, each representing a station with its associated sensor data and metadata.
    """

    datadf = data_parser.get_df()

    # Station creator
    not_an_obstype = ["name", "datetime"]
    stations = []
    for stationname, stationdata in datadf.groupby("name"):
        # initialize a set of sensordata
        all_station_sensor_data = []
        for obstypename in stationdata.columns:
            if obstypename in not_an_obstype:
                continue

            # Get the corresponding obstype
            obstype = known_obstypes[obstypename]

            # Formatting on raw data
            # 1. Skip stations with nan as name (caused by string casting errors)
            if (stationname == str(np.nan)) | (pd.isnull(stationname)):
                logger.warning(
                    "Skipping the records beloging to station with Nan as name. This could be the result from stringcasting the stationnames."
                )
                continue
            # 2. Drop NAT datetimes if present
            stationdata = stationdata.loc[pd.notnull(stationdata["datetime"])]

            # 3. Skip stations if there are less than 2 records (no freq can be estimated)
            if stationdata.shape[0] < 2:
                logger.warning(
                    f"Station {stationname} is skipped because of has only one record."
                )
                continue

            # Get dataseries:
            sensordata = SensorData(
                stationname=stationname,
                datarecords=stationdata[obstype.name].to_numpy(),
                timestamps=stationdata["datetime"].to_numpy(),
                obstype=obstype,
                # timestamps details:
                timezone=timezone,
                freq_estimation_method=freq_estimation_method,
                freq_estimation_simplify_tolerance=freq_estimation_simplify_tolerance,
                origin_simplify_tolerance=origin_simplify_tolerance,
                timestamp_tolerance=timestamp_tolerance,
            )

            # Add Sensordata:
            all_station_sensor_data.append(sensordata)

        # Create a Site object for the station
        if use_metadata:
            stationsite = Site(
                stationname=stationname,
                latitude=metadata_parser.get_station_lat(stationname),
                longitude=metadata_parser.get_station_lon(stationname),
                extradata=metadata_parser.get_station_extra_metadata(stationname),
            )
        else:
            # no metafile is provided
            stationsite = Site(
                stationname=stationname, latitude=np.nan, longitude=np.nan, extradata={}
            )

        # Combine into a Station
        station = Station(
            stationname=stationname,
            site=stationsite,
            all_sensor_data=all_station_sensor_data,
        )

        stations.append(station)

    if use_metadata:
        # Check the stations present in the metadata but not in the data
        missing_in_data = list(
            set(metadata_parser.get_df().index) - set(datadf["name"].unique())
        )
        if bool(missing_in_data):
            logger.warning(
                "The following stations are defined in the metadatafile but no records are found in the data:\n {missing_in_data}"
            )

    return stations


def import_dataset_from_pkl(target_path):
    """Import a Dataset instance from a (pickle) file.

    Parameters
    ----------
    folder_path : str
        The path to the directory where the pickle file is stored.
    filename : str, optional
        The name of the output file. The default is 'saved_dataset.pkl'.

    Returns
    -------
    metobs_toolkit.Dataset
        The Dataset instance.

    See Also
    --------
    Dataset.save_dataset: Save a Dataset as a pickle file.

    Examples
    --------

    As an example, we will import a Dataset that is stored as a pkl file in the
    current working directory (`os.getcwd()`).

    >>> import os
    >>> import metobs_toolkit
    >>>
    >>> dataset=metobs_toolkit.import_dataset(folder_path=os.getcwd(),
    ...                                      filename='your_saved_dataset.pkl')  # doctest: +SKIP
    >>> dataset # doctest: +SKIP
    Instance of Dataset at ...
    """
    # # check if folder_path is known and exists
    # assert os.path.isdir(folder_path), f"{folder_path} is not a directory!"

    # full_path = os.path.join(folder_path, filename)

    # check if file exists

    picklereader = PickleFileReader(file_path=target_path)
    return picklereader.read_as_local_file()
