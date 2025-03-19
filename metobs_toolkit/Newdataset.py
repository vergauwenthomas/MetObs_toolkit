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

from metobs_toolkit.backend_collection.df_helpers import save_concat
from metobs_toolkit.template import Template, update_known_obstype_with_original_data
from metobs_toolkit.station import Station
from metobs_toolkit.io_collection.metadataparser import MetaDataParser
from metobs_toolkit.io_collection.dataparser import DataParser
from metobs_toolkit.io_collection.filereaders import CsvFileReader, PickleFileReader
from metobs_toolkit.site import Site
from metobs_toolkit.sensordata import SensorData
from metobs_toolkit.backend_collection.printing import dataset_string_repr
from metobs_toolkit.backend_collection.argumentcheckers import (
    fmt_timedelta_arg,
    fmt_datetime_arg,
)
from metobs_toolkit.timestampmatcher import simplify_time
from metobs_toolkit.obstypes import tlk_obstypes, ModelObstype
from metobs_toolkit.obstypes import Obstype

import metobs_toolkit.plot_collection as plotting

from metobs_toolkit.qc_collection import toolkit_buddy_check
from metobs_toolkit.backend_collection.docstring_wrapper import copy_doc
from metobs_toolkit.backend_collection.errorclasses import *
from metobs_toolkit.modeltimeseries import ModelTimeSeries
from metobs_toolkit.settings_files.default_formats_settings import label_def

from metobs_toolkit.geedatasetmanagers import (
    GEEStaticDatasetManager,
    GEEDynamicDatasetManager,
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

    def __eq__(self, other):
        if not isinstance(other, Dataset):
            return False
        return (
            self.stations
            == other.stations
            # and self.obstypes == other.obstypes #tested on sensor level
            # self.template == other.template #template is only for creation needed
        )

    def __str__(self):
        """Represent as text."""
        return "Dataset instance"
        # return dataset_string_repr(self)

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
    def template(self):
        return self._template

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
            if stadf.empty:
                continue
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        return combdf

    @property
    def outliersdf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.outliersdf.reset_index()
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @property
    def gapsdf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.gapsdf.reset_index()
            if stadf.empty:
                continue
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "label", "details"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @property
    def modeldatadf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            stadf = sta.modeldatadf.reset_index()
            if stadf.empty:
                continue
            stadf["name"] = sta.name
            concatlist.append(stadf.set_index(["datetime", "obstype", "name"]))

        combdf = save_concat((concatlist))
        combdf.sort_index(inplace=True)
        if combdf.empty:
            combdf = pd.DataFrame(
                columns=["value", "details"],
                index=pd.MultiIndex(
                    levels=[[], [], []],
                    codes=[[], [], []],
                    names=["datetime", "obstype", "name"],
                ),
            )
        return combdf

    @property
    def metadf(self) -> pd.DataFrame:
        concatlist = []
        for sta in self.stations:
            concatlist.append(sta.metadf)
        return save_concat((concatlist)).sort_index()

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
        """Get a Station by the station name.

        Args:
            stationname (str): The stationname

        Raises:
            MetObsStationNotFound: If the station is not found by name.

        Returns:
            Station: metobs_toolkit.Station
        """
        # create lookup by name
        stationlookup = {sta.name: sta for sta in self.stations}
        try:
            station = stationlookup[stationname]
        except KeyError:
            raise MetObsStationNotFound(f"{stationname} is not found in {self}.")

        return station

    def get_info(self, printout=True):
        infostr = ""
        df = self.df
        if df.empty:
            infostr += "Dataset instance without observation records."
        else:
            present_obstypes = list(df.index.get_level_values("obstype").unique())

            infostr += " --- General Info --- \n\n"
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

    def sync_records(
        self,
        target_obstype: str = "temp",
        timestamp_shift_tolerance="2min",
        freq_shift_tolerance="1min",
        fixed_origin=None,
    ):

        # NOTE: cumulative tolerance errors with the settings used to import the dataset !!!

        # format arguments
        fixed_origin = fmt_timedelta_arg(fixed_origin)
        freq_shift_tolerance = fmt_timedelta_arg(freq_shift_tolerance)
        timestamp_shift_tolerance = fmt_timedelta_arg(timestamp_shift_tolerance)

        for sta in self.stations:
            # check if has sensordata
            if target_obstype in sta.obsdata.keys():
                sensor = sta.obsdata[target_obstype]

                freq_target = simplify_time(
                    time=sensor.freq,
                    max_simplify_error=freq_shift_tolerance,
                    zero_protection=True,
                )

                # Note that simplify_time (tries to) simplifies to a target,
                # all targets are natural multipicates of each other --> ensuring syncing
                # over different sensors.

                sensor.resample(
                    target_freq=freq_target,
                    shift_tolerance=timestamp_shift_tolerance,
                    origin=fixed_origin,
                    origin_simplify_tolerance=timestamp_shift_tolerance,
                )
            else:
                logger.warning(
                    f"{sta} does not have {target_obstype} sensordata and is skipped in the synchronisation."
                )

    def resample(
        self,
        target_freq,
        target_obstype: str | None = None,
        shift_tolerance=pd.Timedelta("4min"),
        origin=None,
        origin_simplify_tolerance=pd.Timedelta("4min"),
    ):
        # NOTE: cumulative tolerance errors with the settings used to import the dataset !!!
        # NOTE: if target_obstype is none, all sensors are resampled to the same freq!

        # format arguments
        target_freq = fmt_timedelta_arg(target_freq)
        shift_tolerance = fmt_timedelta_arg(shift_tolerance)

        # apply over all station
        for sta in self.stations:
            sta.resample(
                target_freq=target_freq,
                target_obstype=target_obstype,
                shift_tolerance=shift_tolerance,
                origin=origin,
                origin_simplify_tolerance=origin_simplify_tolerance,
            )

    # ------------------------------------------
    #   Reading/writing data
    # ------------------------------------------

    def import_gee_data_from_file(
        self,
        filepath: str,
        geedynamicdatasetmanager: GEEDynamicDatasetManager,
        force_update=True,
        _force_from_dataframe=None,
    ) -> pd.DataFrame:
        if _force_from_dataframe is None:
            # Reading the data
            reader = CsvFileReader(file_path=filepath)
            data = reader.read_as_local_file()
            force_update = True

            # 1. format the data
            totaldf = geedynamicdatasetmanager._format_gee_df_structure(data)

            # 2. convert units
            totaldf = geedynamicdatasetmanager._convert_units(totaldf)
            # this totaldf will be returned
        else:
            totaldf = _force_from_dataframe

        # 3. Subset to known obstypes
        known_obstypes = list(geedynamicdatasetmanager.modelobstypes.keys())
        cols_to_skip = list(set(totaldf.columns) - set(known_obstypes))
        if bool(cols_to_skip):
            logger.warning(
                f"The folowing columns in the GEE datafile are not present in the known modelobstypes of {geedynamicdatasetmanager}: {cols_to_skip}"
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
                    obstype=geedynamicdatasetmanager.modelobstypes[col],
                    timezone="UTC",
                    modelname=geedynamicdatasetmanager.name,
                    modelvariable=geedynamicdatasetmanager.modelobstypes[
                        col
                    ].model_band,
                )
                # add it to the station
                sta.add_to_modeldata(modeltimeseries, force_update=force_update)

        return totaldf

    def add_new_observationtype(self, obstype: Obstype):
        # Check type
        if not isinstance(obstype, Obstype):
            raise MetObsWrongType(
                f"{obstype} is not an instance of metobs_toolkit.Obstype."
            )
        # Check if the name is unique
        if obstype.name in self.obstypes.keys():
            raise MetObsDataAlreadyPresent(
                f"An Obstype with {obstype.name} as name is already present in the obstypes: {self.obstypes}"
            )

        # add it to the knonw obstypes
        self.obstypes.update({obstype.name: obstype})

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
        template_file: str | Path,
        input_data_file: str | Path = None,
        input_metadata_file: str | Path = None,
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
        if (input_data_file is None) & (template_file is None):
            raise MetObsMissingFile(
                f"No input_data_file or input_metadata_file is provided"
            )

        assert template_file is not None, "No templatefile is specified."

        # Read template
        logger.info(f"Reading the templatefile")
        self.template.read_template_from_file(
            jsonpath=template_file, templatefile_is_url=templatefile_is_url
        )

        # Read Datafile
        if input_data_file is not None:
            use_data = True
            dataparser = DataParser(
                datafilereader=CsvFileReader(file_path=input_data_file),
                template=self.template,
            )
            dataparser.parse(**kwargs_data_read)  # read and parse to a dataframe
        else:
            logger.info("No datafile is provided --> metadata-only mode")
            dataparser = None  # will not be used
            use_data = False

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

        if use_data:
            # Add original columnname and units to the known obstypes
            self.obstypes = update_known_obstype_with_original_data(
                known_obstypes=self.obstypes, template=self.template
            )

        # Special case: in single-station case the name can be define in the template,
        # but a (different) name can be present in the metadata file for the same station.
        # This results in incompatible data-metadata.
        if (use_metadata) & (use_data) & (self.template.data_is_single_station):
            templatename = self.template.single_station_name
            metadataparser._overwrite_name(target_single_name=templatename)

        # Construct Stations
        if use_data:
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
        else:
            # metadata-only mode
            stations = create_metadata_only_stations(metadata_parser=metadataparser)

        # Set stations attribute
        self.stations = stations

    # ------------------------------------------
    #    Plotting
    # ------------------------------------------

    def make_plot_of_modeldata(
        self,
        obstype: str = "temp",
        colormap: dict | None = None,
        ax=None,
        figkwargs: dict = {},
        title: str | None = None,
        linestyle="--",
    ) -> Axes:

        modeldatadf = self.modeldatadf
        if obstype not in modeldatadf.index.get_level_values("obstype"):
            raise MetObsObstypeNotFound(f"There is no modeldata present of {obstype}")

        # Get the modelobstype --> find a Station that holds it
        for sta in self.stations:
            if obstype in sta.modeldata.keys():
                modelobstype = sta.modeldata[obstype].obstype
                modelname = sta.modeldata[obstype].modelname
                modelvar = sta.modeldata[obstype].modelname
                break

        # Create new axes if needed
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        plotdf = (
            modeldatadf.xs(obstype, level="obstype", drop_level=False)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        plotdf = plotdf[["value"]]
        plotdf["label"] = label_def["goodrecord"][
            "label"
        ]  # Just so that they are plotted as lines

        # Define linecolor
        if colormap is None:
            # create color defenitions
            colormap = plotting.create_station_color_map(
                catlist=plotdf.index.get_level_values("name").unique()
            )

        ax = plotting.plot_timeseries_color_by_station(
            plotdf=plotdf,
            colormap=colormap,
            show_outliers=False,  # will not be used,
            show_gaps=False,  # will not be used
            ax=ax,
            linestyle=linestyle,
            legend_prefix=f"{modelname}:{modelvar}@",
        )
        # Styling
        # Set title:
        if title is None:
            plotting.set_title(
                ax, f"{modelobstype.name} data of {modelname} at stations locations."
            )
        else:
            plotting.set_title(ax, title)

        # Set ylabel
        plotting.set_ylabel(ax, modelobstype._get_plot_y_label())

        # Set xlabel
        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax)

        # Add legend
        plotting.set_legend(ax)

        return ax

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

        # Create an axis
        if ax is None:
            ax = plotting.create_axes(**figkwargs)

        # construct plotdf
        plotdf = (
            self.df.xs(obstype, level="obstype", drop_level=False)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        if plotdf.empty:
            raise MetObsObstypeNotFound(
                f"There are no records of {obstype} for plotting."
            )

        colormap = None
        if show_modeldata:
            # create color defenitions
            colormap = plotting.create_station_color_map(
                catlist=plotdf.index.get_level_values("name").unique()
            )
            ax = self.make_plot_of_modeldata(
                obstype=obstype,
                colormap=colormap,
                ax=ax,
                figkwargs=figkwargs,
                title=title,  # will always be overwritten
            )
        # Set colors scheme
        if colorby == "station":
            if colormap is None:  # when no modeldata is added
                # create color defenitions
                colormap = plotting.create_station_color_map(
                    catlist=plotdf.index.get_level_values("name").unique()
                )
            ax = plotting.plot_timeseries_color_by_station(
                plotdf=plotdf,
                colormap=colormap,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
                linestyle="-",
            )

        elif colorby == "label":
            ax = plotting.plot_timeseries_color_by_label(
                plotdf=plotdf,
                # linecolor=linecolor,
                show_outliers=show_outliers,
                show_gaps=show_gaps,
                ax=ax,
            )
        else:
            raise ValueError(
                f'colorby is either "station" or "label" but not {colorby}'
            )

        # Styling
        obstypeinstance = self.obstypes[obstype]

        # Set title:
        if title is None:
            plotting.set_title(ax, f"{obstypeinstance.name} data.")
        else:
            plotting.set_title(ax, title)

        # Set ylabel
        plotting.set_ylabel(ax, obstypeinstance._get_plot_y_label())

        # Set xlabel
        cur_tz = plotdf.index.get_level_values("datetime").tz
        plotting.set_xlabel(ax, f"Timestamps (in {cur_tz})")

        # Format timestamp ticks
        plotting.format_datetime_axes(ax)

        # Add legend
        plotting.set_legend(ax)

        return ax

    def make_gee_plot(
        self,
        geedatasetmanager: (
            GEEStaticDatasetManager | GEEDynamicDatasetManager
        ) = default_datasets["lcz"],
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
        if isinstance(geedatasetmanager, GEEStaticDatasetManager):
            kwargs = dict(
                outputfolder=outputfolder,
                filename=filename,
                save=save,
                vmin=vmin,
                vmax=vmax,
                overwrite=overwrite,
            )

        elif isinstance(geedatasetmanager, GEEDynamicDatasetManager):
            if timeinstance is None:
                raise ValueError(
                    f"Timeinstance is None, but is required for a dynamic dataset like {geedatasetmanager}"
                )
            timeinstance = fmt_datetime_arg(timeinstance, input_tz="UTC")
            if modelobstype not in geedatasetmanager.modelobstypes:
                raise MetObsObstypeNotFound(
                    f"{modelobstype} is not a known modelobstype of {geedatasetmanager}. These are knonw: {geedatasetmanager.modelobstypes} "
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

        return geedatasetmanager.make_gee_plot(metadf=self.metadf, **kwargs)

    # ------------------------------------------
    #    Gee extracting
    # ------------------------------------------

    def get_static_gee_point_data(
        self,
        geestaticdatasetmanager: GEEStaticDatasetManager,
        update_stations: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """simple but slow option is to loop over all stations and get the lcz,
        but this requires N-API calls (N number of stations).

        Faster: construct the metadf with all stations, and get the lcs from one api call
        """
        if not isinstance(geestaticdatasetmanager, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdatasetmanager)}"
            )

        # initialize gee api
        if initialize_gee:
            connect_to_gee()

        geedf = geestaticdatasetmanager.extract_static_point_data(self.metadf)

        # update the station attributes
        varname = geestaticdatasetmanager.name
        if update_stations:
            for staname, geedict in geedf.to_dict(orient="index").items():
                self.get_station(staname).site.set_geedata(varname, geedict[varname])
        return geedf

    def get_static_gee_buffer_fraction_data(
        self,
        geestaticdatasetmanager: GEEStaticDatasetManager,
        buffers=[100],
        aggregate=False,
        update_stations: bool = True,
        initialize_gee: bool = True,
    ) -> pd.DataFrame:
        """simple but slow option is to loop over all stations and get the lcz,
        but this requires N-API calls (N number of stations).

        Faster: construct the metadf with all stations, and get the lcs from one api call
        """
        if not isinstance(geestaticdatasetmanager, GEEStaticDatasetManager):
            raise ValueError(
                f"geestaticdataset should be an isntance of GeeStaticDataset, not {type(geestaticdatasetmanager)}"
            )

        # initialize gee api
        if initialize_gee:
            connect_to_gee()

        dflist = []
        for bufferradius in buffers:
            geedf = geestaticdatasetmanager.extract_static_buffer_frac_data(
                metadf=self.metadf, bufferradius=bufferradius, agg_bool=aggregate
            )
            dflist.append(geedf)

        geedf = save_concat((dflist))

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
            geestaticdatasetmanager=default_datasets["worldcover"],
            buffers=buffers,
            aggregate=aggregate,
            update_stations=update_stations,
            initialize_gee=initialize_gee,
        )

    def get_gee_timeseries_data(
        self,
        geedynamicdatasetmanager: GEEDynamicDatasetManager,
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
        if not isinstance(geedynamicdatasetmanager, GEEDynamicDatasetManager):
            raise ValueError(
                f"geedynamicdataset should be an isntance of GeeDynamicDataset, not {type(geedynamicdatasetmanager)}"
            )

        # Format datetime arguments
        if startdt_utc is None:
            if self.df.empty:
                raise MetObsMissingArgument(
                    "No data is present in the dataset, thus a startdt_utc is required."
                )
            startdt_utc = self.start_datetime.tz_convert("UTC")
        else:
            startdt_utc = fmt_datetime_arg(startdt_utc, input_tz="UTC")

        if enddt_utc is None:
            if self.df.empty:
                raise MetObsMissingArgument(
                    "No data is present in the dataset, thus a enddt_utc is required."
                )
            enddt_utc = self.end_datetime.tz_convert("UTC")
        else:
            enddt_utc = fmt_datetime_arg(enddt_utc, input_tz="UTC")

        # check if target_obstypes are mapped to bands
        for obst in target_obstypes:
            if obst not in geedynamicdatasetmanager.modelobstypes.keys():
                raise MetObsMetadataNotFound(
                    f"{obst} is not a known modelobstype of {geedynamicdatasetmanager}."
                )

        # create specific name for the file that might be written to Drive
        if drive_filename is None:
            drive_filename = f"{geedynamicdatasetmanager.name}_timeseries_data_of_full_dataset_{len(self.stations)}_stations"  # do not include csv extension here

        df = geedynamicdatasetmanager.extract_timeseries_data(
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
            geedynamicdatasetmanager=geedynamicdatasetmanager,
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
        if (lapserate is not None) | (max_alt_diff is not None):
            # test if altitude data is available
            if not all([sta.site.flag_has_altitude() for sta in self.stations]):
                raise MetObsMetadataNotFound(
                    "Not all stations have altitude data, lapserate correction and max_alt_diff filetering could not be applied."
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
        # outlierslist is a list of tuples (stationname, datetime, msg) that are outliers
        # timestamp_map is a dict with keys the stationname and values a series to map the syncronized
        # timestamps to the original timestamps

        # convert to a dataframe
        df = pd.DataFrame(data=outlierslist, columns=["name", "datetime", "detail_msg"])

        # update all the sensordata
        for station in self.stations:
            if target_obstype in station.obsdata.keys():
                # Get the sensordata object
                sensorddata = station.obsdata[target_obstype]

                # get outlier datetimeindex
                outldt = pd.DatetimeIndex(df[df["name"] == station.name]["datetime"])

                if not outldt.empty:
                    # convert to original timestamps
                    dtmap = timestamp_map[station.name]
                    outldt = outldt.map(dtmap)

                # update the sensordata
                sensorddata._update_outliers(
                    qccheckname="buddy_check",
                    outliertimestamps=outldt,
                    check_kwargs=qc_kwargs,
                    extra_columns={
                        "detail_msg": df[df["name"] == station.name][
                            "detail_msg"
                        ].to_numpy()
                    },
                )

    def get_qc_stats(self, target_obstype="temp", make_plot=True):

        freqdf_list = [
            sta.get_qc_stats(target_obstype=target_obstype, make_plot=False)
            for sta in self.stations
        ]

        dfagg = (
            pd.concat(freqdf_list)
            .reset_index()
            .groupby(["qc_check"])
            .sum()
            .drop(columns=["name"])
        )

        if make_plot:
            fig = plotting.qc_overview_pies(df=dfagg)
            fig.suptitle(
                f"QC frequency statistics of {target_obstype} on Dataset level."
            )
            return fig
        else:
            return dfagg

    # ------------------------------------------
    #    Other methods
    # ------------------------------------------

    def rename_stations(self, renamedict: dict):

        if not isinstance(renamedict, dict):
            raise TypeError(f"{renamedict} is not a Dict")

        for origname, trgname in renamedict.items():
            # test if the station exist
            if origname not in [sta.name for sta in self.stations]:
                logger.warning(f"{origname} is not present, skipped.")
                continue

            # test if targetname is unknown
            if trgname in [sta.name for sta in self.stations]:
                logger.warning(
                    f"{trgname} is already present, renaming {origname} --> {trgname} is skipped."
                )
                continue
            # rename
            self.get_station(origname)._rename(targetname=trgname)

    def convert_outliers_to_gaps(self, all_observations=True, obstype="temp"):
        for sta in self.stations:
            sta.convert_outliers_to_gaps(
                all_observations=all_observations, obstype=obstype
            )

    # ------------------------------------------
    #    Gapfilling
    # ------------------------------------------
    @copy_doc(Station.interpolate_gaps)
    def interpolate_gaps(
        self,
        target_obstype: str,
        method: str = "time",
        max_consec_fill: int = 10,
        n_leading_anchors: int = 1,
        n_trailing_anchors: int = 1,
        max_lead_to_gap_distance: pd.Timedelta | None = None,
        max_trail_to_gap_distance: pd.Timedelta | None = None,
        overwrite_fill=False,
        method_kwargs={},
    ):
        # special formatters
        max_lead_to_gap_distance = fmt_timedelta_arg(max_lead_to_gap_distance)
        max_trail_to_gap_distance = fmt_timedelta_arg(max_trail_to_gap_distance)

        for sta in self.stations:
            sta.interpolate_gaps(
                target_obstype=target_obstype,
                method=method,
                max_consec_fill=max_consec_fill,
                n_leading_anchors=n_leading_anchors,
                n_trailing_anchors=n_trailing_anchors,
                max_lead_to_gap_distance=max_lead_to_gap_distance,
                max_trail_to_gap_distance=max_trail_to_gap_distance,
                overwrite_fill=False,
                method_kwargs=method_kwargs,
            )

    @copy_doc(Station.fill_gaps_with_raw_modeldata)
    def fill_gaps_with_raw_modeldata(self, target_obstype: str, overwrite_fill=False):
        for sta in self.stations:
            sta.fill_gaps_with_raw_modeldata(
                target_obstype=target_obstype, overwrite_fill=overwrite_fill
            )

    @copy_doc(Station.fill_gaps_with_debiased_modeldata)
    def fill_gaps_with_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration=pd.Timedelta("24h"),
        min_leading_records_total: int = 60,
        trailing_period_duration=pd.Timedelta("24h"),
        min_trailing_records_total: int = 60,
        overwrite_fill=False,
    ):
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        for sta in self.stations:
            sta.fill_gaps_with_debiased_modeldata(
                target_obstype=target_obstype,
                leading_period_duration=leading_period_duration,
                min_leading_records_total=min_leading_records_total,
                trailing_period_duration=trailing_period_duration,
                min_trailing_records_total=min_trailing_records_total,
                overwrite_fill=overwrite_fill,
            )

    @copy_doc(Station.fill_gaps_with_diurnal_debiased_modeldata)
    def fill_gaps_with_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration=pd.Timedelta("24h"),
        trailing_period_duration=pd.Timedelta("24h"),
        min_debias_sample_size: int = 6,
        overwrite_fill=False,
    ):
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        for sta in self.stations:
            sta.fill_gaps_with_diurnal_debiased_modeldata(
                target_obstype=target_obstype,
                leading_period_duration=leading_period_duration,
                trailing_period_duration=trailing_period_duration,
                min_debias_sample_size=min_debias_sample_size,
                overwrite_fill=overwrite_fill,
            )

    @copy_doc(Station.fill_gaps_with_weighted_diurnal_debiased_modeldata)
    def fill_gaps_with_weighted_diurnal_debiased_modeldata(
        self,
        target_obstype: str,
        leading_period_duration=pd.Timedelta("24h"),
        trailing_period_duration=pd.Timedelta("24h"),
        min_lead_debias_sample_size: int = 2,
        min_trail_debias_sample_size: int = 2,
        overwrite_fill=False,
    ):
        # special formatters
        leading_period_duration = fmt_timedelta_arg(leading_period_duration)
        trailing_period_duration = fmt_timedelta_arg(trailing_period_duration)

        for sta in self.stations:
            sta.fill_gaps_with_weighted_diurnal_debiased_modeldata(
                target_obstype=target_obstype,
                leading_period_duration=leading_period_duration,
                trailing_period_duration=trailing_period_duration,
                min_lead_debias_sample_size=min_lead_debias_sample_size,
                min_trail_debias_sample_size=min_trail_debias_sample_size,
                overwrite_fill=overwrite_fill,
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


def create_metadata_only_stations(metadata_parser: MetaDataParser):

    stations = []

    for stationname in metadata_parser.get_df().index:
        stationsite = Site(
            stationname=stationname,
            latitude=metadata_parser.get_station_lat(stationname),
            longitude=metadata_parser.get_station_lon(stationname),
            extradata=metadata_parser.get_station_extra_metadata(stationname),
        )

        # Combine into a Station
        station = Station(
            stationname=stationname,
            site=stationsite,
            all_sensor_data=[],  # No sensordata!
        )

        stations.append(station)
    return stations


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
                f"The following stations are defined in the metadatafile but no records are found in the data:\n {missing_in_data}"
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
