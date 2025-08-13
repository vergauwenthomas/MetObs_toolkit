"""This module is a collection of get_info functions."""

from typing import Union
from metobs_toolkit.backend_collection.dev_collection import copy_doc
from metobs_toolkit.backend_collection.loggingmodule import log_entry
import metobs_toolkit.backend_collection.printing_collection as printing


@log_entry
def dataset_get_info(dataset, printout: bool = True) -> Union[str, None]:
    """
    Retrieve and optionally print detailed information about the Dataset.

    Parameters
    ----------
    dataset : Dataset
        The dataset instance to retrieve information from.
    printout : bool, optional
        If True, prints the information to the console. If False, returns
        the information as a string. Default is True.

    Returns
    -------
    str or None
        A string containing the dataset information if `printout` is False.
        Otherwise, returns None.
    """
    
    infostr = ""
    infostr += printing.print_fmt_title("General info of Dataset")

    # --- Observational info ---
    infostr += printing.print_fmt_section("Observational info")
    df = dataset.df
    if df.empty:
        infostr += printing.print_fmt_line(
            "Dataset instance without observation records."
        )
    else:
        present_obstypes = list(df.index.get_level_values("obstype").unique())

        infostr += printing.print_fmt_line("Dataset instance with:", 0)
        infostr += printing.print_fmt_line(
            f"{len(dataset.stations)} number of stations"
        )
        infostr += printing.print_fmt_line(f"{len(df.index)} number of records")
        infostr += printing.print_fmt_line(
            f"{len(present_obstypes)} types of sensor data are present."
        )
        infostr += printing.print_fmt_line(
            f"Observations from {dataset.start_datetime} -> {dataset.end_datetime}"
        )

        # -- outlier info --
        outldf = dataset.outliersdf
        infostr += printing.print_fmt_line("Outlier info:")
        if outldf.empty:
            infostr += printing.print_fmt_line("No QC outliers present.", 2)
        else:
            infostr += printing.print_fmt_line(
                f"A total of {outldf.shape[0]} outliers are present.", 2
            )
            infostr += printing.print_fmt_line("label counts:", 3)
            infostr += printing.print_fmt_dict(
                outldf["label"].value_counts().to_dict(), identlvl=4
            )
            infostr += printing.print_fmt_line(
                f"For these obstypes: {list(outldf.index.get_level_values('obstype').unique())}",
                2,
            )
            unique_stations = list(outldf.index.get_level_values("name").unique())
            infostr += printing.print_fmt_line(
                f"For {len(unique_stations)} stations: {unique_stations}", 2
            )

        # -- gap info --
        gapsdf = dataset.gapsdf
        infostr += printing.print_fmt_line("Gaps info:")
        if gapsdf.empty:
            infostr += printing.print_fmt_line("No gaps present.", 2)
        else:
            infostr += printing.print_fmt_line(
                f"A total of {gapsdf.shape[0]} gaps are present.", 2
            )
            infostr += printing.print_fmt_line("label counts: ", 3)
            infostr += printing.print_fmt_dict(
                gapsdf["label"].value_counts().to_dict(), identlvl=4
            )
            infostr += printing.print_fmt_line(
                f"For these obstypes: {list(gapsdf.index.get_level_values('obstype').unique())}",
                2,
            )
            unique_stations = list(gapsdf.index.get_level_values("name").unique())
            infostr += printing.print_fmt_line(
                f"For {len(unique_stations)} stations: {unique_stations}", 2
            )

    # Meta data info
    infostr += printing.print_fmt_section("Metadata info")

    metadf = dataset.metadf
    if metadf.empty:
        infostr += printing.print_fmt_line("Dataset instance without metadata.")
    else:
        infostr += printing.print_fmt_line(
            f"{len(metadf.index)} number of stations"
        )
        infostr += printing.print_fmt_line(
            f"The following metadata is present: {list(metadf.columns)}"
        )

    # Modeldata info
    modeldf = dataset.modeldatadf
    infostr += printing.print_fmt_section("Modeldata info")
    if modeldf.empty:
        infostr += printing.print_fmt_line("Dataset instance without modeldata.")
    else:
        infostr += printing.print_fmt_line(
            f"Modeldata is present for {list(modeldf.index.get_level_values('obstype').unique())}"
        )
        infostr += printing.print_fmt_line(
            f"For period {modeldf.index.get_level_values('datetime').min()} -> {modeldf.index.get_level_values('datetime').max()}"
        )

    if printout:
        print(infostr)
    else:
        return infostr

@log_entry
def station_get_info(station, printout: bool = True) -> Union[str, None]:
    """
    Retrieve and optionally print detailed information about the station.

    Parameters
    ----------
    station : Station
        The station instance to retrieve information from.
    printout : bool, optional
        If True, prints the information to the console. If False, returns
        the information as a string. Default is True.

    Returns
    -------
    str or None
        A string containing the station information if `printout` is False.
        Otherwise, returns None.
    """
   

    infostr = ""
    infostr += printing.print_fmt_title("General info of Station")

    # --- Observational info ---
    infostr += printing.print_fmt_section("Observational info")
    df = station.df
    if df.empty:
        infostr += printing.print_fmt_line(
            "Station instance without observation records."
        )
    else:
        present_obstypes = list(df.index.get_level_values("obstype").unique())

        infostr += printing.print_fmt_line("Station instance with:", 0)
        for obstype in present_obstypes:
            infostr += printing.print_fmt_line(f"{obstype}:", 1)
            infostr += sensordata_get_info(sensordata = station.get_sensor(obstype),
                                           printout=False,
                                           _print_as_sensordata=False,
                                           nident_root=2)

    # Meta data info
    infostr += printing.print_fmt_section("Metadata info")
    infostr += site_get_info(site=station.site,
                             printout=False,
                             _print_as_site=False,
                             nident_root=1)
    
    # Model data info
    infostr += printing.print_fmt_section("Modeldata info")
    if not bool(station.modeldata):
        infostr += printing.print_fmt_line("Station instance without model data.")
    else:
        for obstype, modeldata in station.modeldata.items():
            infostr += printing.print_fmt_line(f"{obstype}:", 1)
            infostr += modeltimeseries_get_info(modeldata=modeldata,
                                                printout=False,
                                                _print_as_modeltimeseries=False,
                                                nident_root=2)

    if printout:
        print(infostr)
    else:
        return infostr

@log_entry
def site_get_info(
    site, printout: bool = True, _print_as_site: bool = True, nident_root: int = 1
) -> Union[str, None]:
    """
    Retrieve and optionally print basic information about the site.

    Parameters
    ----------
    site : Site
        The site instance to retrieve information from.
    printout : bool, optional
        If True, print the information. If False, return as string.
    nident_root : int, optional
        The indentation level for formatting. Default is 1.

    Returns
    -------
    str or None
        Information string if printout is False, else None.
    """
    infostr = ""
    if _print_as_site:
        infostr += printing.print_fmt_title("General Info of Site")
        infostr += printing.print_fmt_line(f"Site of {site.stationname}:", 0)
 
    # Coordinates
    if site.flag_has_coordinates():
        infostr += printing.print_fmt_line(
            f"Coordinates ({site.lat}, {site.lon}) (latitude, longitude)",
            nident_root,
        )
    else:
        infostr += printing.print_fmt_line("Coordinates are unknown", nident_root)

    # Altitude
    if site.flag_has_altitude() & (not site.flag_altitude_from_gee()):
        infostr += printing.print_fmt_line(
            f"Altitude: {site.altitude} (m) (from metadata file)", nident_root
        )
    elif site.flag_has_altitude() & (site.flag_altitude_from_gee()):
        infostr += printing.print_fmt_line(
            f"Altitude: {site.altitude} (m) (from GEE extraction)", nident_root
        )
    else:
        infostr += printing.print_fmt_line("Altitude is unknown", nident_root)

    # LCZ
    if site.flag_has_LCZ() & (site.flag_LCZ_from_gee()):
        infostr += printing.print_fmt_line(
            f"LCZ: {site.LCZ} (from GEE extraction)", nident_root
        )
    elif site.flag_has_LCZ():
        infostr += printing.print_fmt_line(
            f"LCZ: {site.LCZ} (from metadata file)", nident_root
        )
    else:
        infostr += printing.print_fmt_line("LCZ is unknown", nident_root)

    # Buffered fractions
    if site.flag_has_landcoverfractions():
        infostr += printing.print_fmt_line(
            "Land cover fractions are available", nident_root
        )
    else:
        infostr += printing.print_fmt_line(
            "Land cover fractions are unknown", nident_root
        )

    # Extra metadata
    if bool(site.extradata):
        infostr += printing.print_fmt_line(
            "Extra metadata from the metadata file:", nident_root
        )
        for key, value in site.extradata.items():
            infostr += printing.print_fmt_line(f"{key}: {value}", nident_root + 1)
    else:
        infostr += printing.print_fmt_line(
            "No extra metadata available", nident_root
        )

    if printout:
        print(infostr)
    else:
        return infostr

@log_entry
def sensordata_get_info(
    sensordata, printout: bool = True, _print_as_sensordata: bool = True, nident_root: int = 1
) -> Union[str, None]:
    """
    Retrieve and optionally print basic information about the sensor data.

    Parameters
    ----------
    sensordata : SensorData
        The sensor data instance to retrieve information from.
    printout : bool, optional
        If True, the information will be printed to the console. If False,
        the information will be returned as a string. Default is True.
    nident_root : int, optional
        The indentation level for formatting. Default is 1.

    Returns
    -------
    str or None
        If `printout` is False, returns a string containing the information
        about the sensor data. If `printout` is True, returns None.
    """
   

    infostr = ""
    if _print_as_sensordata:
        infostr += printing.print_fmt_title("General info of SensorData")
        infostr += printing.print_fmt_line(
            f"{sensordata.obstype.name} records of {sensordata.stationname}:", 0
        )
    
    infostr += printing.print_fmt_line(
        f"{sensordata.obstype.name} observations in {sensordata.obstype.std_unit}", nident_root
    )
    infostr += printing.print_fmt_line(
        f" from {sensordata.start_datetime} -> {sensordata.end_datetime}", nident_root
    )
    infostr += printing.print_fmt_line(
        f" At a resolution of {sensordata.freq}", nident_root
    )

    # outliers info:
    if sensordata.outliersdf.empty:
        infostr += printing.print_fmt_line("No outliers present.", nident_root)
    else:
        infostr += printing.print_fmt_line(
            f"A total of {sensordata.outliersdf.shape[0]} flagged observations (QC outliers).",
            nident_root,
        )
        infostr += printing.print_fmt_line("label counts: ", nident_root + 1)
        infostr += printing.print_fmt_dict(
            sensordata.outliersdf["label"].value_counts().to_dict(), nident_root + 2
        )

    # gaps info:
    if not sensordata.gaps:
        infostr += printing.print_fmt_line("No gaps present.", nident_root)
    else:
        infostr += printing.print_fmt_line(
            f"{len(sensordata.gaps)} gaps present, a total of {sensordata.gapsdf.shape[0]} missing timestamps.",
            nident_root,
        )
        infostr += printing.print_fmt_line("label counts: ", nident_root + 1)
        infostr += printing.print_fmt_dict(
            sensordata.gapsdf["label"].value_counts().to_dict(), nident_root + 2
        )

    if printout:
        print(infostr)
    else:
        return infostr

@log_entry
def gap_get_info(gap, printout: bool = True) -> Union[str, None]:
    """
    Print or return detailed information about the Gap.

    Parameters
    ----------
    gap : Gap
        The gap instance to retrieve information from.
    printout : bool, optional
        If True, prints the information. If False, returns the information as a string. Default is True.

    Returns
    -------
    str or None
        The gap information as a string if printout is False, otherwise None.
    """

    infostr = ""
    infostr += printing.print_fmt_title("General info of Gap")
    infostr += printing.print_fmt_section("Gap details")

    infostr += printing.print_fmt_line(
        f"Gap of {gap.obstype.name} for station: {gap.stationname}", 0
    )
    infostr += printing.print_fmt_line(
        f"From {gap.start_datetime} -> {gap.end_datetime}", 1
    )
    infostr += printing.print_fmt_line(
        f"Duration gap: {gap.end_datetime - gap.start_datetime}", 1
    )

    infostr += printing.print_fmt_section("Gap filling details")
    infostr += printing.print_fmt_line(f"Gap status: {gap.fillstatus}")
    infostr += printing.print_fmt_line("Gapfill settings used:")
    infostr += printing.print_fmt_dict(d=gap.fillsettings, identlvl=2)

    if printout:
        print(infostr)
    else:
        return infostr

@log_entry
def modeltimeseries_get_info(
    modeltimeseries, printout: bool = True, _print_as_modeltimeseries: bool = True, nident_root: int = 1
) -> Union[str, None]:
    """
    Print or return information about the ModelTimeSeries.

    Parameters
    ----------
    modeltimeseries : ModelTimeSeries
        The model time series instance to retrieve information from.
    printout : bool, optional
        If True, print the information. If False, return as string. Default is True.
    nident_root : int, optional
        The indentation level for formatting. Default is 1.

    Returns
    -------
    str or None
        None if printout is True, otherwise the information string.
    """
   
    infostr = ""
    if _print_as_modeltimeseries:
        infostr += printing.print_fmt_title("General info of ModelTimeSeries")
        infostr += printing.print_fmt_line(
            f"{modeltimeseries.obstype.name} model data at location of {modeltimeseries.stationname}"
        )
    infostr += printing.print_fmt_line(
        f"Origin {modeltimeseries.modelname} -> variable/band: {modeltimeseries.modelvariable}",
        nident_root,
    )
    infostr += printing.print_fmt_line(
        f"From {modeltimeseries.start_datetime} --> {modeltimeseries.end_datetime}", nident_root
    )
    infostr += printing.print_fmt_line(
        f"Assumed frequency: {modeltimeseries.freq}", nident_root
    )
    infostr += printing.print_fmt_line(
        f"Number of records: {modeltimeseries.series.shape[0]}", nident_root
    )
    infostr += printing.print_fmt_line(
        f"Units are converted from {modeltimeseries.obstype.model_unit} --> {modeltimeseries.obstype.std_unit}",
        nident_root,
    )

    if printout:
        print(infostr)
    else:
        return infostr



@log_entry
def dynamicdatasetmanager_get_info(dynamicdatasetmanager, printout: bool = True) -> Union[str, None]:
    """Retrieve and optionally print information about the DynamicDatasetManager."""
    pass

@log_entry
def staticdatasetmanager_get_info(staticdatasetmanager, printout:bool = True)-> Union[str, None]:
    pass

@log_entry
def analysis_get_info(
    analysis, printout: bool = True, possible_time_aggregates: list = []
) -> Union[str, None]:
    """
    Provides information about the Analysis instance, including the number
    of records, observation types, metadata columns, station names, and
    known time derivatives.

    Parameters
    ----------
    analysis : Analysis
        The analysis instance to retrieve information from.
    printout : bool, optional
        If True, prints the information to the console. If False, returns
        the information as a string. Default is True.
    possible_time_aggregates : list, optional
        A list of known time derivatives. Default is an empty list.

    Returns
    -------
    str or None
        Returns None if `printout` is True. Returns the information string
        if `printout` is False.
    """

    infostr = ""
    infostr += printing.print_fmt_title("General info of Analysis")
    infostr += printing.print_fmt_line(f"Number of records: {len(analysis.df)}")
    infostr += printing.print_fmt_line(f"Observation types: {list(analysis._df_cols)}")
    infostr += printing.print_fmt_line(
        f"Available metadata columns: {analysis.metadf.columns.tolist()}"
    )
    infostr += printing.print_fmt_line(
        f"Stations: {analysis.fulldf['name'].unique().tolist()}"
    )
    infostr += printing.print_fmt_line(
        f"All known time-derivatives: {possible_time_aggregates}"
    )

    if printout:
        print(infostr)
    else:
        return infostr

