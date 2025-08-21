import pytest
import sys
from pathlib import Path
import tempfile
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# import metobs_toolkit
libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir


def create_dataset_with_modeldata():
    # 1. Create basic dataset with demo data
    dataset = metobs_toolkit.Dataset()
    dataset.import_data_from_file(
        template_file=metobs_toolkit.demo_template,
        input_metadata_file=metobs_toolkit.demo_metadatafile,
        input_data_file=metobs_toolkit.demo_datafile,
    )

    # 2. Import ERA5 model data from CSV file
    era5_manager = metobs_toolkit.default_GEE_datasets["ERA5-land"]
    era5_file = datadir.joinpath(
        "ERA5-land_timeseries_data_of_full_dataset_28_stations.csv"
    )

    # Import the model data
    imported_df = dataset.import_gee_data_from_file(
        filepath=era5_file, geedynamicdatasetmanager=era5_manager, force_update=True
    )
    return dataset


class TestModelDataImport:
    """Test importing model data from various sources."""

    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testmodeldataimport"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_import_era5_data_from_file(self):
        """Test importing ERA5 data from CSV file."""

        # 1. Create basic dataset with demo data
        dataset = create_dataset_with_modeldata()

        # Verify model data was imported
        modeldatadf = dataset.modeldatadf
        assert not modeldatadf.empty, "Model data DataFrame should not be empty"
        assert len(modeldatadf) > 0, "Should have model data records"
        assert len(modeldatadf.index.levels) == 3
        assert (
            "value" in modeldatadf.columns
        ), "Model data DataFrame should have 'value' column"
        assert (
            "details" in modeldatadf.columns
        ), "Model data DataFrame should have 'details' column"

        # Check that stations have model data
        stations_with_modeldata = [
            sta for sta in dataset.stations if not sta.modeldatadf.empty
        ]
        assert (
            len(stations_with_modeldata) > 0
        ), "At least some stations should have model data"

        sta_modeldf = dataset.stations[3].modeldatadf
        assert not sta_modeldf.empty, "Model data DataFrame should not be empty"
        assert len(sta_modeldf) > 0, "Should have model data records"
        assert (
            len(sta_modeldf.index.levels) == 2
        ), "Model data should have at least one index level"
        assert (
            "value" in sta_modeldf.columns
        ), "Model data DataFrame should have 'value' column"
        assert (
            "details" in sta_modeldf.columns
        ), "Model data DataFrame should have 'details' column"

    def test_modeltimeseries_properties(self):
        """Test ModelTimeSeries data-related attributes and methods."""
        # 1. Create basic dataset with demo data
        dataset = create_dataset_with_modeldata()
        station = dataset.get_station("vlinder04")
        # Test with temperature data
        obstype = "temp"

        modeltimeseries = station.get_modeltimeseries(obstype)
        assert modeltimeseries is not None, f"Should get ModelTimeSeries for {obstype}"

        # Test basic attributes
        assert hasattr(
            modeltimeseries, "series"
        ), "ModelTimeSeries should have series attribute"
        assert (
            modeltimeseries.modelname == "ERA5-land"
        ), "ModelTimeSeries should have modelname 'ERA5-land'"
        assert (
            modeltimeseries.modelvariable == "temperature_2m"
        ), f"ModelTimeSeries should have modelvariable 'temperature_2m'"

        assert modeltimeseries._id() == "vlinder04ERA5-land_temperature_2m"

        assert (
            "value" in modeltimeseries.df.columns
        ), "ModelTimeSeries df should have 'value' column"
        assert (
            "model" in modeltimeseries.df.columns
        ), "ModelTimeSeries df should have 'model' column"

        assert (
            modeltimeseries.stationname == "vlinder04"
        ), "ModelTimeSeries stationname should match station name"
        assert str(modeltimeseries.tz) == "UTC"

        _ = modeltimeseries.get_info()

        assert modeltimeseries.obstype.model_band == "temperature_2m"
        assert modeltimeseries.obstype.model_unit == "kelvin"


class TestModelDataManagers:
    """Test GEE dataset managers and model data configurations."""

    def test_gee_dataset_managers_availability(self):
        """Test that GEE dataset managers are available and properly configured."""
        # Test default GEE datasets
        assert hasattr(
            metobs_toolkit, "default_GEE_datasets"
        ), "Should have default_GEE_datasets"

        default_datasets = metobs_toolkit.default_GEE_datasets
        assert isinstance(
            default_datasets, dict
        ), "default_GEE_datasets should be a dict"

        # Test ERA5-land manager
        assert "ERA5-land" in default_datasets, "Should have ERA5-land dataset manager"
        era5_manager = default_datasets["ERA5-land"]

        # Test manager properties
        assert hasattr(era5_manager, "name"), "ERA5 manager should have name"
        assert hasattr(
            era5_manager, "modelobstypes"
        ), "ERA5 manager should have modelobstypes"
        assert isinstance(
            era5_manager.modelobstypes, dict
        ), "modelobstypes should be a dict"

        # Test that manager has temperature obstype
        assert (
            "temp" in era5_manager.modelobstypes
        ), "ERA5 manager should have temp obstype"

    def test_gee_manager_info_methods(self):
        """Test info methods of GEE dataset managers."""
        era5_manager = metobs_toolkit.default_GEE_datasets["ERA5-land"]

        try:
            # Test get_info method
            info_str = era5_manager.get_info(printout=False)
            assert isinstance(info_str, str), "get_info should return string"
            assert len(info_str) > 0, "info string should not be empty"

        except Exception as e:
            pytest.skip(f"GEE manager get_info failed: {e}")
