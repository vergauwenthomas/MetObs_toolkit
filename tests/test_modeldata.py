import pytest
import sys
from pathlib import Path
import tempfile
import shutil
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit
from metobs_toolkit.backend_collection.errorclasses import (
    MetObsModelDataError,
    MetObsDataAlreadyPresent,
)

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

        assert modeltimeseries.modelobstype.model_band == "temperature_2m"
        assert modeltimeseries.modelobstype.model_unit == "kelvin"


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


class TestStationModelDataMethods:
    """Test Station methods for adding and retrieving model data."""

    # to pass to the solutionfixer
    solkwargs = {
        "testfile": Path(__file__).name,
        "classname": "teststationmodeldatamethods",
    }
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def create_test_station_with_multiple_modeldata(self):
        """Create a station with multiple model data sources for testing."""
        # Create basic dataset
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        # Get first station
        station = dataset.stations[0]

        # Create multiple ModelTimeSeries with same obstype but different modelname/modelvariable
        temp_obstype = dataset.obstypes["temp"]
        temp_modelobstype = metobs_toolkit.ModelObstype(obstype=temp_obstype,
                                                        model_unit="kelvin",
                                                        model_band='fake-band')

        # Create test data
        timestamps = pd.date_range("2022-01-01", periods=24, freq="h")
        data1 = np.array([20.0 + i * 0.1 for i in range(24)])  # ERA5 data
        data2 = np.array([20.5 + i * 0.1 for i in range(24)])  # GFS data
        data3 = np.array([19.8 + i * 0.1 for i in range(24)])  # Different ERA5 variable

        # Create ModelTimeSeries with different modelname/modelvariable combinations
        model_ts1 = metobs_toolkit.ModelTimeSeries(
            site=station.site,
            datarecords=data1,
            timestamps=timestamps,
            modelobstype=temp_modelobstype,
            datadtype=np.float32,
            timezone="UTC",
            modelname="ERA5",
            modelvariable="temperature_2m",
        )

        model_ts2 = metobs_toolkit.ModelTimeSeries(
            site=station.site,
            datarecords=data2,
            timestamps=timestamps,
            modelobstype=temp_modelobstype,
            datadtype=np.float32,
            timezone="UTC",
            modelname="GFS",
            modelvariable="temperature_2m",
        )

        model_ts3 = metobs_toolkit.ModelTimeSeries(
            site=station.site,
            datarecords=data3,
            timestamps=timestamps,
            modelobstype=temp_modelobstype,
            datadtype=np.float32,
            timezone="UTC",
            modelname="ERA5",
            modelvariable="skin_temperature",
        )

        return station, model_ts1, model_ts2, model_ts3

    def test_add_to_modeldata_basic(self):
        """Test basic functionality of add_to_modeldata."""
        station, model_ts1, model_ts2, model_ts3 = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Initially no model data
        assert len(station.modeldata) == 0, "Station should start with no model data"

        # Add first model data
        station.add_to_modeldata(model_ts1, force_update=False)
        assert (
            len(station.modeldata) == 1
        ), "Station should have 1 model data after first addition"

        # Add second model data (different modelname)
        station.add_to_modeldata(model_ts2, force_update=False)
        assert (
            len(station.modeldata) == 2
        ), "Station should have 2 model data after second addition"

        # Add third model data (different modelvariable)
        station.add_to_modeldata(model_ts3, force_update=False)
        assert (
            len(station.modeldata) == 3
        ), "Station should have 3 model data after third addition"

    def test_add_to_modeldata_force_update(self):
        """Test add_to_modeldata with force_update functionality."""
        station, model_ts1, model_ts2, model_ts3 = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Add initial model data
        station.add_to_modeldata(model_ts1, force_update=False)
        assert len(station.modeldata) == 1

        # Try to add same model data again without force_update - should raise error
        with pytest.raises(MetObsDataAlreadyPresent):
            station.add_to_modeldata(model_ts1, force_update=False)

        # Add same model data with force_update=True - should succeed
        station.add_to_modeldata(model_ts1, force_update=True)
        assert (
            len(station.modeldata) == 1
        ), "Should still have only 1 model data after force update"

    def test_add_to_modeldata_wrong_type(self):
        """Test add_to_modeldata with wrong input type."""
        station, _, _, _ = self.create_test_station_with_multiple_modeldata()

        # Try to add non-ModelTimeSeries object
        with pytest.raises(
            TypeError,
            match="new_modeltimeseries must be an instance of ModelTimeSeries",
        ):
            station.add_to_modeldata("not_a_modeltimeseries", force_update=False)

    def test_get_modeltimeseries_by_obstype_only(self):
        """Test get_modeltimeseries with obstype only when unique."""
        station, model_ts1, _, _ = self.create_test_station_with_multiple_modeldata()

        # Add only one model data for temp
        station.add_to_modeldata(model_ts1, force_update=False)

        # Should be able to retrieve by obstype only
        retrieved = station.get_modeltimeseries("temp")
        assert (
            retrieved._id() == model_ts1._id()
        ), "Should retrieve the correct model data"
        assert retrieved.modelname == "ERA5", "Should have correct modelname"
        assert (
            retrieved.modelvariable == "temperature_2m"
        ), "Should have correct modelvariable"

    def test_get_modeltimeseries_multiple_same_obstype(self):
        """Test get_modeltimeseries when multiple model data exist for same obstype."""
        station, model_ts1, model_ts2, model_ts3 = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Add multiple model data for same obstype
        station.add_to_modeldata(model_ts1, force_update=False)
        station.add_to_modeldata(model_ts2, force_update=False)
        station.add_to_modeldata(model_ts3, force_update=False)

        # Should raise error when trying to get by obstype only
        with pytest.raises(MetObsModelDataError, match="Multiple model data found"):
            station.get_modeltimeseries("temp")

    def test_get_modeltimeseries_by_modelname(self):
        """Test get_modeltimeseries filtering by modelname."""
        station, model_ts1, model_ts2, model_ts3 = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Add multiple model data
        station.add_to_modeldata(model_ts1, force_update=False)
        station.add_to_modeldata(model_ts2, force_update=False)
        station.add_to_modeldata(model_ts3, force_update=False)

        # Should still raise error when filtering by modelname only (ERA5 has 2 variables)
        with pytest.raises(MetObsModelDataError, match="Multiple model data found"):
            station.get_modeltimeseries("temp", modelname="ERA5")

        # Should work when filtering by unique modelname
        retrieved = station.get_modeltimeseries("temp", modelname="GFS")
        assert retrieved._id() == model_ts2._id(), "Should retrieve GFS model data"
        assert retrieved.modelname == "GFS", "Should have correct modelname"

    def test_get_modeltimeseries_by_modelvariable(self):
        """Test get_modeltimeseries filtering by modelvariable."""
        station, model_ts1, model_ts2, model_ts3 = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Add multiple model data
        station.add_to_modeldata(model_ts1, force_update=False)
        station.add_to_modeldata(model_ts2, force_update=False)
        station.add_to_modeldata(model_ts3, force_update=False)

        # Should work when filtering by unique modelvariable
        retrieved = station.get_modeltimeseries(
            "temp", modelvariable="skin_temperature"
        )
        assert (
            retrieved._id() == model_ts3._id()
        ), "Should retrieve skin_temperature model data"
        assert (
            retrieved.modelvariable == "skin_temperature"
        ), "Should have correct modelvariable"

    def test_get_modeltimeseries_by_both_filters(self):
        """Test get_modeltimeseries filtering by both modelname and modelvariable."""
        station, model_ts1, model_ts2, model_ts3 = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Add multiple model data
        station.add_to_modeldata(model_ts1, force_update=False)
        station.add_to_modeldata(model_ts2, force_update=False)
        station.add_to_modeldata(model_ts3, force_update=False)

        # Should work when filtering by both modelname and modelvariable
        retrieved = station.get_modeltimeseries(
            "temp", modelname="ERA5", modelvariable="temperature_2m"
        )
        assert (
            retrieved._id() == model_ts1._id()
        ), "Should retrieve ERA5 temperature_2m model data"
        assert retrieved.modelname == "ERA5", "Should have correct modelname"
        assert (
            retrieved.modelvariable == "temperature_2m"
        ), "Should have correct modelvariable"

        # Test another combination
        retrieved2 = station.get_modeltimeseries(
            "temp", modelname="ERA5", modelvariable="skin_temperature"
        )
        assert (
            retrieved2._id() == model_ts3._id()
        ), "Should retrieve ERA5 skin_temperature model data"

    def test_get_modeltimeseries_no_match(self):
        """Test get_modeltimeseries when no model data matches the criteria."""
        station, model_ts1, _, _ = self.create_test_station_with_multiple_modeldata()

        # Add only one model data
        station.add_to_modeldata(model_ts1, force_update=False)

        # Test non-existing obstype
        with pytest.raises(MetObsModelDataError, match="No model data found"):
            station.get_modeltimeseries("humidity")

        # Test non-existing modelname
        with pytest.raises(MetObsModelDataError, match="No model data found"):
            station.get_modeltimeseries("temp", modelname="ECMWF")

        # Test non-existing modelvariable
        with pytest.raises(MetObsModelDataError, match="No model data found"):
            station.get_modeltimeseries("temp", modelvariable="nonexistent_var")

        # Test non-existing combination
        with pytest.raises(MetObsModelDataError, match="No model data found"):
            station.get_modeltimeseries(
                "temp", modelname="ERA5", modelvariable="nonexistent_var"
            )

    def test_get_modeltimeseries_no_modeldata(self):
        """Test get_modeltimeseries when station has no model data."""
        station, _, _, _ = self.create_test_station_with_multiple_modeldata()

        # Don't add any model data
        with pytest.raises(MetObsModelDataError, match="No model data found"):
            station.get_modeltimeseries("temp")

    def test_modeldata_property_list_format(self):
        """Test that modeldata property returns a list."""
        station, model_ts1, model_ts2, _ = (
            self.create_test_station_with_multiple_modeldata()
        )

        # Add model data
        station.add_to_modeldata(model_ts1, force_update=False)
        station.add_to_modeldata(model_ts2, force_update=False)

        # Test modeldata property
        modeldata_list = station.modeldata
        assert isinstance(modeldata_list, list), "modeldata should return a list"
        assert len(modeldata_list) == 2, "Should have 2 model data items"

        # Test that items are ModelTimeSeries
        for item in modeldata_list:
            assert isinstance(
                item, metobs_toolkit.ModelTimeSeries
            ), "Items should be ModelTimeSeries"



if __name__ == "__main__":

    tester = TestStationModelDataMethods()
    tester.test_add_to_modeldata_basic()