"""Test the missing station detection and find_nearest_with_data feature."""

import pytest
import sys
import logging
from pathlib import Path
from unittest.mock import patch, MagicMock

import pandas as pd
import numpy as np

# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))

import metobs_toolkit
from metobs_toolkit.geedatasetmanagers import GEEDynamicDatasetManager


class TestMissingStationWarning:
    """Test that warnings are raised for stations without data."""

    def _get_era5_manager(self):
        """Get the ERA5-land dataset manager."""
        return metobs_toolkit.default_GEE_datasets["ERA5-land"]

    def _get_test_metadf(self):
        """Create a test metadata dataframe with 3 stations."""
        return pd.DataFrame(
            {
                "lat": [51.0, 51.1, 51.2],
                "lon": [3.7, 3.8, 3.9],
            },
            index=pd.Index(["station_a", "station_b", "station_c"], name="name"),
        )

    @patch("metobs_toolkit.geedatasetmanagers.gee_api")
    @patch("metobs_toolkit.geedatasetmanagers.ee")
    @patch("metobs_toolkit.geedatasetmanagers.connect_to_gee")
    def test_warning_raised_for_missing_stations(
        self, mock_connect, mock_ee, mock_gee_api, caplog
    ):
        """Test that a warning is raised when some stations have no data returned."""
        era5 = self._get_era5_manager()
        metadf = self._get_test_metadf()

        # Mock the GEE API responses - only return data for station_a and station_b
        mock_gee_api._df_to_features_point_collection.return_value = MagicMock()
        mock_gee_api._estimate_data_size.return_value = 100  # small enough for direct
        mock_gee_api.get_ee_obj.return_value = MagicMock()
        mock_gee_api.datetime_to_gee_datetime.return_value = MagicMock()
        mock_gee_api._addDate = MagicMock()

        # Simulate GEE returning data only for station_a and station_b (not station_c)
        mock_results = {
            "features": [
                {
                    "properties": {
                        "name": "station_a",
                        "datetime": "20220901120000",
                        "temperature_2m": 293.15,
                    }
                },
                {
                    "properties": {
                        "name": "station_b",
                        "datetime": "20220901120000",
                        "temperature_2m": 294.15,
                    }
                },
            ]
        }

        # Set up the mock chain for raster.filter().map().map().flatten().getInfo()
        mock_raster = MagicMock()
        mock_gee_api.get_ee_obj.return_value = mock_raster
        mock_raster.filter.return_value.map.return_value.map.return_value.flatten.return_value.getInfo.return_value = (
            mock_results
        )

        with caplog.at_level(logging.WARNING):
            df = era5.extract_timeseries_data(
                metadf=metadf,
                startdt_utc=pd.Timestamp("2022-09-01", tz="UTC"),
                enddt_utc=pd.Timestamp("2022-09-01 12:00:00", tz="UTC"),
                obstypes=["temp"],
                get_all_bands=False,
                force_direct_transfer=True,
                initialize_gee=True,
                find_nearest_with_data=False,
            )

        # Check that warning was logged about missing station_c
        warning_messages = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
        assert any("station_c" in msg for msg in warning_messages), (
            f"Expected warning about missing station 'station_c'. Got: {warning_messages}"
        )

    @patch("metobs_toolkit.geedatasetmanagers.gee_api")
    @patch("metobs_toolkit.geedatasetmanagers.ee")
    @patch("metobs_toolkit.geedatasetmanagers.connect_to_gee")
    def test_no_warning_when_all_stations_have_data(
        self, mock_connect, mock_ee, mock_gee_api, caplog
    ):
        """Test that no warning is raised when all stations have data."""
        era5 = self._get_era5_manager()
        metadf = self._get_test_metadf()

        mock_gee_api._df_to_features_point_collection.return_value = MagicMock()
        mock_gee_api._estimate_data_size.return_value = 100
        mock_gee_api.get_ee_obj.return_value = MagicMock()
        mock_gee_api.datetime_to_gee_datetime.return_value = MagicMock()
        mock_gee_api._addDate = MagicMock()

        # Simulate GEE returning data for ALL stations
        mock_results = {
            "features": [
                {
                    "properties": {
                        "name": "station_a",
                        "datetime": "20220901120000",
                        "temperature_2m": 293.15,
                    }
                },
                {
                    "properties": {
                        "name": "station_b",
                        "datetime": "20220901120000",
                        "temperature_2m": 294.15,
                    }
                },
                {
                    "properties": {
                        "name": "station_c",
                        "datetime": "20220901120000",
                        "temperature_2m": 295.15,
                    }
                },
            ]
        }

        mock_raster = MagicMock()
        mock_gee_api.get_ee_obj.return_value = mock_raster
        mock_raster.filter.return_value.map.return_value.map.return_value.flatten.return_value.getInfo.return_value = (
            mock_results
        )

        with caplog.at_level(logging.WARNING):
            df = era5.extract_timeseries_data(
                metadf=metadf,
                startdt_utc=pd.Timestamp("2022-09-01", tz="UTC"),
                enddt_utc=pd.Timestamp("2022-09-01 12:00:00", tz="UTC"),
                obstypes=["temp"],
                get_all_bands=False,
                force_direct_transfer=True,
                initialize_gee=True,
                find_nearest_with_data=False,
            )

        # No warning should be raised about missing stations
        warning_messages = [r.message for r in caplog.records if r.levelno >= logging.WARNING]
        assert not any("no nearest gridpoint" in msg for msg in warning_messages), (
            f"Unexpected warning about missing stations. Got: {warning_messages}"
        )

    @patch("metobs_toolkit.geedatasetmanagers.gee_api")
    @patch("metobs_toolkit.geedatasetmanagers.ee")
    @patch("metobs_toolkit.geedatasetmanagers.connect_to_gee")
    def test_find_nearest_with_data_called_for_missing_stations(
        self, mock_connect, mock_ee, mock_gee_api, caplog
    ):
        """Test that _find_nearest_with_data is called when enabled and stations are missing."""
        era5 = self._get_era5_manager()
        metadf = self._get_test_metadf()

        mock_gee_api._df_to_features_point_collection.return_value = MagicMock()
        mock_gee_api._estimate_data_size.return_value = 100
        mock_gee_api.get_ee_obj.return_value = MagicMock()
        mock_gee_api.datetime_to_gee_datetime.return_value = MagicMock()
        mock_gee_api._addDate = MagicMock()

        # Only return data for 2 out of 3 stations
        mock_results = {
            "features": [
                {
                    "properties": {
                        "name": "station_a",
                        "datetime": "20220901120000",
                        "temperature_2m": 293.15,
                    }
                },
                {
                    "properties": {
                        "name": "station_b",
                        "datetime": "20220901120000",
                        "temperature_2m": 294.15,
                    }
                },
            ]
        }

        mock_raster = MagicMock()
        mock_gee_api.get_ee_obj.return_value = mock_raster
        mock_raster.filter.return_value.map.return_value.map.return_value.flatten.return_value.getInfo.return_value = (
            mock_results
        )

        # Patch _find_nearest_with_data to return None (simulating no recovery)
        with patch.object(
            GEEDynamicDatasetManager, "_find_nearest_with_data", return_value=None
        ) as mock_find:
            with caplog.at_level(logging.WARNING):
                df = era5.extract_timeseries_data(
                    metadf=metadf,
                    startdt_utc=pd.Timestamp("2022-09-01", tz="UTC"),
                    enddt_utc=pd.Timestamp("2022-09-01 12:00:00", tz="UTC"),
                    obstypes=["temp"],
                    get_all_bands=False,
                    force_direct_transfer=True,
                    initialize_gee=True,
                    find_nearest_with_data=True,
                )

            # Verify _find_nearest_with_data was called
            mock_find.assert_called_once()
            # Verify it was called with the missing station's metadf
            call_args = mock_find.call_args
            called_metadf = call_args[1]["metadf"]
            assert "station_c" in called_metadf.index.tolist()

    def test_find_nearest_with_data_parameter_exists_on_station(self):
        """Test that find_nearest_with_data parameter exists on Station.get_gee_timeseries_data."""
        import inspect
        from metobs_toolkit.station import Station

        sig = inspect.signature(Station.get_gee_timeseries_data)
        assert "find_nearest_with_data" in sig.parameters

    def test_find_nearest_with_data_parameter_exists_on_dataset(self):
        """Test that find_nearest_with_data parameter exists on Dataset.get_gee_timeseries_data."""
        import inspect

        sig = inspect.signature(metobs_toolkit.Dataset.get_gee_timeseries_data)
        assert "find_nearest_with_data" in sig.parameters


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
