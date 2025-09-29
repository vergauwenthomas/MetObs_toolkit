#!/usr/bin/env python3
"""
Test for the modeldatadf bug fix
"""

import pytest
import sys
import numpy as np
import pandas as pd
from pathlib import Path

# Add the local source directory to Python path
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))

from metobs_toolkit.modeltimeseries import ModelTimeSeries
from metobs_toolkit.obstypes import ModelObstype, Obstype
from metobs_toolkit.site import Site
from metobs_toolkit.backend_collection.dataframe_constructors import modeltimeseries_df


class TestModelDatadfBugFix:
    """Test the fix for the modeldatadf bug."""

    def test_modeltimeseries_df_normal_case(self):
        """Test that modeltimeseries_df works with a normal ModelTimeSeries instance."""
        # Create components
        temp_obstype = Obstype(
            obsname="temp", std_unit="degree_Celsius", description="Temperature"
        )
        model_obstype = ModelObstype(
            obstype=temp_obstype, model_unit="kelvin", model_band="t2m"
        )
        site = Site(stationname="test_station", latitude=50.0, longitude=4.0)

        # Create test data
        timestamps = pd.date_range("2023-01-01", periods=5, freq="1h")
        data_records = np.array([273.15, 274.15, 275.15, 276.15, 277.15])

        # Create ModelTimeSeries
        model_ts = ModelTimeSeries(
            site=site,
            datarecords=data_records,
            timestamps=timestamps.to_numpy(),
            modelobstype=model_obstype,
            modelname="test_model",
            modelvariable="temperature",
        )

        # Test the df property
        df = model_ts.df
        
        # Verify the structure
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (5, 2)
        assert list(df.columns) == ["value", "model"]
        assert (df["model"] == "test_model").all()
        
        # Verify the values are converted from Kelvin to Celsius
        expected_celsius = np.array([0.0, 1.0, 2.0, 3.0, 4.0])
        np.testing.assert_array_almost_equal(df["value"].values, expected_celsius)

    def test_modeltimeseries_df_missing_modelobstype(self):
        """Test that modeltimeseries_df gives a helpful error for missing modelobstype."""
        # Create a normal ModelTimeSeries first
        temp_obstype = Obstype(
            obsname="temp", std_unit="degree_Celsius", description="Temperature"
        )
        model_obstype = ModelObstype(
            obstype=temp_obstype, model_unit="kelvin", model_band="t2m"
        )
        site = Site(stationname="test_station", latitude=50.0, longitude=4.0)

        timestamps = pd.date_range("2023-01-01", periods=3, freq="1h")
        data_records = np.array([273.15, 274.15, 275.15])

        model_ts = ModelTimeSeries(
            site=site,
            datarecords=data_records,
            timestamps=timestamps.to_numpy(),
            modelobstype=model_obstype,
            modelname="test_model",
            modelvariable="temperature",
        )

        # Remove the modelobstype attribute to simulate the bug
        delattr(model_ts, "modelobstype")

        # Test that accessing df property gives a helpful error message
        with pytest.raises(AttributeError) as exc_info:
            _ = model_ts.df

        assert "ModelTimeSeries instance is missing 'modelobstype' attribute" in str(
            exc_info.value
        )
        assert "incomplete initialization or deserialization" in str(exc_info.value)

    def test_modeltimeseries_df_direct_call(self):
        """Test calling modeltimeseries_df function directly."""
        # Create a normal ModelTimeSeries
        temp_obstype = Obstype(
            obsname="temp", std_unit="degree_Celsius", description="Temperature"
        )
        model_obstype = ModelObstype(
            obstype=temp_obstype, model_unit="kelvin", model_band="t2m"
        )
        site = Site(stationname="test_station", latitude=50.0, longitude=4.0)

        timestamps = pd.date_range("2023-01-01", periods=3, freq="1h")
        data_records = np.array([273.15, 274.15, 275.15])

        model_ts = ModelTimeSeries(
            site=site,
            datarecords=data_records,
            timestamps=timestamps.to_numpy(),
            modelobstype=model_obstype,
            modelname="test_model",
            modelvariable="temperature",
        )

        # Test direct function call
        df = modeltimeseries_df(model_ts)
        
        assert isinstance(df, pd.DataFrame)
        assert df.shape == (3, 2)
        assert list(df.columns) == ["value", "model"]
        assert (df["model"] == "test_model").all()