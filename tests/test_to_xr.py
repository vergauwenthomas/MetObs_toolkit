from pathlib import Path
import pytest
import sys

# import metobs_toolkit
import pandas as pd
import numpy as np
import xarray as xr


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir


class TestDemoData:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_to_xr_on_station(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        station = dataset.get_station("vlinder05")

        station.repetitions_check(max_N_repetitions=200)

        ds = station.to_xr()
        _ = station.get_sensor("temp").to_xr()

        # Check type
        assert isinstance(ds, xr.core.dataset.Dataset)

        # Check dimensions
        assert "datetime" in ds.dims
        assert "kind" in ds.dims
        assert "name" not in ds.dims

        # test datetime
        assert (
            len(ds["datetime"])
            == station.get_sensor("temp").df.index.get_level_values("datetime").shape[0]
        )
        assert ds["datetime"].attrs["standard_name"] == "time"
        assert "timezone" in ds["datetime"].attrs
        assert ds["datetime"].dtype == "datetime64[ns]"

        for obsname in station.obsdata.keys():
            # test observation variables
            assert obsname in ds.data_vars
            var = ds[obsname]
            assert var.data.shape == (
                2,
                station.get_sensor("temp")
                .df.index.get_level_values("datetime")
                .shape[0],
            )
            assert "obstype_name" in var.attrs
            assert "obstype_desc" in var.attrs
            assert "obstype_unit" in var.attrs
            assert "QC checks" in var.attrs
            assert "GF methods" in var.attrs

            _ = var.sel(kind="obs")
            _ = var.sel(kind="label")

        # Test presence of site cooridnates
        assert int(ds["lat"].data) == 51
        assert int(ds["lon"].data) == 3
        assert str(ds["school"].data) == "Sint-Barbara"

        # Test presents of QC labels

        assert "repetitions" in ds["temp"].attrs["QC checks"]
        assert "repetitions" not in ds["humidity"].attrs["QC checks"]
        assert "QC:repetitions.max_N_repetitions" in ds["temp"].attrs
        assert ds["temp"].attrs["QC:repetitions.max_N_repetitions"] == 200

        # test label conversions
        assert ds["temp"].attrs["Label:ok"] == 0
        assert ds["temp"].attrs["Label:repetitions outlier"] == 5
        assert "Label:repetitions outlier" not in ds["humidity"].attrs

        # Test presence of GF labels
        station.convert_outliers_to_gaps(obstype="temp")
        station.interpolate_gaps(
            target_obstype="temp", max_gap_duration_to_fill='2d', overwrite_fill=True
        )
        ds = station.to_xr()

        assert 13 in ds["temp"].sel(kind="label").data
        assert 12 in ds["temp"].sel(kind="label").data  # interpolated
        assert "interpolation" in ds["temp"].attrs["GF methods"]
        assert ds["temp"].attrs["GF:interpolation.method"] == "time"
        assert '2 days' in ds["temp"].attrs["GF:interpolation.max_gap_duration_to_fill"]
        assert (
            ds["temp"].attrs["GF:interpolation.n_leading_anchors"] == 1
        )  # test if default arguments are present

    def test_to_xr_on_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        station = dataset.get_station("vlinder05")

        station.repetitions_check(max_N_repetitions=200)

        ds = dataset.to_xr()

        # Check dimensions
        assert "datetime" in ds.dims
        assert "kind" in ds.dims
        assert "name" in ds.dims

    def test_to_xr_with_modeldata(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        target_era5_csv = datadir.joinpath(
            "ERA5-land_timeseries_data_of_full_dataset_28_stations.csv"
        )
        era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]
        dataset.import_gee_data_from_file(
            filepath=target_era5_csv,
            geedynamicdatasetmanager=era5_model,
            force_update=True,
        )

        station = dataset.get_station("vlinder05")

        station.repetitions_check(max_N_repetitions=200)

        ds = dataset.to_xr()

        # Check dimensions
        assert "datetime" in ds.dims
        assert "kind" in ds.dims
        assert "name" in ds.dims

        var = ds["temp"]

        # Test that 'kind' dimension has 3 elements (obs, label, model)
        assert len(ds.coords["kind"]) == 3
        assert "obs" in ds.coords["kind"].values
        assert "label" in ds.coords["kind"].values
        assert "model" in ds.coords["kind"].values

        # Test var attributes when model data is present
        # obstype related
        assert "obstype_name" in var.attrs
        assert "obstype_desc" in var.attrs
        assert "obstype_unit" in var.attrs
        assert "QC checks" in var.attrs
        assert "GF methods" in var.attrs

        # modelobstype related

        assert "modelobstype_name" in var.attrs
        assert "modelobstype_desc" in var.attrs
        assert "modelobstype_unit" in var.attrs
        assert var.attrs["modelname"] == "ERA5-land"
        assert var.attrs["modelvariable"] == "temperature_2m"

        # Test that observations and labels still exist alongside model data
        obs_temp = ds["temp"].sel(kind="obs")
        label_temp = ds["temp"].sel(kind="label")
        assert obs_temp.data.shape == label_temp.data.shape
        assert obs_temp.data.shape == (28, 4328, 1)

        # Test model coordinate attributes
        assert "models" in ds.coords
        assert "ERA5-land" in ds.coords["models"].values
        assert len(ds.coords["models"]) == 1

        # Test that var has the correct kind dimension
        assert set(var.dims) == set(("name", "kind", "models", "datetime"))

        assert ds["temp"].sel(kind="model").sel(models="ERA5-land").data.shape == (
            28,
            4328,
        )

    def test_station_to_netcdf(self):
        """Test Station.to_netcdf() method."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        station = dataset.get_station("vlinder05")
        station.repetitions_check(max_N_repetitions=200)

        # Test saving to netCDF
        import tempfile
        import os

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_station.nc")
            station.to_netcdf(filepath)

            # Verify file was created
            assert os.path.exists(filepath)

            # Verify file can be read back as xarray Dataset
            ds_from_file = xr.open_dataset(filepath)

            # Basic checks
            assert isinstance(ds_from_file, xr.core.dataset.Dataset)
            assert "datetime" in ds_from_file.dims
            assert "kind" in ds_from_file.dims
            assert "temp" in ds_from_file.data_vars
            assert "humidity" in ds_from_file.data_vars

            ds_from_file.close()

    def test_dataset_to_netcdf(self):
        """Test Dataset.to_netcdf() method."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        station = dataset.get_station("vlinder05")
        station.repetitions_check(max_N_repetitions=200)

        # Test saving to netCDF
        import tempfile
        import os

        with tempfile.TemporaryDirectory() as tmpdir:
            filepath = os.path.join(tmpdir, "test_dataset.nc")
            dataset.to_netcdf(filepath)

            # Verify file was created
            assert os.path.exists(filepath)

            # Verify file can be read back as xarray Dataset
            ds_from_file = xr.open_dataset(filepath)

            # Basic checks
            assert isinstance(ds_from_file, xr.core.dataset.Dataset)
            assert "datetime" in ds_from_file.dims
            assert "kind" in ds_from_file.dims
            assert "name" in ds_from_file.dims  # Dataset should have name dimension
            assert "temp" in ds_from_file.data_vars
            assert "humidity" in ds_from_file.data_vars

            ds_from_file.close()
