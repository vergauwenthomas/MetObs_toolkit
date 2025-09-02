#!/usr/bin/env python3

from pathlib import Path
import tempfile
import os
import pytest
import sys

# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

import xarray as xr

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions") 
from solutionclass import SolutionFixer, assert_equality, datadir


class TestToNetCDF:
    """Test the to_netcdf functionality for Dataset and Station classes."""
    
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testtonetcdf"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_station_to_netcdf_basic(self):
        """Test basic station to_netcdf functionality without model data."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        station = dataset.get_station("vlinder05")

        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
            try:
                # Test saving to netCDF
                station.to_netcdf(tmp.name)
                
                # Verify file was created and is readable
                assert os.path.exists(tmp.name)
                ds_read = xr.open_dataset(tmp.name)
                
                # Check that variables were split into obs and label
                variables = list(ds_read.data_vars.keys())
                obs_vars = [v for v in variables if v.endswith('_obs')]
                label_vars = [v for v in variables if v.endswith('_label')]
                
                assert len(obs_vars) > 0, "Should have observation variables"
                assert len(label_vars) > 0, "Should have label variables"
                
                # Check that attributes are flattened (no nested dicts)
                for var_name in variables:
                    var = ds_read[var_name]
                    for attr_name, attr_val in var.attrs.items():
                        assert not isinstance(attr_val, dict), f"Attribute {attr_name} in {var_name} is still a dict"
                
                # Check that QC attributes are preserved in flattened form
                qc_attrs_found = False
                for var_name in variables:
                    var = ds_read[var_name]
                    for attr_name in var.attrs.keys():
                        if 'QC.' in attr_name:
                            qc_attrs_found = True
                            break
                
                assert qc_attrs_found, "QC attributes should be preserved in flattened form"
                
                # Check timezone handling on datetime coordinate
                if 'datetime' in ds_read.coords:
                    dt_coord = ds_read.coords['datetime']
                    # Should be timezone-naive datetime64[ns]
                    assert 'datetime64[ns]' == str(dt_coord.dtype), f"Datetime should be naive, got {dt_coord.dtype}"
                    # Should have timezone attribute
                    assert 'timezone' in dt_coord.attrs, "Should have timezone attribute"
                
                ds_read.close()
                
            finally:
                if os.path.exists(tmp.name):
                    os.unlink(tmp.name)

    def test_dataset_to_netcdf_basic(self):
        """Test basic dataset to_netcdf functionality."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
            try:
                # Test saving to netCDF
                dataset.to_netcdf(tmp.name)
                
                # Verify file was created and is readable
                assert os.path.exists(tmp.name)
                ds_read = xr.open_dataset(tmp.name)
                
                # Check that we have the expected structure
                variables = list(ds_read.data_vars.keys())
                assert len(variables) > 0, "Should have variables"
                
                # Should have 'name' dimension for multiple stations
                assert 'name' in ds_read.dims, "Should have 'name' dimension for stations"
                
                ds_read.close()
                
            finally:
                if os.path.exists(tmp.name):
                    os.unlink(tmp.name)

    def test_station_to_netcdf_with_model_data(self):
        """Test station to_netcdf functionality with model data if available."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        
        # Try to load model data if available
        target_era5_csv = datadir.joinpath(
            "ERA5-land_timeseries_data_of_full_dataset_28_stations.csv"
        )
        
        has_model_data = False
        if target_era5_csv.exists():
            era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]
            dataset.import_gee_data_from_file(
                filepath=str(target_era5_csv),
                geedynamicdatasetmanager=era5_model,
                force_update=True,
            )
            has_model_data = True

        station = dataset.get_station("vlinder05")

        with tempfile.NamedTemporaryFile(suffix='.nc', delete=False) as tmp:
            try:
                # Test saving to netCDF
                station.to_netcdf(tmp.name)
                
                # Verify file was created and is readable
                assert os.path.exists(tmp.name)
                ds_read = xr.open_dataset(tmp.name)
                
                variables = list(ds_read.data_vars.keys())
                
                if has_model_data:
                    # Check for flattened model attributes
                    model_attrs_found = False
                    for var_name in variables:
                        var = ds_read[var_name]
                        for attr_name in var.attrs.keys():
                            if 'ERA5-land.' in attr_name:
                                model_attrs_found = True
                                break
                    
                    assert model_attrs_found, "Model attributes should be preserved in flattened form"
                
                ds_read.close()
                
            finally:
                if os.path.exists(tmp.name):
                    os.unlink(tmp.name)

    def test_serialization_helpers(self):
        """Test the helper functions for making datasets serializable."""
        from metobs_toolkit.xrconversions import flatten_nested_dict, make_dataset_serializable
        
        # Test nested dict flattening
        nested = {
            "QC": {"gross_value": {"settings": {"threshold": 10}}},
            "GF": {},
            "model": {"ERA5": {"name": "test"}}
        }
        
        flattened = flatten_nested_dict(nested)
        
        assert "QC.gross_value.settings.threshold" in flattened
        assert flattened["QC.gross_value.settings.threshold"] == "10"
        assert "GF" in flattened
        assert flattened["GF"] == "{empty_dict}"
        assert "model.ERA5.name" in flattened
        assert flattened["model.ERA5.name"] == "test"
        
        # Test with empty dict
        empty_nested = {"empty": {}}
        flattened_empty = flatten_nested_dict(empty_nested)
        assert flattened_empty["empty"] == "{empty_dict}"