import pytest
import sys
from pathlib import Path

# import metobs_toolkit
import pandas as pd
import numpy as np
import xarray as xr


libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
# sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir
import shutil
import pytest


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
        
        station = dataset.get_station('vlinder05')
        
        station.repetitions_check(max_N_repetitions=200)
        
        ds = station.to_xr()
        _ = station.get_sensor('temp').to_xr()
        
        # Check type
        assert isinstance(ds, xr.core.dataset.Dataset)
       
        # Check dimensions
        assert "datetime" in ds.dims
        assert 'kind' in ds.dims
        assert 'name' not in ds.dims
        
        assert len(ds['datetime']) == station.get_sensor('temp').df.index.get_level_values('datetime').shape[0]
        
        for obsname in station.obsdata.keys():
            #test observation variables
            assert obsname in ds.data_vars
            var = ds[obsname]
            assert var.data.shape == (2, station.get_sensor('temp').df.index.get_level_values('datetime').shape[0])
            assert 'obstype_name' in var.attrs
            assert 'obstype_desc' in var.attrs
            assert 'obstype_unit' in var.attrs
            assert 'QC' in var.attrs
            assert 'GF' in var.attrs
            
            _ = var.sel(kind='obs')
            _ = var.sel(kind='label')
            

            
        
        #Test presence of site cooridnates
        assert int(ds['lat'].data) == 51
        assert int(ds['lon'].data) == 3
        assert str(ds['school'].data) == 'Sint-Barbara'        
        
        #Test presents of QC labels
       
        assert 'repetitions outlier' in ds['temp'].sel(kind='label').data
        assert 'repetitions' in ds['temp'].attrs['QC']
        assert ds['temp'].attrs['QC']['repetitions']['settings']['max_N_repetitions'] == 200
        
        
        #Test preses of GF labels
        station.convert_outliers_to_gaps(obstype='temp')
        station.interpolate_gaps(target_obstype='temp', max_consec_fill=500, overwrite_fill=True)
        ds = station.to_xr()
        
       
        assert 'failed interpolation' in ds['temp'].sel(kind='label').data
        assert 'interpolation' in ds['temp'].sel(kind='label').data
        assert ds['temp'].attrs['GF']['interpolation']['method'] == 'time'
        assert ds['temp'].attrs['GF']['interpolation']['max_consec_fill'] == 500
        
        
        
    
    def test_to_xr_on_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        
        station = dataset.get_station('vlinder05')
        
        station.repetitions_check(max_N_repetitions=200)
        
        ds = dataset.to_xr()
        
        # Check dimensions
        assert "datetime" in ds.dims
        assert 'kind' in ds.dims
        assert 'name' in ds.dims
        
        
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
        
        station = dataset.get_station('vlinder05')
        
        station.repetitions_check(max_N_repetitions=200)
        
        ds = dataset.to_xr()
        
        # Check dimensions
        assert "datetime" in ds.dims
        assert 'kind' in ds.dims
        assert 'name' in ds.dims
        
        assert ds['temp'].sel(kind='model').sel(models='ERA5-land').data.shape == (28, 4328)
        assert 'ERA5-land' in ds['temp'].attrs
        
        
        
       
            
            
        
        
           
        
 
        

    