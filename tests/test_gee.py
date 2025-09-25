"""This file test the gee interactions using the pytest framework."""

import pytest
import sys
from pathlib import Path
import tempfile

# import metobs_toolkit
import pandas as pd
import numpy as np
import folium
import geemap.foliumap as geemap

# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir


class TestDemoDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodataset"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_import_demo_metadata(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        pass
        # 2. apply a metobs manipulation
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
        )
        data_to_test = dataset

        # test get info on metadata-only dataset
        dataset.get_info(printout=False)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=data_to_test,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # DEBUG
        test = dataset.modeldatadf
        # 5. Construct the equlity tests
        assert_equality(data_to_test, solutionobj)

    def test_gee_connect(self):
        metobs_toolkit.connect_to_gee()

    def test_LCZ_extraction(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )
        LCZ_data = dataset.get_LCZ()

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

        assert_equality(LCZ_data["LCZ"], solutionobj.metadf["LCZ"])
        # calling printoutlc
        _ = dataset.get_station("vlinder18").site.get_info(printout=False)
        assert isinstance(dataset.get_station("vlinder18").site.LCZ, str)

    def test_lcz_seamask_fix(self):
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )

        seastation = dataset.get_station("vlinder15")
        seastation.site._lat = 51.361852
        seastation.site._lon = 3.009151
        # with mask fix
        lcz_return = seastation.get_LCZ(apply_seamask_fix=True, overwrite=True)
        assert lcz_return == "Water (LCZ G)"
        assert seastation.site.LCZ == "Water (LCZ G)"  # LCZ-G is the water LCZ class

        # without mask
        lcz_return = seastation.get_LCZ(apply_seamask_fix=False, overwrite=True)
        assert np.isnan(lcz_return)

        # Now on dataset level
        dataset.stations[5].site._lat = 51.361852
        dataset.stations[5].site._lon = 3.009151

        lczdf = dataset.get_LCZ(apply_seamask_fix=True, overwrite=True)
        assert lczdf["LCZ"].notna().all()
        assert lczdf["LCZ"].eq("Water (LCZ G)").any()
        assert dataset.stations[5].site.LCZ == "Water (LCZ G)"

    def test_altitude_extraction(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )
        alt_data = dataset.get_altitude()

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

        assert_equality(alt_data["altitude"].astype('float64'), solutionobj.metadf["altitude"].astype('float64'))
        # calling printout
        _ = dataset.get_station("vlinder18").site.get_info(printout=False)
        assert isinstance(dataset.get_station("vlinder18").site.altitude, float)

    def test_geemodeldata_getinfo_class(self):
        """Simple tests on GeeModelData instances"""
        for knowngee in metobs_toolkit.default_GEE_datasets.keys():
            model = metobs_toolkit.default_GEE_datasets[knowngee]

            # Try calling all possible methods
            _ = model.get_info(printout=False)

    def test_gee_static_plot(self):
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )
        for geemod in metobs_toolkit.default_GEE_datasets.values():
            if not isinstance(geemod, metobs_toolkit.GEEStaticDatasetManager):
                continue
            mapret = dataset.make_gee_plot(geedatasetmanager=geemod)
            assert type(mapret) == geemap.Map

    def test_gee_dynamic_plot(self):
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )

        geemod = metobs_toolkit.default_GEE_datasets["ERA5-land"]

        with pytest.raises(ValueError):
            mapret = dataset.make_gee_plot(
                geedatasetmanager=geemod,
                modelobstype="temp",
                timeinstance=None,  # must raise error
            )

        timeinstance = pd.Timestamp("2021-01-01 18:12:00")

        from metobs_toolkit.backend_collection.errorclasses import MetObsObstypeNotFound

        with pytest.raises(MetObsObstypeNotFound):
            mapret = dataset.make_gee_plot(
                geedatasetmanager=geemod,
                modelobstype="fakeobstype",  # must raise error
                timeinstance=timeinstance,
            )

        mapret = dataset.make_gee_plot(
            geedatasetmanager=geemod,
            modelobstype="temp",
            timeinstance=timeinstance,
        )

        assert type(mapret) == geemap.Map

    def test_landcover_frac_extraction(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )

        landcover_data = dataset.get_landcover_fractions(
            buffers=[10, 100, 500], aggregate=False
        )

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

        assert isinstance(landcover_data, pd.DataFrame)
        # calling printout
        _ = dataset.get_station("vlinder18").site.get_info(printout=False)
        _ = dataset.metadf
        assert isinstance(
            dataset.get_station("vlinder18").site.buffered_fractions, dict
        )
        assert set(
            dataset.get_station("vlinder18").site.buffered_fractions.keys()
        ) == set([10, 100, 500])

    def test_ERA5_extraction_on_metadata_only(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = (
            "test_ERA5_extraction_on_metadata_only"  # get the name of this method
        )
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )

        era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]

        from metobs_toolkit.backend_collection.errorclasses import MetObsMissingArgument

        with pytest.raises(MetObsMissingArgument):
            era5_data = dataset.get_gee_timeseries_data(
                geedynamicdatasetmanager=era5_model,
                startdt_utc=None,  # raises error in metadata-only case
                enddt_utc=None,
                target_obstypes=["temp"],
                get_all_bands=False,
                drive_filename=None,
                drive_folder="gee_timeseries_data",
                force_direct_transfer=False,
                force_to_drive=False,
            )

        startdt_utc = pd.Timestamp("2021-01-01 16:32:25")
        enddt_utc = pd.Timestamp("2021-01-01 23:16:00")
        era5_data = dataset.get_gee_timeseries_data(
            geedynamicdatasetmanager=era5_model,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            target_obstypes=["temp", "pressure", "wind"],
            get_all_bands=False,
            drive_filename=None,
            drive_folder="gee_timeseries_data",
            force_direct_transfer=False,
            force_to_drive=False,
        )

        assert era5_data.shape == (224, 4)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    def test_get_info_on_modeltimeseries(self):
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )

        # test get info on metadata-only dataset
        for modeltimeseries in dataset.get_station("vlinder18").modeldata:
            _ = modeltimeseries.get_info(printout=False)

    def test_ERA5_extraction(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_ERA5_extraction"  # get the name of this method
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_data_file=metobs_toolkit.demo_datafile,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
        )
        dataset.resample(target_freq="15min")

        era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]

        startdt_utc = pd.Timestamp("2022-08-31 18:32:25")
        enddt_utc = pd.Timestamp("2022-09-01 12:16:00")
        era5_data = dataset.get_gee_timeseries_data(
            geedynamicdatasetmanager=era5_model,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            target_obstypes=["temp"],
            get_all_bands=False,
            drive_filename=None,
            drive_folder="gee_timeseries_data",
            force_direct_transfer=True,
            force_to_drive=False,
        )

        assert era5_data.shape == (532, 1)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    def test_station_timeseries_extraction(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_data_file=metobs_toolkit.demo_datafile,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
        )
        dataset.resample(target_freq="15min")

        era5_manager = metobs_toolkit.default_GEE_datasets["ERA5-land"]
        # Extract the timeseries
        era5_temp = dataset.get_station("vlinder02").get_gee_timeseries_data(
            geedynamicdatasetmanager=era5_manager,  # The datasetmanager to use
            startdt_utc=None,
            enddt_utc=None,
            target_obstypes=[
                "temp"
            ],  # the observationtypes to extract, must be known modelobstypes
            force_direct_transfer=True,
        )

        assert dataset.get_station("vlinder02").modeldatadf.shape == (361, 4)

    def test_pickling(self):
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            # pickle dataset
            dataset.save_dataset_to_pkl(
                target_folder=tmpdir, filename="deleteme.pkl", overwrite=True
            )

            # open datast
            dataset_pkled = metobs_toolkit.import_dataset_from_pkl(
                target_path=tmpdir.joinpath("deleteme.pkl")
            )

        assert_equality(dataset_pkled, dataset)

    def test_era5_modeldata_interactions(self):
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )

        # testing the dataset (with modeldata)
        _ = dataset.get_station("vlinder02").get_info(printout=False)
        _ = dataset.get_info(printout=False)
        _ = (
            dataset.get_station("vlinder02")
            .get_modeltimeseries("temp")
            .get_info(printout=False)
        )

        assert dataset.modeldatadf.shape == (532, 4)
        assert dataset.get_station("vlinder02").modeldatadf.shape == (19, 4)
        assert dataset.get_station("vlinder02").get_modeltimeseries(
            "temp"
        ).series.shape == (19,)

    def test_ERA5_google_drive_interface(self):
        # 1. Test writing to drive file
        # Extract ERA5data and force the storing in drive
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )
        startdt_utc = pd.Timestamp("2021-01-01 16:32:25")
        enddt_utc = pd.Timestamp("2021-01-01 23:16:00")
        era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]
        era5_data = dataset.get_gee_timeseries_data(
            geedynamicdatasetmanager=era5_model,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            target_obstypes=["temp", "pressure", "wind"],
            get_all_bands=False,
            drive_filename=None,
            drive_folder="gee_timeseries_data",
            force_direct_transfer=False,
            force_to_drive=True,
        )  # FORCE DRIVE!

        assert era5_data is None
        # test reading csv file
        target_era5_csv = datadir.joinpath(
            "ERA5-land_timeseries_data_of_full_dataset_28_stations.csv"
        )
        dataset.import_gee_data_from_file(
            filepath=target_era5_csv,
            geedynamicdatasetmanager=era5_model,
            force_update=True,
        )
        # compare with solution of direct import of gee data
        dataset_direct = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs,
            methodname="test_ERA5_extraction_on_metadata_only",
        )

        # test equality
        assert_equality(dataset, dataset_direct)

    def test_modeldata_timeseries_plot(self):
        # 1. get_startpoint data WITH records
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )

        station = dataset.get_station("vlinder05")
        station.make_plot(show_modeldata=True)


if __name__ == "__main__":
    test = TestDemoDataset()
    # test.test_ERA5_extraction_on_metadata_only(overwrite_solution=False)
    # test.test_import_demo_metadata(overwrite_solution=False)
    # test.test_LCZ_extraction(overwrite_solution=False)
    # test.test_altitude_extraction(overwrite_solution=False)
    # test.test_landcover_frac_extraction(overwrite_solution=False)
    # test.test_ERA5_extraction_on_metadata_only(overwrite_solution=False)
    # test.test_ERA5_extraction(overwrite_solution=False)
    # test.test_ERA5_extraction(overwrite_solution=False)
