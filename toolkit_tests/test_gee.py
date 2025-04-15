""" This file test the gee interactions using the pytest framework."""

import pytest
import sys
from pathlib import Path

# import metobs_toolkit
import pandas as pd
import folium
import geemap.foliumap as geemap

libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer
import shutil
import pytest

from gee_service_authenticator import GEE_Authenticator

# authenticate the service account
GEE_Authenticator()

# data folder
datadir = libfolder.joinpath("tests").joinpath("test_data")


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
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # dataset comparison

        # 6. assert the equality
        assert test_expr

    def test_gee_connect(self):
        metobs_toolkit.connect_to_gee()

    def test_lcz_extraction(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_demo_metadata"
        )
        lcz_data = dataset.get_lcz()

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = dataset == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset, solution=solutionobj
            )

        assert dataset == solutionobj
        assert lcz_data["lcz"].equals(solutionobj.metadf["lcz"])
        # calling printoutlc
        _ = dataset.get_station("vlinder18").site.get_info(printout=False)
        assert isinstance(dataset.get_station("vlinder18").site.lcz, str)

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
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = dataset == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset, solution=solutionobj
            )

        assert dataset == solutionobj
        assert alt_data["altitude"].equals(solutionobj.metadf["altitude"])
        # calling printout
        _ = dataset.get_station("vlinder18").site.get_info(printout=False)
        assert isinstance(dataset.get_station("vlinder18").site.altitude, int)

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
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = dataset == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset, solution=solutionobj
            )

        assert dataset == solutionobj

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
        _method_name = sys._getframe().f_code.co_name  # get the name of this method
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
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = dataset == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr

    def test_ERA5_extraction(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method
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
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = dataset == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr

    def test_pickling(self):

        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )

        # pickle dataset
        dataset.save_dataset_to_pkl(
            target_folder=datadir, filename="deleteme.pkl", overwrite=True
        )

        # open datast
        dataset_pkled = metobs_toolkit.import_dataset_from_pkl(
            target_path=datadir.joinpath("deleteme.pkl")
        )

        test_expr = dataset_pkled == dataset
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset_pkled, solution=dataset
            )
        # 6. assert the equality
        assert test_expr

    def test_era5_modeldata_interactions(self):
        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )

        # testing the dataset (with modeldata)
        _ = dataset.get_station("vlinder02").get_info(printout=False)
        _ = dataset.get_info(printout=False)
        _ = dataset.get_station("vlinder02").modeldata["temp"].get_info(printout=False)

        assert dataset.modeldatadf.shape == (532, 2)
        assert dataset.get_station("vlinder02").modeldatadf.shape == (19, 2)

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
            methodname="test_ERA5_extraction_on_metadata_only"
        )

        # test equality
        test_expr = dataset == dataset_direct
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=dataset, solution=dataset_direct
            )
        # 6. assert the equality
        assert test_expr

    def test_modeldata_timeseries_plot(self):
        # 1. get_startpoint data WITH records
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_ERA5_extraction"
        )

        station = dataset.get_station("vlinder05")
        station.make_plot(show_modeldata=True)


# TODO: plot test of modeldata on dataset level


if __name__ == "__main__":

    test = TestDemoDataset()
    # test.test_import_demo_metadata(overwrite_solution=False)
    # test.test_lcz_extraction(overwrite_solution=False)
    # test.test_altitude_extraction(overwrite_solution=False)
    # test.test_landcover_frac_extraction(overwrite_solution=False)
    # test.test_ERA5_extraction_on_metadata_only(overwrite_solution=False)
    # test.test_ERA5_extraction(overwrite_solution=False)
