import pytest
import sys
from pathlib import Path

# import metobs_toolkit
import pandas as pd


libfolder = Path(str(Path(__file__).resolve()).split("MetObs_toolkit")[0]).joinpath(
    "MetObs_toolkit"
)
# point to current version of the toolkit
sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer
import shutil
import pytest

# data folder
datadir = libfolder.joinpath("tests").joinpath("test_data")


class TestDemoData:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_version(self):
        assert isinstance(metobs_toolkit.__version__, str)

    def test_import_demo_data(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        pass

        # 2. apply a metobs manipulation
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        data_to_test = dataset

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solutiondata=data_to_test,
                methodname=_method_name,
                **TestDemoData.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoData.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoData.solutionfixer.create_a_diff(
                to_check=data_to_test, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr

    def test_calling_methods_without_solution_on_dataset(self):
        # 1. get_startpoint data
        dataset = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname="test_import_demo_data"
        )
        # run methods, and see if something breaks

        # testing specials
        _ = dataset.stations
        _ = dataset.obstypes
        _ = dataset.template
        _ = dataset.df
        _ = dataset.outliersdf
        _ = dataset.metadf
        _ = dataset.stations
        _ = dataset.start_datetime
        _ = dataset.end_datetime

        # get info's
        dataset.get_info(printout=False)
        dataset.template.get_info()

        # random collection
        dataset.rename_stations(
            renamedict={"vlinder01": "fakename", "vlider14": "fakename2"}
        )

    def test_calling_methods_without_solution_on_station(self):
        # 1. get_startpoint data
        dataset = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname="test_import_demo_data"
        )
        # run methods, and see if something breaks
        station = dataset.stations[0]

        # testing specials

        _ = station.df
        _ = station.outliersdf
        _ = station.metadf
        _ = station.start_datetime
        _ = station.end_datetime

        # get info's
        station.get_info()
        station.site.get_info()
        station.obsdata["temp"].get_info()

    def test_get_info(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname="test_import_demo_data"
        )
        # 2. apply a metobs manipulation
        data_to_test = dataset.get_info(printout=False)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solutiondata=data_to_test,
                **TestDemoData.solkwargs,
                methodname=_method_name
            )
        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # string comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoData.solutionfixer.create_a_diff(
                to_check=data_to_test, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr

    def test_get_station(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname="test_import_demo_data"
        )
        # 2. apply a metobs manipulation
        data_to_test = dataset.get_station("vlinder05")

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solutiondata=data_to_test,
                **TestDemoData.solkwargs,
                methodname=_method_name
            )

        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # Station comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoData.solutionfixer.create_a_diff(
                to_check=data_to_test, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr

    def test_pickling_dataset(self):
        # 0. Get info of the current check
        # _method_name = sys._getframe().f_code.co_name #get the name of this method

        # 1. get_startpoint data
        dataset = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname="test_import_demo_data"
        )

        # Create a tmp dir
        tmpdir = libfolder.joinpath("tmp")
        tmpdir.mkdir(parents=True, exist_ok=True)
        # pickle dataset
        dataset.save_dataset_to_pkl(target_folder=tmpdir, filename="deleteme")
        # Read in the pickled dataset
        dataset2 = metobs_toolkit.import_dataset_from_pkl(
            target_path=tmpdir.joinpath("deleteme.pkl")
        )
        # Remove the tmp dir
        shutil.rmtree(tmpdir)

        # test if the pickled dataset is equal to the original
        assert dataset == dataset2


class TestWideData:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testwidedata"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("wide_test_data.csv")
    templatefile = datadir.joinpath("wide_test_template.json")

    def test_import_wide_data(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        pass

        # 2. apply a metobs manipulation
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestWideData.templatefile,
            input_metadata_file=None,
            input_data_file=TestWideData.datafile,
            freq_estimation_method="median",
            freq_estimation_simplify_tolerance="2min",
            origin_simplify_tolerance="5min",
            timestamp_tolerance="4min",
        )

        data_to_test = dataset

        # 3. overwrite solution?
        if overwrite_solution:
            TestWideData.solutionfixer.create_solution(
                solutiondata=data_to_test,
                methodname=_method_name,
                **TestWideData.solkwargs
            )

        # 4. Get solution
        solutionobj = TestWideData.solutionfixer.get_solution(
            methodname=_method_name, **TestWideData.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestWideData.solutionfixer.create_a_diff(
                to_check=data_to_test, solution=solutionobj
            )
            print(debug_diff)
        # 6. assert the equality
        assert test_expr

    def test_sync_wide_records(self, overwrite_solution=True):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = TestWideData.solutionfixer.get_solution(
            **TestWideData.solkwargs, methodname="test_import_wide_data"
        )
        # 2. apply a metobs manipulation
        dataset.sync_records(
            timestamp_shift_tolerance="5min2s", freq_shift_tolerance="2min"
        )

        # 3. overwrite solution?
        if overwrite_solution:
            TestWideData.solutionfixer.create_solution(
                solutiondata=dataset, **TestWideData.solkwargs, methodname=_method_name
            )
        # 4. Get solution
        solutionobj = TestWideData.solutionfixer.get_solution(
            **TestWideData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests
        test_expr = dataset == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestWideData.solutionfixer.create_a_diff(
                to_check=dataset, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr
        assert dataset.df.shape == (196, 2)


class TestWideSingleStationData:
    # to pass to the solutionfixer
    solkwargs = {
        "testfile": Path(__file__).name,
        "classname": "testwidesinglestationdata",
    }
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("single_station.csv")
    templatefile = datadir.joinpath("single_station_template.json")
    metadatfile = datadir.joinpath("single_station_metadata.csv")

    def test_import_wide_data(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        pass

        # 2. apply a metobs manipulation
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestWideSingleStationData.templatefile,
            input_metadata_file=TestWideSingleStationData.metadatfile,
            input_data_file=TestWideSingleStationData.datafile,
        )

        data_to_test = dataset

        # 3. overwrite solution?
        if overwrite_solution:
            TestWideSingleStationData.solutionfixer.create_solution(
                solutiondata=data_to_test,
                methodname=_method_name,
                **TestWideSingleStationData.solkwargs
            )

        # 4. Get solution
        solutionobj = TestWideSingleStationData.solutionfixer.get_solution(
            methodname=_method_name, **TestWideSingleStationData.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestWideSingleStationData.solutionfixer.create_a_diff(
                to_check=data_to_test, solution=solutionobj
            )
            print(debug_diff)
        # 6. assert the equality
        assert test_expr


if __name__ == "__main__":
    pytest.main([__file__])

    # widedatatester = TestWideData()
    # widedatatester.test_import_wide_data(False)
    # widedatatester.test_sync_wide_records(False)

    # singletest = TestWideSingleStationData()
    # singletest.test_import_wide_data(False)
