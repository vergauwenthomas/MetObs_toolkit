import pytest
import sys
import logging
from pathlib import Path
import copy
import pandas as pd


libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir

import pytest


class TestDataWithGaps:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdatawithgaps"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_import_data(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        if overwrite_solution:  # GEE INTERACTION !!

            dataset = metobs_toolkit.Dataset()
            dataset.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                input_metadata_file=metobs_toolkit.demo_metadatafile,
                input_data_file=datadir.joinpath("testdata_with_gaps.csv"),
            )
            # To other resolution!!
            dataset.resample(target_freq="15min")

            # extracting modeldata
            era5_manager = metobs_toolkit.default_GEE_datasets["ERA5-land"]
            era5_data = dataset.get_gee_timeseries_data(
                geedynamicdatasetmanager=era5_manager,
                startdt_utc=None,  # raises error in metadata-only case
                enddt_utc=None,
                target_obstypes=["temp"],
                get_all_bands=False,
                drive_filename=None,
                # drive_folder="gee_timeseries_data",
                force_direct_transfer=True,
                force_to_drive=False,
            )

            dataset.save_dataset_to_pkl(
                target_folder=(
                    solutionsdir.joinpath("test_gf_solutions").joinpath(
                        "testdatawithgaps"
                    )
                ),
                filename="test_import_data.pkl",
                overwrite=True,
            )
        else:
            dataset = metobs_toolkit.import_dataset_from_pkl(
                solutionsdir.joinpath("test_gf_solutions")
                .joinpath("testdatawithgaps")
                .joinpath("test_import_data.pkl")
            )
        # 3. overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestDataWithGaps.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDataWithGaps.solutionfixer.get_solution(
            methodname=_method_name, **TestDataWithGaps.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    def test_interpolation_on_station(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_interpolation_on_station"
        #  1. get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        sta = dataset.get_station("vlinder01")

        # test interpolation using higher order cubic spline
        sta.interpolate_gaps(
            target_obstype="temp",
            method="cubicspline",
            max_consec_fill=20,  # is 5 hours at 15min
            n_leading_anchors=3,
            n_trailing_anchors=2,  # only 1 is used
            max_lead_to_gap_distance=pd.Timedelta("3h"),
            max_trail_to_gap_distance=None,
            overwrite_fill=False,
            method_kwargs={"order": 4},
        )

        # regular interpolation iwht overwrite_fill == false -> should only try to fill the failed filled gaps
        sta.interpolate_gaps(
            target_obstype="temp",
            overwrite_fill=False,
        )

        #  3. overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=sta,
                **TestDataWithGaps.solkwargs,
                methodname=_method_name,
            )

        # 4. Get solution
        solutionobj = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests on dataset level
        assert_equality(sta, solutionobj)

        # Test plotting
        _statsdf = sta.make_plot()

    def test_interpolation_on_dataset(self, overwrite_solution=False):

        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # ------------------------------------------
        #   A: Test higher order interpolation on dataset scale
        # ------------------------------------------
        # test interpolation using higher order cubic spline
        dataset.interpolate_gaps(
            target_obstype="temp",
            method="cubicspline",
            max_consec_fill=20,  # is 5 hours at 15min
            n_leading_anchors=3,
            n_trailing_anchors=2,  # only 1 is used
            max_lead_to_gap_distance=pd.Timedelta("3h"),
            max_trail_to_gap_distance=None,
            overwrite_fill=False,
            method_kwargs={"order": 3},
        )

        #  overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}_A",
            )

        # Get solution
        solutionobj_A = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}_A"
        )

        # Construct the equlity tests on dataset level
        assert_equality(dataset, solutionobj_A)

        # ------------------------------------------
        #    B: test overwrite_fill argument
        # ------------------------------------------

        dataset.interpolate_gaps(
            target_obstype="temp",
            method="spline",
            max_consec_fill=20,  # is 5 hours at 15min
            n_leading_anchors=3,
            n_trailing_anchors=2,  # only 1 is used
            max_lead_to_gap_distance=pd.Timedelta("3h"),
            max_trail_to_gap_distance=None,
            overwrite_fill=False,  # This should not do anything, since gaps are already filled
            method_kwargs={"order": 2},
        )
        assert_equality(dataset, solutionobj_A)  # dataset comparison

        # regular interpolation iwht overwrite_fill == True -> should overwrite the data!
        dataset.interpolate_gaps(
            target_obstype="temp",
            overwrite_fill=True,
        )

        assert dataset != solutionobj_A

        #  3. overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}_B",
            )

        # 4. Get solution
        solutionobj_B = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}_B"
        )

        assert_equality(dataset, solutionobj_B)  # dataset comparison

    def test_raw_modeldata_gapfill(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_raw_modeldata_gapfill"
        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # test raw gapfill on dataset
        dataset.fill_gaps_with_raw_modeldata(
            target_obstype="temp", overwrite_fill=False
        )

        #  overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}_A",
            )

        #  Get solution
        solutionobj_A = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}_A"
        )

        assert_equality(dataset, solutionobj_A)  # dataset comparison

        from metobs_toolkit.backend_collection.errorclasses import MetObsModelDataError

        with pytest.raises(MetObsModelDataError):
            dataset.fill_gaps_with_raw_modeldata(
                target_obstype="humidity", overwrite_fill=True
            )

        # test on station and dataset

        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        sta = dataset.get_station("vlinder01")
        sta.fill_gaps_with_raw_modeldata(target_obstype="temp", overwrite_fill=False)

        assert_equality(
            sta, solutionobj_A.get_station("vlinder01")
        )  # station comparison

        # test the overwrite is true option
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        dataset.interpolate_gaps(target_obstype="temp", overwrite_fill=False)
        dataset.fill_gaps_with_raw_modeldata(target_obstype="temp", overwrite_fill=True)

        assert_equality(dataset, solutionobj_A)  # dataset comparison

        # test the plot
        dataset.make_plot()

    def test_debias_modeldata_gapfill(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_debias_modeldata_gapfill"

        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # test debias gapfill on dataset
        dataset.fill_gaps_with_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("6h"),
            min_leading_records_total=5,
            trailing_period_duration=pd.Timedelta("24h"),
            min_trailing_records_total=8,
            overwrite_fill=False,
        )

        #  overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}",
            )

        #  Get solution
        solutionobj = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}"
        )

        # test equality
        assert_equality(to_check=dataset, solution=solutionobj)

        # test on station and dataset
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        sta = dataset.get_station("vlinder01")
        sta.fill_gaps_with_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("6h"),
            min_leading_records_total=5,
            trailing_period_duration=pd.Timedelta("24h"),
            min_trailing_records_total=8,
            overwrite_fill=False,
        )

        assert_equality(sta, solutionobj.get_station("vlinder01"))

    def test_diurnal_debias_modeldata_gapfill(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # test diurnal debias gapfill on dataset
        dataset.fill_gaps_with_diurnal_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_debias_sample_size=2,
            overwrite_fill=False,
        )

        #  overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}_A",
            )

        #  Get solution
        solutionobj_A = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}_A"
        )
        assert_equality(dataset, solutionobj_A)  # dataset comparison

        # test on station and dataset
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        sta = dataset.get_station("vlinder01")
        sta.fill_gaps_with_diurnal_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_debias_sample_size=2,
            overwrite_fill=False,
        )

        assert_equality(
            sta, solutionobj_A.get_station("vlinder01")
        )  # station comparison

    def test_weighted_diurnal_debias_modeldata_gapfill(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # test diurnal debias gapfill on dataset
        dataset.fill_gaps_with_weighted_diurnal_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_lead_debias_sample_size=1,
            min_trail_debias_sample_size=0,  # just testing
            overwrite_fill=False,
        )

        #  overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}_A",
            )

        #  Get solution
        solutionobj_A = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}_A"
        )
        assert_equality(dataset, solutionobj_A)  # dataset comparison

        # test on station and dataset
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        sta = dataset.get_station("vlinder01")
        sta.fill_gaps_with_weighted_diurnal_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_lead_debias_sample_size=1,
            min_trail_debias_sample_size=0,  # just testing
            overwrite_fill=False,
        )

        assert_equality(
            sta, solutionobj_A.get_station("vlinder01")
        )  # station comparison

    # ------------------------------------------
    #    Plotting tests are present in the test_plotting.py
    # ------------------------------------------


if __name__ == "__main__":
    print(
        "To Overwrite the solutions, run: \n pytest test_plotting.py  --mpl --mpl-generate-path=baseline"
    )

    print(
        "To checkout the differences, run: \n pytest test_plotting.py --mpl --mpl-generate-summary=html "
    )

    tester = TestDataWithGaps()
    tester.test_import_data(overwrite_solution=False)
    tester.test_interpolation_on_station(overwrite_solution=False)
    tester.test_interpolation_on_dataset(overwrite_solution=False)
    tester.test_raw_modeldata_gapfill(overwrite_solution=False)
    tester.test_debias_modeldata_gapfill(overwrite_solution=False)
    tester.test_diurnal_debias_modeldata_gapfill(overwrite_solution=False)
    tester.test_weighted_diurnal_debias_modeldata_gapfill(overwrite_solution=False)
