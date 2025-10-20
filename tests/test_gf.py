import pytest
import sys
import logging
from pathlib import Path
import copy
import pandas as pd
import numpy as np

import tempfile


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir


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
            max_gap_duration_to_fill="5h",
            n_leading_anchors=3,
            n_trailing_anchors=2,  # only 1 is used
            max_lead_to_gap_distance=pd.Timedelta("3h"),
            max_trail_to_gap_distance=None,
            overwrite_fill=False,
            method_kwargs={"order": 4},
        )

        # regular interpolation with overwrite_fill == false -> should only try to fill the failed filled gaps
        sta.interpolate_gaps(
            max_gap_duration_to_fill="5h",
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
        _method_name = "test_interpolation_on_dataset"

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
            max_gap_duration_to_fill=pd.Timedelta("5h"),
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
            max_gap_duration_to_fill=pd.Timedelta("5h"),
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

    def test_interpolating_with_station_without_that_obstype(self):
        # goal is to test if metobs is able to interpolate on a dataset,
        # where there is a station without the target obstype.

        df = pd.read_csv(metobs_toolkit.demo_datafile, sep=";")

        trgstation = "vlinder03"
        trg_column = "Temperatuur"
        # all to Nan
        df.loc[df["Vlinder"] == trgstation, trg_column] = np.nan

        # to csv
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            targetfile = tmpdir / "data_with_nans.csv"
            df.to_csv(targetfile, index=False, sep=";")

            dataset = metobs_toolkit.Dataset()
            dataset.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                input_data_file=targetfile,
                input_metadata_file=metobs_toolkit.demo_metadatafile,
            )

        dataset.repetitions_check(max_N_repetitions=8)
        dataset.convert_outliers_to_gaps()
        dataset.interpolate_gaps(target_obstype="temp", method="linear")

    def test_raw_modeldata_gapfill(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_raw_modeldata_gapfill"
        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # test raw gapfill on dataset
        dataset.fill_gaps_with_raw_modeldata(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("6h"),
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
        sta.fill_gaps_with_raw_modeldata(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=False,
        )

        assert_equality(
            sta, solutionobj_A.get_station("vlinder01")
        )  # station comparison

        # test the overwrite is true option
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        dataset.interpolate_gaps(target_obstype="temp", overwrite_fill=False)
        dataset.fill_gaps_with_raw_modeldata(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=True,
        )

        assert_equality(dataset, solutionobj_A)  # dataset comparison

        # test the plot
        dataset.make_plot()

    def test_chaining_gapfill_methods(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_chaining_gapfill_methods"
        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        dataset.interpolate_gaps(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("5h"),
            overwrite_fill=False,
        )

        # test raw gapfill on dataset
        dataset.fill_gaps_with_raw_modeldata(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("10h"),
            overwrite_fill=False,
        )

        dataset.make_plot(colorby="label")

        #  overwrite solution?
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestDataWithGaps.solkwargs,
                methodname=f"{_method_name}",
            )

        #  Get solution
        solutionobj_A = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=f"{_method_name}"
        )

        assert_equality(dataset, solutionobj_A)  # dataset comparison

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
            max_gap_duration_to_fill=pd.Timedelta("12h"),
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
            max_gap_duration_to_fill=pd.Timedelta("12h"),
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
            max_gap_duration_to_fill=pd.Timedelta("12h"),
            min_debias_sample_size=2,
            overwrite_fill=False,
        )

        assert_equality(
            sta, solutionobj_A.get_station("vlinder01")
        )  # station comparison

    def test_get_info_on_objects(self):
        #   get_startpoint data
        dataset_gf = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs,
            methodname="test_diurnal_debias_modeldata_gapfill_A",
        )

        # test on dataset with gapfilled data
        _ = dataset_gf.get_info(printout=False)
        _ = dataset_gf.get_station("vlinder04").get_info(printout=False)
        # test the get_info method on gap
        _ = (
            dataset_gf.get_station("vlinder04")
            .get_sensor("temp")
            .gaps[0]
            .get_info(printout=False)
        )

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

    def test_partially_filled_gaps(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = "test_partially_filled_gaps"

        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # create outliers
        dataset.repetitions_check(max_N_repetitions=8)
        dataset.gross_value_check(upper_threshold=16.1)

        # test diurnal debias gapfill on dataset
        dataset.fill_gaps_with_weighted_diurnal_debiased_modeldata(
            target_obstype="temp",
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_lead_debias_sample_size=1,  # This condition is not always met
            min_trail_debias_sample_size=1,  # This condition is not always met
            overwrite_fill=False,
            max_gap_duration_to_fill=pd.Timedelta("48h"),
        )

        assert (
            "partially successful gapfill" in dataset.gap_overview_df()["label"].values
        )

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

        assert_equality(dataset, solutionobj)  # dataset comparison
        dataset.stations[0].make_plot(colorby="label")

        # Test chaining after partially filled gaps:
        # Choice: when force=False, a partially filled gap will BE filled again when chaining !!

        dataset.stations[0].interpolate_gaps(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("25h"),
            overwrite_fill=False,
        )
        dataset.stations[0].make_plot(colorby="label")

    def test_chaining_on_partially_filled_gaps(self):

        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_partially_filled_gaps"
        )
        assert (
            "partially successful gapfill" in dataset.gap_overview_df()["label"].values
        )

        # Test chaining after partially filled gaps:
        # Choice: when force=False, a partially filled gap will BE filled again when chaining !!

        dataset.stations[0].interpolate_gaps(
            target_obstype="temp",
            max_gap_duration_to_fill=pd.Timedelta("25h"),
            overwrite_fill=False,
        )

        # dataset.stations[0].make_plot(colorby="label")

        assert not (
            "partially successful gapfill" in dataset.gap_overview_df()["label"].values
        )

    def test_add_modeldata_to_station(self):
        #   get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )
        sta = dataset.get_station("vlinder02")

        # create a fake new modeltimesries
        orig_modeltimeseries = sta.get_modeltimeseries("temp")
        fake_data = (orig_modeltimeseries.series + 16.2) / 3.1

        fake_modeltimeseries = metobs_toolkit.ModelTimeSeries(
            site=sta.site,
            datarecords=fake_data.to_numpy(),
            timestamps=fake_data.index.to_numpy(),
            modelobstype=orig_modeltimeseries.modelobstype,
            modelname=orig_modeltimeseries.modelname,
            modelvariable=orig_modeltimeseries.modelvariable,
        )
        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsDataAlreadyPresent,
        )

        with pytest.raises(MetObsDataAlreadyPresent):
            sta.add_to_modeldata(fake_modeltimeseries)

        # change modelname, now modeldata should be accepted
        fake_modeltimeseries.modelname = "fake_model"
        sta.add_to_modeldata(fake_modeltimeseries)

        fake_modeldata2 = copy.deepcopy(fake_modeltimeseries)
        fake_modeldata2.modelvariable = "fake name"
        sta.add_to_modeldata(fake_modeldata2)

        assert len(sta.modeldata) == 3

    def test_gap_status_df(self, overwrite_solution=False):
        """Test gap_overview_df methods on Dataset, Station, and SensorData classes."""
        # 0. Get info of the current check
        _method_name = "test_gap_status_df"

        # 1. Get starting data without gaps/fills
        dataset_original = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # Test 1: Dataset without gaps (original data)
        gap_status_dataset_no_gaps = dataset_original.gap_overview_df()

        # Test 2: Station without gaps (original data)
        station_original = dataset_original.get_station("vlinder01")
        gap_status_station_no_gaps = station_original.gap_overview_df()

        # Test 3: SensorData without gaps (original data)
        sensordata_original = station_original.get_sensor("temp")
        gap_status_sensordata_no_gaps = sensordata_original.gap_overview_df()

        # 2. Get data with gaps (after creating gaps)
        dataset_with_gaps = copy.deepcopy(dataset_original)
        dataset_with_gaps.convert_outliers_to_gaps(all_observations=True)

        # Test 4: Dataset with gaps (no gap filling)
        gap_status_dataset_with_gaps = dataset_with_gaps.gap_overview_df()

        # Test 5: Station with gaps (no gap filling)
        station_with_gaps = dataset_with_gaps.get_station("vlinder01")
        gap_status_station_with_gaps = station_with_gaps.gap_overview_df()

        # Test 6: SensorData with gaps (no gap filling)
        sensordata_with_gaps = station_with_gaps.get_sensor("temp")
        gap_status_sensordata_with_gaps = sensordata_with_gaps.gap_overview_df()

        # 3. Get data with gap filling applied
        dataset_filled = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_interpolation_on_dataset_A"
        )

        # Test 7: Dataset with gap filling
        gap_status_dataset_filled = dataset_filled.gap_overview_df()

        # Test 8: Station with gap filling
        station_filled = dataset_filled.get_station("vlinder01")
        gap_status_station_filled = station_filled.gap_overview_df()

        # Test 9: SensorData with gap filling
        sensordata_filled = station_filled.get_sensor("temp")
        gap_status_sensordata_filled = sensordata_filled.gap_overview_df()

        # 4. Combine all results for solution comparison
        test_results = {
            "dataset_no_gaps": gap_status_dataset_no_gaps,
            "station_no_gaps": gap_status_station_no_gaps,
            "sensordata_no_gaps": gap_status_sensordata_no_gaps,
            "dataset_with_gaps": gap_status_dataset_with_gaps,
            "station_with_gaps": gap_status_station_with_gaps,
            "sensordata_with_gaps": gap_status_sensordata_with_gaps,
            "dataset_filled": gap_status_dataset_filled,
            "station_filled": gap_status_station_filled,
            "sensordata_filled": gap_status_sensordata_filled,
        }

        # 5. Overwrite solution if requested
        if overwrite_solution:
            TestDataWithGaps.solutionfixer.create_solution(
                solutiondata=test_results,
                **TestDataWithGaps.solkwargs,
                methodname=_method_name,
            )

        # 6. Get solution and test equality
        solutionobj = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname=_method_name
        )

        # 7. Assert equality for all test cases
        for key in test_results:
            assert_equality(test_results[key], solutionobj[key])

    def test_min_max_value_clipping(self):
        """Test that min_value and max_value parameters work for all model-based gap filling methods."""
        # Get test data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            **TestDataWithGaps.solkwargs, methodname="test_import_data"
        )

        # Test parameters
        target_obstype = "temp"
        min_threshold = 11.0
        max_threshold = 16.0

        # Test 1: Raw model data gap fill with min/max clipping
        dataset_raw = copy.deepcopy(dataset)
        dataset_raw.fill_gaps_with_raw_modeldata(
            target_obstype=target_obstype,
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=False,
            min_value=min_threshold,
            max_value=max_threshold,
        )

        # Verify that filled values are within bounds
        assert (
            dataset_raw.gapsdf.xs("temp", level="obstype")["value"].min()
            == min_threshold
        )
        assert (
            dataset_raw.gapsdf.xs("temp", level="obstype")["value"].max()
            == max_threshold
        )

        # Test 2: Debiased model data gap fill with min/max clipping
        dataset_debiased = copy.deepcopy(dataset)
        dataset_debiased.fill_gaps_with_debiased_modeldata(
            target_obstype=target_obstype,
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_leading_records_total=10,
            min_trailing_records_total=10,
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=False,
            min_value=None,
            max_value=max_threshold,
        )

        # Verify that filled values are within bounds
        assert (
            dataset_debiased.gapsdf.xs("temp", level="obstype")["value"].max()
            == max_threshold
        )

        # Test 3: Diurnal debiased model data gap fill with min/max clipping
        dataset_diurnal = copy.deepcopy(dataset)
        dataset_diurnal.fill_gaps_with_diurnal_debiased_modeldata(
            target_obstype=target_obstype,
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_debias_sample_size=2,
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=False,
            min_value=min_threshold,
            max_value=None,
        )
        assert (
            dataset_diurnal.gapsdf.xs("temp", level="obstype")["value"].min()
            == min_threshold
        )

        # Test 4: Weighted diurnal debiased model data gap fill with min/max clipping
        dataset_weighted = copy.deepcopy(dataset)
        dataset_weighted.fill_gaps_with_weighted_diurnal_debiased_modeldata(
            target_obstype=target_obstype,
            leading_period_duration=pd.Timedelta("24h"),
            trailing_period_duration=pd.Timedelta("24h"),
            min_lead_debias_sample_size=1,
            min_trail_debias_sample_size=1,
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=False,
            min_value=min_threshold,
            max_value=max_threshold,
        )

        assert (
            dataset_weighted.gapsdf.xs("temp", level="obstype")["value"].min()
            == min_threshold
        )
        assert (
            dataset_weighted.gapsdf.xs("temp", level="obstype")["value"].max()
            == max_threshold
        )

        # Test 5: Test on individual Station (not Dataset)
        station = dataset.get_station("vlinder01")
        station.fill_gaps_with_raw_modeldata(
            target_obstype=target_obstype,
            max_gap_duration_to_fill=pd.Timedelta("6h"),
            overwrite_fill=False,
            min_value=12.90,
            max_value=max_threshold,
        )

        assert station.gapsdf.xs("temp", level="obstype")["value"].min() == 12.90
        assert (
            station.gapsdf.xs("temp", level="obstype")["value"].max() == max_threshold
        )

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
    # tester.test_import_data(overwrite_solution=False)
    # tester.test_interpolation_on_station(overwrite_solution=False)
    # tester.test_interpolation_on_dataset(overwrite_solution=False)
    # tester.test_raw_modeldata_gapfill(overwrite_solution=False)
    # tester.test_partially_filled_gaps(overwrite_solution=False)
    # tester.test_chaining_gapfill_methods(overwrite_solution=False)
    # tester.test_debias_modeldata_gapfill(overwrite_solution=False)
    # tester.test_diurnal_debias_modeldata_gapfill(overwrite_solution=False)
    # tester.test_weighted_diurnal_debias_modeldata_gapfill(overwrite_solution=False)
    # tester.test_gap_status_df(overwrite_solution=False)
    # tester.test_min_max_value_clipping()
