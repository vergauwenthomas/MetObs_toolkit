import sys
from pathlib import Path

import copy


# import metobs_toolkit
import pandas as pd


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir

import pytest


class TestBreakingDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testbreakingdata"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("testdata_breaking.csv")
    templatefile = datadir.joinpath("template_breaking.json")

    def test_import_data(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        pass

        # 2. apply a metobs manipulation
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestBreakingDataset.templatefile,
            #   input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=TestBreakingDataset.datafile,
        )

        # 3. overwrite solution?
        if overwrite_solution:
            TestBreakingDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestBreakingDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestBreakingDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestBreakingDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    def test_apply_qc(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        dataset = TestBreakingDataset.solutionfixer.get_solution(
            **TestBreakingDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        # apply QC
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=-15.0,
            upper_threshold=29.0,
            use_mp=False,
        )
        # fake check on humidity to see if this does not affect temp records
        dataset.gross_value_check(
            obstype="humidity",
            lower_threshold=5.0,
            upper_threshold=10.0,
            use_mp=False,
        )

        dataset.persistence_check(
            obstype="temp",
            timewindow="1h",
            min_records_per_window=3,
            use_mp=False,
        )

        dataset.repetitions_check(obstype="temp", max_N_repetitions=5, use_mp=False)
        dataset.step_check(
            obstype="temp",
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
            use_mp=False,
        )
        dataset.window_variation_check(
            obstype="temp",
            timewindow="1h",
            min_records_per_window=3,
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
            # use_mp=True,
            use_mp=False,
        )

        #  3. overwrite solution?
        if overwrite_solution:
            TestBreakingDataset.solutionfixer.create_solution(
                solutiondata=dataset,
                **TestBreakingDataset.solkwargs,
                methodname=_method_name,
            )
        # 4. Get solution
        solutionobj = TestBreakingDataset.solutionfixer.get_solution(
            **TestBreakingDataset.solkwargs, methodname=_method_name
        )

        # reading manual labels
        man_df = pd.read_csv(TestBreakingDataset.datafile)
        # create datetime coumn
        man_df["datetime"] = man_df["date"] + " " + man_df["time"]
        man_df["datetime"] = pd.to_datetime(man_df["datetime"])
        man_df["datetime"] = man_df["datetime"].dt.tz_localize("Europe/Berlin")
        # drop the rows with an irregular timestamp
        man_df = man_df[man_df["irregular stamps"] != "irregular timestamp"]
        # create double index
        man_df["name"] = man_df["station"]
        man_df = man_df.set_index(["name", "datetime"]).sort_index()
        # format the label column
        man_df = man_df.rename(columns={"qc_flags": "label_manual"})
        # check if labels are the same
        compare = man_df[["label_manual"]]

        checkdf = (
            dataset.df.xs("temp", level="obstype")
            .reset_index()
            .set_index(["name", "datetime"])
        )

        compare = man_df.merge(checkdf, how="left", left_index=True, right_index=True)

        compare["correct_labeled"] = compare["label_manual"].eq(compare["label"])
        diff = compare.loc[compare["correct_labeled"] == False, :]
        assert diff.empty, f"Incorrect labels compared to the manual labels: \n {diff}"

        # 5. Construct the equlity tests on dataset level

        assert_equality(dataset, solutionobj)  # dataset comparison

    def test_qc_statistics(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name
        #  1. get_startpoint data
        dataset = TestBreakingDataset.solutionfixer.get_solution(
            **TestBreakingDataset.solkwargs, methodname="test_apply_qc"
        )

        #  2. apply a metobs manipulation
        # apply QC
        statsdf = dataset.get_qc_stats(obstype="temp", make_plot=False)
        #  3. overwrite solution?
        if overwrite_solution:
            TestBreakingDataset.solutionfixer.create_solution(
                solutiondata=statsdf,
                **TestBreakingDataset.solkwargs,
                methodname=_method_name,
            )
        # 4. Get solution
        solutionobj = TestBreakingDataset.solutionfixer.get_solution(
            **TestBreakingDataset.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests on dataset level
        assert_equality(statsdf, solutionobj)

        # Test plotting
        _statsdf = dataset.get_qc_stats(obstype="temp", make_plot=True)

    def test_get_info(self):
        #  1. get_startpoint data
        dataset = TestBreakingDataset.solutionfixer.get_solution(
            **TestBreakingDataset.solkwargs, methodname="test_apply_qc"
        )
        # call get info on dataset, station and sensor level
        _ = dataset.get_info(printout=True)
        _ = dataset.get_station("Fictional").get_info(printout=True)
        _ = dataset.get_station("Fictional").get_sensor("temp").get_info(printout=True)


class TestDemoDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_import_data(self, overwrite_solution=False):
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
        # To hourly !!
        dataset.resample(target_freq="1h")

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

    def test_qc_when_some_stations_missing_obs(self):
        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )
        obstype = "temp"

        # Drop the temp of 2 random stations
        orig_count = len(dataset.stations)
        del dataset.stations[6].obsdata["temp"]
        del dataset.stations[13].obsdata["temp"]

        # See if qc pipeline still runs

        # Run all QC checks with default settings
        dataset.gross_value_check(obstype=obstype)
        dataset.persistence_check(obstype=obstype)
        dataset.repetitions_check(obstype=obstype)
        dataset.step_check(obstype=obstype)
        dataset.window_variation_check(obstype=obstype)
        dataset.buddy_check(obstype=obstype)

        assert orig_count == len(dataset.stations)

    def test_buddy_check(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        # Test if errors ar raised

        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsMetadataNotFound,
        )

        with pytest.raises(MetObsMetadataNotFound):
            # Should raise error because no altitude info is available and lapsrate is not none
            dataset.buddy_check(
                obstype="temp",
                spatial_buddy_radius=17000,
                min_sample_size=3,
                max_alt_diff=150,
                min_std=1.0,
                spatial_z_threshold=2.4,
                N_iter=1,
                instantaneous_tolerance=pd.Timedelta("4min"),
                lapserate=-0.0065,  # -0.0065
                use_mp=False,
            )

        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsMetadataNotFound,
        )

        with pytest.raises(MetObsMetadataNotFound):
            # Should raise error because no altitude info is available and max_alt_diffis not none
            dataset.buddy_check(
                obstype="temp",
                spatial_buddy_radius=17000,
                min_sample_size=3,
                max_alt_diff=150,
                min_std=1.0,
                spatial_z_threshold=2.4,
                N_iter=1,
                instantaneous_tolerance=pd.Timedelta("4min"),
                lapserate=None,  # -0.0065
                use_mp=False,
            )

        # Test that buddy check runs, with settings that does not create outliers
        dataset.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            spatial_z_threshold=5.9,  # this does noet create outliers
            N_iter=1,
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        assert dataset.outliersdf.empty

        # Now create outliers with the buddy check

        dataset1 = copy.deepcopy(dataset)  # used to test 1 iteration
        dataset2 = copy.deepcopy(dataset)  # use to test 2 iterations

        # test one iteration
        dataset1.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,  # one iteration test
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        # outliersdf_1_iter = dataset1.outliersdf

        # test two iteration
        dataset2.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            spatial_z_threshold=2.1,
            N_iter=2,  # one iteration test
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        # overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset1,
                **TestDemoDataset.solkwargs,
                methodname=_method_name + "_1_iter",
            )
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset2,
                **TestDemoDataset.solkwargs,
                methodname=_method_name + "_2_iter",
            )

        # 4. Get solution
        solutionobj_1iter = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname=_method_name + "_1_iter"
        )
        solutionobj_2iter = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname=_method_name + "_2_iter"
        )

        # validate expression
        assert_equality(dataset1, solutionobj_1iter)  # dataset comparison

        assert_equality(dataset2, solutionobj_2iter)  # dataset comparison

    def test_buddy_check_with_big_radius(self):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # Tricky thing is that with big radii, a station can appear in multiple
        # buddy groups, which can lead to edge cases. Here we test that the code
        # runs without errors

        dataset.buddy_check(
            obstype="temp",
            spatial_buddy_radius=50000,  # Large radius
            min_sample_size=2,
            spatial_z_threshold=1.8,  # Lower threshold
            N_iter=1,  # Multiple iterations increases chance of edge cases
            instantaneous_tolerance=pd.Timedelta("5min"),
            use_mp=False,  # Deterministic behavior
        )

    def test_buddy_check_with_safety_nets(self, overwrite_solution=False):
        """Test the generalized buddy_check_with_safety_nets method."""
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        # 1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # Ensure LCZ data is present for all stations
        if not all(sta.site.flag_has_LCZ() for sta in dataset.stations):
            dataset.get_LCZ()

        # Test 1: Using safety_net_configs with LCZ should match the LCZ safety net method
        dataset1 = copy.deepcopy(dataset)
        dataset1.buddy_check_with_safetynets(
            obstype="temp",
            spatial_buddy_radius=25000,
            safety_net_configs=[
                {
                    "category": "LCZ",
                    "buddy_radius": 100000,
                    "z_threshold": 1.4,
                    "min_sample_size": 2,
                }
            ],
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,
            use_mp=False,
        )
        assert (
            dataset1.outliersdf.shape[0] == 74
        ), f"Expected 74 outliers, got {dataset1.outliersdf.shape[0]}"

        # overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=dataset1,
                **TestDemoDataset.solkwargs,
                methodname=_method_name,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname=_method_name
        )

        # validate expression
        assert_equality(dataset1, solutionobj)

    def test_buddy_check_with_safety_nets_missing_min_sample_size(self):
        """Test that an error is raised when min_sample_size is missing from safety_net_configs."""
        # Get dataset
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # Ensure LCZ data is present
        if not all(sta.site.flag_has_LCZ() for sta in dataset.stations):
            dataset.get_LCZ()

        # Test that missing min_sample_size raises an error
        with pytest.raises(KeyError):
            dataset.buddy_check_with_safetynets(
                obstype="temp",
                spatial_buddy_radius=25000,
                safety_net_configs=[
                    {
                        "category": "LCZ",
                        "buddy_radius": 100000,
                        "z_threshold": 1.4,
                        # "min_sample_size" is intentionally missing
                    }
                ],
                min_sample_size=3,
                spatial_z_threshold=2.1,
                use_mp=False,
            )


class TestWhiteRecords:
    """Test white_records functionality for all QC checks on both Dataset and Station level."""

    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testwhiterecords"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_whiterecords_reprs(self):
        """Test the WhiteSet and SensorWhiteSet classes."""
        # Create a WhiteSet with mixed levels
        white_records = pd.date_range(
            start="2023-01-01 00:00",
            periods=10,
            freq="h",
            tz="UTC",
        )
        multi_index = pd.MultiIndex.from_product(
            [
                ["vlinder05", "vlinder06", "vlinder07"],
                white_records,
                ["temp", "humidity"],
            ],
            names=["name", "datetime", "obstype"],
        )
        whiteset = metobs_toolkit.WhiteSet(multi_index)

        # Check repr and str
        repr_str = repr(whiteset)
        str_str = str(whiteset)

        assert "WhiteSet" in repr_str
        assert "n_records=" in repr_str
        assert "levels=" in repr_str
        assert str_str == repr_str

        # Create a SensorWhiteSet
        sensor_whiteset = whiteset.create_sensorwhitelist(
            stationname="vlinder06", obstype="temp"
        )
        # Check repr and str
        repr_sensor = repr(sensor_whiteset)
        str_sensor = str(sensor_whiteset)

        assert "SensorWhiteSet" in repr_sensor
        assert "n_timestamps=" in repr_sensor
        assert str_sensor == repr_sensor

    def test_import_data(self, overwrite_solution=False):
        """Import demo dataset for white_records testing."""
        _method_name = sys._getframe().f_code.co_name

        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        # Resample to hourly for consistent testing
        dataset.resample(target_freq="1h")

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj)

    def test_white_records_input_combinations(self, overwrite_solution=False):
        """Test white_records with gross_value_check on Dataset level."""
        _method_name = sys._getframe().f_code.co_name

        dataset = TestWhiteRecords.solutionfixer.get_solution(
            **TestWhiteRecords.solkwargs, methodname="test_import_data"
        )

        # Get some timestamps that would be flagged as outliers
        # First run without white_records to identify outliers
        test_dataset = copy.deepcopy(dataset)
        test_dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            use_mp=False,
        )

        # Get outliers to use as white_records
        outliers = test_dataset.outliersdf
        if outliers.empty:
            pytest.skip("No outliers found for white_records testing")

        # Test 1a: Index with only datetimes
        dataset1a = copy.deepcopy(dataset)
        white_dt_only = pd.Index(
            outliers.reset_index()["datetime"].head(20), name="datetime"
        )
        dataset1a.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
            use_mp=False,
        )
        outliers1a = dataset1a.outliersdf

        # Test 1b: Index with only name
        dataset1b = copy.deepcopy(dataset)
        white_name_only = pd.Index(
            data=["vlinder05", "vlinder05", "vlinder06", "fake"], name="name"
        )
        dataset1b.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_only),
            use_mp=False,
        )

        outliers1b = dataset1b.outliersdf

        # test on station object
        copy.deepcopy(dataset).get_station("vlinder05").gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_only),
        )

        # Test 2: MultiIndex with name and datetime
        dataset2 = copy.deepcopy(dataset)
        white_name_dt = (
            outliers.head(20)
            .reset_index()[["name", "datetime"]]
            .set_index(["name", "datetime"])
            .index
        )
        dataset2.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
            use_mp=False,
        )
        outliers2 = dataset2.outliersdf

        # Test 3: MultiIndex with obstype, name, and datetime
        dataset3 = copy.deepcopy(dataset)
        white_full = outliers.head(25).index
        dataset3.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_full),
            use_mp=False,
        )
        outliers3 = dataset3.outliersdf

        # Quick Verify that white-listed records are not in the outliers
        for white_record in white_dt_only:
            assert not any(
                outliers1a.reset_index()["datetime"] == white_record
            ), f"White-listed record {white_record} found in outliers (datetime only)"

        # Store results
        results = {
            "outliers_no_white": outliers,
            "outliers_dt_only": outliers1a,
            "outliers_name_only": outliers1b,
            "outliers_name_dt": outliers2,
            "outliers_full": outliers3,
        }

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solutiondata=results,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        for key in results:
            # Drop 'details' column if present for comparison
            df_to_compare = results[key].drop(columns=["details"], errors="ignore")
            sol_to_compare = solutionobj[key].drop(columns=["details"], errors="ignore")
            assert_equality(df_to_compare, sol_to_compare)

    def test_white_records_buddy_check_dataset(self):
        """Test white_records with buddy_check on Dataset level."""

        dataset = TestWhiteRecords.solutionfixer.get_solution(
            **TestWhiteRecords.solkwargs, methodname="test_import_data"
        )

        # First run without white_records
        # test_dataset = copy.deepcopy(dataset)
        # test_dataset.buddy_check(
        #     obstype="temp",
        #     spatial_buddy_radius=25000,
        #     min_sample_size=3,
        #     spatial_z_threshold=1.8,
        #     N_iter=2,
        #     use_mp=False,
        # )

        # outliers = test_dataset.outliersdf
        # white_dt_only = pd.Index(
        #     outliers.reset_index()["datetime"].sample(n=33, random_state=42),
        #     name="datetime",
        # )
        white_dt_only = pd.DatetimeIndex(['2022-09-11 16:00:00+00:00', '2022-09-09 01:00:00+00:00',
               '2022-09-01 15:00:00+00:00', '2022-09-05 03:00:00+00:00',
               '2022-09-12 16:00:00+00:00', '2022-09-15 04:00:00+00:00',
               '2022-09-04 18:00:00+00:00', '2022-09-04 05:00:00+00:00',
               '2022-09-01 19:00:00+00:00', '2022-09-14 13:00:00+00:00',
               '2022-09-12 18:00:00+00:00', '2022-09-01 14:00:00+00:00',
               '2022-09-14 05:00:00+00:00', '2022-09-04 10:00:00+00:00',
               '2022-09-09 03:00:00+00:00', '2022-09-09 02:00:00+00:00',
               '2022-09-12 04:00:00+00:00', '2022-09-01 05:00:00+00:00',
               '2022-09-05 02:00:00+00:00', '2022-09-15 06:00:00+00:00',
               '2022-09-09 00:00:00+00:00', '2022-09-06 06:00:00+00:00',
               '2022-09-11 16:00:00+00:00', '2022-09-03 09:00:00+00:00',
               '2022-09-05 07:00:00+00:00', '2022-09-04 17:00:00+00:00',
               '2022-09-11 00:00:00+00:00', '2022-09-14 17:00:00+00:00',
               '2022-09-02 01:00:00+00:00', '2022-09-04 18:00:00+00:00',
               '2022-09-05 05:00:00+00:00', '2022-09-09 07:00:00+00:00',
               '2022-09-15 02:00:00+00:00',
               '2023-09-15 02:00:00+00:00'],  #not in range
              dtype='datetime64[ns, UTC]', name='datetime', freq=None)
        
        # Test with different structures
        dataset1 = copy.deepcopy(dataset)
        dataset1.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            spatial_z_threshold=2.1,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
            use_mp=False,
        )
        outlier_timestamps = dataset1.outliersdf.index.get_level_values('datetime').unique()
        intersect = outlier_timestamps.intersection(white_dt_only)
        #Check if the white dt only timestamps are not in the outliers
        assert intersect.empty , "outlier timestamps found in white dt only"
        assert not outlier_timestamps.empty , "not outliers found"

        dataset2 = copy.deepcopy(dataset)
        
        # white_name_dt = (
        #     outliers.sample(n=18, random_state=42)
        #     .reset_index()[["name", "datetime"]]
        #     .set_index(["name", "datetime"])
        #     .index
        # )
        
        white_name_dt = pd.MultiIndex.from_arrays(
            [['vlinder27', 'vlinder05', 'vlinder27', 'vlinder07', 'vlinder05',
       'vlinder05', 'vlinder09', 'vlinder06', 'vlinder06', 'vlinder05',
       'vlinder05', 'vlinder27', 'vlinder05', 'vlinder05', 'vlinder05',
       'vlinder05', 'vlinder06', 'vlinder06', 'dummy_station', 'vlinder05',],
            pd.DatetimeIndex(['2022-09-11 16:00:00+00:00', '2022-09-09 01:00:00+00:00',
               '2022-09-01 15:00:00+00:00', '2022-09-05 03:00:00+00:00',
               '2022-09-12 16:00:00+00:00', '2022-09-15 04:00:00+00:00',
               '2022-09-04 18:00:00+00:00', '2022-09-04 05:00:00+00:00',
               '2022-09-01 19:00:00+00:00', '2022-09-14 13:00:00+00:00',
               '2022-09-12 18:00:00+00:00', '2022-09-01 14:00:00+00:00',
               '2022-09-14 05:00:00+00:00', '2022-09-04 10:00:00+00:00',
               '2022-09-09 03:00:00+00:00', '2022-09-09 02:00:00+00:00',
               '2022-09-12 04:00:00+00:00', '2022-09-01 05:00:00+00:00',
               '2022-09-12 04:00:00+00:00', '2023-09-01 05:00:00+00:00',], #fake records
              dtype='datetime64[ns, UTC]', name='datetime', freq=None)],
            names=('name', 'datetime'))
        

        dataset2.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            spatial_z_threshold=2.1,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
            use_mp=False,
        )
        
        outlier_idxs = dataset2.outliersdf.index.get_level_values('datetime').unique()
        
        intersect = outlier_idxs.intersection(white_name_dt)
        #Check if the white dt only timestamps are not in the outliers
        assert intersect.empty , "outlier timestamps found in white name dt"
        assert not outlier_timestamps.empty , "not outliers found"
        assert dataset2.outliersdf.shape[0] > dataset1.outliersdf.shape[0], 'something wrong'

       

    def test_white_records_buddy_check_with_safety_nets_dataset(
        self
    ):
        """Test white_records with buddy_check_with_safety_nets on Dataset level."""
       

        dataset = TestWhiteRecords.solutionfixer.get_solution(
            **TestWhiteRecords.solkwargs, methodname="test_import_data"
        )
        
        
        # Ensure LCZ data is present
        if not all(sta.site.flag_has_LCZ() for sta in dataset.stations):
            dataset.get_LCZ()

        # First run without whiteset
        # test_dataset = copy.deepcopy(dataset)
        # test_dataset.buddy_check_with_safetynets(
        #     obstype="temp",
        #     spatial_buddy_radius=25000,
        #     safety_net_configs=[
        #         {
        #             "category": "LCZ",
        #             "buddy_radius": 40000,
        #             "z_threshold": 1.8,
        #             "min_sample_size": 3,
        #         }
        #     ],
        #     min_sample_size=3,
        #     spatial_z_threshold=1.8,
        #     N_iter=2,
        #     use_mp=False,
        # )

        # outliers = test_dataset.outliersdf


        # Test with different structures
    
        # white_dt_only = pd.Index(
        #     outliers.reset_index()["datetime"].sample(n=33, random_state=42),
        #     name="datetime",
        # )
        white_dt_only = pd.DatetimeIndex([
            '2022-09-11 17:00:00+00:00', '2022-09-01 15:00:00+00:00',
               '2022-09-13 01:00:00+00:00', '2022-09-06 13:00:00+00:00',
               '2022-09-04 18:00:00+00:00', '2022-09-12 17:00:00+00:00',
               '2022-09-04 05:00:00+00:00', '2022-09-01 14:00:00+00:00',
               '2022-09-09 02:00:00+00:00', '2022-09-09 01:00:00+00:00',
               '2022-09-15 05:00:00+00:00', '2022-09-01 19:00:00+00:00',
               '2022-09-05 05:00:00+00:00', '2022-09-01 05:00:00+00:00',
               '2022-09-14 10:00:00+00:00', '2022-09-09 07:00:00+00:00',
               '2022-09-04 18:00:00+00:00', '2022-09-09 03:00:00+00:00',
               '2022-09-15 21:00:00+00:00', '2022-09-09 09:00:00+00:00',
               '2022-09-10 13:00:00+00:00', '2022-09-12 07:00:00+00:00',
               '2022-09-04 16:00:00+00:00', '2022-09-03 09:00:00+00:00',
               '2022-09-14 04:00:00+00:00', '2022-09-09 05:00:00+00:00',
               '2022-09-04 17:00:00+00:00', '2022-09-06 06:00:00+00:00',
               '2022-09-02 01:00:00+00:00', '2022-09-07 04:00:00+00:00',
               '2022-09-09 08:00:00+00:00', '2022-09-05 03:00:00+00:00',
               '2022-09-01 22:00:00+00:00'],
              dtype='datetime64[ns, UTC]', name='datetime', freq=None)
        
    
        dataset1 = copy.deepcopy(dataset)
        dataset1.buddy_check_with_safetynets(
            obstype="temp",
            spatial_buddy_radius=25000,
            safety_net_configs=[
                {
                    "category": "LCZ",
                    "buddy_radius": 40000,
                    "z_threshold": 1.8,
                    "min_sample_size": 3,
                }
            ],
            min_sample_size=3,
            spatial_z_threshold=1.8,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
            use_mp=False,
        )
        
        outlier_timestamps = dataset1.outliersdf.index.get_level_values('datetime').unique()
        intersect = outlier_timestamps.intersection(white_dt_only)
        #Check if the white dt only timestamps are not in the outliers
        assert intersect.empty , "outlier timestamps found in white dt only"
        assert not outlier_timestamps.empty , "not outliers found"
        assert dataset1.outliersdf.shape[0] == 150, 'something wrong with outlier count'
        
        # white_name_dt = (
        #     outliers.sample(n=21, random_state=42)
        #     .reset_index()[["name", "datetime"]]
        #     .set_index(["name", "datetime"])
        #     .index
        # )
        white_name_dt =pd.MultiIndex.from_arrays(
            [['vlinder05', 'vlinder27', 'vlinder06', 'vlinder12', 'vlinder05',
                'vlinder05', 'vlinder06', 'vlinder27', 'vlinder05', 'vlinder05',
                'vlinder05', 'vlinder06', 'vlinder28', 'vlinder06', 'vlinder05',
                'vlinder05', 'vlinder09', 'vlinder05', 'vlinder05', 'vlinder05',
                'vlinder05', 'dummy_station', 'vlinder09'],
            pd.DatetimeIndex(['2022-09-11 17:00:00+00:00', '2022-09-01 15:00:00+00:00',
               '2022-09-13 01:00:00+00:00', '2022-09-06 13:00:00+00:00',
               '2022-09-04 18:00:00+00:00', '2022-09-12 17:00:00+00:00',
               '2022-09-04 05:00:00+00:00', '2022-09-01 14:00:00+00:00',
               '2022-09-09 02:00:00+00:00', '2022-09-09 01:00:00+00:00',
               '2022-09-15 05:00:00+00:00', '2022-09-01 19:00:00+00:00',
               '2022-09-05 05:00:00+00:00', '2022-09-01 05:00:00+00:00',
               '2022-09-14 10:00:00+00:00', '2022-09-09 07:00:00+00:00',
               '2022-09-04 18:00:00+00:00', '2022-09-09 03:00:00+00:00',
               '2022-09-15 21:00:00+00:00', '2022-09-09 09:00:00+00:00',
               '2022-09-10 13:00:00+00:00',
               '2022-09-10 13:00:00+00:00', '2023-09-10 13:00:00+00:00'], #fake records
              dtype='datetime64[ns, UTC]', name='datetime', freq=None)],
            names=('name', 'datetime'))
        

        dataset2 = copy.deepcopy(dataset)
        dataset2.buddy_check_with_safetynets(
            obstype="temp",
            spatial_buddy_radius=25000,
            safety_net_configs=[
                {
                    "category": "LCZ",
                    "buddy_radius": 40000,
                    "z_threshold": 1.8,
                    "min_sample_size": 3,
                }
            ],
            min_sample_size=3,
            spatial_z_threshold=1.8,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
            use_mp=False,
        )
        
        outlier_idxs = dataset2.outliersdf.index.get_level_values('datetime').unique()
        
        intersect = outlier_idxs.intersection(white_name_dt)
        #Check if the white dt only timestamps are not in the outliers
        assert intersect.empty , "outlier timestamps found in white name dt"
        assert not outlier_timestamps.empty , "not outliers found"
        assert dataset2.outliersdf.shape[0] > dataset1.outliersdf.shape[0], 'something wrong'
        #extra sanity check to ensure the test is valid
        assert dataset2.outliersdf.shape[0] == 171, 'something wrong with outlier count'


    def test_whiterecords_get_info(self):
        """Test the WhiteSet and SensorWhiteSet classes."""
        # Create a WhiteSet with mixed levels
        white_records = pd.date_range(
            start="2023-01-01 00:00",
            periods=10,
            freq="h",
            tz="UTC",
        )
        multi_index = pd.MultiIndex.from_product(
            [
                ["vlinder05", "vlinder06", "vlinder07"],
                white_records,
                ["temp", "humidity"],
            ],
            names=["name", "datetime", "obstype"],
        )
        whiteset = metobs_toolkit.WhiteSet(multi_index)
        _ = whiteset.get_info(printout=False)

    def test_import_data(self, overwrite_solution=False):
        """Import demo dataset for white_records testing."""
        _method_name = sys._getframe().f_code.co_name

        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        # Resample to hourly for consistent testing
        dataset.resample(target_freq="1h")

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solutiondata=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj)

    def test_white_records_input_combinations(self, overwrite_solution=False):
        """Test white_records with gross_value_check on Dataset level."""
        _method_name = sys._getframe().f_code.co_name

        dataset = TestWhiteRecords.solutionfixer.get_solution(
            **TestWhiteRecords.solkwargs, methodname="test_import_data"
        )

        # Get some timestamps that would be flagged as outliers
        # First run without white_records to identify outliers
        test_dataset = copy.deepcopy(dataset)
        test_dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            use_mp=False,
        )

        # Get outliers to use as white_records
        outliers = test_dataset.outliersdf
        if outliers.empty:
            pytest.skip("No outliers found for white_records testing")

        # Test 1a: Index with only datetimes
        dataset1a = copy.deepcopy(dataset)
        white_dt_only = pd.Index(
            outliers.reset_index()["datetime"].head(20), name="datetime"
        )
        dataset1a.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
            use_mp=False,
        )
        outliers1a = dataset1a.outliersdf

        # Test 1b: Index with only name
        dataset1b = copy.deepcopy(dataset)
        white_name_only = pd.Index(
            data=["vlinder05", "vlinder05", "vlinder06", "fake"], name="name"
        )
        dataset1b.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_only),
            use_mp=False,
        )

        outliers1b = dataset1b.outliersdf

        # test on station object
        copy.deepcopy(dataset).get_station("vlinder05").gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_only),
        )

        # Test 2: MultiIndex with name and datetime
        dataset2 = copy.deepcopy(dataset)
        white_name_dt = (
            outliers.head(20)
            .reset_index()[["name", "datetime"]]
            .set_index(["name", "datetime"])
            .index
        )
        dataset2.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
            use_mp=False,
        )
        outliers2 = dataset2.outliersdf

        # Test 3: MultiIndex with obstype, name, and datetime
        dataset3 = copy.deepcopy(dataset)
        white_full = outliers.head(25).index
        dataset3.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_full),
            use_mp=False,
        )
        outliers3 = dataset3.outliersdf

        # Quick Verify that white-listed records are not in the outliers
        for white_record in white_dt_only:
            assert not any(
                outliers1a.reset_index()["datetime"] == white_record
            ), f"White-listed record {white_record} found in outliers (datetime only)"

        # Store results
        results = {
            "outliers_no_white": outliers,
            "outliers_dt_only": outliers1a,
            "outliers_name_only": outliers1b,
            "outliers_name_dt": outliers2,
            "outliers_full": outliers3,
        }

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solutiondata=results,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        for key in results:
            assert_equality(results[key], solutionobj[key])

    def test_all_qc_methods_with_whiteset(self):
        """Test all QC methods on Dataset and Station with non-default whiteset.

        This test ensures that all QC methods can accept and work with a WhiteSet
        parameter without raising errors. It doesn't validate specific behavior,
        only that the methods execute successfully.
        """
        # Get dataset
        dataset = TestWhiteRecords.solutionfixer.get_solution(
            **TestWhiteRecords.solkwargs, methodname="test_import_data"
        )

        # Create a non-default whiteset with various index structures
        # Create timestamps to whitelist
        sample_times = pd.date_range("2022-09-01", periods=10, freq="1h")

        # Test 1: WhiteSet with datetime only
        whiteset_dt = metobs_toolkit.WhiteSet(pd.Index(sample_times, name="datetime"))

        # Test 2: WhiteSet with name and datetime
        white_records_multiindex = pd.MultiIndex.from_arrays(
            [["vlinder05"] * 5 + ["vlinder06"] * 5, sample_times],
            names=["name", "datetime"],
        )
        whiteset_multi = metobs_toolkit.WhiteSet(white_records_multiindex)

        # Test 3: WhiteSet with name, obstype, and datetime
        white_records_full = pd.MultiIndex.from_arrays(
            [["vlinder05"] * 5 + ["vlinder06"] * 5, ["temp"] * 10, sample_times],
            names=["name", "obstype", "datetime"],
        )
        whiteset_full = metobs_toolkit.WhiteSet(white_records_full)

        # Test all Dataset-level QC methods with different whitesets
        try:
            # gross_value_check
            ds1 = copy.deepcopy(dataset)
            ds1.gross_value_check(
                obstype="temp",
                lower_threshold=10.0,
                upper_threshold=25.0,
                whiteset=whiteset_dt,
                use_mp=True,
            )

            # persistence_check
            ds2 = copy.deepcopy(dataset)
            ds2.persistence_check(
                obstype="temp",
                timewindow="4h",
                min_records_per_window=3,
                whiteset=whiteset_multi,
                use_mp=True,
            )

            # repetitions_check
            ds3 = copy.deepcopy(dataset)
            ds3.repetitions_check(
                obstype="temp",
                max_N_repetitions=4,
                whiteset=whiteset_full,
                use_mp=True,
            )

            # step_check
            ds4 = copy.deepcopy(dataset)
            ds4.step_check(
                obstype="temp",
                max_increase_per_second=3 / 3600,
                max_decrease_per_second=-1.0 * 3 / 3600,
                whiteset=whiteset_dt,
                use_mp=True,
            )

            # window_variation_check
            ds5 = copy.deepcopy(dataset)
            ds5.window_variation_check(
                obstype="temp",
                max_increase_per_second=3 / 3600,
                max_decrease_per_second=-1.0 * 3 / 3600,
                timewindow="3h",
                whiteset=whiteset_multi,
                use_mp=True,
            )

        except Exception as e:
            pytest.fail(f"Dataset-level QC method failed with whiteset: {e}")

        # Test all Station-level QC methods with whiteset
        try:
            station = copy.deepcopy(dataset.stations[0])

            # gross_value_check
            st1 = copy.deepcopy(station)
            st1.gross_value_check(
                obstype="temp",
                lower_threshold=10.0,
                upper_threshold=25.0,
                whiteset=whiteset_dt,
            )

            # persistence_check
            st2 = copy.deepcopy(station)
            st2.persistence_check(
                obstype="temp",
                timewindow="4h",
                min_records_per_window=3,
                whiteset=whiteset_multi,
            )

            # repetitions_check
            st3 = copy.deepcopy(station)
            st3.repetitions_check(
                obstype="temp",
                max_N_repetitions=4,
                whiteset=whiteset_full,
            )

            # step_check
            st4 = copy.deepcopy(station)
            st4.step_check(
                obstype="temp",
                max_increase_per_second=3 / 3600,
                max_decrease_per_second=-1.0 * 3 / 3600,
                whiteset=whiteset_dt,
            )

            # window_variation_check
            st5 = copy.deepcopy(station)
            st5.window_variation_check(
                obstype="temp",
                max_increase_per_second=3 / 3600,
                max_decrease_per_second=-1.0 * 3 / 3600,
                timewindow="3h",
                whiteset=whiteset_multi,
            )

        except Exception as e:
            pytest.fail(f"Station-level QC method failed with whiteset: {e}")

        # If we get here, all methods executed successfully
        assert True


# if __name__ == "__main__":
# pytest.main([__file__])
# Run all methods with overwrite_solution=False
# test_breaking_dataset = TestBreakingDataset()
# test_breaking_dataset.test_import_data(overwrite_solution=False)
# test_breaking_dataset.test_apply_qc(overwrite_solution=False)
# test_breaking_dataset.test_qc_statistics(overwrite_solution=False)

# test_demo_dataset = TestDemoDataset()
# test_demo_dataset.test_import_data(overwrite_solution=True)
# test_demo_dataset.test_buddy_check(overwrite_solution=False)
# test_demo_dataset.test_buddy_check_with_safety_nets(overwrite_solution=False)
# test_demo_dataset.test_buddy_check_with_safety_nets(overwrite_solution=False)
# test_demo_dataset.test_buddy_check_with_LCZ_safety_net(overwrite_solution=False)

# Run white_records tests
# test_white_records = TestWhiteRecords()
# test_white_records.test_import_data(overwrite_solution=False)
# test_white_records.test_white_records_input_combinations(overwrite_solution=False)
# test_white_records.test_white_records_buddy_check_dataset(overwrite_solution=True)
# test_white_records.test_white_records_buddy_check_with_safety_nets_dataset(overwrite_solution=False)
# test_white_records.test_white_records_buddy_check_with_LCZ_safety_net_dataset(overwrite_solution=True)
