import sys
from pathlib import Path

import copy

import numpy as np

# import metobs_toolkit
import pandas as pd


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer2, assert_equality, datadir

import pytest


class TestBreakingDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testbreakingdata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("testdata_breaking.csv")
    templatefile = datadir.joinpath("template_breaking.json")

    @pytest.fixture(scope="class")
    def import_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestBreakingDataset.templatefile,
            #   input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=TestBreakingDataset.datafile,
        )
        return dataset

    @pytest.fixture(scope="class")
    def regular_qc_on_dataset(self, import_dataset):
        dataset = copy.deepcopy(import_dataset)

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
        return dataset

    def test_qc_labels(self, regular_qc_on_dataset):
        dataset = copy.deepcopy(regular_qc_on_dataset)

        # reading manual labels
        man_df = pd.read_csv(TestBreakingDataset.datafile)
        # create datetime column
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

    def test_qc_with_solution(self, regular_qc_on_dataset, overwrite_solution=False):
        method_name = "test_qc_with_solution"

        dataset = copy.deepcopy(regular_qc_on_dataset)

        if overwrite_solution:
            TestBreakingDataset.solutionfixer.create_solution(
                solution=dataset,
                methodname=method_name,
                **TestBreakingDataset.solkwargs,
            )

        solutionobj = TestBreakingDataset.solutionfixer.get_solution(
            **TestBreakingDataset.solkwargs, methodname=method_name
        )

        assert_equality(dataset, solutionobj, exclude_columns=['details'] )  # dataset comparison

    def test_qc_stats_check(self, regular_qc_on_dataset, overwrite_solution=False):
        method_name = "test_qc_stats_check"
        dataset = copy.deepcopy(regular_qc_on_dataset)

        countdicts = dataset.get_qc_stats(obstype="temp", make_plot=False)

        if overwrite_solution:
            TestBreakingDataset.solutionfixer.create_solution(
                solution=countdicts,
                methodname=method_name,
                **TestBreakingDataset.solkwargs,
            )

        solutiondf = TestBreakingDataset.solutionfixer.get_solution(
            methodname=method_name, **TestBreakingDataset.solkwargs
        )
        
        for key, val in countdicts.items():
            assert_equality(val, solutiondf[key])

    @pytest.mark.mpl_image_compare
    def test_make_plot_by_label_with_outliers(self, regular_qc_on_dataset):
        dataset = copy.deepcopy(regular_qc_on_dataset)
        ax = dataset.make_plot(colorby="label", obstype="temp")
        fig = ax.get_figure()
        fig.set_size_inches(15, 5)
        return fig

    def test_get_info(self, regular_qc_on_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(regular_qc_on_dataset)
        # call get info on dataset, station and sensor level
        _ = dataset.get_info(printout=True)
        _ = dataset.get_station("Fictional").get_info(printout=True)
        _ = dataset.get_station("Fictional").get_sensor("temp").get_info(printout=True)


class TestBuddyCheck:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testbuddycheck"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("random_belgian_temp_data.csv")
    metadatafile = datadir.joinpath("random_belgian_temp_metadata.csv")
    templatefile = datadir.joinpath("template_random_belgium.json")
    
    @pytest.fixture(scope="class")
    def import_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestBuddyCheck.templatefile,
            input_metadata_file=TestBuddyCheck.metadatafile,
            input_data_file=TestBuddyCheck.datafile,
        )
        # To hourly !!
        # dataset.resample(target_freq="1h")

        dataset.stations= dataset.stations[:400] #to speedup the tests
        # dataset.get_LCZ()
        # dataset.get_altitude()
        return dataset
    
    def test_import_data(self, import_dataset, overwrite_solution=False):
        """Import demo dataset for QC testing."""
        _method_name = "test_import_data"
        dataset = copy.deepcopy(import_dataset)

        if overwrite_solution:
            TestBuddyCheck.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestBuddyCheck.solkwargs,
            )

        solutionobj = TestBuddyCheck.solutionfixer.get_solution(
            methodname=_method_name, **TestBuddyCheck.solkwargs
        )

        assert_equality(dataset, solutionobj)


    

    def test_buddy_check_one_iteration(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name
        dataset = copy.deepcopy(import_dataset)

        # test one iteration
        dataset.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,  # one iteration test
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        # overwrite solution?
        if overwrite_solution:
            TestBuddyCheck.solutionfixer.create_solution(
                solution=dataset,
                **TestBuddyCheck.solkwargs,
                methodname=_method_name,
            )

        solutionobj = TestBuddyCheck.solutionfixer.get_solution(
            **TestBuddyCheck.solkwargs, methodname=_method_name
        )

        # validate expression
        assert_equality(
            dataset, solutionobj, exclude_columns=["details"]
        )  # dataset comparison

    def test_buddy_check_more_iterations(
        self, import_dataset, overwrite_solution=False
    ):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name
        

        # test one iteration
        dataset_1iter = copy.deepcopy(import_dataset)
        dataset_1iter.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,  # one iteration test
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=True,
        )
        
        #Test 2 iterations
        dataset_2iter = copy.deepcopy(import_dataset)
        dataset_2iter.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=2,  # two iteration test
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=True,
        )
        
        #apply relative tests
        outl2 = dataset_2iter.outliersdf
        outl1 = dataset_1iter.outliersdf

        #test if all oult indexes are in outl2
        assert outl1.index.isin(outl2.index).all()
        assert outl2.shape[0] > outl1.shape[0]
        
        #absolute testings

        # overwrite solution?
        if overwrite_solution:
            TestBuddyCheck.solutionfixer.create_solution(
                solution=dataset_2iter,
                **TestBuddyCheck.solkwargs,
                methodname=_method_name,
            )

        solutionobj = TestBuddyCheck.solutionfixer.get_solution(
            **TestBuddyCheck.solkwargs, methodname=_method_name
        )

        # validate expression
        assert_equality(
            dataset_2iter, solutionobj, exclude_columns=["details"]
        )  # dataset comparison

    def test_buddy_check_no_outliers(self, import_dataset):

        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # Test that buddy check runs, with settings that does not create outliers
        dataset.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_sample_spread=1.0,
            spatial_z_threshold=12.8,  # this does noet create outliers
            use_z_robust_method=True,
            N_iter=1,
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        assert dataset.outliersdf.empty
        
        #Extra test, same settings with non-robust z-method will make outliers
        dataset = copy.deepcopy(import_dataset)
        dataset.buddy_check(
                    obstype="temp",
                    spatial_buddy_radius=25000,
                    min_sample_size=3,
                    max_alt_diff=None,
                    min_sample_spread=1.0,
                    spatial_z_threshold=7.4,  # this does noet create outliers
                    N_iter=1,
                    instantaneous_tolerance=pd.Timedelta("4min"),
                    use_z_robust_method=False,
                    lapserate=None,  # -0.0065
                    use_mp=False,
                )

        assert not dataset.outliersdf.empty

    def test_buddy_check_with_big_radius(
        self, import_dataset, overwrite_solution=False
    ):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

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

        if overwrite_solution:
            TestBuddyCheck.solutionfixer.create_solution(
                solution=dataset,
                **TestBuddyCheck.solkwargs,
                methodname=_method_name,
            )
        solutionobj = TestBuddyCheck.solutionfixer.get_solution(
            **TestBuddyCheck.solkwargs, methodname=_method_name
        )
        assert_equality(dataset, solutionobj, exclude_columns=["details"])

    def test_buddy_check_with_safety_nets(
        self, import_dataset, overwrite_solution=False
    ):
        """Test the generalized buddy_check_with_safety_nets method."""
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        # 1. get_startpoint data
        dataset_with_saftynet = copy.deepcopy(import_dataset)
        # Test 1: Using safety_net_configs with LCZ should match the LCZ safety net method
        dataset_with_saftynet = copy.deepcopy(dataset_with_saftynet)
        dataset_with_saftynet.buddy_check_with_safetynets(
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
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,
            use_mp=True,
        )
        
        #Use same settings without saftyenet to make a relative comparison
        dataset_without_saftynet = copy.deepcopy(import_dataset)
        dataset_without_saftynet.buddy_check_with_safetynets(
            obstype="temp",
            spatial_buddy_radius=25000,
            safety_net_configs=[], # No safety nets
            min_sample_size=3,
            max_alt_diff=None,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            instantaneous_tolerance=pd.Timedelta("4min"),
            lapserate=None,
            use_mp=True,
        )
        
        #Relative tests
        oult_saftynet = dataset_with_saftynet.outliersdf 
        oult_without_saftynet = dataset_without_saftynet.outliersdf


        assert oult_saftynet.index.isin(oult_without_saftynet.index).all()
        assert oult_without_saftynet.shape[0] > oult_saftynet.shape[0]
        
        

        # overwrite solution?
        if overwrite_solution:
            TestBuddyCheck.solutionfixer.create_solution(
                solution=dataset_with_saftynet,
                **TestBuddyCheck.solkwargs,
                methodname=_method_name,
            )

        # 4. Get solution
        solutionobj = TestBuddyCheck.solutionfixer.get_solution(
            **TestBuddyCheck.solkwargs, methodname=_method_name
        )

        # validate expression
        assert_equality(dataset_with_saftynet, solutionobj, exclude_columns=["details"])

    def test_buddy_check_with_safety_nets_missing_min_sample_size(self, import_dataset):
        """Test that an error is raised when min_sample_size is missing from safety_net_configs."""
        # Get dataset
        dataset = copy.deepcopy(import_dataset)

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

    def test_buddy_check_filter_args(self, import_dataset):
        """Test buddy_check with max_alt_diff, max_sample_size and min_sample_size."""
        # Baseline: default run without filters
        dataset_baseline = copy.deepcopy(import_dataset)
        dataset_baseline.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_sample_size=None,
            max_alt_diff=None,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_mp=False,
        )
        baseline_outliers = dataset_baseline.outliersdf
        assert not baseline_outliers.empty, "Baseline run must produce outliers"

        # Test 1: max_alt_diff filters buddies by altitude difference.
        # With a very small max_alt_diff many buddy groups become too small,
        # so fewer outliers are expected.
        dataset_alt = copy.deepcopy(import_dataset)
        dataset_alt.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=5,  # very strict altitude filter
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_mp=False,
        )
        alt_outliers = dataset_alt.outliersdf
        # Strict altitude filter should produce fewer (or equal) outliers
        assert len(alt_outliers) <= len(baseline_outliers)

        # Test 2: max_sample_size limits the number of buddies per station.
        dataset_maxsample = copy.deepcopy(import_dataset)
        dataset_maxsample.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=2,
            max_sample_size=3,  # limit to 3 nearest buddies
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_mp=False,
        )
        # Should run without error and produce outliers
        assert not dataset_maxsample.outliersdf.empty
        assert not dataset_maxsample.outliersdf.index.equals(baseline_outliers.index)

        # Test 3: Raising min_sample_size high enough should reduce or eliminate outliers
        # because buddy groups become too small.
        dataset_highmin = copy.deepcopy(import_dataset)
        dataset_highmin.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=50,  # very high: most groups won't meet this
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_mp=False,
        )
        highmin_outliers = dataset_highmin.outliersdf
        assert len(highmin_outliers) <= len(baseline_outliers)

    def test_buddy_check_with_lapserate(self, import_dataset):
        """Test buddy_check with lapserate altitude correction."""
        # The fixture already has altitude via get_altitude(), so lapserate can be used.

        # Run without lapserate
        dataset_no_lapse = copy.deepcopy(import_dataset)
        dataset_no_lapse.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            lapserate=None,
            use_mp=False,
        )
        outliers_no_lapse = dataset_no_lapse.outliersdf
        assert not outliers_no_lapse.empty

        # Run with standard temperature lapserate (-6.5°C/km)
        dataset_lapse = copy.deepcopy(import_dataset)
        dataset_lapse.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            lapserate=-0.0065,
            use_mp=False,
        )
        outliers_lapse = dataset_lapse.outliersdf
        assert not outliers_lapse.empty

        # Lapserate correction adjusts for altitude, so the outlier sets
        # should differ (at least in number or composition).
        differ_in_count = len(outliers_lapse) != len(outliers_no_lapse)
        differ_in_index = not outliers_lapse.index.equals(outliers_no_lapse.index)
        assert differ_in_count or differ_in_index, (
            "Lapserate correction should change the outlier set"
        )

    def test_buddy_check_only_if_previous_has_no_buddies(self, import_dataset):
        """Test safety_net_configs with only_if_previous_had_no_buddies=True."""

        # Setup: Two cascading safety nets. The second net with
        # only_if_previous_had_no_buddies=True should only test records that
        # had insufficient buddies in the first safety net.
        dataset_cascade = copy.deepcopy(import_dataset)
        dataset_cascade.buddy_check_with_safetynets(
            obstype="temp",
            spatial_buddy_radius=25000,
            safety_net_configs=[
                {
                    "category": "LCZ",
                    "buddy_radius": 25000,  # small radius: some stations won't have LCZ buddies
                    "z_threshold": 1.4,
                    "min_sample_size": 3,
                },
                {
                    "category": "LCZ",
                    "buddy_radius": 100000,  # larger radius fallback
                    "z_threshold": 1.8,
                    "min_sample_size": 2,
                    "only_if_previous_had_no_buddies": True,
                },
            ],
            min_sample_size=3,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_mp=False,
        )
        cascade_outliers = dataset_cascade.outliersdf

        # Compare against a single safety net (no cascade)
        dataset_single = copy.deepcopy(import_dataset)
        dataset_single.buddy_check_with_safetynets(
            obstype="temp",
            spatial_buddy_radius=25000,
            safety_net_configs=[
                {
                    "category": "LCZ",
                    "buddy_radius": 25000,
                    "z_threshold": 1.4,
                    "min_sample_size": 3,
                },
            ],
            min_sample_size=3,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_mp=False,
        )
        single_outliers = dataset_single.outliersdf

        # The cascade fallback can only save additional outliers (or same),
        # so cascade should have <= outliers than the single net.
        assert len(cascade_outliers) <= len(single_outliers)

    def test_buddy_check_use_z_robust_method(self, import_dataset):
        """Test that use_z_robust_method toggles between robust (MAD) and classic (std) z-scores."""
        # Run with robust method (default)
        dataset_robust = copy.deepcopy(import_dataset)
        dataset_robust.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_z_robust_method=True,
            use_mp=False,
        )
        outliers_robust = dataset_robust.outliersdf
        assert not outliers_robust.empty

        # Run with classic std method
        dataset_classic = copy.deepcopy(import_dataset)
        dataset_classic.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            min_sample_spread=1.0,
            spatial_z_threshold=2.1,
            N_iter=1,
            use_z_robust_method=False,
            use_mp=False,
        )
        outliers_classic = dataset_classic.outliersdf
        assert not outliers_classic.empty

        # The two methods use different spread estimators (MAD vs std),
        # so the outlier sets should differ in count or composition.
        differ_in_count = len(outliers_robust) != len(outliers_classic)
        differ_in_index = not outliers_robust.index.equals(outliers_classic.index)
        assert differ_in_count or differ_in_index, (
            "Robust (MAD) and classic (std) z-score methods should produce different outlier sets"
        )


class TestDemoDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    @pytest.fixture(scope="class")
    def import_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        # To hourly !!
        dataset.resample(target_freq="1h")

        dataset.get_LCZ()
        return dataset

    def test_import_data(self, import_dataset, overwrite_solution=False):
        """Import demo dataset for QC testing."""
        _method_name = "test_import_data"
        dataset = copy.deepcopy(import_dataset)

        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        assert_equality(dataset, solutionobj)
    
    def test_buddy_check_raise_errors(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsMetadataNotFound,
        )
        
        with pytest.raises(DeprecationWarning):
            # Should raise error because no altitude info is available and lapsrate is not none
            dataset.buddy_check(
                obstype="temp",
                spatial_buddy_radius=17000,
                min_sample_size=3,
                max_alt_diff=150,
                min_std=1.0, #cause a deprecation warning
                spatial_z_threshold=2.4,
                N_iter=1,
                instantaneous_tolerance=pd.Timedelta("4min"),
                lapserate=-0.0065,  # -0.0065
                use_mp=False,
            )


        with pytest.raises(MetObsMetadataNotFound):
            # Should raise error because no altitude info is available and lapsrate is not none
            dataset.buddy_check(
                obstype="temp",
                spatial_buddy_radius=17000,
                min_sample_size=3,
                max_alt_diff=150,
                min_sample_spread=1.0,
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
                min_sample_spread=1.0,
                spatial_z_threshold=2.4,
                N_iter=1,
                instantaneous_tolerance=pd.Timedelta("4min"),
                lapserate=None,  # -0.0065
                use_mp=False,
            )
        
    def test_buddy_check_with_paralelism(self, import_dataset):
        dataset = copy.deepcopy(import_dataset)
        obstype = "temp"
        
        dataset.buddy_check(obstype=obstype, use_mp=True)
        assert not dataset.outliersdf.empty
        
    def test_if_qc_can_be_chained(self, import_dataset):
        
        obstype = "temp"
        
        #Create solution
        dataset = copy.deepcopy(import_dataset)
        dataset.gross_value_check(obstype=obstype,
                                    lower_threshold=14.6,
                                    upper_threshold=26.9,
                                   use_mp=False)
        
        
        # Chain the same check
        dataset_chain = copy.deepcopy(import_dataset)
        dataset_chain.gross_value_check(obstype=obstype,
                                    lower_threshold=11.6,
                                    upper_threshold=28.9,
                                   use_mp=False)
        assert dataset.outliersdf.shape[0] > dataset_chain.outliersdf.shape[0]
        
        dataset_chain.gross_value_check(obstype=obstype,
                                    lower_threshold=11.8,
                                    upper_threshold=26.9,
                                   use_mp=False)
        dataset_chain.gross_value_check(obstype=obstype,
                                    lower_threshold=14.6,
                                    upper_threshold=31.9,
                                   use_mp=False)
        
        assert_equality(dataset, dataset_chain, exclude_columns=['details']) 
        
        #test gap related df constructions
        _ = dataset_chain.get_qc_stats(obstype=obstype, make_plot=False)
        _ = dataset_chain.qc_overview_df()        

    def test_qc_when_some_stations_missing_obs(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
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
        
        #test if test_qc works if some stations are missing obstype
        dataset.get_qc_stats(obstype=obstype)

    def test_persistence_check_unmet_window(self, import_dataset):
        """Details mention unmet window condition when min_records_per_window is too large.

        When min_records_per_window exceeds the number of records that can fit
        inside the timewindow (given the dataset frequency), the persistence check
        sets all records to 'condition_unmet' instead of 'flagged'. These records
        do NOT appear in outliersdf (which contains only flagged outliers), but the
        unmet-condition details are visible via qc_overview_df.
        """
        dataset = copy.deepcopy(import_dataset)
        # Dataset is resampled to 1h. A min_records_per_window of 100 in a 1h
        # window is impossible to satisfy, forcing the unmet-condition path.
        dataset.persistence_check(
            obstype="temp",
            timewindow="1h",
            min_records_per_window=100,
            use_mp=False,
        )

        # Persistence produced no flagged outliers because condition was unmet
        if not dataset.outliersdf.empty:
            assert "persistence" not in dataset.outliersdf["label"].values, (
                "Persistence check should not flag any outliers when window "
                "condition is not met"
            )

        # The unmet-condition details ARE accessible via qc_overview_df
        overview = dataset.qc_overview_df(subset_obstypes=["temp"])
        assert not overview.empty

        # Extract details for the persistence checkname
        persistence_details = overview[("details", "persistence")]
        unique_details = persistence_details.dropna().unique()
        # Every record should carry the "not met" detail message
        assert len(unique_details) == 1, (
            "Expected a single uniform detail message when window condition is unmet"
        )
        assert "not met" in unique_details[0].lower(), (
            f"Expected 'not met' in details, got: {unique_details[0]!r}"
        )

    def test_window_variation_check_unmet_window(self, import_dataset):
        """Details mention unmet window condition when min_records_per_window is too large.

        When min_records_per_window exceeds the number of records that can fit
        inside the timewindow (given the dataset frequency), the window_variation
        check sets all records to 'condition_unmet' instead of 'flagged'. These
        records do NOT appear in outliersdf, but the unmet-condition details are
        visible via qc_overview_df.
        """
        dataset = copy.deepcopy(import_dataset)
        # Dataset is resampled to 1h. A min_records_per_window of 100 in a 1h
        # window is impossible to satisfy, forcing the unmet-condition path.
        dataset.window_variation_check(
            obstype="temp",
            timewindow="1h",
            min_records_per_window=100,
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
            use_mp=False,
        )

        # window_variation produced no flagged outliers because condition was unmet
        if not dataset.outliersdf.empty:
            assert "window_variation" not in dataset.outliersdf["label"].values, (
                "window_variation check should not flag any outliers when window "
                "condition is not met"
            )

        # The unmet-condition details ARE accessible via qc_overview_df
        overview = dataset.qc_overview_df(subset_obstypes=["temp"])
        assert not overview.empty

        # Extract details for the window_variation checkname
        wv_details = overview[("details", "window_variation")]
        unique_details = wv_details.dropna().unique()
        # Every record should carry the "not met" detail message
        assert len(unique_details) == 1, (
            "Expected a single uniform detail message when window condition is unmet"
        )
        assert "not met" in unique_details[0].lower(), (
            f"Expected 'not met' in details, got: {unique_details[0]!r}"
        )
    


class TestQcOverviewDf:
    """Test qc_overview_df functions at sensor, station and dataset level."""

    # ------------------------------------------------------------------ #
    # fixtures
    # ------------------------------------------------------------------ #

    @pytest.fixture(scope="class")
    def breaking_import(self):
        """Import the breaking testdata and return raw dataset."""
        datafile = datadir.joinpath("testdata_breaking.csv")
        templatefile = datadir.joinpath("template_breaking.json")
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=templatefile,
            input_data_file=datafile,
        )
        return dataset

    @pytest.fixture(scope="class")
    def breaking_with_qc(self, breaking_import):
        """Run a typical QC pipeline on the breaking dataset."""
        dataset = copy.deepcopy(breaking_import)
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=-15.0,
            upper_threshold=29.0,
            use_mp=False,
        )
        dataset.gross_value_check(
            obstype="humidity",
            lower_threshold=5.0,
            upper_threshold=10.0,
            use_mp=False,
        )
        dataset.persistence_check(
            obstype="temp", timewindow="1h", min_records_per_window=3, use_mp=False,
        )
        dataset.repetitions_check(obstype="temp", max_N_repetitions=5, use_mp=False)
        dataset.step_check(
            obstype="temp",
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
            use_mp=False,
        )
        dataset.window_variation_check(
            obstype="temp", timewindow="1h", min_records_per_window=3,
            max_increase_per_second=8.0 / 3600.0, max_decrease_per_second=-10.0 / 3600.0,
            use_mp=False,
        )
        return dataset

    # ------------------------------------------------------------------ #
    # Sensordata level
    # ------------------------------------------------------------------ #

    def test_sensor_overview_after_qc(self, breaking_with_qc):
        """Sensor-level overview contains expected columns and index after QC."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")
        sensor = station.get_sensor("temp")

        df = sensor.qc_overview_df()
        assert not df.empty, "Expected non-empty overview for checked sensor"
        assert df.index.name == "datetime"
        assert "value" in df.columns.get_level_values(0)
        assert "label" in df.columns.get_level_values(0)
        assert "details" in df.columns.get_level_values(0)

    def test_sensor_overview_no_explicit_qc(self, breaking_import):
        """Sensor-level overview after import contains automatic import-time checks."""
        dataset = copy.deepcopy(breaking_import)
        station = dataset.get_station("Fictional")
        sensor = station.get_sensor("temp")

        df = sensor.qc_overview_df()
        # Import-time checks (like duplicated_timestamp) make the overview non-empty
        assert not df.empty
        assert df.index.name == "datetime"
        check_names = df.columns.get_level_values("checkname").unique().dropna()
        assert "duplicated_timestamp" in check_names

    def test_sensor_overview_humidity(self, breaking_with_qc):
        """Sensor-level overview works for non-temp obstypes that had QC."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")
        sensor = station.get_sensor("humidity")

        df = sensor.qc_overview_df()
        # gross_value_check was applied to humidity
        assert not df.empty, "Expected non-empty overview for humidity sensor"

    def test_sensor_overview_unchecked_obstype(self, breaking_with_qc):
        """Sensor with no explicit QC still has import-time check results."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")
        sensor = station.get_sensor("wind_speed")

        df = sensor.qc_overview_df()
        # Import-time checks run for every obstype
        assert not df.empty
        check_names = df.columns.get_level_values("checkname").unique().dropna()
        assert "duplicated_timestamp" in check_names
        # Explicit QC checks (like gross_value) should NOT be present
        assert "gross_value" not in check_names

    # ------------------------------------------------------------------ #
    # Station level
    # ------------------------------------------------------------------ #

    def test_station_overview_after_qc(self, breaking_with_qc):
        """Station-level overview has a (datetime, obstype) MultiIndex after QC."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")

        df = station.qc_overview_df()
        assert not df.empty
        assert df.index.names == ["datetime", "obstype"]
        assert "value" in df.columns.get_level_values(0)
        assert "label" in df.columns.get_level_values(0)
        assert "details" in df.columns.get_level_values(0)

        # Should contain at least temp and humidity (both had QC)
        obstypes_present = df.index.get_level_values("obstype").unique()
        assert "temp" in obstypes_present
        assert "humidity" in obstypes_present

    def test_station_overview_no_explicit_qc(self, breaking_import):
        """Station-level overview after import contains import-time check results."""
        dataset = copy.deepcopy(breaking_import)
        station = dataset.get_station("Fictional")

        df = station.qc_overview_df()
        assert not df.empty
        assert list(df.index.names) == ["datetime", "obstype"]
        check_names = df.columns.get_level_values("checkname").unique().dropna()
        assert "duplicated_timestamp" in check_names

    def test_station_overview_subset_obstypes(self, breaking_with_qc):
        """Station-level overview with subset_obstypes returns only requested types."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")

        df = station.qc_overview_df(subset_obstypes=["temp"])
        assert not df.empty
        obstypes_present = df.index.get_level_values("obstype").unique()
        assert list(obstypes_present) == ["temp"]

    def test_station_overview_subset_obstypes_unknown(self, breaking_with_qc):
        """Station-level overview silently ignores unknown obstype names."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")

        df = station.qc_overview_df(subset_obstypes=["nonexistent_obstype"])
        assert df.empty

    def test_station_overview_subset_obstypes_mixed(self, breaking_with_qc):
        """Subset with valid + unknown names returns only the valid one."""
        dataset = copy.deepcopy(breaking_with_qc)
        station = dataset.get_station("Fictional")

        df = station.qc_overview_df(subset_obstypes=["temp", "fake_obstype"])
        assert not df.empty
        obstypes_present = df.index.get_level_values("obstype").unique()
        assert list(obstypes_present) == ["temp"]

    # ------------------------------------------------------------------ #
    # Dataset level
    # ------------------------------------------------------------------ #

    def test_dataset_overview_after_qc(self, breaking_with_qc):
        """Dataset-level overview has (datetime, obstype, name) MultiIndex."""
        dataset = copy.deepcopy(breaking_with_qc)

        df = dataset.qc_overview_df()
        assert not df.empty
        assert df.index.names == ["datetime", "obstype", "name"]
        assert "value" in df.columns.get_level_values(0)
        assert "label" in df.columns.get_level_values(0)
        assert "details" in df.columns.get_level_values(0)

    def test_dataset_overview_no_explicit_qc(self, breaking_import):
        """Dataset-level overview after import contains import-time check results."""
        dataset = copy.deepcopy(breaking_import)

        df = dataset.qc_overview_df()
        assert not df.empty
        assert list(df.index.names) == ["datetime", "obstype", "name"]
        check_names = df.columns.get_level_values("checkname").unique().dropna()
        assert "duplicated_timestamp" in check_names

    def test_dataset_overview_subset_stations(self, breaking_with_qc):
        """Dataset-level overview respects subset_stations filter."""
        dataset = copy.deepcopy(breaking_with_qc)

        df = dataset.qc_overview_df(subset_stations=["Fictional"])
        assert not df.empty
        names = df.index.get_level_values("name").unique()
        assert list(names) == ["Fictional"]

    def test_dataset_overview_subset_obstypes(self, breaking_with_qc):
        """Dataset-level overview respects subset_obstypes filter."""
        dataset = copy.deepcopy(breaking_with_qc)

        df = dataset.qc_overview_df(subset_obstypes=["temp"])
        assert not df.empty
        obstypes_present = df.index.get_level_values("obstype").unique()
        assert list(obstypes_present) == ["temp"]

    def test_dataset_overview_both_subsets(self, breaking_with_qc):
        """Dataset-level overview works with both subsets applied."""
        dataset = copy.deepcopy(breaking_with_qc)

        df = dataset.qc_overview_df(
            subset_stations=["Fictional"], subset_obstypes=["temp"]
        )
        assert not df.empty
        assert list(df.index.get_level_values("name").unique()) == ["Fictional"]
        assert list(df.index.get_level_values("obstype").unique()) == ["temp"]

    def test_dataset_overview_matches_station_overview(self, breaking_with_qc):
        """Dataset-level overview for one station matches station-level overview."""
        dataset = copy.deepcopy(breaking_with_qc)

        ds_df = dataset.qc_overview_df(
            subset_stations=["Fictional"], subset_obstypes=["temp"]
        )
        sta_df = dataset.get_station("Fictional").qc_overview_df(
            subset_obstypes=["temp"]
        )

        # After dropping the extra 'name' level from the dataset result,
        # the two frames should have identical shapes and datetimes.
        ds_df_dropped = ds_df.droplevel("name")
        assert ds_df_dropped.shape == sta_df.shape
        assert ds_df_dropped.index.equals(sta_df.index)


class TestWhiteRecords:
    """Test white_records functionality for all QC checks on both Dataset and Station level."""

    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testwhiterecords"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    @pytest.fixture(scope="class")
    def import_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        # Resample to hourly for consistent testing
        dataset.resample(target_freq="1h")

        dataset.get_LCZ()
        return dataset

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

    def test_import_data(self, import_dataset, overwrite_solution=False):
        """Import demo dataset for white_records testing."""
        _method_name = sys._getframe().f_code.co_name

        dataset = copy.deepcopy(import_dataset)
        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj)

    @pytest.fixture(scope="class")
    def dataset_with_outliers(self, import_dataset):
        """Run gross_value_check to identify outliers for white_records testing."""
        dataset = copy.deepcopy(import_dataset)
        test_dataset = copy.deepcopy(dataset)
        test_dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            use_mp=False,
        )
        outliers = test_dataset.outliersdf
        if outliers.empty:
            pytest.skip("No outliers found for white_records testing")
        return {"dataset": dataset, "outliers": outliers}

    def test_whiteset_datetime_only(
        self, dataset_with_outliers, overwrite_solution=False
    ):
        """Test white_records with Index containing only datetimes."""
        _method_name = sys._getframe().f_code.co_name

        dataset = copy.deepcopy(dataset_with_outliers["dataset"])
        outliers = copy.deepcopy(dataset_with_outliers["outliers"])

        white_dt_only = pd.Index(
            outliers.reset_index()["datetime"].head(20), name="datetime"
        )
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
            use_mp=False,
        )
        outliers_result = dataset.outliersdf

        # Verify that white-listed records are not in the outliers
        for white_record in white_dt_only:
            assert not any(
                outliers_result.reset_index()["datetime"] == white_record
            ), f"White-listed record {white_record} found in outliers (datetime only)"

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj, exclude_columns=['details'] )

    def test_whiteset_name_only(self, dataset_with_outliers, overwrite_solution=False):
        """Test white_records with Index containing only station names."""
        _method_name = sys._getframe().f_code.co_name

        dataset = copy.deepcopy(dataset_with_outliers["dataset"])

        white_name_only = pd.Index(
            data=["vlinder05", "vlinder05", "vlinder06", "fake"], name="name"
        )
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_only),
            use_mp=False,
        )

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj, exclude_columns=['details'])

    def test_whiteset_name_only_on_station(self, dataset_with_outliers):
        """Test white_records with name-only Index on Station object."""
        dataset = copy.deepcopy(dataset_with_outliers["dataset"])

        white_name_only = pd.Index(
            data=["vlinder05", "vlinder05", "vlinder06", "fake"], name="name"
        )
        # Should not raise - test on station object
        dataset.get_station("vlinder05").gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_only),
        )

    def test_whiteset_name_and_datetime(
        self, dataset_with_outliers, overwrite_solution=False
    ):
        """Test white_records with MultiIndex containing name and datetime."""
        _method_name = sys._getframe().f_code.co_name

        dataset = copy.deepcopy(dataset_with_outliers["dataset"])
        outliers = copy.deepcopy(dataset_with_outliers["outliers"])

        white_name_dt = (
            outliers.head(20)
            .reset_index()[["name", "datetime"]]
            .set_index(["name", "datetime"])
            .index
        )
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
            use_mp=False,
        )

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj, exclude_columns=['details'])

    def test_whiteset_full_multiindex(
        self, dataset_with_outliers, overwrite_solution=False
    ):
        """Test white_records with full MultiIndex (obstype, name, datetime)."""
        _method_name = sys._getframe().f_code.co_name

        dataset = copy.deepcopy(dataset_with_outliers["dataset"])
        outliers = copy.deepcopy(dataset_with_outliers["outliers"])

        white_full = outliers.head(25).index
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=10.0,
            upper_threshold=20.0,
            whiteset=metobs_toolkit.WhiteSet(white_full),
            use_mp=False,
        )

        if overwrite_solution:
            TestWhiteRecords.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWhiteRecords.solkwargs,
            )

        solutionobj = TestWhiteRecords.solutionfixer.get_solution(
            methodname=_method_name, **TestWhiteRecords.solkwargs
        )

        assert_equality(dataset, solutionobj, exclude_columns=['details'])

    def test_white_dt_only_records_buddy_check(
        self, import_dataset,
    ):
        
        white_dt_only = pd.DatetimeIndex(
            [
                "2022-09-11 16:00:00+00:00",
                "2022-09-09 01:00:00+00:00",
                "2022-09-01 15:00:00+00:00",
                "2022-09-05 03:00:00+00:00",
                "2022-09-12 16:00:00+00:00",
                "2022-09-15 04:00:00+00:00",
                "2022-09-04 18:00:00+00:00",
                "2022-09-04 05:00:00+00:00",
                "2022-09-01 19:00:00+00:00",
                "2022-09-14 13:00:00+00:00",
                "2022-09-12 18:00:00+00:00",
                "2022-09-01 14:00:00+00:00",
                "2022-09-14 05:00:00+00:00",
                "2022-09-04 10:00:00+00:00",
                "2022-09-09 03:00:00+00:00",
                "2022-09-09 02:00:00+00:00",
                "2022-09-12 04:00:00+00:00",
                "2022-09-01 05:00:00+00:00",
                "2022-09-05 02:00:00+00:00",
                "2022-09-15 06:00:00+00:00",
                "2022-09-09 00:00:00+00:00",
                "2022-09-06 06:00:00+00:00",
                "2022-09-11 16:00:00+00:00",
                "2022-09-03 09:00:00+00:00",
                "2022-09-05 07:00:00+00:00",
                "2022-09-04 17:00:00+00:00",
                "2022-09-11 00:00:00+00:00",
                "2022-09-14 17:00:00+00:00",
                "2022-09-02 01:00:00+00:00",
                "2022-09-04 18:00:00+00:00",
                "2022-09-05 05:00:00+00:00",
                "2022-09-09 07:00:00+00:00",
                "2022-09-15 02:00:00+00:00",
                "2023-09-15 02:00:00+00:00",
            ],  # not in range
            dtype="datetime64[ns, UTC]",
            name="datetime",
            freq=None,
        )
        # Test 1: see if there is overlap with the whites, if no whites are given
        dataset1 = copy.deepcopy(import_dataset)
        dataset1.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            spatial_z_threshold=2.1,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(), #empty
            use_mp=False,
        )
        outlier_timestamps = dataset1.outliersdf.index.get_level_values(
            "datetime"
        ).unique()
        intersect = outlier_timestamps.intersection(white_dt_only)
        # Check if the white dt only timestamps are not in the outliers
        assert not intersect.empty, "There must be overlap, else the test is bad designed"
        assert not outlier_timestamps.empty, "not outliers found"
        
        
        # Test 2: test that the white dt only records are not in the outliers after buddy check
        dataset2 = copy.deepcopy(import_dataset)
        dataset2.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            spatial_z_threshold=2.1,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
            use_mp=False,
        )
        outlier_timestamps = dataset2.outliersdf.index.get_level_values(
            "datetime"
        ).unique()
        intersect = outlier_timestamps.intersection(white_dt_only)
        # Check if the white dt only timestamps are not in the outliers
        assert intersect.empty, "outlier timestamps found in white dt only"
        assert not outlier_timestamps.empty, "not outliers found"

       
    def test_white_multi_idx_records_buddy_check(
        self, import_dataset,
    ):
        """Test white_records with buddy_check on Dataset level."""
    
        white_name_dt = pd.MultiIndex.from_arrays(
            [
                [
                    "vlinder27",
                    "vlinder05",
                    "vlinder27",
                    "vlinder07",
                    "vlinder05",
                    "vlinder05",
                    "vlinder09",
                    "vlinder06",
                    "vlinder06",
                    "vlinder05",
                    "vlinder05",
                    "vlinder27",
                    "vlinder05",
                    "vlinder05",
                    "vlinder05",
                    "vlinder05",
                    "vlinder06",
                    "vlinder06",
                    "dummy_station",
                    "vlinder05",
                ],
                pd.DatetimeIndex(
                    [
                        "2022-09-11 16:00:00+00:00",
                        "2022-09-09 01:00:00+00:00",
                        "2022-09-01 15:00:00+00:00",
                        "2022-09-05 03:00:00+00:00",
                        "2022-09-12 16:00:00+00:00",
                        "2022-09-15 04:00:00+00:00",
                        "2022-09-04 18:00:00+00:00",
                        "2022-09-04 05:00:00+00:00",
                        "2022-09-01 19:00:00+00:00",
                        "2022-09-14 13:00:00+00:00",
                        "2022-09-12 18:00:00+00:00",
                        "2022-09-01 14:00:00+00:00",
                        "2022-09-14 05:00:00+00:00",
                        "2022-09-04 10:00:00+00:00",
                        "2022-09-09 03:00:00+00:00",
                        "2022-09-09 02:00:00+00:00",
                        "2022-09-12 04:00:00+00:00",
                        "2022-09-01 05:00:00+00:00",
                        "2022-09-12 04:00:00+00:00",
                        "2023-09-01 05:00:00+00:00",
                    ],  # fake records
                    dtype="datetime64[ns, UTC]",
                    name="datetime",
                    freq=None,
                ),
            ],
            names=("name", "datetime"),
        )

        # Test 1: see if there is overlap with the whites, if no whites are given
        dataset1 = copy.deepcopy(import_dataset)
        dataset1.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            spatial_z_threshold=2.1,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(), #empty
            use_mp=False,
        )
        
        #There must be overlap
        assert not dataset1.outliersdf.empty, "No outliers found, test is bad designed"
        outlier_name_dt = dataset1.outliersdf[
            dataset1.outliersdf.index.get_level_values("obstype") == "temp"
        ].index.droplevel("obstype")

        outlier_name_dt = outlier_name_dt.reorder_levels(["name", "datetime"])
        intersect = outlier_name_dt.intersection(white_name_dt.reorder_levels(["name", "datetime"]))
        assert not intersect.empty, "There must be overlap, else the test is bad designed"


        
        
        # Test 2: test that the white dt only records are not in the outliers after buddy check
        dataset2 = copy.deepcopy(import_dataset)
        dataset2.buddy_check(
            obstype="temp",
            spatial_buddy_radius=25000,
            min_sample_size=3,
            spatial_z_threshold=2.1,
            N_iter=2,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
            use_mp=False,
        )
        
        
        assert not dataset2.outliersdf.empty, "No outliers found, test is bad designed"
        outlier_name_dt = dataset2.outliersdf[
            dataset2.outliersdf.index.get_level_values("obstype") == "temp"
        ].index.droplevel("obstype")
        outlier_name_dt = outlier_name_dt.reorder_levels(["name", "datetime"])
        intersect = outlier_name_dt.intersection(white_name_dt.reorder_levels(["name", "datetime"]))
        assert intersect.empty, "outlier timestamps found in white name dt"
        
        

    def test_white_dt_only_records_buddy_check_with_safety_nets(
        self, import_dataset,
    ):
        """Test white_records with buddy_check_with_safety_nets on Dataset level."""
    
        white_dt_only = pd.DatetimeIndex(
            [
                "2022-09-11 17:00:00+00:00",
                "2022-09-01 15:00:00+00:00",
                "2022-09-13 01:00:00+00:00",
                "2022-09-06 13:00:00+00:00",
                "2022-09-04 18:00:00+00:00",
                "2022-09-12 17:00:00+00:00",
                "2022-09-04 05:00:00+00:00",
                "2022-09-01 14:00:00+00:00",
                "2022-09-09 02:00:00+00:00",
                "2022-09-09 01:00:00+00:00",
                "2022-09-15 05:00:00+00:00",
                "2022-09-01 19:00:00+00:00",
                "2022-09-05 05:00:00+00:00",
                "2022-09-01 05:00:00+00:00",
                "2022-09-14 10:00:00+00:00",
                "2022-09-09 07:00:00+00:00",
                "2022-09-04 18:00:00+00:00",
                "2022-09-09 03:00:00+00:00",
                "2022-09-15 21:00:00+00:00",
                "2022-09-09 09:00:00+00:00",
                "2022-09-10 13:00:00+00:00",
                "2022-09-12 07:00:00+00:00",
                "2022-09-04 16:00:00+00:00",
                "2022-09-03 09:00:00+00:00",
                "2022-09-14 04:00:00+00:00",
                "2022-09-09 05:00:00+00:00",
                "2022-09-04 17:00:00+00:00",
                "2022-09-06 06:00:00+00:00",
                "2022-09-02 01:00:00+00:00",
                "2022-09-07 04:00:00+00:00",
                "2022-09-09 08:00:00+00:00",
                "2022-09-05 03:00:00+00:00",
                "2022-09-01 22:00:00+00:00",
            ],
            dtype="datetime64[ns, UTC]",
            name="datetime",
            freq=None,
        )

        args = {
            "obstype": "temp",
            "spatial_buddy_radius": 25000,
            "safety_net_configs": [
                {
                    "category": "LCZ",
                    "buddy_radius": 40000,
                    "z_threshold": 1.8,
                    "min_sample_size": 3,
                }
            ],
            "min_sample_size": 3,
            "spatial_z_threshold": 1.8,
            "N_iter": 2,
            # "whiteset": metobs_toolkit.WhiteSet(white_dt_only),
            "use_mp": False,
        }

        # Test 1: see if there is overlap with the whites, if no whites are given
        dataset1 = copy.deepcopy(import_dataset)
        dataset1.buddy_check_with_safetynets(
            **args,
            whiteset=metobs_toolkit.WhiteSet(), #empty
        )
        outlier_timestamps = dataset1.outliersdf.index.get_level_values(
            "datetime"
        ).unique()
        intersect = outlier_timestamps.intersection(white_dt_only)
        # Check if the white dt only timestamps are not in the outliers
        assert not intersect.empty, "There must be overlap, else the test is bad designed"
        assert not outlier_timestamps.empty, "not outliers found"
        
        
        # Test 2: test that the white dt only records are not in the outliers after buddy check
        dataset2 = copy.deepcopy(import_dataset)
        dataset2.buddy_check_with_safetynets(
            **args,
            whiteset=metobs_toolkit.WhiteSet(white_dt_only),
        )
        outlier_timestamps = dataset2.outliersdf.index.get_level_values(
            "datetime"
        ).unique()
        intersect = outlier_timestamps.intersection(white_dt_only)
        # Check if the white dt only timestamps are not in the outliers
        assert intersect.empty, "outlier timestamps found in white dt only"
        assert not outlier_timestamps.empty, "not outliers found"


     

    def test_white_multi_idx_records_buddy_check_with_safety_nets(
        self, import_dataset,
    ):
        """Test white_records with buddy_check_with_safety_nets on Dataset level."""

        white_name_dt = pd.MultiIndex.from_arrays(
            [
                [
                    "vlinder05",
                    "vlinder27",
                    "vlinder06",
                    "vlinder12",
                    "vlinder05",
                    "vlinder05",
                    "vlinder06",
                    "vlinder27",
                    "vlinder05",
                    "vlinder05",
                    "vlinder05",
                    "vlinder06",
                    "vlinder28",
                    "vlinder06",
                    "vlinder05",
                    "vlinder05",
                    "vlinder09",
                    "vlinder05",
                    "vlinder05",
                    "vlinder05",
                    "vlinder05",
                    "dummy_station",
                    "vlinder09",
                ],
                pd.DatetimeIndex(
                    [
                        "2022-09-11 17:00:00+00:00",
                        "2022-09-01 15:00:00+00:00",
                        "2022-09-13 01:00:00+00:00",
                        "2022-09-06 13:00:00+00:00",
                        "2022-09-04 18:00:00+00:00",
                        "2022-09-12 17:00:00+00:00",
                        "2022-09-04 05:00:00+00:00",
                        "2022-09-01 14:00:00+00:00",
                        "2022-09-09 02:00:00+00:00",
                        "2022-09-09 01:00:00+00:00",
                        "2022-09-15 05:00:00+00:00",
                        "2022-09-01 19:00:00+00:00",
                        "2022-09-05 05:00:00+00:00",
                        "2022-09-01 05:00:00+00:00",
                        "2022-09-14 10:00:00+00:00",
                        "2022-09-09 07:00:00+00:00",
                        "2022-09-04 18:00:00+00:00",
                        "2022-09-09 03:00:00+00:00",
                        "2022-09-15 21:00:00+00:00",
                        "2022-09-09 09:00:00+00:00",
                        "2022-09-10 13:00:00+00:00",
                        "2022-09-10 13:00:00+00:00",
                        "2023-09-10 13:00:00+00:00",
                    ],  # fake records
                    dtype="datetime64[ns, UTC]",
                    name="datetime",
                    freq=None,
                ),
            ],
            names=("name", "datetime"),
        )

        args = {
            "obstype": "temp",
            "spatial_buddy_radius": 25000,
            "safety_net_configs": [
                {
                    "category": "LCZ",
                    "buddy_radius": 40000,
                    "z_threshold": 1.8,
                    "min_sample_size": 3,
                }
            ],
            "min_sample_size": 3,
            "spatial_z_threshold": 1.8,
            "N_iter": 2,
            # "whiteset": metobs_toolkit.WhiteSet(white_dt_only),
            "use_mp": False,
        }

        # Test 1: see if there is overlap with the whites, if no whites are given
        dataset1 = copy.deepcopy(import_dataset)
        dataset1.buddy_check_with_safetynets(
            **args,
            whiteset=metobs_toolkit.WhiteSet(), #empty
        )
        
        #There must be overlap
        assert not dataset1.outliersdf.empty, "No outliers found, test is bad designed"
        outlier_name_dt = dataset1.outliersdf[
            dataset1.outliersdf.index.get_level_values("obstype") == "temp"
        ].index.droplevel("obstype")

        outlier_name_dt = outlier_name_dt.reorder_levels(["name", "datetime"])
        intersect = outlier_name_dt.intersection(white_name_dt.reorder_levels(["name", "datetime"]))
        assert not intersect.empty, "There must be overlap, else the test is bad designed"


        
        
        # Test 2: test that the white dt only records are not in the outliers after buddy check
        dataset2 = copy.deepcopy(import_dataset)
        dataset2.buddy_check_with_safetynets(
            **args,
            whiteset=metobs_toolkit.WhiteSet(white_name_dt),
        )
        
        assert not dataset2.outliersdf.empty, "No outliers found, test is bad designed"
        outlier_name_dt = dataset2.outliersdf[
            dataset2.outliersdf.index.get_level_values("obstype") == "temp"
        ].index.droplevel("obstype")
        outlier_name_dt = outlier_name_dt.reorder_levels(["name", "datetime"])
        intersect = outlier_name_dt.intersection(white_name_dt.reorder_levels(["name", "datetime"]))
        assert intersect.empty, "outlier timestamps found in white name dt"



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

    def test_all_qc_methods_with_whiteset(self, import_dataset):
        """Test all QC methods on Dataset and Station with non-default whiteset.

        This test ensures that all QC methods can accept and work with a WhiteSet
        parameter without raising errors. It doesn't validate specific behavior,
        only that the methods execute successfully.
        """
        # Get dataset
        dataset = copy.deepcopy(import_dataset)

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


class TestQCresult:
    """Tests for the QCresult class methods."""

    def _make_qcresult(self):
        """Helper to create a simple QCresult for testing."""
        from metobs_toolkit.qcresult import QCresult

        dt_index = pd.date_range("2022-09-01", periods=5, freq="h", tz="UTC")
        flags = pd.Series(
            ["passed", "flagged", "passed", "flagged", "passed"], index=dt_index
        )
        return QCresult(
            checkname="gross_value",
            checksettings={"lower": -15, "upper": 35},
            flags=flags,
            detail="test detail",
        )

    def test_repr(self):
        """Test __repr__ returns expected string."""
        qcr = self._make_qcresult()
        assert "gross_value" in repr(qcr)

    def test_get_outlier_timestamps(self):
        """Test that get_outlier_timestamps returns flagged timestamps."""
        qcr = self._make_qcresult()
        outlier_ts = qcr.get_outlier_timestamps()
        assert len(outlier_ts) == 2  # two flagged entries
        assert isinstance(outlier_ts, pd.DatetimeIndex)

    def test_add_details_by_series(self):
        """Test updating details by series."""
        qcr = self._make_qcresult()
        dt_index = qcr.flags.index
        detail_update = pd.Series(
            ["updated detail 1", "updated detail 2"],
            index=dt_index[:2],
        )
        qcr.add_details_by_series(detail_update)
        assert qcr.details.iloc[0] == "updated detail 1"
        assert qcr.details.iloc[1] == "updated detail 2"
        # Unchanged entries still have the original detail
        assert qcr.details.iloc[2] == "test detail"

    def test_remap_timestamps(self):
        """Test remap_timestamps remaps and drops unmapped entries."""
        qcr = self._make_qcresult()
        old_index = qcr.flags.index
        new_ts = old_index + pd.Timedelta("30min")
        # Only map the first 3 timestamps, the last 2 should be dropped
        mapping = {old_index[i]: new_ts[i] for i in range(3)}

        qcr.remap_timestamps(mapping)

        assert len(qcr.flags) == 3
        assert len(qcr.details) == 3
        assert qcr.flags.index[0] == new_ts[0]

    def test_create_outliersdf_subset(self):
        """Test create_outliersdf with subset_to_outliers=True."""
        qcr = self._make_qcresult()
        df = qcr.create_outliersdf(subset_to_outliers=True)
        assert len(df) == 2  # only the flagged entries
        assert "label" in df.columns
        assert "details" in df.columns

    def test_create_outliersdf_all(self):
        """Test create_outliersdf with subset_to_outliers=False."""
        qcr = self._make_qcresult()
        df = qcr.create_outliersdf(subset_to_outliers=False)
        assert len(df) == 5  # all entries

    def test_create_outliersdf_empty(self):
        """Test create_outliersdf returns empty df when no outliers."""
        from metobs_toolkit.qcresult import QCresult

        dt_index = pd.date_range("2022-09-01", periods=3, freq="h", tz="UTC")
        flags = pd.Series(["passed", "passed", "passed"], index=dt_index)
        qcr = QCresult(
            checkname="gross_value",
            checksettings={},
            flags=flags,
        )
        df = qcr.create_outliersdf(subset_to_outliers=True)
        assert df.empty

    def test_init_invalid_index_raises(self):
        """Test that non-DatetimeIndex raises TypeError."""
        from metobs_toolkit.qcresult import QCresult

        flags = pd.Series(["passed", "flagged"], index=[0, 1])
        with pytest.raises(TypeError, match="DatetimeIndex"):
            QCresult(checkname="test", checksettings={}, flags=flags)


class TestConvertOutliersToGaps:
    """Test convert_outliers_to_gaps at Dataset and Station level."""

    @pytest.fixture(scope="class")
    def dataset_with_outliers(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        dataset.resample(target_freq="1h")
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=-15.0,
            upper_threshold=29.0,
            use_mp=False,
        )
        return dataset

    def test_dataset_has_outliers_before(self, dataset_with_outliers):
        """Verify fixture has outliers."""
        dataset = copy.deepcopy(dataset_with_outliers)
        assert not dataset.outliersdf.empty

    def test_convert_outliers_to_gaps_dataset(self, dataset_with_outliers):
        """Test converting outliers to gaps removes outliers and creates gaps."""
        dataset = copy.deepcopy(dataset_with_outliers)
        n_outliers_before = len(dataset.outliersdf)
        assert n_outliers_before > 0

        dataset.convert_outliers_to_gaps(obstype="temp")

        # After conversion, outliers for temp should be empty
        if not dataset.outliersdf.empty:
            # If non-empty, temp should not be present
            assert "temp" not in dataset.outliersdf.index.get_level_values(
                "obstype"
            ).unique()

    def test_convert_outliers_to_gaps_station(self, dataset_with_outliers):
        """Test converting outliers at station level."""
        dataset = copy.deepcopy(dataset_with_outliers)
        station = dataset.stations[0]
        # Ensure station has outliers
        sensor = station.get_sensor("temp")
        if not sensor.outliersdf.empty:
            station.convert_outliers_to_gaps(obstype="temp")
            assert sensor.outliersdf.empty


class TestQCStats:
    """Test get_qc_stats and get_qc_freq_statistics."""

    @pytest.fixture(scope="class")
    def dataset_after_qc(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        dataset.resample(target_freq="1h")
        dataset.gross_value_check(
            obstype="temp",
            lower_threshold=-15.0,
            upper_threshold=29.0,
            use_mp=False,
        )
        dataset.persistence_check(obstype="temp", use_mp=False)
        return dataset

    def test_get_qc_stats_dict(self, dataset_after_qc):
        """Test get_qc_stats returns dict with expected keys when make_plot=False."""
        dataset = copy.deepcopy(dataset_after_qc)
        result = dataset.get_qc_stats(obstype="temp", make_plot=False)
        assert isinstance(result, dict)
        assert "all_labels" in result
        assert "outlier_labels" in result
        assert "per_check_labels" in result
        assert not result["all_labels"].empty
        assert not result["per_check_labels"].empty

    def test_get_qc_freq_statistics_sensor(self, dataset_after_qc):
        """Test SensorData.get_qc_freq_statistics returns expected structure."""
        dataset = copy.deepcopy(dataset_after_qc)
        station = dataset.stations[0]
        sensor = station.get_sensor("temp")
        stats = sensor.get_qc_freq_statistics()
        assert isinstance(stats, pd.Series)
        assert stats.name == "counts"
        # MultiIndex with checkname and flag levels
        assert stats.index.names == ["checkname", "flag"]
        # Should have entries for gross_value and persistence
        checknames = stats.index.get_level_values("checkname").unique()
        assert "gross_value" in checknames
        assert "persistence" in checknames

    def test_get_qc_stats_missing_obstype(self, dataset_after_qc):
        """Test get_qc_stats returns None for nonexistent obstype."""
        dataset = copy.deepcopy(dataset_after_qc)
        result = dataset.get_qc_stats(obstype="nonexistent_obs", make_plot=False)
        assert result is None




if __name__ == "__main__":
    # When running outside pytest
    OVERWRITE = False
    # test_breaking_dataset = TestBreakingDataset()
    # Manually call fixtures and pass results to tests
    # Access the original unwrapped function via __wrapped__
    # imported_dataset = test_breaking_dataset.import_dataset.__wrapped__(test_breaking_dataset)
    # qc_dataset = test_breaking_dataset.regular_qc_on_dataset.__wrapped__(
    #     test_breaking_dataset, imported_dataset
    # )

    # test_breaking_dataset.test_qc_labels(qc_dataset)
    # test_breaking_dataset.test_qc_with_solution(qc_dataset, overwrite_solution=False)
    # test_breaking_dataset.test_qc_stats_check(qc_dataset, overwrite_solution=OVERWRITE)
    # test_breaking_dataset.test_make_plot_by_label_with_outliers(qc_dataset)
    # test_breaking_dataset.test_get_info(qc_dataset)

    test_buddy_check = TestBuddyCheck()
    buddycheckdataset = test_buddy_check.import_dataset.__wrapped__(
        test_buddy_check
    )
    # test_buddy_check.test_import_data(buddycheckdataset, overwrite_solution=OVERWRITE)
   
    # test_buddy_check.test_buddy_check_one_iteration(buddycheckdataset, overwrite_solution=OVERWRITE)
    # test_buddy_check.test_buddy_check_more_iterations(buddycheckdataset, overwrite_solution=OVERWRITE)
    # test_buddy_check.test_buddy_check_no_outliers(buddycheckdataset)
    # test_buddy_check.test_buddy_check_with_big_radius(
    #     buddycheckdataset, overwrite_solution=OVERWRITE
    # )
    # test_buddy_check.test_buddy_check_with_safety_nets(buddycheckdataset, overwrite_solution=OVERWRITE)
    # test_buddy_check.test_buddy_check_with_safety_nets_missing_min_sample_size(buddycheckdataset)


    # test_demo_dataset = TestDemoDataset()
    # imported_demo_dataset = test_demo_dataset.import_dataset.__wrapped__(
    #     test_demo_dataset
    # )
    # test_demo_dataset.test_import_data(imported_demo_dataset, overwrite_solution=OVERWRITE)
    # test_demo_dataset.test_qc_when_some_stations_missing_obs(imported_demo_dataset)
    
    # Run white_records tests
    # test_white_records = TestWhiteRecords()
    # imported_wr_dataset = test_white_records.import_dataset.__wrapped__(test_white_records)
    # dataset_with_outliers = test_white_records.dataset_with_outliers.__wrapped__(test_white_records, imported_wr_dataset)

    # test_white_records.test_whiterecords_reprs()
    # test_white_records.test_import_data(imported_wr_dataset, overwrite_solution=OVERWRITE)
    # test_white_records.test_whiteset_datetime_only(dataset_with_outliers, overwrite_solution=OVERWRITE)
    # test_white_records.test_whiteset_name_only(dataset_with_outliers, overwrite_solution=OVERWRITE)
    # test_white_records.test_whiteset_name_only_on_station(dataset_with_outliers)
    # test_white_records.test_whiteset_name_and_datetime(dataset_with_outliers, overwrite_solution=OVERWRITE)
    # test_white_records.test_whiteset_full_multiindex(dataset_with_outliers, overwrite_solution=OVERWRITE)
    # test_white_records.test_white_dt_only_records_buddy_check(imported_wr_dataset, overwrite_solution=OVERWRITE)
    # test_white_records.test_white_multi_idx_records_buddy_check(imported_wr_dataset, overwrite_solution=OVERWRITE)
    # test_white_records.test_white_dt_only_records_buddy_check_with_safety_nets(imported_wr_dataset, overwrite_solution=OVERWRITE)
    # test_white_records.test_white_multi_idx_records_buddy_check_with_safety_nets(imported_wr_dataset, overwrite_solution=OVERWRITE)
    # test_white_records.test_whiterecords_get_info()
    # test_white_records.test_all_qc_methods_with_whiteset(imported_wr_dataset)
