import pytest
import sys
from pathlib import Path
import copy

# import metobs_toolkit
import pandas as pd


libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
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
            target_obstype="temp", lower_threshold=-15.0, upper_threshold=29.0
        )
        # fake check on humidity to see if this does not affect temp records
        dataset.gross_value_check(
            target_obstype="humidity", lower_threshold=5.0, upper_threshold=10.0
        )

        dataset.persistence_check(
            target_obstype="temp",
            timewindow="1h",
            min_records_per_window=3,
            use_mp=True,
        )

        dataset.repetitions_check(
            target_obstype="temp", max_N_repetitions=5, use_mp=True
        )
        dataset.step_check(
            target_obstype="temp",
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
        )
        dataset.window_variation_check(
            target_obstype="temp",
            timewindow="1h",
            min_records_per_window=3,
            max_increase_per_second=8.0 / 3600.0,
            max_decrease_per_second=-10.0 / 3600.0,
            use_mp=True,
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
        statsdf = dataset.get_qc_stats(target_obstype="temp", make_plot=False)
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
        _statsdf = dataset.get_qc_stats(target_obstype="temp", make_plot=True)


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
                target_obstype="temp",
                buddy_radius=17000,
                min_sample_size=3,
                max_alt_diff=150,
                min_std=1.0,
                std_threshold=2.4,
                N_iter=1,
                instantanious_tolerance=pd.Timedelta("4min"),
                lapserate=-0.0065,  # -0.0065
                use_mp=False,
            )

        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsMetadataNotFound,
        )

        with pytest.raises(MetObsMetadataNotFound):
            # Should raise error because no altitude info is available and max_alt_diffis not none
            dataset.buddy_check(
                target_obstype="temp",
                buddy_radius=17000,
                min_sample_size=3,
                max_alt_diff=150,
                min_std=1.0,
                std_threshold=2.4,
                N_iter=1,
                instantanious_tolerance=pd.Timedelta("4min"),
                lapserate=None,  # -0.0065
                use_mp=False,
            )

        # Test that buddy check runs, with settings that does not create outliers
        dataset.buddy_check(
            target_obstype="temp",
            buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            std_threshold=5.9,  # this does noet create outliers
            N_iter=1,
            instantanious_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        assert dataset.outliersdf.empty

        # Now create outliers with the buddy check

        dataset1 = copy.deepcopy(dataset)  # used to test 1 iteration
        dataset2 = copy.deepcopy(dataset)  # use to test 2 iterations

        # test one iteration
        dataset1.buddy_check(
            target_obstype="temp",
            buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            std_threshold=2.1,
            N_iter=1,  # one iteration test
            instantanious_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        outliersdf_1_iter = dataset1.outliersdf

        # test two iteration
        dataset2.buddy_check(
            target_obstype="temp",
            buddy_radius=25000,
            min_sample_size=3,
            max_alt_diff=None,
            min_std=1.0,
            std_threshold=2.1,
            N_iter=2,  # one iteration test
            instantanious_tolerance=pd.Timedelta("4min"),
            lapserate=None,  # -0.0065
            use_mp=False,
        )

        outliersdf_2_iter = dataset2.outliersdf

        assert not outliersdf_1_iter.equals(
            outliersdf_2_iter
        )  # else this check is not relevant

        # overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=outliersdf_1_iter,
                **TestDemoDataset.solkwargs,
                methodname=_method_name + "_1_iter",
            )
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=outliersdf_2_iter,
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
        assert_equality(outliersdf_1_iter, solutionobj_1iter)  # dataset comparison

        assert_equality(outliersdf_2_iter, solutionobj_2iter)  # dataset comparison


if __name__ == "__main__":
    # pytest.main([__file__])
    # Run all methods with overwrite_solution=True
    test_breaking_dataset = TestBreakingDataset()
    test_breaking_dataset.test_import_data(overwrite_solution=False)
    test_breaking_dataset.test_apply_qc(overwrite_solution=False)
    test_breaking_dataset.test_qc_statistics(overwrite_solution=False)

    test_demo_dataset = TestDemoDataset()
    test_demo_dataset.test_import_data(overwrite_solution=False)
    test_demo_dataset.test_buddy_check(overwrite_solution=False)
