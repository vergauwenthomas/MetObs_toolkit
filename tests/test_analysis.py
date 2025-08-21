import pytest
import sys
from pathlib import Path
import pandas as pd



# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir


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
        dataset.get_LCZ()  # for aggregation
        # create analysis
        dataset.resample(target_freq="30min")
        ana = metobs_toolkit.Analysis(Dataholder=dataset)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=ana, methodname=_method_name, **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana, solutionobj)  # dataset comparison

    def test_if_analysis_can_be_created_from_station(self, overwrite_solution=False):
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
        # create analysis
        dataset.resample(target_freq="15min")
        station = dataset.get_station("vlinder06")
        ana = metobs_toolkit.Analysis(Dataholder=station)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=ana.df,  # test dataframe
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana.df, solutionobj)  # dataframe comparison

    def test_basic_analysis_method_calls(self):
        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # calls without value test
        ana.get_info()
        ana.df

    def test_filtering(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # filter data test
        ana.apply_filter_on_records("(wind_speed <= 2.5) & (hour > 12) & (hour < 20)")
        ana.apply_filter_on_records('season=="autumn" | season=="winter"')
        ana.apply_filter_on_metadata("LCZ == 'Large lowrise'")

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=ana.fulldf,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana.fulldf, solutionobj)  # dataframe comparison

    def test_subsetting_time(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # filter data test
        startstr = pd.Timestamp("2022-09-04 16:29:36")
        endstr = pd.Timestamp("2022-09-13 07:00:18")
        ana.subset_period(startdt=startstr, enddt=endstr)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=ana.fulldf,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana.fulldf, solutionobj)  # dataframe comparison

    def test_aggregate_df_method(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        aggdf = ana.aggregate_df(trgobstype="humidity", agg=["LCZ", "season", "hour"])

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solutiondata=aggdf, methodname=_method_name, **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(aggdf, solutionobj)  # dataframe comparison

    @pytest.mark.mpl_image_compare
    def test_diurnal_cycle_plot(self):
        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = ana.plot_diurnal_cycle(trgobstype="temp", colorby="LCZ")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_diurnal_cycle_plot_with_reference(self):
        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = ana.plot_diurnal_cycle_with_reference_station(
            ref_station="vlinder02", trgobstype="temp", colorby="LCZ"
        )
        fig = ax.get_figure()
        return fig


if __name__ == "__main__":
    print(
        "To Overwrite the solutions, run: \n pytest test_plotting.py  --mpl --mpl-generate-path=baseline"
    )

    print(
        "To checkout the differences, run: \n pytest test_plotting.py --mpl --mpl-generate-summary=html "
    )

    test = TestDemoDataset()
    # test.test_import_data(overwrite_solution=False)
    # test.test_if_analysis_can_be_created_from_station(overwrite_solution=False)
    # test.test_aggregate_df_method(overwrite_solution=False)
    # test.test_basic_analysis_method_calls()
    # test.test_filtering(overwrite_solution=False)
    # test.test_subsetting_time(overwrite_solution=False)
