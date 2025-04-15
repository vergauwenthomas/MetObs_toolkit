import pytest
import sys
from pathlib import Path

# import metobs_toolkit
libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("toolkit_tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer
import pytest

# data folder
datadir = libfolder.joinpath("tests").joinpath("test_data")


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
        dataset.get_lcz()  # for aggregation
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
        test_expr = ana == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=ana, solution=solutionobj
            )
        # 6. assert the equality
        # IF THIS FAILS, it can be an issue with the dataset.resample method
        assert test_expr

    def test_if_analysis_can_be_created_from_dataset(self, overwrite_solution=False):
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
                **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = ana.df.equals(solutionobj)  # dataframe comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=ana.df, solution=solutionobj
            )
        # 6. assert the equality
        # IF THIS FAILS, it can be an issue with the dataset.resample method
        assert test_expr

    def test_aggregate_df_method(self, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )
        ana.get_info()
        aggdf = ana.aggregate_df(trgobstype="humidity", agg=["lcz", "season", "hour"])

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
        test_expr = aggdf.equals(solutionobj)  # dataframe comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestDemoDataset.solutionfixer.create_a_diff(
                to_check=aggdf, solution=solutionobj
            )
        # 6. assert the equality
        # IF THIS FAILS, it can be an issue with the dataset.resample method
        assert test_expr

    @pytest.mark.mpl_image_compare
    def test_diurnal_cycle_plot(self):
        #  1. get_startpoint data
        ana = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = ana.plot_diurnal_cycle(trgobstype="temp", colorby="lcz")
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
    # test.test_if_analysis_can_be_created_from_dataset(overwrite_solution=False)
    # test.test_aggregate_df_method(overwrite_solution=False)
    # test.test_import_data(overwrite_solution=False)
