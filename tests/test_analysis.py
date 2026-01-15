import pytest
import sys
import copy
from pathlib import Path
import pandas as pd


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer2, assert_equality, datadir


class TestDemoDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    @pytest.fixture(scope='class')
    def import_analysis(self):
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
        
        return ana
          
    def test_import_data(self, import_analysis, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        ana = copy.deepcopy(import_analysis)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=ana, methodname=_method_name, **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana, solutionobj)  # analysis comparison

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
                solution=ana,  
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana, solutionobj)  

    def test_basic_analysis_method_calls(self, import_analysis, overwrite_solution=False):
        #  1. get_startpoint data
        ana = copy.deepcopy(import_analysis)

        # calls without value test
        ana.get_info()
        df = ana.df

        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=df,
                methodname=sys._getframe().f_code.co_name,
                **TestDemoDataset.solkwargs,
            )
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=sys._getframe().f_code.co_name, **TestDemoDataset.solkwargs
        )   
        assert_equality(df, solutionobj)  # dataframe comparison

    def test_filtering(self, import_analysis, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana =  copy.deepcopy(import_analysis)
        # filter data test
        ana.apply_filter_on_records("(wind_speed <= 2.5) & (hour > 12) & (hour < 20)")
        ana.apply_filter_on_records('season=="autumn" | season=="winter"')
        ana.apply_filter_on_metadata("LCZ == 'Large lowrise'")

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=ana,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana, solutionobj)  # dataframe comparison
        
        with pytest.raises(AssertionError):
            assert_equality(ana, copy.deepcopy(import_analysis))

    def test_subsetting_time(self, import_analysis, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana = copy.deepcopy(import_analysis)

        # filter data test
        startstr = pd.Timestamp("2022-09-04 16:29:36")
        endstr = pd.Timestamp("2022-09-13 07:00:18")
        ana.subset_period(startdt=startstr, enddt=endstr)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=ana,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(ana, solutionobj)  # dataframe comparison
        
        with pytest.raises(AssertionError):
            assert_equality(ana, copy.deepcopy(import_analysis))
            
    def test_aggregate_df_method(self, import_analysis, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name

        #  1. get_startpoint data
        ana = copy.deepcopy(import_analysis)

        aggdf = ana.aggregate_df(obstype="humidity", agg=["LCZ", "season", "hour"])

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=aggdf, methodname=_method_name, **TestDemoDataset.solkwargs
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(aggdf, solutionobj)  # dataframe comparison
        assert_equality(ana, copy.deepcopy(import_analysis))

    @pytest.mark.mpl_image_compare
    def test_diurnal_cycle_plot(self, import_analysis):
        #  1. get_startpoint data
        ana = copy.deepcopy(import_analysis)
        # 2. apply a metobs manipulation
        ax = ana.plot_diurnal_cycle(obstype="temp", colorby="LCZ")
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.mpl_image_compare
    def test_diurnal_cycle_plot_with_reference(self, import_analysis):
        #  1. get_startpoint data
        ana = copy.deepcopy(import_analysis)

        # 2. apply a metobs manipulation
        ax = ana.plot_diurnal_cycle_with_reference_station(
            ref_station="vlinder02", obstype="temp", colorby="LCZ"
        )
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig


if __name__ == "__main__":
    print(
        "To Overwrite the solutions, run: \n pytest test_plotting.py  --mpl --mpl-generate-path=baseline"
    )

    print(
        "To checkout the differences, run: \n pytest test_plotting.py --mpl --mpl-generate-summary=html "
    )

    OVERWRITE_SOLUTION = False
    test = TestDemoDataset()
    
    analysis = test.import_analysis.__wrapped__(test)
    
    test.test_import_data(analysis, overwrite_solution=OVERWRITE_SOLUTION)
    test.test_if_analysis_can_be_created_from_station(overwrite_solution=OVERWRITE_SOLUTION)
    test.test_aggregate_df_method(analysis, overwrite_solution=OVERWRITE_SOLUTION)
    test.test_basic_analysis_method_calls(analysis, overwrite_solution=OVERWRITE_SOLUTION)
    test.test_filtering(analysis, overwrite_solution=OVERWRITE_SOLUTION)
    test.test_subsetting_time(analysis, overwrite_solution=OVERWRITE_SOLUTION)
    
    
    
    test.test_diurnal_cycle_plot(analysis)
    test.test_diurnal_cycle_plot_with_reference(analysis)
   
