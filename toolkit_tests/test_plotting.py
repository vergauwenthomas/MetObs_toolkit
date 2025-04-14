import pytest
import sys
from pathlib import Path
import copy
import matplotlib.pyplot as plt

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


print(
    "To Overwrite the solutions, run: \n pytest test_plotting.py --mpl-generate-path=baseline"
)


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
        # IF THIS FAILS, it can be an issue with the dataset.resample method
        assert test_expr

    @pytest.mark.mpl_image_compare
    def test_dataset_timeseries_plotting_by_label(self):

        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="label")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_timeseries_plotting_by_station(self):

        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="station")
        fig = ax.get_figure()
        return fig

    # @pytest.mark.mpl_image_compare
    def test_station_timeseries_plotting_existing_ax(self):

        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        sta = dataset.get_station("vlinder02")
        fig, ax = plt.subplots(figsize=(8, 4))

        # 2. apply a metobs manipulation
        ax = sta.make_plot(colorby="label", show_outliers=False, ax=ax)
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_station_timeseries_with_modeldata(self):
        #  1. get_startpoint data
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot(show_modeldata=True)
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_station_modeldata_timeseries(self):
        #  1. get_startpoint data
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot_of_modeldata(obstype="temp")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_modeldata_timeseries_plot(self):
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        ax = dataset_with_era.make_plot_of_modeldata()
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_color_by_station_and_modeldata_timeseries_plot(self):
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        ax = dataset_with_era.make_plot(colorby="station", show_modeldata=True)
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_color_by_label_and_modeldata_timeseries_plot(self):
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        ax = dataset_with_era.make_plot(colorby="label", show_modeldata=True)
        fig = ax.get_figure()
        return fig


# ------------------------------------------
#    test plotting timeseries with GF labels
# ------------------------------------------


class TestDataWithGaps:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdatawithgaps"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    @pytest.mark.mpl_image_compare
    def test_interpolated_timeseries_plot(self):
        dataset_with_gf = TestDataWithGaps.solutionfixer.get_solution(
            testfile="test_gf",  # OTHER TEST FILE!
            classname="TestDataWithGaps",
            methodname="test_interpolation_on_dataset_A",
        )
        ax = dataset_with_gf.make_plot(colorby="label")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_raw_modeldata_gf_timeseries_plot(self):
        dataset_with_gf = TestDataWithGaps.solutionfixer.get_solution(
            testfile="test_gf",  # OTHER TEST FILE!
            classname="TestDataWithGaps",
            methodname="test_raw_modeldata_gapfill_A",
        )
        ax = dataset_with_gf.make_plot(colorby="label")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_debias_modeldata_gf_timeseries_plot(self):
        dataset_with_gf = TestDataWithGaps.solutionfixer.get_solution(
            testfile="test_gf",  # OTHER TEST FILE!
            classname="TestDataWithGaps",
            methodname="test_debias_modeldata_gapfill_A",
        )
        ax = dataset_with_gf.make_plot(colorby="label")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_diurnal_debias_modeldata_gf_timeseries_plot(self):
        dataset_with_gf = TestDataWithGaps.solutionfixer.get_solution(
            testfile="test_gf",  # OTHER TEST FILE!
            classname="TestDataWithGaps",
            methodname="test_diurnal_debias_modeldata_gapfill_A",
        )
        ax = dataset_with_gf.make_plot(colorby="label")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_weighted_diurnal_debias_modeldata_gf_timeseries_plot(self):
        dataset_with_gf = TestDataWithGaps.solutionfixer.get_solution(
            testfile="test_gf",  # OTHER TEST FILE!
            classname="TestDataWithGaps",
            methodname="test_weighted_diurnal_debias_modeldata_gapfill_A",
        )
        ax = dataset_with_gf.make_plot(colorby="label")
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
    test.test_import_data(overwrite_solution=False)
