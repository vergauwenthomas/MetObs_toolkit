import pytest
import sys
from pathlib import Path
import copy
import matplotlib.pyplot as plt

# import metobs_toolkit
import pandas as pd


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir
import shutil


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
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    @pytest.mark.mpl_image_compare
    def test_dataset_timeseries_plotting_by_label(self):
        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="label", figkwargs={"figsize": (10, 5)})
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_timeseries_plotting_by_station(self):
        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="station", figkwargs={"figsize": (10, 5)})
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_test_show_outliers_labelby_station(self):
        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        dataset.repetitions_check(max_N_repetitions=8)
        ax = dataset.make_plot(colorby="station", show_outliers=False)
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_test_show_outliers_labelby_labels(self):
        #  1. get_startpoint data
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # 2. apply a metobs manipulation
        dataset.repetitions_check(max_N_repetitions=8)
        ax = dataset.make_plot(colorby="label", show_outliers=False)
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
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
    def test_station_plot_of_modeldata_with_modelname(self):
        """Test station.make_plot_of_modeldata with modelname argument."""
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot_of_modeldata(obstype="temp", modelname="ERA5-land")
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_plot_of_modeldata_with_modelname(self):
        """Test station.make_plot_of_modeldata with modelname argument."""
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        ax = dataset_with_era.make_plot_of_modeldata(
            obstype="temp", modelname="ERA5-land"
        )
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_station_plot_humidity_with_temp_modeldata(self):
        """Test station.make_plot with humidity obstype and temp modeldata."""
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot(
            obstype="humidity",
            modelobstype="temp",
            show_modeldata=True,
            modeldata_kwargs={
                "modelname": "ERA5-land",
                "modelvariable": "temperature_2m",
            },
        )
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_station_plot_temp_with_modeldata_kwargs(self):
        """Test station.make_plot with temp obstype and modeldata kwargs."""
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot(
            obstype="temp",
            show_modeldata=True,
            modeldata_kwargs={
                "modelname": "ERA5-land",
                "modelvariable": "temperature_2m",
            },
        )
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_dataset_plot_humidity_with_modelvariable(self):
        """Test dataset.make_plot with humidity obstype and modelvariable."""
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        ax = dataset_with_era.make_plot(
            obstype="humidity",
            show_modeldata=True,
            modelobstype="temp",
            modeldata_kwargs={
                "modelvariable": "temperature_2m",
            },
        )
        fig = ax.get_figure()
        return fig

    @pytest.mark.mpl_image_compare
    def test_modeldatatimeseries_timeseries(self):
        #  1. get_startpoint data
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        modelseries = dataset_with_era.get_station("vlinder05").get_modeltimeseries(
            "temp"
        )
        ax = modelseries.make_plot()
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

    def test_linekwargs_in_make_plot_of_modeldata(self):
        """Test that linekwargs works properly in make_plot_of_modeldata."""
        # Get a dataset with model data
        dataset_with_era = TestDemoDataset.solutionfixer.get_solution(
            testfile="test_gee",  # OTHER TEST FILE!
            classname="TestDemoDataset",
            methodname="test_ERA5_extraction",
        )

        station = dataset_with_era.get_station("vlinder05")

        # Test 1: Station.make_plot_of_modeldata with linekwargs
        fig, ax = plt.subplots()
        try:
            ax = station.make_plot_of_modeldata(
                ax=ax,
                obstype="temp",
                modelname="ERA5-land",
                linekwargs={"ls": ":", "color": "red"},
            )
            # Check that the line was actually added
            lines = ax.get_lines()
            assert len(lines) > 0, "No lines were plotted"
            # Check that at least one line has the specified properties
            has_red_dotted = False
            for line in lines:
                if line.get_color() == "red" and line.get_linestyle() == ":":
                    has_red_dotted = True
                    break
            assert has_red_dotted, "No line with red color and dotted style found"
        finally:
            plt.close(fig)

        # Test 2: Dataset.make_plot_of_modeldata with linekwargs
        fig, ax = plt.subplots()
        try:
            ax = dataset_with_era.make_plot_of_modeldata(
                ax=ax,
                obstype="temp",
                modelname="ERA5-land",
                linekwargs={"ls": "--", "color": "blue"},
            )
            # Check that lines were added
            lines = ax.get_lines()
            assert len(lines) > 0, "No lines were plotted"
            # Check that at least one line has the specified properties
            has_blue_dashed = False
            for line in lines:
                if line.get_color() == "blue" and line.get_linestyle() == "--":
                    has_blue_dashed = True
                    break
            assert has_blue_dashed, "No line with blue color and dashed style found"
        finally:
            plt.close(fig)

        # Test 3: Multiple models with different styles on same axes
        fig, ax = plt.subplots()
        try:
            # This is the use case from the original issue
            ax = station.make_plot(obstype="temp", show_outliers=False, ax=ax)
            
            # Note: Since we only have ERA5-land in the test data, we can't test 
            # multiple different models, but we can verify the mechanism works
            ax = station.make_plot_of_modeldata(
                ax=ax,
                obstype="temp",
                modelname="ERA5-land",
                linekwargs={"ls": ":", "color": "red"},
            )
            
            lines = ax.get_lines()
            assert len(lines) > 0, "No lines were plotted"
        finally:
            plt.close(fig)

        # Test 4: Verify that invalid figkwargs in make_plot doesn't cause issues
        # This used to cause an AttributeError
        try:
            # figkwargs should only be used for figure creation, not passed down
            ax = station.make_plot(
                obstype="temp",
                show_modeldata=False,
                figkwargs={"figsize": (10, 4)},
            )
            plt.close(ax.get_figure())
        except AttributeError as e:
            pytest.fail(f"make_plot raised AttributeError with figkwargs: {e}")


# ------------------------------------------
#    test plotting timeseries with GF labels
# ------------------------------------------


class TestDataWithGaps:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdatawithgaps"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    @pytest.mark.mpl_image_compare
    def test_dataset_test_show_gaps_labelby_labels(self):
        #  1. get_startpoint data
        dataset = TestDataWithGaps.solutionfixer.get_solution(
            testfile="test_gf",  # OTHER TEST FILE!
            classname="TestDataWithGaps",
            methodname="test_weighted_diurnal_debias_modeldata_gapfill_A",
        )

        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="label", show_gaps=False)
        fig = ax.get_figure()
        return fig

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
            methodname="test_debias_modeldata_gapfill",
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
    # test.test_import_data(overwrite_solution=True)
    # test.test_import_data(overwrite_solution=True)
    # test.test_dataset_timeseries_plotting_by_label()
    # test.test_dataset_timeseries_plotting_by_station()
    # test.test_station_timeseries_plotting_existing_ax()
    # test.test_station_timeseries_with_modeldata()
    # test.test_station_modeldata_timeseries()
    # test.test_dataset_modeldata_timeseries_plot()
    # test.test_dataset_color_by_station_and_modeldata_timeseries_plot()
    # test.test_dataset_color_by_label_and_modeldata_timeseries_plot()

    # test_gaps = TestDataWithGaps()
    # test_gaps.test_interpolated_timeseries_plot()
    # test_gaps.test_raw_modeldata_gf_timeseries_plot()
    # test_gaps.test_debias_modeldata_gf_timeseries_plot()
    # test_gaps.test_diurnal_debias_modeldata_gf_timeseries_plot()
    # test_gaps.test_weighted_diurnal_debias_modeldata_gf_timeseries_plot()
