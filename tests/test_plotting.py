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

    def test_show_outliers_false_behavior(self):
        """Test that show_outliers=False actually removes outliers from the plot data."""
        # Get a dataset and apply some quality control to create outliers
        dataset = TestDemoDataset.solutionfixer.get_solution(
            **TestDemoDataset.solkwargs, methodname="test_import_data"
        )

        # Apply a gross value check to create some outliers
        dataset.gross_value_check(
            target_obstype="temp", lower_threshold=-10, upper_threshold=30
        )

        # Check that we have some outliers
        outliers = dataset.outliersdf
        assert not outliers.empty, "Expected some outliers to be created"

        # Create plot data for testing
        plotdf = (
            dataset.df.xs("temp", level="obstype", drop_level=False)
            .reset_index()
            .set_index(["name", "obstype", "datetime"])
            .sort_index()
        )

        # Get available stations for color mapping
        stations = plotdf.index.get_level_values("name").unique().tolist()

        # Test colorby="label" behavior
        from metobs_toolkit.plot_collection import timeseries_plotting
        from metobs_toolkit.plot_collection import general_functions
        import matplotlib.pyplot as plt

        # Create axes for testing
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

        # Plot with show_outliers=True (should include outliers)
        ax1 = timeseries_plotting.plot_timeseries_color_by_label(
            plotdf=plotdf.copy(),
            show_outliers=True,
            show_gaps=True,
            ax=ax1,
        )
        general_functions.set_legend(ax1)

        # Plot with show_outliers=False (should exclude outliers)
        ax2 = timeseries_plotting.plot_timeseries_color_by_label(
            plotdf=plotdf.copy(),
            show_outliers=False,
            show_gaps=True,
            ax=ax2,
        )
        general_functions.set_legend(ax2)

        # Check that the legends are different (fewer items when outliers are hidden)
        legend1_labels = [t.get_text() for t in ax1.get_legend().get_texts()]
        legend2_labels = [t.get_text() for t in ax2.get_legend().get_texts()]

        # The plot with show_outliers=False should have fewer legend items
        assert len(legend2_labels) < len(
            legend1_labels
        ), "Expected fewer legend items when outliers are hidden"

        plt.close(fig)

        # Test colorby="station" behavior
        fig, (ax3, ax4) = plt.subplots(1, 2, figsize=(12, 5))

        # Create colormap for station plotting using actual station names
        colormap = {station: f"C{i}" for i, station in enumerate(stations)}

        # Plot with show_outliers=True (should include outliers)
        ax3 = timeseries_plotting.plot_timeseries_color_by_station(
            plotdf=plotdf.copy(),
            colormap=colormap,
            show_outliers=True,
            show_gaps=True,
            ax=ax3,
        )

        # Plot with show_outliers=False (should exclude outliers by setting them to NaN)
        ax4 = timeseries_plotting.plot_timeseries_color_by_station(
            plotdf=plotdf.copy(),
            colormap=colormap,
            show_outliers=False,
            show_gaps=True,
            ax=ax4,
        )

        # For station plotting, outliers should be set to NaN, so we should see fewer data points
        # Get all the lines from both plots and compare data point counts
        lines1 = ax3.get_lines()
        lines2 = ax4.get_lines()

        # Count non-NaN data points in each plot
        total_points_with_outliers = sum(
            len([x for x in line.get_ydata() if not pd.isna(x)]) for line in lines1
        )
        total_points_without_outliers = sum(
            len([x for x in line.get_ydata() if not pd.isna(x)]) for line in lines2
        )

        # The plot without outliers should have fewer data points
        assert (
            total_points_without_outliers < total_points_with_outliers
        ), f"Expected fewer data points when outliers are hidden. With: {total_points_with_outliers}, Without: {total_points_without_outliers}"

        plt.close(fig)


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
    test.test_dataset_timeseries_plotting_by_label()
    test.test_dataset_timeseries_plotting_by_station()
    test.test_station_timeseries_plotting_existing_ax()
    test.test_station_timeseries_with_modeldata()
    test.test_station_modeldata_timeseries()
    test.test_dataset_modeldata_timeseries_plot()
    test.test_dataset_color_by_station_and_modeldata_timeseries_plot()
    test.test_dataset_color_by_label_and_modeldata_timeseries_plot()

    test_gaps = TestDataWithGaps()
    test_gaps.test_interpolated_timeseries_plot()
    test_gaps.test_raw_modeldata_gf_timeseries_plot()
    test_gaps.test_debias_modeldata_gf_timeseries_plot()
    test_gaps.test_diurnal_debias_modeldata_gf_timeseries_plot()
    test_gaps.test_weighted_diurnal_debias_modeldata_gf_timeseries_plot()
