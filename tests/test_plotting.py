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
from solutionclass import SolutionFixer2, assert_equality, datadir
import shutil


print(
    "To Overwrite the solutions, run: \n pytest test_plotting.py --mpl-generate-path=baseline"
)


# Fixture to ensure matplotlib figures are cleaned up after each test
@pytest.fixture(autouse=True)
def cleanup_figures():
    """Close all matplotlib figures before and after each test to prevent state leakage."""
    # Close any figures that might exist before the test
    plt.close("all")
    yield
    # Close all figures created during the test
    plt.close("all")


class TestDemoDataset:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)
    
    @pytest.fixture(scope='class')
    def import_dataset(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        # To hourly !!
        dataset.resample(target_freq="1h")
        return dataset
    
    @pytest.fixture(autouse=True)
    def import_dataset_with_era5(self):
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_data_file=metobs_toolkit.demo_datafile,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
        )
        dataset.resample(target_freq="15min")

        era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]

        startdt_utc = pd.Timestamp("2022-08-31 18:32:25")
        enddt_utc = pd.Timestamp("2022-09-01 12:16:00")
        era5_data = dataset.get_gee_timeseries_data(
            gee_dynamic_manager=era5_model,
            startdt_utc=startdt_utc,
            enddt_utc=enddt_utc,
            obstypes=["temp"],
            get_all_bands=False,
            drive_filename=None,
            drive_folder="gee_timeseries_data",
            force_direct_transfer=True,
            force_to_drive=False,
        )
        
        assert era5_data.shape == (532, 1)
        return dataset

    @pytest.mark.dependency()
    def test_import_data(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison
    
    @pytest.mark.dependency()
    def test_import_data_with_era5(self, import_dataset_with_era5, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset_with_era5)
        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoDataset.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestDemoDataset.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoDataset.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoDataset.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison
        
    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data"])
    @pytest.mark.mpl_image_compare
    def test_dataset_timeseries_plotting_by_label(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="label", figkwargs={"figsize": (10, 5)})
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data"])
    @pytest.mark.mpl_image_compare
    def test_dataset_timeseries_plotting_by_station(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

        # 2. apply a metobs manipulation
        ax = dataset.make_plot(colorby="station", figkwargs={"figsize": (10, 5)})
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data"])
    @pytest.mark.mpl_image_compare
    def test_dataset_test_show_outliers_labelby_station(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # 2. apply a metobs manipulation
        dataset.repetitions_check(max_N_repetitions=8)
        ax = dataset.make_plot(colorby="station", show_outliers=False)
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data"])
    @pytest.mark.mpl_image_compare
    def test_dataset_test_show_outliers_labelby_labels(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

        # 2. apply a metobs manipulation
        dataset.repetitions_check(max_N_repetitions=8)
        ax = dataset.make_plot(colorby="label", show_outliers=False)
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data"])
    @pytest.mark.mpl_image_compare
    def test_station_timeseries_plotting_existing_ax(self, import_dataset):
        #  1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        sta = dataset.get_station("vlinder02")
        fig, ax = plt.subplots(figsize=(8, 4))

        # 2. apply a metobs manipulation
        ax = sta.make_plot(colorby="label", show_outliers=False, ax=ax)
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_station_timeseries_with_modeldata(self, import_dataset_with_era5):
        #  1. get_startpoint data
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot(show_modeldata=True)
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_station_modeldata_timeseries(self, import_dataset_with_era5):
        #  1. get_startpoint data
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot_of_modeldata(obstype="temp")
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_station_plot_of_modeldata_with_modelname(self, import_dataset_with_era5):
        """Test station.make_plot_of_modeldata with modelname argument."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        station = dataset_with_era.get_station("vlinder05")
        ax = station.make_plot_of_modeldata(obstype="temp", modelname="ERA5-land")
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_dataset_plot_of_modeldata_with_modelname(self, import_dataset_with_era5):
        """Test Dataset.make_plot_of_modeldata with modelname argument."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        ax = dataset_with_era.make_plot_of_modeldata(
            obstype="temp", modelname="ERA5-land"
        )
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig
    
    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_station_plot_humidity_with_temp_modeldata(self, import_dataset_with_era5):
        """Test station.make_plot with humidity obstype and temp modeldata."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

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
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_station_plot_temp_with_modeldata_kwargs(self, import_dataset_with_era5):
        """Test station.make_plot with temp obstype and modeldata kwargs."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)


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
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_dataset_plot_humidity_with_modelvariable(self, import_dataset_with_era5):
        """Test dataset.make_plot with humidity obstype and modelvariable."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        ax = dataset_with_era.make_plot(
            obstype="humidity",
            show_modeldata=True,
            modelobstype="temp",
            modeldata_kwargs={
                "modelvariable": "temperature_2m",
            },
        )
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_modeldatatimeseries_timeseries(self, import_dataset_with_era5):
        #  1. get_startpoint data
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        modelseries = dataset_with_era.get_station("vlinder05").get_modeltimeseries(
            "temp"
        )
        ax = modelseries.make_plot()
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_dataset_modeldata_timeseries_plot(self, import_dataset_with_era5):
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        ax = dataset_with_era.make_plot_of_modeldata()
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_dataset_color_by_station_and_modeldata_timeseries_plot(self, import_dataset_with_era5):
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        ax = dataset_with_era.make_plot(colorby="station", show_modeldata=True)
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_dataset_color_by_label_and_modeldata_timeseries_plot(self, import_dataset_with_era5):
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        ax = dataset_with_era.make_plot(colorby="label", show_modeldata=True)
        fig = ax.get_figure()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_sensordata_pd_plot(self, import_dataset_with_era5):
        """Test SensorData.pd_plot() method with specific pandas plot kwargs."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        station = dataset_with_era.get_station("vlinder01")
        sensordata = station.get_sensor("temp")

        # Test with specific pandas plot kwargs
        ax = sensordata.pd_plot(
            show_labels=["ok"],
            color="red",
            linewidth=2,
            linestyle="-.",
            alpha=0.8,
            figsize=(10, 6),
        )
        ax.legend()
        fig = plt.gcf()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_sensordata_pd_plot_with_filters(self, import_dataset_with_era5):
        """Test SensorData.pd_plot() method with specific pandas plot kwargs."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        station = dataset_with_era.get_station("vlinder01")
        station.repetitions_check()
        sensordata = station.get_sensor("temp")

        # Test with specific pandas plot kwargs
        ax = sensordata.pd_plot(
            show_labels=["ok"],  # outliers are not present
            color="red",
            linewidth=2,
            linestyle="-.",
            alpha=0.8,
            figsize=(10, 6),
        )
        ax.legend()
        fig = plt.gcf()
        fig.set_size_inches(15, 5) 
        return fig

    @pytest.mark.dependency(depends=["TestDemoDataset::test_import_data_with_era5"])
    @pytest.mark.mpl_image_compare
    def test_modeltimeseries_pd_plot(self, import_dataset_with_era5):
        """Test ModelTimeSeries.pd_plot() method with specific pandas plot kwargs."""
        dataset_with_era = copy.deepcopy(import_dataset_with_era5)

        station = dataset_with_era.get_station("vlinder05")
        modeltimeseries = station.get_modeltimeseries("temp")

        # Test with specific pandas plot kwargs
        ax = modeltimeseries.pd_plot(
            color="blue",
            linewidth=3,
            linestyle="--",
            marker="o",
            markersize=4,
            alpha=0.7,
            figsize=(12, 8),
        )
        ax.legend()
        fig = plt.gcf()
        fig.set_size_inches(15, 5) 
        return fig



if __name__ == "__main__":
    print(
        "To Overwrite the solutions, run: \n pytest test_plotting.py --mpl --mpl-generate-path=baseline"
    )

    print(
        "To checkout the differences, run: \n pytest test_plotting.py --mpl --mpl-generate-summary=html"
    )

    OVERWRITE_SOLUTION = False

    tester = TestDemoDataset()
    
    # Prepare fixtures
    import_dataset = tester.import_dataset.__wrapped__(tester)
    import_dataset_with_era5 = tester.import_dataset_with_era5.__wrapped__(tester)
    
    # Run tests with solutions
    tester.test_import_data(import_dataset, overwrite_solution=OVERWRITE_SOLUTION)
    tester.test_import_data_with_era5(import_dataset_with_era5, overwrite_solution=OVERWRITE_SOLUTION)
    
    # # Run plotting tests (using import_dataset)
    # tester.test_dataset_timeseries_plotting_by_label(import_dataset)
    # tester.test_dataset_timeseries_plotting_by_station(import_dataset)
    # tester.test_dataset_test_show_outliers_labelby_station(import_dataset)
    # tester.test_dataset_test_show_outliers_labelby_labels(import_dataset)
    # tester.test_station_timeseries_plotting_existing_ax(import_dataset)
    
    # # Run plotting tests (using import_dataset_with_era5)
    # tester.test_station_timeseries_with_modeldata(import_dataset_with_era5)
    # tester.test_station_modeldata_timeseries(import_dataset_with_era5)
    # tester.test_station_plot_of_modeldata_with_modelname(import_dataset_with_era5)
    # tester.test_dataset_plot_of_modeldata_with_modelname(import_dataset_with_era5)
    # tester.test_station_plot_humidity_with_temp_modeldata(import_dataset_with_era5)
    # tester.test_station_plot_temp_with_modeldata_kwargs(import_dataset_with_era5)
    # tester.test_dataset_plot_humidity_with_modelvariable(import_dataset_with_era5)
    # tester.test_modeldatatimeseries_timeseries(import_dataset_with_era5)
    # tester.test_dataset_modeldata_timeseries_plot(import_dataset_with_era5)
    # tester.test_dataset_color_by_station_and_modeldata_timeseries_plot(import_dataset_with_era5)
    # tester.test_dataset_color_by_label_and_modeldata_timeseries_plot(import_dataset_with_era5)
    # tester.test_sensordata_pd_plot(import_dataset_with_era5)
    # tester.test_sensordata_pd_plot_with_filters(import_dataset_with_era5)
    # tester.test_modeltimeseries_pd_plot(import_dataset_with_era5)
