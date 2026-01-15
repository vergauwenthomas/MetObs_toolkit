import pytest
import sys
from pathlib import Path
import copy
import tempfile

# import metobs_toolkit
import pandas as pd
import numpy as np


# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit
from metobs_toolkit.backend_collection import errorclasses as err

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer2, assert_equality, datadir


class TestDemoData:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testdemodata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    @pytest.fixture(scope="class")
    def import_dataset(self):
        """Import demo dataset fixture."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        return dataset

    def test_version(self):
        # check if the local version is used
        initpath = libfolder.joinpath(
            "src", "metobs_toolkit", "settings_collection", "version.py"
        )
        with open(initpath, "r") as f:
            content = f.read()
        version_line = [line for line in content.splitlines() if "__version__" in line][
            0
        ]
        local_version = version_line.split("=")[1].strip().strip('"').strip("'").strip()
        assert metobs_toolkit.__version__ == local_version
        assert isinstance(metobs_toolkit.__version__, str)

    @pytest.mark.dependency()
    def test_import_demo_data(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestDemoData.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            methodname=_method_name, **TestDemoData.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_calling_methods_without_solution_on_dataset(self, import_dataset):
        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # run methods, and see if something breaks

        # testing specials
        _ = dataset.stations
        _ = dataset.obstypes
        _ = dataset.template
        _ = dataset.df
        _ = dataset.outliersdf
        _ = dataset.metadf
        _ = dataset.stations
        _ = dataset.start_datetime
        _ = dataset.end_datetime
        _ = dataset.present_observations

        # get info's
        dataset.get_info(printout=False)
        dataset.template.get_info()

        # random collection
        dataset.rename_stations(
            renamedict={"vlinder01": "fakename", "vlider14": "fakename2"}
        )

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_calling_methods_without_solution_on_station(self, import_dataset):
        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # run methods, and see if something breaks
        station = dataset.stations[0]

        # testing specials

        _ = station.df
        _ = station.outliersdf
        _ = station.metadf
        _ = station.start_datetime
        _ = station.end_datetime
        _ = station.present_observations

        # get info's
        station.get_info()
        station.site.get_info()
        station.get_sensor("temp").get_info()

        # make plot
        station.make_plot()

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_subset_by_stations(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

        # 2. apply a metobs manipulation
        # Subset by valid stations
        data_to_test = dataset.subset_by_stations(
            stationnames=[" blabla", "vlinder01", "vlinder02"]
        )

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solution=data_to_test,
                **TestDemoData.solkwargs,
                methodname=_method_name,
            )

        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equality tests
        assert_equality(data_to_test, solutionobj)  # Dataset comparison

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_subset_by_stations_invalid(self, import_dataset, caplog):
        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        #  Test invalid input IDs
        with pytest.raises(ValueError):
            dataset.subset_by_stations(stationnames="vlinder01")

        #  Test invalid input IDs
        with pytest.raises(ValueError):
            dataset.subset_by_stations(stationnames=["vlinder01"])

        #  Test if a warning is logged for invalid station names
        with caplog.at_level("WARNING"):
            dataset.subset_by_stations(stationnames=["a", "b"])
        assert "No stations matched the provided station names" in caplog.text

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_get_info(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # 2. apply a metobs manipulation
        data_to_test = dataset.get_info(printout=False)

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solution=data_to_test,
                **TestDemoData.solkwargs,
                methodname=_method_name,
            )
        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests
        assert_equality(data_to_test, solutionobj)  # string comparison

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_get_station(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # 2. apply a metobs manipulation
        data_to_test = dataset.get_station("vlinder05")

        # 3. overwrite solution?
        if overwrite_solution:
            TestDemoData.solutionfixer.create_solution(
                solution=data_to_test,
                **TestDemoData.solkwargs,
                methodname=_method_name,
            )

        # 4. Get solution
        solutionobj = TestDemoData.solutionfixer.get_solution(
            **TestDemoData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests
        assert_equality(data_to_test, solutionobj)  # Station comparison

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_pickling_dataset(self, import_dataset):
        # 0. Get info of the current check
        # _method_name = sys._getframe().f_code.co_name #get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)

        # Create a tmp dir
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            trgfile = tmpdir.joinpath("deleteme.pkl")
            # pickle dataset
            dataset.save_dataset_to_pkl(filepath=trgfile)
            # Read in the pickled dataset
            dataset2 = metobs_toolkit.import_dataset_from_pkl(target_path=trgfile)

        # test if the pickled dataset is equal to the original
        assert_equality(dataset, dataset2)

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_dataset_to_parquet(self, import_dataset):
        """Test Dataset.to_parquet method"""
        # 1. get dataset data
        dataset = copy.deepcopy(import_dataset)
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            # Save to parquet
            parquet_file = tmpdir / "test_dataset.parquet"
            dataset.to_parquet(parquet_file)

            # Read back and compare
            df_original = dataset.df
            df_read = pd.read_parquet(parquet_file)

        # Test if dataframes are equal
        assert_equality(df_original, df_read)

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_dataset_to_csv(self, import_dataset):
        """Test Dataset.to_csv method"""
        # 1. get dataset data
        dataset = copy.deepcopy(import_dataset)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            # Save to CSV
            csv_file = tmpdir / "test_dataset.csv"
            dataset.to_csv(csv_file)

            # Read back and compare
            df_original = dataset.df
            df_read = pd.read_csv(csv_file, index_col=[0, 1, 2])  # Multi-index

            # ----Typecasting for compatibility ----
            # Convert datetime index level to datetime format to match original
            df_read.index = df_read.index.set_levels(
                pd.to_datetime(df_read.index.levels[0]), level=0
            )

            # Convert 'value' column to float32 to match original
            df_read["value"] = df_read["value"].astype("float32")

        # Test if dataframes are equal
        assert_equality(df_original, df_read)

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_station_to_parquet(self, import_dataset):
        """Test Station.to_parquet method"""
        # 1. get dataset data
        dataset = copy.deepcopy(import_dataset)

        # Get a station
        station = dataset.get_station("vlinder05")

        # Create a tmp dir
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            # Save to parquet
            parquet_file = tmpdir / "test_station.parquet"
            station.to_parquet(parquet_file)

            # Read back and compare
            df_original = station.df
            df_read = pd.read_parquet(parquet_file)

        # Test if dataframes are equal
        assert_equality(df_original, df_read)

    @pytest.mark.dependency(depends=["TestDemoData::test_import_demo_data"])
    def test_station_to_csv(self, import_dataset):
        """Test Station.to_csv method"""
        # 1. get dataset data
        dataset = copy.deepcopy(import_dataset)

        # Get a station
        station = dataset.get_station("vlinder05")

        # Create a tmp dir
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            # Save to CSV
            csv_file = tmpdir / "test_station.csv"
            station.to_csv(csv_file)

            # Read back and compare
            df_original = station.df
            df_read = pd.read_csv(csv_file, index_col=[0, 1])  # Multi-index

            # ----Typecasting for compatibility ----
            # Convert datetime index level to datetime format to match original
            df_read.index = df_read.index.set_levels(
                pd.to_datetime(df_read.index.levels[0]), level=0
            )

            # Convert 'value' column to float32 to match original
            df_read["value"] = df_read["value"].astype("float32")

        # Test if dataframes are equal
        assert_equality(df_original, df_read)

    def test_importing_data_with_nans_for_single_station(self):
        # goal is to test if metobs is able to import a datafile,
        # that has nans for a specific obstype for a specific station.

        df = pd.read_csv(metobs_toolkit.demo_datafile, sep=";")

        trgstation = "vlinder03"
        trg_column = "Vochtigheid"
        # all to Nan
        df.loc[df["Vlinder"] == trgstation, trg_column] = np.nan

        # to csv
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            targetfile = tmpdir / "data_with_nans.csv"
            df.to_csv(targetfile, index=False, sep=";")

            dataset = metobs_toolkit.Dataset()
            dataset.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                input_data_file=targetfile,
                input_metadata_file=metobs_toolkit.demo_metadatafile,
            )

        sta = dataset.get_station(trgstation)
        assert len(sta.obsdata) == 3


class TestWideData:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testwidedata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("wide_test_data.csv")
    templatefile = datadir.joinpath("wide_test_template.json")

    @pytest.fixture(scope="class")
    def import_wide_dataset(self):
        """Import wide format dataset fixture."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestWideData.templatefile,
            input_metadata_file=None,
            input_data_file=TestWideData.datafile,
            freq_estimation_method="median",
            freq_estimation_simplify_tolerance="2min",
            origin_simplify_tolerance="5min",
            timestamp_tolerance="4min",
        )
        return dataset

    @pytest.mark.dependency()
    def test_import_wide_data(self, import_wide_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_wide_dataset)

        # 3. overwrite solution?
        if overwrite_solution:
            TestWideData.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWideData.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestWideData.solutionfixer.get_solution(
            methodname=_method_name, **TestWideData.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    @pytest.mark.dependency(depends=["TestWideData::test_import_wide_data"])
    def test_sync_wide_records(self, import_wide_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_wide_dataset)
        # 2. apply a metobs manipulation
        dataset.sync_records(
            timestamp_shift_tolerance="5min2s", freq_shift_tolerance="2min"
        )

        # 3. overwrite solution?
        if overwrite_solution:
            TestWideData.solutionfixer.create_solution(
                solution=dataset, **TestWideData.solkwargs, methodname=_method_name
            )
        # 4. Get solution
        solutionobj = TestWideData.solutionfixer.get_solution(
            **TestWideData.solkwargs, methodname=_method_name
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison
        assert dataset.df.shape == (196, 2)


class TestWideSingleStationData:
    # to pass to the solutionfixer
    solkwargs = {
        "testfile": Path(__file__).name,
        "classname": "testwidesinglestationdata",
    }
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    # paths to data
    datafile = datadir.joinpath("single_station.csv")
    templatefile = datadir.joinpath("single_station_template.json")
    metadatfile = datadir.joinpath("single_station_metadata.csv")

    @pytest.fixture(scope="class")
    def import_single_station_dataset(self):
        """Import single station dataset fixture."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=TestWideSingleStationData.templatefile,
            input_metadata_file=TestWideSingleStationData.metadatfile,
            input_data_file=TestWideSingleStationData.datafile,
        )
        return dataset

    @pytest.mark.dependency()
    def test_import_wide_data(
        self, import_single_station_dataset, overwrite_solution=False
    ):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. get_startpoint data
        dataset = copy.deepcopy(import_single_station_dataset)

        # 3. overwrite solution?
        if overwrite_solution:
            TestWideSingleStationData.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestWideSingleStationData.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestWideSingleStationData.solutionfixer.get_solution(
            methodname=_method_name, **TestWideSingleStationData.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison


class TestStationAddMethods:
    @pytest.fixture(scope="class")
    def import_dataset(self):
        """Import demo dataset fixture."""
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        return dataset

    def test_station_add_sensordata_and_modeldata(self, import_dataset):
        # Get a dataset and a station
        dataset = copy.deepcopy(import_dataset)
        station = dataset.get_station("vlinder01")

        # --- Test add_to_sensordata ---
        orig_sensor = station.get_sensor("humidity")
        new_sensor = orig_sensor.copy(deep=True)
        # Change obstype name to avoid collision
        new_sensor.obstype.name = "new_obstype"
        # Change values to distinguish
        new_sensor.series = new_sensor.series + 100
        # Add new SensorData
        station.add_to_sensordata(new_sensor)
        assert "new_obstype" in station.sensordata

        # Try to add again without force_update, should raise
        import metobs_toolkit.backend_collection.errorclasses as err

        with pytest.raises(err.MetObsDataAlreadyPresent):
            station.add_to_sensordata(new_sensor)

        # Add with force_update, should succeed
        station.add_to_sensordata(new_sensor, force_update=True)


class TestParquetData:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testparquetdata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    def test_import_single_station(self):
        # 1. csv to parquet
        csv_file = datadir.joinpath("single_station.csv")
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            parquet_file = tmpdir / "single_station.parquet"
            df = pd.read_csv(csv_file, sep=",")
            df.to_parquet(parquet_file)

            # 2. metobs dataset from csv files
            dataset_a = metobs_toolkit.Dataset()
            dataset_a.import_data_from_file(
                template_file=datadir.joinpath("single_station_template.json"),
                # input_metadata_file=self.metadatfile,
                input_data_file=csv_file,
            )
            # 3. metobs dataset from parquet files
            dataset_parq = metobs_toolkit.Dataset()
            dataset_parq.import_data_from_file(
                template_file=datadir.joinpath("single_station_template.json"),
                # input_metadata_file=self.metadatfile,
                input_data_file=parquet_file,
            )

        assert_equality(dataset_a, dataset_parq)

    def test_import_demo_data(self):
        # 1. csv to parquet
        csv_file = Path(metobs_toolkit.demo_datafile)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pqfile = tmpdir / "demo_data.parquet"
            df = pd.read_csv(csv_file, sep=";")
            df.to_parquet(pqfile)

            # 2. metobs dataset from csv files
            dataset_a = metobs_toolkit.Dataset()
            dataset_a.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                # input_metadata_file=self.metadatfile,
                input_data_file=csv_file,
            )
            # 3. metobs dataset from parquet files
            dataset_parq = metobs_toolkit.Dataset()
            dataset_parq.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                # input_metadata_file=self.metadatfile,
                input_data_file=pqfile,
            )
        assert_equality(dataset_a, dataset_parq)

    def test_import_demo_metadata_only(self):
        # 1. csv to parquet
        csv_file = Path(metobs_toolkit.demo_metadatafile)

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pqfile = tmpdir / "demo_metadata.parquet"
            metadf = pd.read_csv(csv_file, sep=",")
            metadf.to_parquet(pqfile)

            # 2. metobs dataset from csv files
            dataset_a = metobs_toolkit.Dataset()
            dataset_a.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                input_metadata_file=metobs_toolkit.demo_metadatafile,
                # input_data_file=csv_file,
            )
            # 3. metobs dataset from parquet files
            dataset_parq = metobs_toolkit.Dataset()
            dataset_parq.import_data_from_file(
                template_file=metobs_toolkit.demo_template,
                input_metadata_file=pqfile,
                # input_data_file=parquet_file,
            )
        assert_equality(dataset_a, dataset_parq)

    def test_import_wide_data(self):
        # 1. csv to parquet
        datafile = datadir.joinpath("wide_test_data.csv")
        templatefile = datadir.joinpath("wide_test_template.json")

        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            pqfile = tmpdir / "demo_wide_data.parquet"
            df = pd.read_csv(datafile, sep=",")
            df.to_parquet(pqfile)

            # 2. metobs dataset from csv files
            dataset_a = metobs_toolkit.Dataset()
            dataset_a.import_data_from_file(
                template_file=templatefile,
                # input_metadata_file=self.metadatfile,
                input_data_file=datafile,
            )
            # 3. metobs dataset from parquet files
            dataset_parq = metobs_toolkit.Dataset()
            dataset_parq.import_data_from_file(
                template_file=templatefile,
                # input_metadata_file=self.metadatfile,
                input_data_file=pqfile,
            )

        assert_equality(dataset_a, dataset_parq)

    def test_import_with_timezone(self):
        inputdatafile = metobs_toolkit.demo_datafile
        # read the csv
        rawdf = pd.read_csv(inputdatafile, sep=";")
        rawdf = rawdf[:300]  # reduce storage
        # format datetime and something else
        rawdf["datetime"] = pd.to_datetime(
            rawdf["Datum"] + rawdf["Tijd (UTC)"], format="%Y-%m-%d%H:%M:%S"
        )
        rawdf.drop(columns=["Datum", "Tijd (UTC)"], inplace=True)

        rawdf_utc = copy.deepcopy(rawdf)
        rawdf_utc["datetime"] = rawdf_utc["datetime"].dt.tz_localize("UTC")

        # to parquet
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            target_parquetfile_utc = tmpdir.joinpath(
                "demo_data_with_timezone_utc.parquet"
            )
            rawdf_utc.to_parquet(target_parquetfile_utc)

            template_file = datadir.joinpath("demo_template_for_parquet.json")

            dataset = metobs_toolkit.Dataset()
            dataset.import_data_from_file(
                template_file=template_file, input_data_file=target_parquetfile_utc
            )

        # Now make sure the tz is mismatched between the parquet file and template
        rawdf_paris = copy.deepcopy(rawdf)
        rawdf_paris["datetime"] = rawdf_paris["datetime"].dt.tz_localize("Europe/Paris")

        # to parquet
        # to parquet
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            target_parquetfile_paris = tmpdir.joinpath(
                "demo_data_with_timezone_paris.parquet"
            )

            rawdf_paris.to_parquet(target_parquetfile_paris)

            dataset = metobs_toolkit.Dataset()

            with pytest.raises(err.MetObsTemplateError):
                dataset.import_data_from_file(
                    template_file=template_file,
                    input_data_file=target_parquetfile_paris,
                )


if __name__ == "__main__":
    print("To Overwrite the solutions, run: \n pytest test_importing.py -v")

    OVERWRITE_SOLUTION = True

    # TestDemoData
    demo_tester = TestDemoData()
    demo_dataset = demo_tester.import_dataset.__wrapped__(demo_tester)

    demo_tester.test_version()
    demo_tester.test_import_demo_data(
        demo_dataset, overwrite_solution=OVERWRITE_SOLUTION
    )
    demo_tester.test_calling_methods_without_solution_on_dataset(demo_dataset)
    demo_tester.test_calling_methods_without_solution_on_station(demo_dataset)
    demo_tester.test_subset_by_stations(
        demo_dataset, overwrite_solution=OVERWRITE_SOLUTION
    )
    # demo_tester.test_subset_by_stations_invalid(demo_dataset, caplog)  # requires caplog fixture
    demo_tester.test_get_info(demo_dataset, overwrite_solution=OVERWRITE_SOLUTION)
    demo_tester.test_get_station(demo_dataset, overwrite_solution=OVERWRITE_SOLUTION)
    demo_tester.test_pickling_dataset(demo_dataset)
    demo_tester.test_dataset_to_parquet(demo_dataset)
    demo_tester.test_dataset_to_csv(demo_dataset)
    demo_tester.test_station_to_parquet(demo_dataset)
    demo_tester.test_station_to_csv(demo_dataset)
    demo_tester.test_importing_data_with_nans_for_single_station()

    # TestWideData
    wide_data_tester = TestWideData()
    wide_dataset = wide_data_tester.import_wide_dataset.__wrapped__(wide_data_tester)

    wide_data_tester.test_import_wide_data(
        wide_dataset, overwrite_solution=OVERWRITE_SOLUTION
    )
    wide_data_tester.test_sync_wide_records(
        wide_dataset, overwrite_solution=OVERWRITE_SOLUTION
    )

    # TestWideSingleStationData
    single_station_tester = TestWideSingleStationData()
    single_station_dataset = (
        single_station_tester.import_single_station_dataset.__wrapped__(
            single_station_tester
        )
    )

    single_station_tester.test_import_wide_data(
        single_station_dataset, overwrite_solution=OVERWRITE_SOLUTION
    )

    # TestStationAddMethods
    station_add_tester = TestStationAddMethods()
    station_add_dataset = station_add_tester.import_dataset.__wrapped__(
        station_add_tester
    )

    station_add_tester.test_station_add_sensordata_and_modeldata(station_add_dataset)

    # TestParquetData
    parquet_tester = TestParquetData()
    parquet_tester.test_import_single_station()
    parquet_tester.test_import_demo_data()
    parquet_tester.test_import_demo_metadata_only()
    parquet_tester.test_import_wide_data()
    parquet_tester.test_import_with_timezone()
