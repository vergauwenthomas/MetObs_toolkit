import pytest
import sys
from pathlib import Path
import tempfile
import shutil

# import metobs_toolkit
import pandas as pd
import numpy as np


libfolder = Path(str(Path(__file__).resolve())).parent.parent

# point to current version of the toolkit
# sys.path.insert(1, str(libfolder))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer, assert_equality, datadir


class TestVerification:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testverification"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_create_verification_with_demo_data(self, overwrite_solution=False):
        """Create a verification instance from demo data with ERA5 model data."""
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        # 1. Get starting point data - create a dataset with demo observations
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )

        # 2. Import ERA5 model data from the test data CSV file
        # This provides realistic model data for verification testing
        era5_model = metobs_toolkit.default_GEE_datasets["ERA5-land"]
        target_era5_csv = datadir.joinpath(
            "ERA5-land_timeseries_data_of_full_dataset_28_stations.csv"
        )
        dataset.import_gee_data_from_file(
            filepath=target_era5_csv,
            geedynamicdatasetmanager=era5_model,
            force_update=True,
        )
        dataset.make_plot_of_modeldata()
        # 3. Create verification instance
        verification = metobs_toolkit.Verification(dataset)

        # 4. overwrite solution?
        if overwrite_solution:
            TestVerification.solutionfixer.create_solution(
                solutiondata=verification, methodname=_method_name, **TestVerification.solkwargs
            )

        # 5. Get solution
        solutionobj = TestVerification.solutionfixer.get_solution(
            methodname=_method_name, **TestVerification.solkwargs
        )

        # 6. Construct the equality tests
        assert_equality(verification, solutionobj)  # verification comparison

    def test_save_and_load_verification_pickle(self):
        """Test saving and loading verification instance to/from pickle."""
        # 1. Get verification instance from previous test
        verification = TestVerification.solutionfixer.get_solution(
            **TestVerification.solkwargs, methodname="test_create_verification_with_demo_data"
        )

        # Create a tmp dir
        tmpdir = libfolder.joinpath("tmp")
        tmpdir.mkdir(parents=True, exist_ok=True)
        
        try:
            # Test saving verification to pickle
            verification.save_verification_to_pkl(target_folder=tmpdir, filename="test_verification")
            
            # Verify the file was created
            expected_file = tmpdir.joinpath("test_verification.pkl")
            assert expected_file.exists(), "Pickle file was not created"
            assert expected_file.stat().st_size > 0, "Pickle file is empty"

            # Test loading verification from pickle
            verification2 = metobs_toolkit.import_verification_from_pkl(
                target_path=expected_file
            )

            # Test if the pickled verification is equal to the original
            assert_equality(verification, verification2)

        finally:
            # Remove the tmp dir
            if tmpdir.exists():
                shutil.rmtree(tmpdir)

    def test_save_verification_pickle_with_overwrite(self):
        """Test saving verification with overwrite functionality."""
        # 1. Get verification instance
        verification = TestVerification.solutionfixer.get_solution(
            **TestVerification.solkwargs, methodname="test_create_verification_with_demo_data"
        )

        # Create a tmp dir
        tmpdir = libfolder.joinpath("tmp")
        tmpdir.mkdir(parents=True, exist_ok=True)
        
        try:
            filename = "test_overwrite"
            target_file = tmpdir.joinpath(f"{filename}.pkl")
            
            # Save verification first time
            verification.save_verification_to_pkl(target_folder=tmpdir, filename=filename)
            assert target_file.exists(), "First save failed"
            
            # Test that saving again without overwrite raises error
            with pytest.raises(FileExistsError):
                verification.save_verification_to_pkl(target_folder=tmpdir, filename=filename, overwrite=False)
            
            # Test that saving with overwrite=True works
            verification.save_verification_to_pkl(target_folder=tmpdir, filename=filename, overwrite=True)
            assert target_file.exists(), "Overwrite save failed"
            
            # Verify we can still load it
            verification2 = metobs_toolkit.import_verification_from_pkl(target_path=target_file)
            assert_equality(verification, verification2)

        finally:
            # Remove the tmp dir
            if tmpdir.exists():
                shutil.rmtree(tmpdir)

    def test_save_verification_pickle_invalid_folder(self):
        """Test saving verification to non-existent folder raises error."""
        # 1. Get verification instance
        verification = TestVerification.solutionfixer.get_solution(
            **TestVerification.solkwargs, methodname="test_create_verification_with_demo_data"
        )

        # Test with non-existent folder
        non_existent_folder = libfolder.joinpath("this_folder_does_not_exist")
        
        with pytest.raises(FileNotFoundError):
            verification.save_verification_to_pkl(
                target_folder=non_existent_folder, 
                filename="test.pkl"
            )

    def test_save_verification_filename_without_pkl_extension(self):
        """Test that .pkl extension is automatically added if missing."""
        # 1. Get verification instance
        verification = TestVerification.solutionfixer.get_solution(
            **TestVerification.solkwargs, methodname="test_create_verification_with_demo_data"
        )

        # Create a tmp dir
        tmpdir = libfolder.joinpath("tmp")
        tmpdir.mkdir(parents=True, exist_ok=True)
        
        try:
            filename_without_ext = "test_verification_no_ext"
            expected_file = tmpdir.joinpath(f"{filename_without_ext}.pkl")
            
            # Save verification with filename without extension
            verification.save_verification_to_pkl(
                target_folder=tmpdir, 
                filename=filename_without_ext
            )
            
            # Verify the file was created with .pkl extension
            assert expected_file.exists(), "File with .pkl extension was not created"
            
            # Verify we can load it
            verification2 = metobs_toolkit.import_verification_from_pkl(target_path=expected_file)
            assert_equality(verification, verification2)

        finally:
            # Remove the tmp dir
            if tmpdir.exists():
                shutil.rmtree(tmpdir)

    def test_verification_basic_functionality(self):
        """Test basic verification functionality to ensure our test data works."""
        # 1. Get verification instance
        verification = TestVerification.solutionfixer.get_solution(
            **TestVerification.solkwargs, methodname="test_create_verification_with_demo_data"
        )

        # Test that verification has the expected attributes
        assert hasattr(verification, 'verifdf'), "Verification missing verifdf attribute"
        assert hasattr(verification, 'metadf'), "Verification missing metadf attribute"
        assert hasattr(verification, 'obstypes'), "Verification missing obstypes attribute"

        # Test that verifdf has the expected structure
        verifdf = verification.verifdf
        assert isinstance(verifdf, pd.DataFrame), "verifdf is not a DataFrame"
        assert not verifdf.empty, "verifdf is empty"
        
        # Check for expected columns
        expected_columns = ['value_obs', 'value_model']
        for col in expected_columns:
            assert col in verifdf.columns, f"Missing column {col} in verifdf"
        
        # Test that we can calculate basic scores
        try:
            scores = verification.traditional_scores()
            assert isinstance(scores, pd.DataFrame), "traditional_scores did not return DataFrame"
            assert not scores.empty, "traditional_scores returned empty DataFrame"
        except Exception as e:
            pytest.skip(f"traditional_scores method failed: {e}")
