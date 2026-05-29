"""This file test the Obstype class and its methods using the pytest framework."""

import pytest
import sys
import copy
from pathlib import Path

# import metobs_toolkit
import pandas as pd

# Get the path to the MetObs_toolkit directory
# Add the local source directory to Python path for development
libfolder = Path(str(Path(__file__).resolve())).parent.parent
sys.path.insert(0, str(libfolder / "src"))
import metobs_toolkit

# solutionfolder
solutionsdir = libfolder.joinpath("tests").joinpath("pkled_solutions")
from solutionclass import SolutionFixer2, assert_equality, datadir

import pint


class TestObstype:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testobstypedata"}
    solutionfixer = SolutionFixer2(solutiondir=solutionsdir)

    @pytest.fixture(scope="class")
    def import_dataset(self):
        # 2. apply a metobs manipulation
        dataset = metobs_toolkit.Dataset()
        dataset.import_data_from_file(
            template_file=metobs_toolkit.demo_template,
            input_metadata_file=metobs_toolkit.demo_metadatafile,
            input_data_file=metobs_toolkit.demo_datafile,
        )
        return dataset

    @pytest.mark.dependency()
    def test_import_demo_data(self, import_dataset, overwrite_solution=False):
        # 0. Get info of the current check
        _method_name = sys._getframe().f_code.co_name  # get the name of this method

        dataset = copy.deepcopy(import_dataset)

        # 3. overwrite solution?
        if overwrite_solution:
            TestObstype.solutionfixer.create_solution(
                solution=dataset,
                methodname=_method_name,
                **TestObstype.solkwargs,
            )

        # 4. Get solution
        solutionobj = TestObstype.solutionfixer.get_solution(
            methodname=_method_name, **TestObstype.solkwargs
        )

        # 5. Construct the equlity tests
        assert_equality(dataset, solutionobj)  # dataset comparison

    @pytest.mark.dependency(depends=["TestObstype::test_import_demo_data"])
    def test_calling_methods_without_solution_on_obstypes(self, import_dataset):
        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # run methods, and see if something breaks
        temp = dataset.obstypes["temp"]
        # testing specials
        _ = temp.name
        _ = temp.std_unit
        _ = temp.description
        _ = temp.original_name
        _ = temp.original_unit

        # test setting methods
        temp.original_name = "dummy orig name"

        temp.description = "dummy description"

        # get info's
        _ = temp.get_info(printout=False)

        # Method calls on modelobstypes
        era5_manager = metobs_toolkit.default_GEE_datasets["ERA5-land"]
        _ = era5_manager.modelobstypes["temp"].get_info(printout=False)

        # Method calls on modelobstypes_vectorfields
        _ = era5_manager.modelobstypes["wind"].get_info(printout=False)

    @pytest.mark.dependency(depends=["TestObstype::test_import_demo_data"])
    def test_units_io(self, import_dataset):
        # 1. get_startpoint data
        dataset = copy.deepcopy(import_dataset)
        # run methods, and see if something breaks
        temp = dataset.obstypes["temp"]

        # Test if error is trown for unknow units
        from metobs_toolkit.obstypes import MetObsUnitUnknown

        with pytest.raises(MetObsUnitUnknown):
            temp.original_unit = "dummy orig unit"

        # test str return types
        assert isinstance(temp.std_unit, str)
        assert isinstance(temp.original_unit, str)
        temp.description = None
        assert isinstance(temp.description, str)

    @pytest.mark.dependency(depends=["TestObstype::test_import_demo_data"])
    def test_creating_new_obstypes(self, import_dataset):
        # 1. Create a new Obstype instance
        new_obstype = metobs_toolkit.obstypes.Obstype(
            name="new_temp", std_unit="°F", description="New temperature observation"
        )

        # Test if error is trown for incompatible units
        from metobs_toolkit.obstypes import MetObsUnitsIncompatible

        with pytest.raises(MetObsUnitsIncompatible):
            new_obstype.original_unit = "kilometer/hour"

        assert (
            len(new_obstype.get_compatible_units()) == 7
        )  # 7 base units (Pint) for temperatures

        new_obstype.original_name = 168
        assert new_obstype.original_name == "168"

        # test units conversion
        new_obstype.original_unit = "°C"
        origdata = pd.Series([60, 20, 30])  # in °C

        in_farenheit = new_obstype.convert_to_standard_units(
            input_data=origdata, input_unit=new_obstype.original_unit
        )
        assert (in_farenheit.astype(int) == pd.Series([139, 67, 85])).all()

        # Test adding units to the dataset
        dataset = copy.deepcopy(import_dataset)

        dataset.add_new_observationtype(obstype=new_obstype)
        assert "new_temp" in dataset.obstypes.keys()

        # Test if an error is trown when the obstype is already present
        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsDataAlreadyPresent,
        )

        with pytest.raises(MetObsDataAlreadyPresent):
            dataset.add_new_observationtype(obstype=dataset.obstypes["humidity"])

    def test_calling_methods_on_modelobstypes(self):
        # 1. get_startpoint data
        modeltemp = metobs_toolkit.default_GEE_datasets["ERA5-land"].modelobstypes[
            "temp"
        ]

        # testing specials
        _ = modeltemp.name
        _ = modeltemp.std_unit
        _ = modeltemp.description
        _ = modeltemp.original_name
        _ = modeltemp.original_unit
        # specific to modelobstype
        _ = modeltemp.model_band
        _ = modeltemp.model_unit
        _ = modeltemp.get_info(printout=False)

    def test_model_unit_setter(self):
        """Test the model_unit setter of ModelObstype class."""
        # Create a base temperature obstype
        temp_obstype = metobs_toolkit.obstypes.Obstype(
            name="temp", std_unit="°C", description="Temperature"
        )

        # Create a ModelObstype with initial model_unit
        model_obstype = metobs_toolkit.obstypes.ModelObstype(
            obstype=temp_obstype, model_unit="kelvin", model_band="temperature_2m"
        )

        # Test initial value
        assert model_obstype.model_unit == "kelvin"

        # Test setting compatible unit (string)
        model_obstype.model_unit = "fahrenheit"
        assert model_obstype.model_unit == "degree_Fahrenheit"

        # Test setting compatible unit (pint.Unit)
        ureg = pint.UnitRegistry()
        model_obstype.model_unit = ureg.degC
        assert model_obstype.model_unit == "degree_Celsius"

        # Test setting incompatible unit should raise error
        from metobs_toolkit.backend_collection.errorclasses import (
            MetObsUnitsIncompatible,
        )

        with pytest.raises(MetObsUnitsIncompatible):
            model_obstype.model_unit = (
                "meter"  # Length unit incompatible with temperature
            )

    def test_model_band_setter(self):
        """Test the model_band setter of ModelObstype class."""
        # Create a base temperature obstype
        temp_obstype = metobs_toolkit.obstypes.Obstype(
            name="temp", std_unit="°C", description="Temperature"
        )

        # Create a ModelObstype with initial model_band
        model_obstype = metobs_toolkit.obstypes.ModelObstype(
            obstype=temp_obstype, model_unit="kelvin", model_band="temperature_2m"
        )

        # Test initial value
        assert model_obstype.model_band == "temperature_2m"

        # Test setting new string value
        model_obstype.model_band = "skin_temperature"
        assert model_obstype.model_band == "skin_temperature"


if __name__ == "__main__":
    # pytest.main([__file__])
    OVERWRITE_SOLUTION = False

    testobstype = TestObstype()

    demodata = testobstype.import_dataset.__wrapped__(testobstype)
    testobstype.test_import_demo_data(
        overwrite_solution=OVERWRITE_SOLUTION, import_dataset=demodata
    )
    testobstype.test_calling_methods_without_solution_on_obstypes(
        import_dataset=demodata
    )
    testobstype.test_units_io(import_dataset=demodata)
    testobstype.test_creating_new_obstypes(import_dataset=demodata)
    testobstype.test_calling_methods_on_modelobstypes()
    testobstype.test_model_unit_setter()
    testobstype.test_model_band_setter()
