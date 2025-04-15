""" This file test the Obstype class and its methods using the pytest framework."""

import pytest
import sys
from pathlib import Path

# import metobs_toolkit
import pandas as pd

# Get the path to the MetObs_toolkit directory

libfolder = Path(str(Path(__file__).resolve())).parent.parent

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


class TestObstype:
    # to pass to the solutionfixer
    solkwargs = {"testfile": Path(__file__).name, "classname": "testobstypedata"}
    solutionfixer = SolutionFixer(solutiondir=solutionsdir)

    def test_import_demo_data(self, overwrite_solution=False):
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
        data_to_test = dataset

        # 3. overwrite solution?
        if overwrite_solution:
            TestObstype.solutionfixer.create_solution(
                solutiondata=data_to_test,
                methodname=_method_name,
                **TestObstype.solkwargs
            )

        # 4. Get solution
        solutionobj = TestObstype.solutionfixer.get_solution(
            methodname=_method_name, **TestObstype.solkwargs
        )

        # 5. Construct the equlity tests
        test_expr = data_to_test == solutionobj  # dataset comparison

        # 5. save comparison, create difference (only used when debugging, so no termina output)
        if not test_expr:
            debug_diff = TestObstype.solutionfixer.create_a_diff(
                to_check=data_to_test, solution=solutionobj
            )
        # 6. assert the equality
        assert test_expr

    def test_calling_methods_without_solution_on_obstypes(self):
        # 1. get_startpoint data
        dataset = TestObstype.solutionfixer.get_solution(
            **TestObstype.solkwargs, methodname="test_import_demo_data"
        )
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
        temp.get_info(printout=False)

    def test_units_io(self):
        # 1. get_startpoint data
        dataset = TestObstype.solutionfixer.get_solution(
            **TestObstype.solkwargs, methodname="test_import_demo_data"
        )
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

    def test_creating_new_obstypes(self):

        # 1. Create a new Obstype instance
        new_obstype = metobs_toolkit.obstypes.Obstype(
            obsname="new_temp", std_unit="°F", description="New temperature observation"
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
        dataset = TestObstype.solutionfixer.get_solution(
            **TestObstype.solkwargs, methodname="test_import_demo_data"
        )

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


if __name__ == "__main__":
    # pytest.main([__file__])

    testobstype = TestObstype()
    testobstype.test_import_demo_data(overwrite_solution=False)
