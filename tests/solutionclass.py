from pathlib import Path
import pickle

import pandas as pd
import pandas.testing

libfolder = Path(str(Path(__file__).resolve())).parent.parent
# testdatadir
datadir = libfolder.joinpath("tests").joinpath("data")


class SolutionFixer:
    def __init__(self, solutiondir):
        self.solutiondir = Path(solutiondir)

    def get_solution(self, testfile, classname, methodname):
        basefoldername = f"{str(testfile).strip('.py')}_solutions"
        subfoldername = f"{classname.lower()}"
        trgfilename = f"{str(methodname)}.pkl"

        solutionfile = (
            self.solutiondir.joinpath(basefoldername)
            .joinpath(subfoldername)
            .joinpath(trgfilename)
        )
        if not solutionfile.exists():
            raise SolutionNotExisting(
                f"The solution of {testfile} -->{classname}:{methodname} does not exist (at {solutionfile})!"
            )

        with open(solutionfile, "rb") as file:
            solution_data = pickle.load(file)
        return solution_data

    def create_solution(self, solutiondata, testfile, classname, methodname):
        # construct path
        basefoldername = f"{str(testfile).strip('.py')}_solutions"
        subfoldername = f"{classname.lower()}"
        trgfilename = f"{str(methodname)}.pkl"

        solutionfile = (
            self.solutiondir.joinpath(basefoldername)
            .joinpath(subfoldername)
            .joinpath(trgfilename)
        )
        # clear previous solutions
        if solutionfile.exists():
            print(
                f"!! OVERWRITING SOLUTION FOR  {testfile} --> {classname}:{methodname} !!! "
            )
            # delete file
            solutionfile.unlink()

        # Create directory if it does not exist
        solutionfile.parent.mkdir(parents=True, exist_ok=True)

        # Pickle data object
        with open(solutionfile, "wb") as file:
            pickle.dump(solutiondata, file)


def assert_equality(to_check, solution):
    """Returns some debug help when the to_check is not equalt to the solution

    This output is used for developping to point in more details why a test
    is failing.

    The object retured can have differnt types, and is desined for the developper
    to help in the debugging process.

    """

    # type equal test
    if type(to_check) != type(solution):
        retstr = f"DIFF: to_check type is {type(to_check)}, while solution is of type {type(solution)} "
        raise AssertionError(retstr)

    # float test
    elif isinstance(to_check, float):
        if to_check != solution:
            retstr = f"DIFF: to_check is not equal to solution {to_check} !== {solution} (float comparison)"
            raise AssertionError(retstr)

    # int test
    elif isinstance(to_check, int):
        if to_check != solution:
            retstr = f"DIFF: to_check is not equal to solution {to_check} !== {solution} (int comparison)"
            raise AssertionError(retstr)

    # sting test
    elif isinstance(to_check, str):
        if to_check != solution:
            retstr = f"DIFF: to_check is not equal to solution (string comparison)"
            retstr += f"to check: \n =================== \n{to_check}"
            retstr += f"solution: \n =================== \n{solution}"
            raise AssertionError(retstr)

    # tuple comparison
    elif isinstance(to_check, tuple):
        if to_check != solution:
            retstr = f"DIFF: to_check is not equal to solution {to_check} !== {solution} (tuple comparison)"
            raise AssertionError(retstr)

    # list test
    elif isinstance(to_check, list):
        if to_check != solution:
            # length
            if len(to_check) != len(solution):
                retstr = f"DIFF: to_check length is {len(to_check)} !== {len(solution)} (list comparison)"
                raise AssertionError(retstr)
            # set check
            if set(to_check) == set(solution):
                retstr = f"DIFF: to_check list {len(to_check)}!== {len(solution)} (list comparison), BUT if converted to sets they are identical !! "
                raise AssertionError(retstr)
            else:
                retstr = f"DIFF: to_check list {len(to_check)}!== {len(solution)} (list comparison), no hints found, debug further."
                raise AssertionError(retstr)
    # set test
    elif isinstance(to_check, set):
        if to_check != solution:
            # length
            if len(to_check) != len(solution):
                retstr = f"DIFF: to_check length is {len(to_check)} !== {len(solution)} (set comparison)"
                raise AssertionError(retstr)
            else:
                retstr = f"DIFF: to_check list {len(to_check)}!== {len(solution)} (set comparison), no hints found, debug further."
                raise AssertionError(retstr)

    # pandas seriers
    elif isinstance(to_check, pd.Series):
        pd.testing.assert_series_equal(
            left=to_check,
            right=solution,
            check_exact=False,
        )

    # pandas dataframes
    elif isinstance(to_check, pd.DataFrame):
        pd.testing.assert_frame_equal(
            left=to_check,
            right=solution,
            check_exact=False,
        )

    # metobs_toolkit.Dataset test
    elif to_check.__class__.__name__ == "Dataset":
        compare_df_attr(to_check, solution, "metadf")
        compare_df_attr(to_check, solution, "gapsdf")
        compare_df_attr(to_check, solution, "modeldatadf")
        compare_df_attr(to_check, solution, "outliersdf")
        compare_df_attr(to_check, solution, "df")
        assert (
            to_check.obstypes == solution.obstypes
        ), "There is a mismatch in obstypes with the solution!"

    # metobs_toolkit.Station test
    elif to_check.__class__.__name__ == "Station":
        compare_df_attr(to_check, solution, "metadf")
        compare_df_attr(to_check, solution, "gapsdf")
        compare_df_attr(to_check, solution, "modeldatadf")
        compare_df_attr(to_check, solution, "outliersdf")
        compare_df_attr(to_check, solution, "df")
    # metobs_toolkit.Station test
    elif to_check.__class__.__name__ == "Analysis":
        compare_df_attr(to_check, solution, "metadf")
        compare_df_attr(to_check, solution, "df")

    # Else
    else:
        retstr = f"DIFF: to_check list {to_check}!== {solution} (NotImplemented type comparison), no hints found, debug further."
        raise AssertionError(retstr)


def compare_df_attr(testobj, solutionobj, attr):
    try:
        pd.testing.assert_frame_equal(
            left=getattr(testobj, attr),
            right=getattr(solutionobj, attr),
            check_exact=False,
        )
    except AssertionError as e:
        raise AssertionError(f"DIFF in {attr}-attribute:\n " + str(e))


class UnforseenDifference(Exception):
    """Raise when encountering an unforseen difference case"""


class SolutionNotExisting(Exception):
    """Raise when the solutionfile does not exist"""
