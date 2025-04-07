from pathlib import Path
import pickle

import pandas as pd


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
                f"The solution of {testfile} -->{classname}:{methodname} does not exist!"
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

    def create_a_diff(self, to_check, solution):
        """Returns some debug help when the to_check is not equalt to the solution

        This output is used for developping to point in more details why a test
        is failing.

        The object retured can have differnt types, and is desined for the developper
        to help in the debugging process.

        """

        # type equal test
        if type(to_check) != type(solution):
            retstr = f"DIFF: to_check type is {type(to_check)}, while solution is of type {type(solution)} "
            print(retstr)
            return retstr

        # float test
        elif isinstance(to_check, float):
            retstr = f"DIFF: to_check is not equal to solution {to_check} !== {solution} (float comparison)"
            print(retstr)
            return retstr

        # int test
        elif isinstance(to_check, int):
            retstr = f"DIFF: to_check is not equal to solution {to_check} !== {solution} (int comparison)"
            print(retstr)
            return retstr

        # sting test
        elif isinstance(to_check, str):
            retstr = f"DIFF: to_check is not equal to solution (string comparison)"
            retstr += f"to check: \n =================== \n{to_check}"
            retstr += f"solution: \n =================== \n{solution}"
            print(retstr)
            return retstr

        # tuple comparison
        elif isinstance(to_check, tuple):
            retstr = f"DIFF: to_check is not equal to solution {to_check} !== {solution} (tuple comparison)"
            print(retstr)
            return retstr

        # list test
        elif isinstance(to_check, list):
            # length
            if len(to_check) != len(solution):
                retstr = f"DIFF: to_check length is {len(to_check)} !== {len(solution)} (list comparison)"
                print(retstr)
                return retstr
            # set check
            if set(to_check) == set(solution):
                retstr = f"DIFF: to_check list {len(to_check)}!== {len(solution)} (list comparison), BUT if converted to sets they are identical !! "
                print(retstr)
                return retstr
            else:
                retstr = f"DIFF: to_check list {len(to_check)}!== {len(solution)} (list comparison), no hints found, debug further."
                print(retstr)
                return retstr
        # set test
        elif isinstance(to_check, set):
            # length
            if len(to_check) != len(solution):
                retstr = f"DIFF: to_check length is {len(to_check)} !== {len(solution)} (set comparison)"
                print(retstr)
                return retstr
            else:
                retstr = f"DIFF: to_check list {len(to_check)}!== {len(solution)} (set comparison), no hints found, debug further."
                print(retstr)
                return retstr

        # pandas seriers
        elif isinstance(to_check, pd.Series):
            return compare_series_are_equal(
                testseries=to_check, solutionseries=solution
            )

        # pandas dataframes
        elif isinstance(to_check, pd.DataFrame):
            return compare_df_are_equal(testdf=to_check, solutiondf=solution)

        # metobs_toolkit.Dataset test
        elif to_check.__class__.__name__ == "Dataset":
            return compare_dataset_are_equal(testdataset=to_check, solution=solution)

        # metobs_toolkit.Station test
        elif to_check.__class__.__name__ == "Station":
            return compare_dataset_are_equal(
                testdataset=to_check, solution=solution  # also works for stations
            )

        # Else
        else:
            retstr = f"DIFF: to_check list {to_check}!== {solution} (NotImplemented type comparison), no hints found, debug further."
            print(retstr)
            return retstr


def compare_dataset_are_equal(testdataset, solution):

    # test if the metadf attribute is equal
    if not testdataset.metadf.equals(solution.metadf):
        # diff can be spotted in metadf
        retstr = "DIFF: A difference between two Datasets can be seen by inspecting the metadf!! \n "
        retstr += "\n ................................\n"
        print(retstr)
        return compare_df_are_equal(
            testdf=testdataset.metadf, solutiondf=solution.metadf
        )

    # test if the gaps are different
    if not testdataset.gapsdf.equals(solution.gapsdf):
        # diff can be spotted in metadf
        retstr = "DIFF: A difference between two Datasets can be seen by inspecting the gapsdf!! \n "
        retstr += "\n ................................\n"
        print(retstr)
        return compare_df_are_equal(
            testdf=testdataset.gapsdf, solutiondf=solution.gapsdf
        )

    # test if the outliers are different
    if not testdataset.outliersdf.equals(solution.outliersdf):
        # diff can be spotted in metadf
        retstr = "DIFF: A difference between two Datasets can be seen by inspecting the outliersdf!! \n "
        retstr += "\n ................................\n"
        print(retstr)
        return compare_df_are_equal(
            testdf=testdataset.outliersdf, solutiondf=solution.outliersdf
        )

    # test if the df is different
    if not testdataset.df.equals(solution.df):
        # diff can be spotted in metadf
        retstr = "DIFF: A difference between two Datasets can be seen by inspecting the df!! \n "
        retstr += "\n ................................\n"
        print(retstr)
        return compare_df_are_equal(testdf=testdataset.df, solutiondf=solution.df)

    # else deep investigation needed
    retstr = "DIFF: A difference between two Datasets is detected, but no hints are found where they differ. \n "
    retstr += "\n ................................\n"
    retstr += "Deep investigation is required !! "
    return retstr


def compare_series_are_equal(testseries, solutionseries):
    retstr = "DIFF: "

    # check for duplicates in index
    if (not solutionseries.index.duplicated().any()) & (
        testseries.index.duplicated().any()
    ):
        retstr += "The tested series contains duplicated indexes!: \n"
        retstr += f"{testseries.loc[testseries.index.duplicated()]} \n"
        print(retstr)
        return testseries.loc[testseries.index.duplicated()]  # gives more insight

    added_data = testseries[~testseries.index.isin(solutionseries.index)]
    missing_data = solutionseries[~solutionseries.index.isin(testseries.index)]
    # situation 1: testseries has additive rows wrt solution
    if (not added_data.empty) & (missing_data.empty):
        retstr += "These rows are found IN the to_test but NOT IN the solution: \n"
        retstr += f"{added_data}"
        print(retstr)
        return testseries

    # situation 2 tes # situation 1: testseries has additive rows wrt solution
    if (added_data.empty) & (not missing_data.empty):
        retstr += "These rows are missing in the to_test, that are in the solution : \n"
        retstr += f"{missing_data}"
        print(retstr)
        return missing_data

    # situation 3 some rows are lacking, and some rows are added
    if (not added_data.empty) & (not missing_data.empty):
        retstr += f"The following rows are IN the to_test but NOT in the solution: \n {added_data}\n"
        retstr += "\n ================== \n"
        retstr += f"The following rows missing in the to_test: \n {missing_data}\n"
        return retstr

    # test index levels
    if list(testseries.index.names) != list(solutionseries.index.names):
        retstr += (
            "The index structure is not the same between the test and the solution \n"
        )

        retstr += f"Index structure of the test:\n {testseries.index.names} \n"
        retstr += "\n ================== \n"
        retstr += f"Index structure of the solution:\n {solutionseries.index.names}"
        return retstr

    # test is indexes are the same
    if not (testseries.index == solutionseries.index).all():
        retstr += "The index of the test and solutionseries are not the same! \n"
        index_diff = testseries.index.difference(solutionseries.index)
        retstr += f"The following indexes are in the testseries but not in the solutionseries: {index_diff}\n"
        retstr += "\n ================= \n"
        index_diff = solutionseries.index.difference(testseries.index)
        retstr += f"The following indexes are in the solutionseries but not in the testseries: {index_diff}\n"
        return retstr

    are_equal = testseries.equals(solutionseries)
    if not are_equal:

        diffdf = testseries.compare(
            solutionseries, keep_shape=False, result_names=("TO_TEST", "SOLUTION")
        )

        retstr += "The stucture of series is equal, but the content differs!! \n"
        retstr == f"See the diff_df for more details:\n {diffdf}\n"
        print(retstr)
        return diffdf  # gives more insight

    raise UnforseenDifference(
        "Encountering a unforseen differnce case for dataframe comparison."
    )


def compare_df_are_equal(testdf, solutiondf):
    retstr = "DIFF: "
    # test is columns are the same
    if not set(testdf.columns) == set(solutiondf.columns):
        retstr += "The columns of the to_check and solution are not the same\n"
        retstr += f"Columns in the testdf :\n {testdf.columns} \n"
        retstr += "\n ================\n"
        retstr += f"Columns in the solutionsdf :\n {solutiondf.columns}"
        return retstr

    # order columns alphabetically
    testdf = testdf.reindex(sorted(testdf.columns), axis=1)
    solutiondf = solutiondf.reindex(sorted(solutiondf.columns), axis=1)
    # common error/bug signal checks

    # check for duplicates in index
    if (not solutiondf.index.duplicated().any()) & (testdf.index.duplicated().any()):
        retstr += "The tested dataframe contains duplicated indexes!: \n"
        retstr += f"{testdf.loc[testdf.index.duplicated()]} \n"
        print(retstr)
        return testdf.loc[testdf.index.duplicated()]  # gives more insight

    # overlap_data = testdf[testdf.index.isin(solutiondf.index)]
    added_data = testdf[~testdf.index.isin(solutiondf.index)]
    missing_data = solutiondf[~solutiondf.index.isin(testdf.index)]
    # situation 1: testdf has additive rows wrt solution
    if (not added_data.empty) & (missing_data.empty):
        retstr += "These rows are found IN the to_test but NOT IN the solution: \n"
        retstr += f"{added_data}"
        print(retstr)
        return added_data  # gives more insight

    # situation 2 tes # situation 1: testdf has additive rows wrt solution
    if (added_data.empty) & (not missing_data.empty):
        retstr += "These rows are missing in the to_test, that are in the solution : \n"
        retstr += f"{missing_data}"
        print(retstr)
        return missing_data  # gives more insight

    # situation 3 some rows are lacking, and some rows are added
    if (not added_data.empty) & (not missing_data.empty):
        retstr += f"The following rows are IN the to_test but NOT in the solution: \n {added_data}\n"
        retstr += "\n ================== \n"
        retstr += f"The following rows missing in the to_test: \n {missing_data}\n"
        return retstr

    # test index levels
    if list(testdf.index.names) != list(solutiondf.index.names):
        retstr += (
            "The index structure is not the same between the test and the solution \n"
        )

        retstr += f"Index structure of the test:\n {testdf.index.names} \n"
        retstr += "\n ================== \n"
        retstr += f"Index structure of the solution:\n {solutiondf.index.names}"
        return retstr

    # test is indexes are the same
    if not (testdf.index == solutiondf.index).all():
        retstr += "The index of the test and solutiondf are not the same! \n"
        index_diff = testdf.index.difference(solutiondf.index)
        retstr += f"The following indexes are in the testdf but not in the solutiondf: {index_diff}\n"
        retstr += "\n ================= \n"
        index_diff = solutiondf.index.difference(testdf.index)
        retstr += f"The following indexes are in the solutiondf but not in the testdf: {index_diff}\n"
        return retstr

    # Convert numeric columns to float32
    # numeric_cols = testdf.select_dtypes(include=['number']).columns
    # testdf[numeric_cols] = testdf[numeric_cols].astype('float32')
    # solutiondf[numeric_cols] = solutiondf[numeric_cols].astype('float32')
    are_equal = testdf.equals(solutiondf)
    if not are_equal:

        diffdf = testdf.compare(
            solutiondf, keep_shape=False, result_names=("TO_TEST", "SOLUTION")
        )

        retstr += "The stucture of datasets is equal, but the content differs!! \n"
        retstr == f"See the diff_df for more details:\n {diffdf}\n"
        print(retstr)
        return diffdf  # gives more insight

    raise UnforseenDifference(
        "Encountering a unforseen differnce case for dataframe comparison."
    )


class UnforseenDifference(Exception):
    """Raise when encountering an unforseen difference case"""


class SolutionNotExisting(Exception):
    """Raise when the solutionfile does not exist"""
