from pathlib import Path
import pickle

import pandas as pd
import geopandas as gpd
import json
import pandas.testing

libfolder = Path(str(Path(__file__).resolve())).parent.parent
# testdatadir
datadir = libfolder.joinpath("tests").joinpath("data")


# class SolutionFixer:
#     def __init__(self, solutiondir):
#         self.solutiondir = Path(solutiondir)

#     def get_solution(self, testfile, classname, methodname):
#         basefoldername = f"{str(testfile).strip('.py')}_solutions"
#         subfoldername = f"{classname.lower()}"
#         trgfilename = f"{str(methodname)}.pkl"

#         solutionfile = (
#             self.solutiondir.joinpath(basefoldername)
#             .joinpath(subfoldername)
#             .joinpath(trgfilename)
#         )
#         if not solutionfile.exists():
#             raise SolutionNotExisting(
#                 f"The solution of {testfile} -->{classname}:{methodname} does not exist (at {solutionfile})!"
#             )

#         with open(solutionfile, "rb") as file:
#             solution_data = pickle.load(file)
#         return solution_data

#     def create_solution(self, solutiondata, testfile, classname, methodname):
#         # construct path
#         basefoldername = f"{str(testfile).strip('.py')}_solutions"
#         subfoldername = f"{classname.lower()}"
#         trgfilename = f"{str(methodname)}.pkl"

#         solutionfile = (
#             self.solutiondir.joinpath(basefoldername)
#             .joinpath(subfoldername)
#             .joinpath(trgfilename)
#         )
#         # clear previous solutions
#         if solutionfile.exists():
#             print(
#                 f"!! OVERWRITING SOLUTION FOR  {testfile} --> {classname}:{methodname} !!! "
#             )
#             # delete file
#             solutionfile.unlink()

#         # Create directory if it does not exist
#         solutionfile.parent.mkdir(parents=True, exist_ok=True)

#         # Pickle data object
#         with open(solutionfile, "wb") as file:
#             pickle.dump(solutiondata, file)


class SolutionFixer2:
    """Store solutions as data files, not pickled objects."""
    
    def __init__(self, solutiondir: Path):
        self.solutiondir = solutiondir
        
    
    def get_solution_dir(self, methodname: str, testfile: str, classname: str) -> Path:

        solutiondir = (
            self.solutiondir.joinpath(f"{str(testfile).strip('.py')}_solutions")
            .joinpath(f"{classname.lower()}")
            .joinpath(f"{methodname}")
        )
        
        if not solutiondir.exists():
            solutiondir.mkdir(parents=True, exist_ok=True)
        return solutiondir
    
    def create_solution(self, solution, methodname: str, testfile: str, classname: str):
        """Extract and save only the data attributes."""
        base_path = self.get_solution_dir(methodname, testfile, classname)
        base_path.mkdir(exist_ok=True)
        
       
        print(
            f"!! OVERWRITING SOLUTION FOR  {testfile} --> {classname}:{methodname} !!! "
        )
        #if solution is Metobs-toolkit.Dataset
        if solution.__class__.__name__ == "Dataset":
            store_dataset(solution, base_path)
        elif solution.__class__.__name__ == "Station":
            store_station(solution, base_path)
        elif isinstance(solution, pd.DataFrame):
            store_pandasdf(solution, base_path)
        elif isinstance(solution, dict):
            store_dict(solution, base_path)
        elif isinstance(solution, str):
            store_string(solution, base_path)
        else:
            raise NotImplementedError(
                "SolutionFixer2.create_solution only supports Dataset, Station, and Analysis objects.")
        
           
    def get_solution(self, methodname: str, testfile: str, classname: str) -> dict:
        """Load solution data as a dict of DataFrames."""
        base_path = self.get_solution_dir(methodname, testfile, classname)
        if not base_path.exists():
            raise SolutionNotExisting(
                f"The solution at {base_path} does not exist!"
            )

        with open(base_path / "datatype.json", "r") as f:
            metobsobj = json.load(f)
        
        solutionclass = metobsobj.get("class")
        

        if solutionclass == 'Dataset':
            return read_dataset(base_path)
        elif solutionclass == 'Station':
            return read_station(base_path)
        elif solutionclass == 'Dict':
            return read_dict(base_path)
        elif solutionclass == 'DataFrame':
            return read_pandasdf(base_path)
        elif solutionclass == 'String':
            return read_string(base_path)
        
        else:
            raise NotImplementedError(
                "SolutionFixer2.get_solution only supports Dataset, Station, and DataFrame objects.")
            
def store_string(data_str: str, dir: Path):
    """Store the dataset as a json-serializable dict."""
    
    dir.mkdir(parents=True, exist_ok=True)
    with open(dir / "datatype.json", "w") as f:
            json.dump({"class": "String"}, f)
            
    with open(dir / "solution_string.txt", "w") as f:
        f.write(data_str)
        
def read_string(dir: Path) -> str:
    with open(dir / "solution_string.txt", "r") as f:
        data_str = f.read()
    return data_str

def store_dict(data_dict: dict, dir: Path):
    """Store the dataset as a json-serializable dict."""
    
    dir.mkdir(parents=True, exist_ok=True)
    with open(dir / "datatype.json", "w") as f:
            json.dump({"class": "Dict"}, f)
            
    for key, val in data_dict.items():
        keydir = dir / f"{key}"
        if val.__class__.__name__ == "Dataset":
            store_dataset(val, keydir)
            
        elif val.__class__.__name__ == "Station":  
            store_station(val, keydir)
        
        elif isinstance(val, pd.DataFrame):
            store_pandasdf(df=val, dir=keydir)
        else:
            raise NotImplementedError(
                "store_dict only supports Dataset, Station, and DataFrame objects.")
def read_dict(dir: Path) -> dict:
    returndict = {}
    for subdir in dir.iterdir():
        if subdir.is_dir():
            with open(subdir / "datatype.json", "r") as f:
                metobsobj = json.load(f)
            solutionclass = metobsobj.get("class")
            if solutionclass == 'Dataset':
                returndict[subdir.stem] = read_dataset(subdir)
            elif solutionclass == 'Station':
                returndict[subdir.stem] = read_station(subdir)
            elif solutionclass == 'DataFrame':
                returndict[subdir.stem] = read_pandasdf(subdir)
            else:
                raise NotImplementedError(
                    "read_dict only supports Dataset, Station, and DataFrame objects.")
    return returndict
    
    
    
def read_dataset(dir: Path):
    def file_to_attr(fpath: Path) -> str:
        return fpath.stem.replace("solution_", "")

    result = {}
    for f in dir.glob("solution_*.parquet"):
        result[file_to_attr(f)] = pd.read_parquet(f)
    return SerializedDataset(result)

def read_station(dir: Path):
    def file_to_attr(fpath: Path) -> str:
        return fpath.stem.replace("solution_", "")

    result = {}
    for f in dir.glob("solution_*.parquet"):
        result[file_to_attr(f)] = pd.read_parquet(f)
    return SerializedStation(result)

def read_pandasdf(dir: Path) -> pd.DataFrame:
    df = pd.read_parquet(dir / "solution_df.parquet")
    return df


        
def store_dataset(dataset, dir: Path):
    """Store the dataset as a json-serializable dict."""
    
    dir.mkdir(parents=True, exist_ok=True)
    with open(dir / "datatype.json", "w") as f:
            json.dump({"class": "Dataset"}, f)
    dataset.df.to_parquet(dir / "solution_df.parquet")
    dataset.metadf.to_parquet(dir / "solution_metadf.parquet")
    dataset.gapsdf.to_parquet(dir / "solution_gapsdf.parquet")
    dataset.modeldatadf.to_parquet(dir / "solution_modeldatadf.parquet")
    dataset.outliersdf.to_parquet(dir / "solution_outliersdf.parquet")
    
    #Specific attributes
    # solution.obstypes.to_parquet(dir / "solution_obstypes.parquet")
    
def store_station(station, dir: Path):
    """Store the dataset as a json-serializable dict."""
    
    dir.mkdir(parents=True, exist_ok=True)
    with open(dir / "datatype.json", "w") as f:
            json.dump({"class": "Station"}, f)
    station.df.to_parquet(dir / "solution_df.parquet")
    station.metadf.to_parquet(dir / "solution_metadf.parquet")
    station.gapsdf.to_parquet(dir / "solution_gapsdf.parquet")
    station.modeldatadf.to_parquet(dir / "solution_modeldatadf.parquet")
    station.outliersdf.to_parquet(dir / "solution_outliersdf.parquet")

def store_pandasdf(df: pd.DataFrame, dir: Path):
    """Store the dataset as a json-serializable dict."""
    
    dir.mkdir(parents=True, exist_ok=True)
    with open(dir / "datatype.json", "w") as f:
            json.dump({"class": "DataFrame"}, f)
    df.to_parquet(dir / "solution_df.parquet")
    
class SerializedDataset():
    """A class to hold the deserialized data of a Dataset for comparison."""
    
    def __init__(self, data: dict):
        self.df = data.get("df")
        self.metadf = format_metadata(data.get("metadf"))
        self.gapsdf = data.get("gapsdf")
        self.modeldatadf = data.get("modeldatadf")
        self.outliersdf = data.get("outliersdf")
        
        
class SerializedStation():
    """A class to hold the deserialized data of a Dataset for comparison."""
    
    def __init__(self, data: dict):
        self.df = data.get("df")
        self.metadf = format_metadata(data.get("metadf"))
        self.gapsdf = data.get("gapsdf")
        self.modeldatadf = data.get("modeldatadf")
        self.outliersdf = data.get("outliersdf")
       
def format_metadata(metadf):
    
    #geoseries not serializable to json, so recreate
    if 'geometry' in metadf.columns:
        metadf = gpd.GeoDataFrame(metadf, geometry=gpd.points_from_xy(metadf['lon'],
                                                                      metadf['lat']))
    return metadf






def assert_equality(to_check, solution, **kwargs):
    """Returns some debug help when the to_check is not equalt to the solution

    This output is used for developping to point in more details why a test
    is failing.

    The object retured can have differnt types, and is desined for the developper
    to help in the debugging process.

    """

    if ((to_check.__class__.__name__ == "Dataset") and 
        (solution.__class__.__name__ == "SerializedDataset")):
        seriealize_comparison(to_check, solution, **kwargs)
    
    elif ((to_check.__class__.__name__ == "Dataset") and 
        (solution.__class__.__name__ == "Dataset")):
        
        compare_df_attr(to_check, solution, "metadf", **kwargs)
        compare_df_attr(to_check, solution, "gapsdf", **kwargs)
        compare_df_attr(to_check, solution, "modeldatadf", **kwargs)
        compare_df_attr(to_check, solution, "outliersdf", **kwargs)
        compare_df_attr(to_check, solution, "df", **kwargs)
        assert (
            to_check.obstypes == solution.obstypes
        ), "There is a mismatch in obstypes with the solution!"
        
        
        
    elif ((to_check.__class__.__name__ == "Station") and 
        (solution.__class__.__name__ == "SerializedStation")):
        seriealize_comparison(to_check, solution, **kwargs)
        
        
    elif ((to_check.__class__.__name__ == "Station") and 
        (solution.__class__.__name__ == "Station")):
        
        compare_df_attr(to_check, solution, "metadf", **kwargs)
        compare_df_attr(to_check, solution, "gapsdf", **kwargs)
        compare_df_attr(to_check, solution, "modeldatadf", **kwargs)
        compare_df_attr(to_check, solution, "outliersdf", **kwargs)
        compare_df_attr(to_check, solution, "df", **kwargs)
        
         
    
    # type equal test
    elif type(to_check) != type(solution):
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
            rtol=0.001,
            
        )

    # pandas dataframes
    elif isinstance(to_check, pd.DataFrame):
        pd.testing.assert_frame_equal(
            left=to_check,
            right=solution,
            check_exact=False,
            rtol=0.001,
        )

    # # metobs_toolkit.Dataset test
    # elif to_check.__class__.__name__ == "Dataset":
    #     

    # # metobs_toolkit.Station test
    # elif to_check.__class__.__name__ == "Station":
    #     compare_df_attr(to_check, solution, "metadf")
    #     compare_df_attr(to_check, solution, "gapsdf")
    #     compare_df_attr(to_check, solution, "modeldatadf")
    #     compare_df_attr(to_check, solution, "outliersdf", exclude_columns="details")
    #     compare_df_attr(to_check, solution, "df")
    # # metobs_toolkit.Station test
    # elif to_check.__class__.__name__ == "Analysis":
    #     compare_df_attr(to_check, solution, "metadf")
    #     compare_df_attr(to_check, solution, "df")

    # Else
    else:
        retstr = f"DIFF: to_check list {to_check}!== {solution} (NotImplemented type comparison), no hints found, debug further."
        raise AssertionError(retstr)


def compare_df_attr(testobj, solutionobj, attr, exclude_columns=None):
    try:
        left_df = getattr(testobj, attr)
        right_df = getattr(solutionobj, attr)

        # Exclude specified columns if provided
        if exclude_columns:
            if isinstance(exclude_columns, str):
                exclude_columns = [exclude_columns]
            left_df = left_df.drop(
                columns=[c for c in exclude_columns if c in left_df.columns]
            )
            right_df = right_df.drop(
                columns=[c for c in exclude_columns if c in right_df.columns]
            )

        pd.testing.assert_frame_equal(
            left=left_df,
            right=right_df,
            check_exact=False,
        )
    except AssertionError as e:
        raise AssertionError(f"DIFF in {attr}-attribute:\n " + str(e))

def seriealize_comparison(to_check, solution, **kwargs):
    solutiondict = solution.__dict__
    for key in solutiondict.keys():
        compare_df_attr(testobj=to_check, solutionobj=solution, attr=key, **kwargs)
        
    
    



class UnforseenDifference(Exception):
    """Raise when encountering an unforseen difference case"""


class SolutionNotExisting(Exception):
    """Raise when the solutionfile does not exist"""
