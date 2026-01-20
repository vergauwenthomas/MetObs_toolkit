"""Solution management for test fixtures.

This module provides utilities to store and retrieve test solutions as parquet files,
enabling reproducible testing without relying on pickled objects.
"""

from pathlib import Path
from typing import Union

import geopandas as gpd
import json
import pandas as pd


# Project paths
libfolder = Path(__file__).resolve().parent.parent
datadir = libfolder / "tests" / "data"


# =============================================================================
# Exceptions
# =============================================================================


class SolutionNotExisting(Exception):
    """Raised when the solution file does not exist."""


class UnforeseenDifference(Exception):
    """Raised when encountering an unforeseen difference case."""


# =============================================================================
# Serialized Data Classes
# =============================================================================


class SerializedDataset:
    """Container for deserialized Dataset data for comparison."""

    def __init__(self, data: dict):
        self.df = data.get("df")
        self.metadf = _format_metadata(data.get("metadf"))
        self.gapsdf = data.get("gapsdf")
        self.modeldatadf = data.get("modeldatadf")
        self.outliersdf = data.get("outliersdf")


class SerializedStation:
    """Container for deserialized Station data for comparison."""

    def __init__(self, data: dict):
        self.df = data.get("df")
        self.metadf = _format_metadata(data.get("metadf"))
        self.gapsdf = data.get("gapsdf")
        self.modeldatadf = data.get("modeldatadf")
        self.outliersdf = data.get("outliersdf")


class SerializedAnalysis:
    """Container for deserialized Analysis data for comparison."""

    def __init__(self, data: dict):
        self.fulldf = data.get("fulldf")
        self.metadf = _format_metadata(data.get("metadf"))


# =============================================================================
# Solution Fixer Class
# =============================================================================


class SolutionFixer2:
    """Store and retrieve test solutions as parquet files."""

    def __init__(self, solutiondir: Path):
        self.solutiondir = solutiondir

    def get_solution_dir(self, methodname: str, testfile: str, classname: str) -> Path:
        """Get the directory path for a specific test solution."""
        testfile_stem = str(testfile).rstrip(".py")
        solutiondir = (
            self.solutiondir
            / f"{testfile_stem}_solutions"
            / classname.lower()
            / methodname
        )
        solutiondir.mkdir(parents=True, exist_ok=True)
        return solutiondir

    def create_solution(
        self, solution, methodname: str, testfile: str, classname: str
    ) -> None:
        """Save a solution object to disk."""
        base_path = self.get_solution_dir(methodname, testfile, classname)

        print(
            f"!! OVERWRITING SOLUTION FOR {testfile} --> {classname}:{methodname} !!!"
        )

        classname_map = {
            "Dataset": _store_dataset,
            "Station": _store_station,
            "Analysis": _store_analysis,
        }

        obj_classname = solution.__class__.__name__
        if obj_classname in classname_map:
            classname_map[obj_classname](solution, base_path)
        elif isinstance(solution, pd.DataFrame):
            _store_dataframe(solution, base_path)
        elif isinstance(solution, dict):
            _store_dict(solution, base_path)
        elif isinstance(solution, str):
            _store_string(solution, base_path)
        else:
            raise NotImplementedError(
                f"create_solution does not support {obj_classname} objects."
            )

    def get_solution(self, methodname: str, testfile: str, classname: str) -> Union[
        SerializedDataset,
        SerializedStation,
        SerializedAnalysis,
        pd.DataFrame,
        dict,
        str,
    ]:
        """Load a solution from disk."""
        base_path = self.get_solution_dir(methodname, testfile, classname)

        if not base_path.exists():
            raise SolutionNotExisting(f"The solution at {base_path} does not exist!")

        with open(base_path / "datatype.json", "r") as f:
            metadata = json.load(f)

        solution_class = metadata.get("class")

        reader_map = {
            "Dataset": _read_dataset,
            "Station": _read_station,
            "Analysis": _read_analysis,
            "Dict": _read_dict,
            "DataFrame": _read_dataframe,
            "String": _read_string,
        }

        if solution_class not in reader_map:
            raise NotImplementedError(
                f"get_solution does not support {solution_class} objects."
            )

        return reader_map[solution_class](base_path)


# =============================================================================
# Storage Functions
# =============================================================================


def _write_datatype(dir_path: Path, class_name: str) -> None:
    """Write the datatype metadata file."""
    dir_path.mkdir(parents=True, exist_ok=True)
    with open(dir_path / "datatype.json", "w") as f:
        json.dump({"class": class_name}, f)


def _store_string(data_str: str, dir_path: Path) -> None:
    """Store a string solution."""
    _write_datatype(dir_path, "String")
    with open(dir_path / "solution_string.txt", "w") as f:
        f.write(data_str)


def _store_dataframe(df: pd.DataFrame, dir_path: Path) -> None:
    """Store a DataFrame solution."""
    _write_datatype(dir_path, "DataFrame")
    df.to_parquet(dir_path / "solution_df.parquet")


def _store_dataset(dataset, dir_path: Path) -> None:
    """Store a Dataset solution."""
    _write_datatype(dir_path, "Dataset")
    dataset.df.to_parquet(dir_path / "solution_df.parquet")
    dataset.metadf.to_parquet(dir_path / "solution_metadf.parquet")
    dataset.gapsdf.to_parquet(dir_path / "solution_gapsdf.parquet")
    dataset.modeldatadf.to_parquet(dir_path / "solution_modeldatadf.parquet")
    dataset.outliersdf.to_parquet(dir_path / "solution_outliersdf.parquet")


def _store_station(station, dir_path: Path) -> None:
    """Store a Station solution."""
    _write_datatype(dir_path, "Station")
    station.df.to_parquet(dir_path / "solution_df.parquet")
    station.metadf.to_parquet(dir_path / "solution_metadf.parquet")
    station.gapsdf.to_parquet(dir_path / "solution_gapsdf.parquet")
    station.modeldatadf.to_parquet(dir_path / "solution_modeldatadf.parquet")
    station.outliersdf.to_parquet(dir_path / "solution_outliersdf.parquet")


def _store_analysis(analysis, dir_path: Path) -> None:
    """Store an Analysis solution."""
    _write_datatype(dir_path, "Analysis")
    analysis.df.to_parquet(dir_path / "solution_df.parquet")
    analysis.fulldf.to_parquet(dir_path / "solution_fulldf.parquet")
    analysis.metadf.to_parquet(dir_path / "solution_metadf.parquet")


def _store_dict(data_dict: dict, dir_path: Path) -> None:
    """Store a dict solution containing Dataset, Station, or DataFrame values."""
    _write_datatype(dir_path, "Dict")

    for key, val in data_dict.items():
        key_dir = dir_path / str(key)
        obj_classname = val.__class__.__name__

        if obj_classname == "Dataset":
            _store_dataset(val, key_dir)
        elif obj_classname == "Station":
            _store_station(val, key_dir)
        elif isinstance(val, pd.DataFrame):
            _store_dataframe(val, key_dir)
        else:
            raise NotImplementedError(
                f"store_dict does not support {obj_classname} objects."
            )


# =============================================================================
# Reading Functions
# =============================================================================


def _read_parquet_files(dir_path: Path) -> dict:
    """Read all solution parquet files from a directory."""
    result = {}
    for f in dir_path.glob("solution_*.parquet"):
        attr_name = f.stem.replace("solution_", "")
        result[attr_name] = pd.read_parquet(f)
    return result


def _read_string(dir_path: Path) -> str:
    """Read a string solution."""
    with open(dir_path / "solution_string.txt", "r") as f:
        return f.read()


def _read_dataframe(dir_path: Path) -> pd.DataFrame:
    """Read a DataFrame solution."""
    return pd.read_parquet(dir_path / "solution_df.parquet")


def _read_dataset(dir_path: Path) -> SerializedDataset:
    """Read a Dataset solution."""
    return SerializedDataset(_read_parquet_files(dir_path))


def _read_station(dir_path: Path) -> SerializedStation:
    """Read a Station solution."""
    return SerializedStation(_read_parquet_files(dir_path))


def _read_analysis(dir_path: Path) -> SerializedAnalysis:
    """Read an Analysis solution."""
    return SerializedAnalysis(_read_parquet_files(dir_path))


def _read_dict(dir_path: Path) -> dict:
    """Read a dict solution."""
    result = {}
    for subdir in dir_path.iterdir():
        if not subdir.is_dir():
            continue

        with open(subdir / "datatype.json", "r") as f:
            metadata = json.load(f)

        solution_class = metadata.get("class")

        if solution_class == "Dataset":
            result[subdir.stem] = _read_dataset(subdir)
        elif solution_class == "Station":
            result[subdir.stem] = _read_station(subdir)
        elif solution_class == "DataFrame":
            result[subdir.stem] = _read_dataframe(subdir)
        else:
            raise NotImplementedError(
                f"read_dict does not support {solution_class} objects."
            )

    return result


# =============================================================================
# Utility Functions
# =============================================================================


def _format_metadata(metadf: pd.DataFrame) -> pd.DataFrame:
    """Convert metadata DataFrame to GeoDataFrame if geometry column exists."""
    if metadf is not None and "geometry" in metadf.columns:
        metadf = gpd.GeoDataFrame(
            metadf, geometry=gpd.points_from_xy(metadf["lon"], metadf["lat"])
        )
    return metadf


# =============================================================================
# Comparison Functions
# =============================================================================


def _compare_df_attr(
    testobj, solutionobj, attr: str, exclude_columns: Union[str, list, None] = None
) -> None:
    """Compare a DataFrame attribute between test and solution objects."""
    try:
        left_df = getattr(testobj, attr)
        right_df = getattr(solutionobj, attr)

        if exclude_columns:
            if isinstance(exclude_columns, str):
                exclude_columns = [exclude_columns]
            left_df = left_df.drop(
                columns=[c for c in exclude_columns if c in left_df.columns]
            )
            right_df = right_df.drop(
                columns=[c for c in exclude_columns if c in right_df.columns]
            )

        pd.testing.assert_frame_equal(left=left_df, right=right_df, check_exact=False)
    except AssertionError as e:
        raise AssertionError(f"DIFF in {attr}-attribute:\n {e}")


def _serialized_comparison(to_check, solution, **kwargs) -> None:
    """Compare all attributes from a serialized solution."""
    for key in solution.__dict__.keys():
        _compare_df_attr(testobj=to_check, solutionobj=solution, attr=key, **kwargs)


def assert_equality(to_check, solution, **kwargs) -> None:
    """Assert equality between a test result and a solution.

    Provides detailed debug information when differences are found.
    Supports Dataset, Station, Analysis, DataFrame, Series, and primitive types.
    """
    to_check_class = to_check.__class__.__name__
    solution_class = solution.__class__.__name__

    # Dataset comparisons
    if to_check_class == "Dataset" and solution_class == "SerializedDataset":
        _serialized_comparison(to_check, solution, **kwargs)
        return

    if to_check_class == "Dataset" and solution_class == "Dataset":
        for attr in ["metadf", "gapsdf", "modeldatadf", "outliersdf", "df"]:
            _compare_df_attr(to_check, solution, attr, **kwargs)
        assert to_check.obstypes == solution.obstypes, "Mismatch in obstypes!"
        return

    # Station comparisons
    if to_check_class == "Station" and solution_class == "SerializedStation":
        _serialized_comparison(to_check, solution, **kwargs)
        return

    if to_check_class == "Station" and solution_class == "Station":
        for attr in ["metadf", "gapsdf", "modeldatadf", "outliersdf", "df"]:
            _compare_df_attr(to_check, solution, attr, **kwargs)
        return

    # Analysis comparisons
    if to_check_class == "Analysis" and solution_class == "SerializedAnalysis":
        _serialized_comparison(to_check, solution, **kwargs)
        return

    if to_check_class == "Analysis" and solution_class == "Analysis":
        for attr in ["df", "fulldf", "metadf"]:
            _compare_df_attr(to_check, solution, attr, **kwargs)
        return

    # Type mismatch
    if type(to_check) != type(solution):
        raise AssertionError(
            f"DIFF: to_check type is {type(to_check)}, "
            f"while solution is of type {type(solution)}"
        )

    # Primitive type comparisons
    if isinstance(to_check, (float, int)):
        if to_check != solution:
            raise AssertionError(
                f"DIFF: {to_check} !== {solution} ({type(to_check).__name__} comparison)"
            )
        return

    if isinstance(to_check, str):
        if to_check != solution:
            raise AssertionError(
                f"DIFF: string mismatch\n"
                f"to_check:\n===================\n{to_check}\n"
                f"solution:\n===================\n{solution}"
            )
        return

    if isinstance(to_check, tuple):
        if to_check != solution:
            raise AssertionError(f"DIFF: {to_check} !== {solution} (tuple comparison)")
        return

    # Collection comparisons
    if isinstance(to_check, list):
        if to_check != solution:
            if len(to_check) != len(solution):
                raise AssertionError(
                    f"DIFF: list length {len(to_check)} !== {len(solution)}"
                )
            if set(to_check) == set(solution):
                raise AssertionError(
                    f"DIFF: lists differ but sets are identical (order mismatch?)"
                )
            raise AssertionError("DIFF: lists differ, debug further.")
        return

    if isinstance(to_check, set):
        if to_check != solution:
            if len(to_check) != len(solution):
                raise AssertionError(
                    f"DIFF: set length {len(to_check)} !== {len(solution)}"
                )
            raise AssertionError("DIFF: sets differ, debug further.")
        return

    # Pandas comparisons
    if isinstance(to_check, pd.Series):
        pd.testing.assert_series_equal(
            left=to_check, right=solution, check_exact=False, rtol=0.001
        )
        return

    if isinstance(to_check, pd.DataFrame):
        pd.testing.assert_frame_equal(
            left=to_check, right=solution, check_exact=False, rtol=0.001
        )
        return

    # Fallback
    raise AssertionError(
        f"DIFF: {to_check} !== {solution} (unsupported type comparison)"
    )
