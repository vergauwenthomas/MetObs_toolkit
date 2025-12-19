import logging
import json
from os import PathLike
from pathlib import Path

from metobs_toolkit.backend_collection.decorators import log_entry

logger = logging.getLogger("<metobs_toolkit>")


@log_entry
def write_dict_to_json(dictionary: dict, trgfile: str) -> None:
    """
    Write a dictionary to a JSON file.

    Parameters
    ----------
    templdict : dict
        The template dictionary to write.
    trgfile : str
        The target file path for the JSON output.

    Returns
    -------
    None
    """
    j = json.dumps(dictionary, indent=4)
    with open(trgfile, "w") as f:
        print(j, file=f)
    return None


@log_entry
def fmt_output_filepath(
    filepath: str | PathLike | None,
    default_filename: str,
    overwrite: bool = False,
    suffix: str | None = None,
) -> Path:
    """
    Validate and prepare an output file path for writing.

    This function checks if the target file path is valid and writable.
    If the file already exists, it can optionally be removed. If the
    parent directory does not exist, it will be created.

    Parameters
    ----------
    filepath : str or path-like or None
        The file path to validate and prepare. If None, defaults to
        the current working directory with the default_filename.
    default_filename : str
        The default filename to use if filepath is None. This will be
        joined with the current working directory.
    overwrite : bool, optional
        If True, an existing file at the path will be removed.
        If False, a FileExistsError is raised when the file exists.
        Default is False.
    suffix : str or None, optional
        If provided, ensures the file path has this suffix (e.g., '.html', '.csv').
        If the filepath does not end with the specified suffix, it will be appended.
        Default is None.

    Returns
    -------
    Path
        The validated file path as a Path object.

    Raises
    ------
    FileExistsError
        If the file already exists and overwrite is False.
    IOError
        If the parent directory cannot be created.
    """
    # use default filepath if none provided
    if filepath is None:
        filepath = Path.cwd() / default_filename

    # Convert to Path object
    if not isinstance(filepath, Path):
        filepath = Path(filepath)

    # Check suffix
    if suffix:
        if filepath.suffix != suffix:
            logger.info(f"Appending suffix '{suffix}' to file path.")
            filepath = filepath.with_suffix(suffix)

    if filepath.exists():
        if overwrite:
            logger.info(f"File {filepath} exists but will be overwritten.")
            filepath.unlink()
        else:
            raise FileExistsError(
                f"File {filepath} already exists and overwrite is set to False."
            )

    parent_dir = filepath.parent
    if not parent_dir.exists():
        try:
            parent_dir.mkdir(parents=True, exist_ok=True)
            logger.info(f"Created directory {parent_dir}.")
        except Exception as e:
            raise IOError(f"Cannot create directory {parent_dir}: {e}")

    return filepath
