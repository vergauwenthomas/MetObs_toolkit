import logging
import json

from metobs_toolkit.backend_collection.loggingmodule import log_entry

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
    return
