import logging
import json

logger = logging.getLogger("<metobs_toolkit>")


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
    logger.debug("Entering write_dict_to_json()")
    j = json.dumps(dictionary, indent=4)
    with open(trgfile, "w") as f:
        print(j, file=f)
    return
