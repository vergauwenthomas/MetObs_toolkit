from typing import Iterable, List, Any
import logging

from metobs_toolkit.backend_collection.loggingmodule import log_entry

logger = logging.getLogger("<metobs_toolkit>")

"""Set comparison logic for MetObs functions.

The uniqueness is assumed if the ._id() of the instances differs.
"""

# Note: another design would be to implement the __eq__ and __ne__ methods in
# the metobs classes, so that set() operations could be used. The issue with that
# is that np.nan != np.nan. So 'identical gaps' are then assumed to be non-equal.


@log_entry
def metobs_intersection(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:
    """
    Get the intersection of two collections based on their '_id()' method.

    Parameters
    ----------
    A : Iterable[Any]
        First collection of elements with a '_id()' method.
    B : Iterable[Any]
        Second collection of elements with a '_id()' method.

    Returns
    -------
    List[Any]
        List of elements from A whose '_id()' is also present in B.
    """
    # NOTE: elements returned originate from A
    # transform to id-space (these are strings)
    A_ids = set([a._id() for a in A])
    B_ids = set([b._id() for b in B])
    # Intersection of IDs
    intersection_ids = A_ids.intersection(B_ids)
    # Revert to the object-space
    intersection = [a for a in A if a._id() in intersection_ids]
    return intersection


@log_entry
def metobs_B_not_in_A(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:
    """
    Get elements from B whose '_id()' is not present in A.

    Parameters
    ----------
    A : Iterable[Any]
        First collection of elements with a '_id()' method.
    B : Iterable[Any]
        Second collection of elements with a '_id()' method.

    Returns
    -------
    List[Any]
        List of elements from B whose '_id()' is not present in A.
    """
    # transform to id-space (these are strings)
    A_ids = set([a._id() for a in A])
    B_ids = set([b._id() for b in B])
    # Intersection of IDs
    b_not_in_a_ids = B_ids - A_ids
    # Revert to the object-space
    return [b for b in B if b._id() in b_not_in_a_ids]


@log_entry
def metobs_A_not_in_B(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:
    """
    Get elements from A whose '_id()' is not present in B.

    Parameters
    ----------
    A : Iterable[Any]
        First collection of elements with a '_id()' method.
    B : Iterable[Any]
        Second collection of elements with a '_id()' method.

    Returns
    -------
    List[Any]
        List of elements from A whose '_id()' is not present in B.
    """
    # transform to id-space (these are strings)
    return metobs_B_not_in_A(A=B, B=A)


@log_entry
def get_by_id(A: Iterable[Any], id: str) -> Any:
    """
    Retrieve an element from a collection by its '_id()' value.

    Parameters
    ----------
    A : Iterable[Any]
        Collection of elements with a '_id()' method.
    id : str
        The id to search for.

    Returns
    -------
    Any
        The element with the matching id.

    Raises
    ------
    IndexError
        If no element or more than one element is found with the given id.
    """
    candidates = [a for a in A if a._id() == id]
    if len(candidates) > 1:
        raise IndexError(
            f"More than one element found with id={id} in the collection: {A}"
        )
    if len(candidates) == 0:
        raise IndexError(f"No element found with id={id} in the collection: {A}")
    return candidates[0]


@log_entry
def join_collections(col_A: Iterable[Any], col_B: Iterable[Any]) -> List[Any]:
    """
    Join two collections, applying in-depth addition for elements with the same id.

    Parameters
    ----------
    col_A : Iterable[Any]
        Collection of elements with a '_id()' method.
    col_B : Iterable[Any]
        Collection of elements with a '_id()' method.

    Returns
    -------
    List[Any]
        The joined collection. This collection holds elements without duplicated ids.
    """
    joined_col = []
    # Triage elements using metobs logic
    same_id_elem = metobs_intersection(A=list(col_A), B=list(col_B))
    # Note: same_id_elem originates from col_A!

    unique_A_elem = metobs_A_not_in_B(A=list(col_A), B=list(col_B))
    joined_col.extend(unique_A_elem)  # simple join

    unique_B_elem = metobs_B_not_in_A(A=list(col_A), B=list(col_B))
    joined_col.extend(unique_B_elem)  # simple join

    # Rely on the __add__ methods of the elements for in-depth join
    for elem in same_id_elem:
        joined_col.append(elem + get_by_id(col_B, id=elem._id()))

    return joined_col
