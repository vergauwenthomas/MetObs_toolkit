from typing import Iterable, List, Any

""" Set comparison logic for Metobs functions.

    The uniqueness is assumed if the ._id() of the instances differs.
"""


#Note: another design would be to implement the __eq__ and __ne__ methods in 
#the metobs classes, so that set() operations could be used. The issue with that 
#is that np.nan != np.nan. So 'identical gaps' are then assumed to be non-equal.


def metobs_intersection(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:

    #NOTE elements returned originates from A
    # transform to id-space (these are strings)
    A_ids = set([a._id() for a in A])
    B_ids = set([b._id() for b in B])
    # Intersection of ID's
    intersection_ids = A_ids.intersection(B_ids)
    # Revert to the object-space
    intersection = [a for a in A if a._id() in intersection_ids]
    return intersection

def metobs_B_not_in_A(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:
    #transform to id-space (these are strings)
    A_ids = set([a._id() for a in A])
    B_ids = set([b._id() for b in B])
    #Intersection of ID's
    b_not_in_a_ids= B_ids - A_ids
    #Revert to the object-space
    return [b for b in B if b._id() in b_not_in_a_ids]

def metobs_A_not_in_B(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:
    #transform to id-space (these are strings)
    return metobs_B_not_in_A(A=B, B=A)
   
def get_by_id(A: Iterable[Any], id: str) -> Any:
    candidates = [a for a in A if a._id() == id]
    if len(candidates) > 1:
        raise IndexError(f'More than one element found with id={id} in the collection: {A}')
    if len(candidates) == 0:
        raise IndexError(f'No element found with id={id} in the collection: {A}')
    return candidates[0]


def join_collections(col_A: Iterable[any], col_B:Iterable[any]) -> list[any]:
    """ join two collections, in dept addition applied for same-id-elements.

    Parameters
    ----------
    col_A : Iterable[any]
        Collection of elemlents with a '_id()' method.
    col_B : Iterable[any]
        Collection of elemlents with a '_id()' method.

    Returns
    -------
    list[any]
        The joint collection. This collection hold element without duplicated
        id's.
    """

    joined_col = []
    #trigage elements using metobs-logics
    same_id_elem = metobs_intersection(
            A=list(col_A),              
            B=list(col_B))
    #Note: same_id_elem originates form col_A!
    
    unique_A_elem = metobs_A_not_in_B(
            A=list(col_A),              
            B=list(col_B))
    joined_col.extend(unique_A_elem) #simple join
    
    unique_B_elem = metobs_B_not_in_A(
            A=list(col_A),              
            B=list(col_B))
    joined_col.extend(unique_B_elem) #simple join

    #Rely on the __add__ methods of the elements for indept-join
    for elem in same_id_elem:
        joined_col.append(elem + get_by_id(col_B, id=elem._id()))
    
    return joined_col
        