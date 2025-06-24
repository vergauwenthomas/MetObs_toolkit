from typing import Iterable, List, Any

""" Set comparison logic for Metobs functions.

    The uniqueness is assumed if the ._id() of the instances differs.
"""


#Note: another design would be to implement the __eq__ and __ne__ methods in 
#the metobs classes, so that set() operations could be used. The issue with that 
#is that np.nan != np.nan. So 'identical gaps' are then assumed to be non-equal.


def metobs_intersection(A: Iterable[Any], B: Iterable[Any]) -> List[Any]:
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
    A_ids = set([a._id() for a in A])
    B_ids = set([b._id() for b in B])
    #Intersection of ID's
    a_not_in_b_ids= A_ids - B_ids
    #Revert to the object-space
    return [a for a in A if a._id() in a_not_in_b_ids]

   
def get_by_id(A: Iterable[Any], id: str) -> Any:
    candidates = [a for a in A if a._id() == id]
    if len(candidates) > 1:
        raise IndexError(f'More than one element found with id={id} in the collection: {A}')
    if len(candidates) == 0:
        raise IndexError(f'No element found with id={id} in the collection: {A}')
    return candidates[0]