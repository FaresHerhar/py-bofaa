from typing import List
from math import factorial, floor

from models.Fragment import Fragment


def max_possible_samples(n: int, k: int) -> int:
    """This function is for calculating the maximum number of possible
    combination with no repetition, while keeping different order
    (look for the law of combinations/permutation, in Probabilistic).

    ...

    Parameters
    ----------
    n: int
        The size of the set.
    k: int
        The wanted size, for the generated combinations.

    Rturns
    ------
    int
        The maximum number of possible combination.
    """

    if n < k:
        raise ValueError(
            "There is in error in your input, the size of the sample k,should not exeed the size of the space n.")

    return int(factorial(n)/factorial(n-k))


def read_fragments(file_name: str) -> List[Fragment]:
    """This function is for the initial step of the application,
    we read all the fragments, from the DNA sample file found.

    ...

    Parameters
    ----------
    file_name: str
        This is the DNA file's full path.

    Rturns
    ------
    list
        A list of Fragments.
    """

    fragments = list()

    file = open(file_name, "r")
    data = file.read().split("\n")
    file.close()

    # To understand better, check the benchmarks/test.data
    # to habe an idea about the file's format.
    # We need to parse the list by a step of 2,
    # to only get the sequences.
    count = 0
    for frag in data[1::2]:
        fragments.append(Fragment(frag, len(frag), count))
        count += 1

    return fragments


def kthperm(list_, k) -> List[int]:
    """This function is for calculating the combinaison of the 
    index number K in the lexecographic order of the vector S
    it is a non recursive version with complexité O(n)
    ...

    Parameters
    ----------
    list_: list
        The table of set to be permuted to get the combinaison of the k-th index.
    k: int
        The index from the lexicographie order of the combinations wanted.

    Rturns
    ------
    list
        The combination of the index k.
    """

    P = []
    while list_ != []:
        f = factorial(len(list_)-1)
        i = floor(k//f)
        if i > len(list_)-1:
            i = len(list_)-1
        x = list_[i]
        k = k % f
        P.append(x)
        list_ = list_[:i] + list_[i+1:]
    return P
