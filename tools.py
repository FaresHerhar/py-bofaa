from typing import List
from math import factorial

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
