from typing import List

from models.Fragment import Fragment


def waterman_algorithm(str_1: str, str_2: str, match_score: int, mismatch_score: int, gap_cost: int) -> float:
    """Take a look at The smith waterman algorithm.

    ...

    Parameters
    ----------
    str_1: str
        The first sequence(string)
    str_1: str
        The second sequence(string
    match_score: int
        A positive int, we add in case to caracters match.
    mismatch_score: int
        A negative int, we add in case to caracters don't match.
    gap_cost: int
        A negative int, for the gap.

    Rturns
    ------
    float
        An float, represeting the lenght of The Longest Common Substring.
    """

    len_1 = len(str_1)
    len_2 = len(str_2)
    # The scoring result
    result = 0

    # Make a mtrix, to store the matches,
    # and mismatches, and to calculate the final result.
    # Both the row, and column number 0, are skiped.
    # Check The Longest Common Subsequence problem.
    H = [[0 for k in range(len_2 + 1)] for l in range(len_1 + 1)]

    for i in range(len_1 + 1):
        for j in range(len_2 + 1):
            # skip this one because it doesn't represent any thing.
            if (i == 0) or (j == 0):
                H[i][j] = 0
            else:
                match = H[i - 1][j - 1] + \
                    (match_score if str_1[i - 1] ==
                     str_2[j - 1] else + mismatch_score)
                delete = H[i - 1][j] + gap_cost
                insert = H[i][j - 1] + gap_cost
                result = H[i][j] = max(match, delete, insert, 0)

    return result


def overlap(frag_1: Fragment, frag_2: Fragment, match_score: int, mismatch_score: int, gap_cost: int) -> float:
    """This function calculates the overlap of two fragments,
    in our case means the longest common nucleotides sequence,
    based on the same method of the waterman_algorithm.
    Take a look at The smith waterman algorithm.

    ...

    Parameters
    ----------
    frag_1: Fragment
        The first fragment.
    frag_2: Fragment
        The second fragment.
    match_score: int
        A positive int, we add in case to caracters match.
    mismatch_score: int
        A negative int, we add in case to caracters don't match.
    gap_cost: int
        A negative int, for the gap.

    Returns
    -------
    float
       The value of the walterman algorith applied on the fragments sequences.
    """
    return waterman_algorithm(frag_1.sequence, frag_2.sequence, match_score, mismatch_score, gap_cost)


def overlap_scores(fragments: List[Fragment], match_score: int, mismatch_score: int, gap_cost: int) -> List[List[float]]:
    """This function calculates the overlap scores between each fragment
    and another, then the scores are all stored in a matrix.

    PS:
    an overlap score between a fragment and itself is set to -1,
    because such a thing doesn't exist.

    ...


    Parameters
    ----------
    fragments: list
        A list of fragments.
    match_score: int
        A positive int, we add in case to caracters match.
    mismatch_score: int
        A negative int, we add in case to caracters don't match.
    gap_cost: int
        A negative int, for the gap.
    Returns
    -------
    list
        A list of lists(matrix) of float, that contains the overlaping scores.
    """

    len_frag = len(fragments)

    # Creating the matrix
    scores = [[-1.0 for i in range(len_frag)] for j in range(len_frag)]

    for i in range(0, len_frag - 1):
        for j in range(i + 1, len_frag):
            # Since waterman_algorithm(a, b) == waterman_algorithm(b, a)
            # we don't have to calculate twice, therefore we do
            # the calculations once, and we assign twice.
            scores[i][j] = scores[j][i] = overlap(
                fragments[i], fragments[j], match_score, mismatch_score, gap_cost)

    return scores
