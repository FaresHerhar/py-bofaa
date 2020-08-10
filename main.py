from typing import List
from random import sample
from math import factorial

from models.Fragment import Fragment
from models.Solution import Solution


def max_possible_samples(n, k):
    """
    the max number of the solutions with no repetition and no order switching
    """
    if n < k:
        raise "There is in error in your input, the size of the sample k,should not exeed the size of the space n."
    return factorial(n)/factorial(n-k)


def read_fragments(file_name: str) -> List[Fragment]:
    fragments = list()

    file = open(file_name, "r")
    data = file.read().split("\n")
    file.close()

    count = 0
    for frag in data[1::2]:
        fragments.append(Fragment(frag, len(frag), count))
        count += 1

    return fragments


def longest_common(str_1: str, str_2: str) -> int:
    len_1 = len(str_1)
    len_2 = len(str_2)
    result = 0

    H = [[0 for k in range(len_2 + 1)] for l in range(len_1 + 1)]

    for i in range(len_1 + 1):
        for j in range(len_2 + 1):
            if (i == 0) or (j == 0):
                H[i][j] = 0
            elif str_1[i - 1] == str_2[j - 1]:
                H[i][j] = H[i - 1][j - 1] + 1
                result = max(result, H[i][j])
            else:
                H[i][j] = 0

    return result


def overlap(frag_1: Fragment, frag_2: Fragment) -> int:
    return longest_common(frag_1.sequence, frag_2.sequence)


def overlap_scores(fragments: List[Fragment]) -> List[List[int]]:
    len_frag = len(fragments)
    scores = [[-1 for i in range(len_frag)] for j in range(len_frag)]
    score = 0

    for i in range(0, len_frag - 1):
        for j in range(i + 1, len_frag):
            score = overlap(fragments[i], fragments[j])
            scores[i][j] = scores[j][i] = score

    return scores


def init_population(fragments_number: int, population_size: int) -> List[List[int]]:
    l = [i for i in range(fragments_number)]
    if population_size > max_possible_samples(fragments_number,
                                              fragments_number):
        raise "The number of the population is bigger than the maximum possible number of samples"

    hash_values = set()
    solutions = list()

    count = 0
    while count != population_size:
        sol = sample(l, fragments_number)
        hash_val = hash(tuple(sol))

        if hash_val not in hash_values:
            hash_values.add(hash_val)
            solutions.append(Solution(sol))
            count += 1

    return solutions


def oaf(sol: Solution, scores: List[List[int]]) -> int:
    oaf_value = 0
    sol_len = sol.genome_size
    temp_genome = sol.genome.copy()

    for i in range(sol_len - 1):
        oaf_value += scores[temp_genome[i]][temp_genome[i + 1]] * 2
    del temp_genome

    return oaf_value


def odf(sol: Solution, scores: List[List[int]]) -> int:
    odf_value = 0
    sol_len = sol.genome_size
    temp_genome = sol.genome.copy()

    for i in range(sol_len - 2):
        p = i
        for j in range(i + 2, sol_len):
            odf_value += ((j - p) * scores[temp_genome[i]][temp_genome[j]]) * 2
    del temp_genome

    return odf_value


if __name__ == "__main__":
    fragments = read_fragments("benchmarks/test")
    scores = overlap_scores(fragments)
    sol = Solution([5, 2, 6, 7, 4, 0, 1, 3, 8, 9])
    print(oaf(sol, scores))
    print(odf(sol, scores))
