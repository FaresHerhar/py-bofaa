from typing import List, Tuple, Set
from random import sample, randint
from math import factorial

from models.Fragment import Fragment
from models.Solution import Solution


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


def longest_common(str_1: str, str_2: str) -> int:
    """This function Calculate, The Longest Common Subsequence
    between two sequences(in our case strings).

    ...

    Parameters
    ----------
    str_1: str
        The first sequence(string)
    str_1: str
        The second sequence(string)

    Rturns
    ------
    int
        An int, represeting the lenght of The Longest Common Substring.
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
            # In case of a match update, the result variable.
            elif str_1[i - 1] == str_2[j - 1]:
                H[i][j] = H[i - 1][j - 1] + 1
                result = max(result, H[i][j])
            # In case of a mismatch, just put zero.
            else:
                H[i][j] = 0

    return result


def overlap(frag_1: Fragment, frag_2: Fragment) -> int:
    """This function calculates the overlap of two fragments,
    in our case means the longest common nucleotides sequence,
    based on the same method of the longest common string
    of two string, since we're working on DNA sequences that are
    Strings based.

    ...

    Parameters
    ----------
    frag_1: Fragment
        The first fragment.
    frag_2: Fragment
        The second fragment.

    Returns
    -------
    int
       The longest common nucleotides sequence.
    """

    return longest_common(frag_1.sequence, frag_2.sequence)


def overlap_scores(fragments: List[Fragment]) -> List[List[int]]:
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
    Returns
    -------
    list
        A list of lists(matrix) of int, that contains the overlaping scores.
    """

    len_frag = len(fragments)

    # Creating the matrix
    scores = [[-1 for i in range(len_frag)] for j in range(len_frag)]

    for i in range(0, len_frag - 1):
        for j in range(i + 1, len_frag):
            # Since longest_common(a, b) == longest_common(b, a)
            # we don't have to calculate twice, therefore we do
            # the calculations once, and we assign twice.
            scores[i][j] = scores[j][i] = overlap(fragments[i], fragments[j])

    return scores


def init_population(fragments_number: int, population_size: int) -> Tuple[List[Solution], Set[int]]:
    """This function create the initial population for the NSGA-II algorithm.

    ...

    Parameters
    ----------
    fragments_number: int
        The number of the fragments, since our solution must be a combination
        of all the fragments, therefore the solution will be only a combination
        of indexes for these fragments, stored in another variable.
    population_size: int
        The size of the wanted initial population.

    Returns
    -------
    list
        A list of Solution.
    """

    # Check if the wanted population, is larger than
    # Maximum possible number of samples, see max_possible_samples()
    if population_size > max_possible_samples(fragments_number,
                                              fragments_number):
        raise ValueError(
            "The number of the population is bigger than the maximum possible number of samples")

    l = [i for i in range(fragments_number)]  # Our fragments, indexes

    # The hash values of each solution to avoid redundancy.
    hash_values = set()
    solutions = list()  # The list of solutions.

    count = 0  # To check if we've reached the number of the wanted population
    while count != population_size:
        sol = sample(l, fragments_number)
        # We can't calculate the hash value of mutable objects.
        hash_val = hash(tuple(sol))

        # Check if the solution already exists.
        if hash_val not in hash_values:
            hash_values.add(hash_val)
            solutions.append(Solution(sol, generation=1))
            count += 1

    return solutions, hash_values


def oaf(sol: Solution, scores: List[List[int]]) -> int:
    """It is  the first objective function, Overlaping Adjacent Fragments.

    ...

    Parameters
    ----------
    sol: Solution
        The solution that we want to calculate its first objective function.
    scores: list
        A list of list(matrix) of int, that contains the overlaping scores.

    Returns
    -------
    int
        The Overlaping Adjacent Fragments values, of the given solution.
    """

    oaf_value = 0
    sol_len = sol.genome_size
    temp_genome = sol.genome.copy()  # To avoid the over memotry access

    # Can't explain, take a look at the research paper(/papers)
    for i in range(sol_len - 1):
        oaf_value += scores[temp_genome[i]][temp_genome[i + 1]] * 2
    del temp_genome

    return oaf_value


def odf(sol: Solution, scores: List[List[int]]) -> int:
    """It is  the second objective function, Overlaping Distant Fragments.

    ...

    Parameters
    ----------
    sol: Solution
        The solution that we want to calculate its second objective function.
    scores: list
        A list of list(matrix) of int, that contains the overlaping scores.

    Returns
    -------
    int
        The Overlaping Distant Fragments values, of the given solution.
    """
    odf_value = 0
    sol_len = sol.genome_size
    temp_genome = sol.genome.copy()  # To avoid the over memotry access

    # Can't explain, take a look at the research paper(/papers)
    for i in range(sol_len - 2):
        p = i
        for j in range(i + 2, sol_len):
            odf_value += ((j - p) * scores[temp_genome[i]][temp_genome[j]]) * 2
    del temp_genome

    return odf_value


def domination(sol_1: Solution, sol_2: Solution) -> int:
    """This function test the dominance, between two solutions
    based on the two objective function listed above.
    It is for The Non-Domination Sorting Algorithm.

    ...

    Parameters
    ----------
    sol_1: Solution
        The first solution.
    sol_2: Solution
        The second solution.

    Returns
    -------
    int
        An int value,
        if -1 sol_1 dominates sol_2.
        if 1 sol_2 dominates sol_1.
        if 0 neither dominates the other.
    """
    if (sol_1.oaf >= sol_2.oaf) and (sol_1.odf < sol_2.odf):
        return -1
    elif (sol_2.oaf >= sol_1.oaf) and (sol_2.odf < sol_1.odf):
        return 1
    else:
        return 0


def non_dominate_sorting(population: List[Solution]) -> List[List[int]]:
    """This fonction is for the non dominate sorting for the NSGA-II Algorithm,
    The how the function works wont be explained here, look at the full Algorithm
    online, or take a look at the research paper.
    It returns the fonts, where each font is a list of integers, that represents
    the indexes of the solutions, in the population list.

    ...

    Parameters
    ----------
    population: list
        A list of solutions.


    Returns
    -------
    list
        A list of lists of integers.
    """
    # It contains all our fonts, in which each solutoin belongs to only one.
    fonts = [[]]
    # Each solution my dominate other solutions, in that case they are stores here.
    dominated = [[] for i in range(len(population))]
    # In case a solution gets dominated by another solution, we store how many times.
    dominate_number = [0 for i in range(len(population))]

    # STEP-1: Calculation the first Font.
    for p in range(len(population)):
        for q in range(len(population)):
            # Go to the domination function.
            d = domination(population[p], population[q])
            if d == -1:
                dominated[p].append(q)
            if d == 1:
                dominate_number[p] += 1
        if dominate_number[p] == 0:
            fonts[0].append(p)
            # set the rank to one, since it belongs to the first font.
            population[p].rank = 1

    # STEP-2: Calculating the rest of the Fonts.
    fonts_counter = 1

    while fonts_counter <= len(fonts):
        for p in fonts[fonts_counter - 1]:
            for q in dominated[p]:
                dominate_number[q] -= 1

                # if the solution is no longer dominated, add it to the next Front.
                if dominate_number[q] == 0:
                    # set the rank to one, since it belongs to the next font.
                    population[q].rank = fonts_counter + 1
                    if len(fonts) == fonts_counter:
                        fonts.append([q])
                    else:
                        fonts[-1].append(q)
        fonts_counter += 1

    return fonts


def crowding_distance(population: List[Solution], fonts: List[List[int]]) -> List[float]:
    """This function is for calculation the crowding distance of each solution.
    at first we sort the solution, for each front based on oaf values,
    and odf values separately, while odf ascendant dort, and oaf is descendant sort.
    After that, we calculate the crownding distance for each solution,
    with oaf and odf sorts, and exluding the limits of each sorted font(first and last element),
    we set them to -1. Than, we chack if the values match the excluded and sum them.
    It returns a list of floats, that represents the crowding distances for each solution.

    ...

    Parameters
    ----------
    population: list
        A list of solutions.
    fonts: list
        A list of Fonts(list of int).


    Returns
    -------
    list
        A list of integers, that represents the crowding distance of each solution.
    """
    # Each fonts, will be sorted first by oaf, then odf
    # while maximize oaf, and minimize odf
    oaf_sorted_fonts = []
    odf_sorted_fonts = []

    # The crowding will be calculated, by oaf, odf separately, than summed
    oaf_cwrowding = [0 for i in range(len(population))]
    odf_cwrowding = [0 for i in range(len(population))]
    crowding = [0.0 for i in range(len(population))]

    # STEP-0 sorting
    for p in fonts:
        if len(p) == 1:
            oaf_sorted_fonts.append(p)
            odf_sorted_fonts.append(p)
        if len(p) > 1:
            temp = [population[q] for q in p]

            # Sort by oaf, reversed for acendant
            temp.sort(key=lambda sol: sol.oaf, reverse=True)
            oaf_sorted_fonts.append([population.index(q) for q in temp])
            # Sort by odf, not reversed for desandant
            temp.sort(key=lambda sol: sol.odf)
            odf_sorted_fonts.append([population.index(q) for q in temp])

    # STEP-1 calculate crowding distances
    # STEP-1.1 for oaf
    for p in oaf_sorted_fonts:
        oaf_cwrowding[p[0]] = -1
        # The fraction of max - min
        kill = population[p[0]].oaf - \
            population[p[-1]].oaf
        if len(p) > 1:
            temp_oaf = 0
            oaf_cwrowding[p[-1]] = -1
            # calculating the oaf crowding for each element of the font, excluding limits
            for q in range(1, len(p) - 1):
                temp_oaf += (population[p[q - 1]].oaf -
                             population[p[q + 1]].oaf) / kill
                oaf_cwrowding[p[q]] = temp_oaf

    # STEP-1.2 for odf
    for p in odf_sorted_fonts:
        odf_cwrowding[p[0]] = -1
        # The fraction of (max - min)
        kill = population[p[-1]].odf - \
            population[p[0]].odf
        if len(p) > 1:
            temp_odf = 0
            odf_cwrowding[p[-1]] = -1
            # calculating the oaf crowding for each element of the font, excluding limits
            for q in range(1, len(p) - 1):
                temp_odf += (population[p[q + 1]].odf -
                             population[p[q - 1]].odf) / kill
                odf_cwrowding[p[q]] = temp_odf

    # STEP-2 assembling crowding distances
    # since we sort twice, for odf and oaf, so some solution get exluded
    # in one sort, and not in another, there for we will check
    for index in range(len(population)):
        if oaf_cwrowding[index] == -1 or odf_cwrowding[index] == -1:
            crowding[index] = float(-1)
            population[index].crowding_distance = -float(-1)
        else:
            crowding[index] = float(
                oaf_cwrowding[index] + odf_cwrowding[index])
            population[index].crowding_distance = float(oaf_cwrowding[index] +
                                                        odf_cwrowding[index])

    return crowding


def select_cross_solutions(population: List[Solution], cross_over_propability: int) -> List[int]:
    """This function use the binary tournament selection method to select the layouts
    for crossover and mutation i.e to generate the parent population PP.
    It uses a tousize of 2, i.e selecting two solution randomly, comparing then picking,
    untill the number of the solutions it equals to the poolsize, where the
    poolsize is (cross_over_propability * the size of the population.)

    ...

    Parameters
    ----------
    cross_over_probability: int
        The probability of cross over.
    population: list
        A list of solutions.


    Returns
    -------
    list
        A list of integers, that represents the selected solution for mutation.
    """
    population_size = len(population)
    pool_size = round(population_size * cross_over_propability)
    selection_counter = 0
    selection = list()

    while selection_counter != pool_size:
        # Randomly select two solution from the population, and
        # since we're using indexes, its easier to use integers.
        first_selection = randint(0, population_size - 1)
        second_selection = randint(0, population_size - 1)

        # if the rank is not the same take the one with the less rank.
        if population[first_selection].rank != population[second_selection].rank:
            if population[first_selection].rank < population[second_selection].rank:
                selection.append(first_selection)

            if population[first_selection].rank > population[second_selection].rank:
                selection.append(second_selection)

            selection_counter += 1
        # add the one with gratter crowding distance
        else:
            if population[first_selection].crowding_distance <= population[second_selection].crowding_distance:
                selection.append(first_selection)

            if population[first_selection].crowding_distance >= population[second_selection].crowding_distance:
                selection.append(second_selection)

            selection_counter += 1

    return selection
