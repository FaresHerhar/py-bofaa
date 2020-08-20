from tools import *
from time import time


if __name__ == "__main__":
    CROSS_OVER_PROBABILITY = 0.8
    MUTATION_PROBABILITY = None
    GENERATIONS_NUMBER = None
    POPULATION_SIZE = 10
    OVECTIVE_FUNCTIONS_NUMBER = 2

    # STEP 0, reading fragments from file
    fragments = read_fragments("benchmarks/test.dat")
    # STEP 1, compute pair wise overlap
    scores = overlap_scores(fragments)
    # STEP 2, generate initial population
    population = init_population(len(fragments), POPULATION_SIZE)

    while True:
        # STEP 3, compute ODF and OAF fitness
        for sol in population:
            sol.oaf = oaf(sol, scores)
            sol.odf = odf(sol, scores)

        # STEP 4, calculate the fonts
        fonts = non_dominate_sorting(population)
        print(fonts)

        # STEP 5, calculate the crowding distances
        crownding = crowding_distance(population, fonts)

        # STEP 6, select solution for mutation
        selection = select_cross_solutions(population, CROSS_OVER_PROBABILITY)

        for i in population:
            print(i)
            print("---------------------------")
        break
