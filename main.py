from tools import *
from time import time


if __name__ == "__main__":
    # STEP 0, reading fragments from file
    fragments = read_fragments("benchmarks/test.dat")
    # STEP 1, compute pair wise overlap
    scores = overlap_scores(fragments)
    # STEP 2, generate initial population
    population = init_population(len(fragments), 20)

    while True:
        # STEP 3, compute ODF and OAF fitness
        for sol in population:
            sol.oaf = oaf(sol, scores)
            sol.odf = odf(sol, scores)

        # STEP 4, calculate the fonts
        fonts = non_dominate_sorting(population)

        # STEP 5, calculate the crowding distances
        crownding = crowding_distance(population, fonts)
        print(crownding)
        
        break
