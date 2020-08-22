from tools import *
from time import time


if __name__ == "__main__":
    CROSS_OVER_PROBABILITY = 0.8
    MUTATION_PROBABILITY = None
    GENERATIONS_NUMBER = None
    POPULATION_SIZE = 100
    OVECTIVE_FUNCTIONS_NUMBER = 2

    start = time()
    # STEP 0, reading fragments from file
    print("STEP-0 :: READING FRAGMENTS.")
    fragments = read_fragments("benchmarks/dna-instances/x60189_4.dat")
    
    # STEP 1, compute pair wise overlap
    print("STEP-1 :: CALCULATING THE OVERLAP SCORES.")
    scores = overlap_scores(fragments)
    print(time() - start)

    # STEP 2, generate initial population, and retreving the set of the solutions
    print("STEP-2 :: GENERATING SOLUTIONS (INITIAL POPULATION).")
    population = init_population(len(fragments), POPULATION_SIZE)
    population, hash_values = population[0], population[-1]
    

    while True:
        # STEP 3, compute ODF and OAF fitness
        for sol in population:
            sol.oaf = oaf(sol, scores)
            sol.odf = odf(sol, scores)

        # STEP 4, calculate the fonts
        fonts = non_dominate_sorting(population)

        # STEP 5, calculate the crowding distances
        crownding = crowding_distance(population, fonts)

        # STEP 6, select solution for mutation
        selection = select_cross_solutions(population, CROSS_OVER_PROBABILITY)

        childs = cross_over(population, selection, hash_values)
        childs, hash_values = childs[0], childs[-1]

        for i in childs:
            print(i)
            print("---------------------------")
        break

    print("EXECTION TIME:: {} Seconds.".format(time() - start))