from time import time

from tools import read_fragments
from MultiObjective import MultiObjective as mo
from NsGa2 import NsGa2 as nsga2
from config import *
from scoring import *


def run_nsga2(benchmark_file: str) -> None:
    print("USING THE NSGA-II Algorithm.")
    # Counting the number of generations
    generation_counter = 1

    start = time()
    # STEP 0, reading fragments from file
    print("STEP-0 :: READING FRAGMENTS FROM FILE --> {}".format(benchmark_file))
    fragments = read_fragments(benchmark_file)

    # STEP 1, compute pair wise overlap
    print("STEP-1 :: CALCULATING THE OVERLAP SCORES.")
    scores = overlap_scores(fragments, MATCH_SCORE, MISMATCH_SCORE, GAP_COST)

    # STEP 2, generate initial population, and retreving the set of the solutions
    print("STEP-2 :: GENERATING SOLUTIONS (INITIAL POPULATION).")
    population, hash_values = mo.init_population(
        len(fragments), POPULATION_SIZE)

    print("STEP-3 :: CALCULATING OBJECTIVE FUNCTIONS.")
    # STEP 3, compute ODF and OAF fitness
    for sol in population:
        sol.oaf_objective(scores)
        sol.odf_objective(scores)

    print("STEP-4 :: CALCULATING AND ATTRIBUTING FONTS.")
    # STEP 4, calculate the fonts
    fonts = mo.non_dominate_sorting(population)

    print("STEP-5 :: CALCULATING CROWDING DISTANCES.")
    # STEP 5, calculate the crowding distances
    crownding = mo.crowding_distance(population, fonts)

    while generation_counter <= GENERATIONS_NUMBER:
        print("GENERATION :: {}".format(generation_counter))

        print("\tG-{} --> STEP-6 :: SELECTING SOLUTIONS POOL.".format(generation_counter))
        # STEP 6, select solutions for pool
        selection = nsga2.select_cross_solutions(
            population, CROSS_OVER_PROBABILITY)

        # STEP 7, crossover and mutation
        print("\tG-{} --> STEP-7.1 :: OPERATING CROSSOVER.".format(generation_counter))
        # STEP 7.1, crossover
        childs = nsga2.crossover(population, selection,
                                 hash_values, generation_counter)

        print("\tG-{} --> STEP-7.2 :: OPERATING MUTATION.".format(generation_counter))
        # STEP 7.2, mutation
        childs.extend(nsga2.mutation(population, selection,
                                     hash_values, MUTATION_PROBABILITY, generation_counter))

        # STEP 8, offsoring
        print("\tG-{} --> STEP-8.1 :: CALCULATING OBJECTIVE FUNCTIONS FOR CHILDS.".format(generation_counter))
        # STEP 8.1, calculate oaf, odf to the childs
        for child in childs:
            child.oaf_objective(scores)
            child.odf_objective(scores)

        # STEP 8.2, merge child with current population
        print(
            "\tG-{} --> STEP-8.2 :: CREATING OFFSPRING POPULATION.".format(generation_counter))
        # to create offspring population
        population.extend(childs)

        print("\tG-{} --> STEP-8.3 :: CALCULATING AND ATTRIBUTING FONTS FOR THE OFFSPRING POPULATION.".format(generation_counter))
        # STEP 8.3, recalculate the fonts for the offspring population
        fonts = mo.non_dominate_sorting(population)

        print("\tG-{} --> STEP-8.4 :: CALCULATING CROWDING DISTANCES FOR THE OFFSPRING POPULATION.".format(generation_counter))
        # STEP 8.4, recalculate the crowding distances for the offspring population
        crownding = mo.crowding_distance(population, fonts)

        print("\tG-{} --> STEP-9 :: PASSING THE FIRST {} OFFSSPRING SOLUTION THE NEXT GENERATION POPULATION.".format(
            generation_counter, GENERATIONS_NUMBER))
        # STEP 9, passing the next first POPULATION_SIZE solutions
        temp = [index for indexes in fonts for index in indexes]
        population = [population[p] for p in temp[:POPULATION_SIZE]]

        generation_counter += 1

    # Gtting the somution
    print("\nSOLUTIONS::\n")
    # Take only the elemnts of the first fonts i.e rnak=1
    population = [p for p in population if p.rank == 1]

    # Calculate the number of contigs
    for p in population:
        p.contigs_number(scores)

    # Sort population by the number of contigs
    population.sort(key=lambda x: x.contigs)

    out = "* Genome:: {}\n* Genome size:: {}\n* OAF::{}\n* ODF:: {}\n* Rank:: {}\n* Crowding distance:: {}\n* Contigs number:: {}\n* Generation:: {}"
    
    # Print the solution
    for p in population:
        file = open("/home/fares/lol.txt", "a+")
        file.write(out.format(p.genome, p.genome_size, p.oaf, p.odf,
                              p.rank, p.crowding_distance, p.contigs, p.generation))
        file.close()
        break

        print(p)
        print("------------")

    print("DONE.\n")
    print("EXECTION TIME:: {} Seconds.".format(round(time() - start)))


if __name__ == "__main__":
    run_nsga2(BECHMARK_FILE)
