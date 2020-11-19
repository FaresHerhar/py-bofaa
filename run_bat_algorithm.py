from time import time

from config import *
from use.tools import read_fragments
from use.scoring import *
from algorithm.BatAlgorithm import *


if __name__ == "__main__":
    # Counting the number of generations
    generation_counter = 1

    start = time()
    # STEP 0, reading fragments from file
    print("STEP-0 :: READING FRAGMENTS FROM FILE --> {}".format(BECHMARK_FILE))
    fragments = read_fragments(BECHMARK_FILE)
    print(len(fragments))
    print("--------------------------------")

    # STEP 1, compute pair wise overlap
    print("STEP-1 :: CALCULATING THE OVERLAP SCORES.")
    scores = overlap_scores(fragments, MATCH_SCORE, MISMATCH_SCORE, GAP_COST)
    # print(scores)
    Algorithm = BatAlgorithm(DIMENTION_NUMBER, MOBA_POPULATION_SIZE, GENERATIONS_NUMBER, len(fragments),
                             LOUDNESS, RATE_PLUSSE, ALPHA, GAMA, MINIMUM_FREQUANCY, MAXIMUM_FREQUANCY, scores)
    Algorithm.move_bat()

    print("DONE.\n")
    print("EXECTION TIME:: {} Seconds.".format(round(time() - start)))
