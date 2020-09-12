benchmarks = [
    "benchmarks/test.dat",
    "benchmarks/dna-instances/x60189_4.dat",
    "benchmarks/dna-instances/x60189_5.dat",
    "benchmarks/dna-instances/x60189_6.dat",
    "benchmarks/dna-instances/x60189_7.dat",
    "benchmarks/dna-instances/m15421_5.dat",
    "benchmarks/dna-instances/m15421_6.dat",
    "benchmarks/dna-instances/m15421_7.dat",
    "benchmarks/dna-instances/38524243_4.dat",
    "benchmarks/dna-instances/38524243_7.dat",
    "benchmarks/dna-instances/j02459_7.dat"
]

# Variables for the NSGA-II Algorithm
CROSS_OVER_PROBABILITY = 0.8
MUTATION_PROBABILITY = 0.8
GENERATIONS_NUMBER = 2000
NSGA_POPULATION_SIZE = 100
OVECTIVE_FUNCTIONS_NUMBER = 2
MATCH_SCORE = 1
MISMATCH_SCORE = -1
GAP_COST = -1.33

# Variables for the MOBA Algorithm
DIMENTION_NUMBER = 20
MOBA_POPULATION_SIZE = 10
LOUDNESS = 9
RATE_PLUSSE = 2
ALPHA = 0.9
GAMA = 0.9
MINIMUM_FREQUANCY = 0
MAXIMUM_FREQUANCY = 15

BECHMARK_FILE = benchmarks[1]
