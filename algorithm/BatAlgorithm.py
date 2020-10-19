from random import randint, gauss, uniform
from use.tools import kthperm
from algorithm.MultiObjective import MultiObjective as mo
from math import factorial, floor, exp
from typing import List, Tuple, Set
from models.Solution import Solution


class BatAlgorithm():
    def __init__(self, D, NP, N_Gen, NF, A, r, Alpha, Gama, Qmin, Qmax, scores):
        self.D = D  # number of Bats for each individual in the population
        self.NP = NP  # population size
        self.N_Gen = N_Gen  # generations number
        self.NF = NF  # fragments number
        # loudness of Bats for each solutions
        self.A = [A for i in range(self.NP)]
        # pulse rate of Bats for each solutions
        self.r = [r for i in range(self.NP)]
        self.A0 = A  # initial loudness
        self.r0 = r  # initial pulse
        self.Alpha = Alpha  # parametre of loudness calculation
        self.Gama = Gama  # parametre of pulse calculation

        self.Qmin = Qmin  # frequency min
        self.Qmax = Qmax  # frequency max
        self.Q = [0] * self.NP  # frequency of Bats for each solutions

        #self.scores = scores
        self.scores = [i for i in scores]
        self.l = [i for i in range(self.NF)]  # fragments index sequance

        self.min_index = 0  # the minimum index in lexecographie ordre
        # the maximum index in lexecographie ordre
        self.max_index = factorial(self.NF)-1

        self.x_best_pos = 0  # the position of the best individual in our solution
        self.v = [[0 for i in range(self.D)] for j in range(
            self.NP)]  # velocity of each Bats
        self.Sol = [[0 for i in range(self.D)]
                    for j in range(self.NP)]  # index of each Bats

        # intermediate population used in the non dominated sorting
        self.inter_Population = List[Solution]
        self.inter_Population = [Solution(kthperm(self.l, 0), generation=0)]
        self.Population, self.Positions = self.init_bat_population(
            NF, NP)  # the initial population

        for P in self.Population:
            P.oaf_objective(self.scores)
            P.odf_objective(self.scores)
        self.x_best = Solution("", generation=0)  # the best solution
        self.x_best = self.Update_solution(self.x_best, self.Population[0], 0)

    @staticmethod
    def init_bat_population(fragments_number: int, population_size: int) -> Tuple[List[Solution], List[int]]:
        """This function create the initial population for the Bat algorithm.

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
        list
            A list of position of the Solutions in the lexicographie ordre.
        """

        # STEP 2, generate initial population, and retreving the set of the solutions
        print("STEP-2 :: GENERATING SOLUTIONS (INITIAL POPULATION).")

        l = [i for i in range(fragments_number)]  # Our fragments, indexes

        solutions = list()  # The list of solutions.
        position = list()   # The list of positions.
        # we generate a uniforme distributed set of solution in a uniforme destributed solution space
        # the solution n is selected randomly from the interval [n*fragments_number!/population_size (n+1)*fragments_number!/population_size]
        equal_intervale = factorial(fragments_number)//population_size
        for i in range(population_size):
            combinaison = randint(
                (i)*equal_intervale, (i+1)*equal_intervale-1)
            sol = kthperm(l, combinaison)
            init_sol = Solution(sol, generation=0)
            solutions.append(init_sol)
            position.append(combinaison)
        return solutions, position

    def correct(self, x: int) -> int:
        """This function correct the index of the Bat to avoide 
        wrong index (out of range index).

        ...

        Parameters
        ----------
        x: int
            the initial index 

        Returns
        -------
        int
            the realisable index.

        """

        if x > self.max_index or x < self.min_index:
            return randint(0, self.max_index)
        return x

    @staticmethod
    def Update_solution(sol_1: Solution, sol_2: Solution, Generation: int) -> Solution:
        """This function update all attribute of the Solution sol_1
        by affecting value of attribute of sol_2.

        ...

        Parameters
        ----------
        sol_1: Solution
            the solution we want to update it.
        sol_2: Solution
            the new solution.      
        Generation: int
            generation of changes.

        Returns
        -------
        Solution
            the solution updated.

        """
        sol_1.genome = sol_2.genome
        sol_1.genome_size = sol_2.genome_size
        sol_1.generation = Generation
        sol_1.odf = sol_2.odf
        sol_1.oaf = sol_2.oaf
        return sol_1

    def best_bat(self, Generation: int):
        """This function search for the global optimum in the
        new solution by comparing Overlap Adjacent Fragment and
        Overlap Distant Fragment and number of contigue between 
        X_best and every individual in the new solution and 
        update X_best if an individual dominate it.

        ...

        Parameters
        ----------
        Generation: int
            the generation when X_best is updated.

        Returns
        -------


        """
        for i in range(self.NP):
            cond = mo.domination(self.Population[i], self.x_best)
            if cond == -1:
                self.x_best.contigs_number(self.scores)
                self.Population[i].contigs_number(self.scores)
                X_nc = self.x_best.contigs
                P_nc = self.Population[i].contigs
                if P_nc < X_nc:
                    self.x_best = self.Update_solution(
                        self.x_best, self.Population[i], Generation)
                    self.x_best_pos = self.Positions[i]
            elif cond == 0:
                self.x_best.contigs_number(self.scores)
                self.Population[i].contigs_number(self.scores)
                if self.Population[i].contigs < self.x_best.contigs:
                    self.x_best = self.Update_solution(
                        self.x_best, self.Population[i], Generation)
                    self.x_best_pos = self.Positions[i]

    def init_bat(self):
        """This function initialise the D Bats for every individual 
        in the initial solution and get the best NP individual after
         using the non dominate sorting algorithme.


        ...

        Parameters
        ----------

        Returns
        -------


        """
        # we generate uniforme distributed set of Bats around every solution in the population
        equal_intervale_sol = self.max_index//self.NP
        equal_intervale_bat = equal_intervale_sol//self.D
        equal_intervale = factorial(self.NF)//self.NP
        self.inter_Population = list()
        for i in range(self.NP):
            self.Q[i] = 0
            for j in range(self.D):
                rnd = uniform(0, 1)
                self.v[i][j] = 0.0
                self.Sol[i][j] = randint(
                    i*equal_intervale + j*equal_intervale_bat, i*equal_intervale + (j+1)*equal_intervale_bat-1)
                self.Sol[i][j] = self.correct(self.Sol[i][j])
                x = Solution(kthperm(self.l, self.Sol[i][j]), generation=0)
                x.oaf_objective(self.scores)
                x.odf_objective(self.scores)
                self.inter_Population.append(x)
        # we get NP first solutions from the K first front
        inter_population = mo.non_dominate_sorting(self.inter_Population)
        i = 0
        end = False
        for ip in inter_population:
            if end:
                break
            for f in ip:
                if i >= self.NP:
                    end = True
                    break
                Sol_i = f//self.D
                Sol_j = f % self.D
                self.Population[i] = self.Update_solution(
                    self.Population[i], self.inter_Population[f], 0)
                self.Positions[i] = self.Sol[Sol_i][Sol_j]
                i += 1

    def move_bat(self):
        """ We apply the Bat Algorithme to solve the DNA FAP

        ...

        Parameters
        ----------

        Returns
        -------


        """

        # STEP 3, generate initial bats, and retreving the set of the solutions
        print("STEP-3 :: GENERATING D BATS ARROUNG EACH SOLUTIONS (INITIAL BATS).")

        # we generate a uniforme distributed set of bats in equal intervals in the search space and compute ODF and OAF fitness
        self.init_bat()

        for t in range(self.N_Gen):

            # STEP 4, update Qi,Vi and Xi then move bats to generate a new local solution
            print("GENERATION :: {}".format(t))
            print(
                "\tG-{} --> STEP-4 :: GENERATING NEW SOLUTION AND UPDATING  Qi,Vi AND Xi PARAMETRES.".format(t))
            self.inter_Population = list()
            for i in range(self.NP):
                rnd = uniform(-1, 1)
                self.Q[i] = int(self.Qmin + (self.Qmax - self.Qmin) * rnd)
                for j in range(self.D):
                    self.v[i][j] = int(self.v[i][j]) + (self.Sol[i][j] -
                                                        self.x_best_pos) * self.Q[i]
                    self.Sol[i][j] = self.Sol[i][j] + int(self.v[i][j])
                    self.Sol[i][j] = self.correct(self.Sol[i][j])
                    x = Solution(kthperm(self.l, self.Sol[i][j]), generation=t)
                    x.oaf_objective(self.scores)
                    x.odf_objective(self.scores)
                    self.inter_Population.append(x)

            # STEP 5.1, compute ODF and OAF fitness and apply non dominated sorting
            print("\tG-{} --> STEP-5.1 :: APPLY NON DOMINATED SORTING TO GET BEST NP INDIVIDUAL FROM THE LOCAL SOLUTION.".format(t))
            inter_population = mo.non_dominate_sorting(self.inter_Population)

            # STEP 5.2, get the first NP solution as our new best population
            print("\tG-{} --> STEP-5.2 :: SELECT FIRST NP INDIVIDUAL FROM THE HIGHER FRONT TO BE OUR POPULATION .".format(t))
            i = 0
            end = False
            for ip in inter_population:
                if end:
                    break
                for f in ip:
                    if i >= self.NP:
                        end = True
                        break
                    Sol_i = f//self.D
                    Sol_j = f % self.D
                    self.Population[i] = self.Update_solution(
                        self.Population[i], self.inter_Population[f], t)
                    self.Positions[i] = self.Sol[Sol_i][Sol_j]
                    i += 1

            # STEP 6, select the global optimum in the solution to do a local search arround it
            print("\tG-{} --> STEP-6 :: SELECTING THE GLOBAL BEST POSITION.".format(t))
            self.best_bat(t)

            # STEP 7.1, generate a local solution arround the global optimum
            print("\tG-{} --> STEP-7.1 :: GENERATE A RANDOM NUMBER AND CREATE A LOCAL SOLUTION ARROND THE BEST SOLUTION.".format(t))
            for i in range(self.NP):
                rnd = uniform(0, 1)
                if rnd > self.r[i]:
                    for j in range(self.D):
                        new_pos = self.x_best_pos + \
                            self.A[i] // (gauss(-1, 1)**-(1))
                        new_pos = self.correct(new_pos)
                        x = Solution(kthperm(self.l, new_pos), generation=t)
                        x.oaf_objective(self.scores)
                        x.odf_objective(self.scores)

                        if mo.domination(self.Population[i], x) == -1:
                            self.Population[i] = self.Update_solution(
                                self.Population[i], x, t)
                            self.Positions[i] = new_pos

                    rnd = int(uniform(0, 100))
                    new_pos = self.Positions[i] + int(rnd)
                    new_pos = self.correct(new_pos)
                    x = Solution(kthperm(self.l, new_pos), generation=t)
                    x.oaf_objective(self.scores)
                    x.odf_objective(self.scores)

                    # STEP 7.2, if the random number generated < Ai we update Ai and ri
                    print(
                        "\tG-{} --> STEP-7.2 :: GENERATE A RANDOM NUMBER AND UPDATE Ai AND ri if it's < Ai.".format(t))
                    rnd = uniform(0, 1)
                    if rnd < self.A[i] and mo.domination(x, self.x_best) == -1:
                        self.Population[i] = self.Update_solution(
                            self.Population[i], x, t)
                        self.Positions[i] = new_pos
                        self.A[i] = self.A[i]*self.Alpha
                        self.r[i] = self.r0*(1-exp(-self.Gama*t))

        # STEP 8.1, select the global optimum in the final solution
        print("\tG-{} --> STEP-8.1 :: SELECTING THE GLOBAL BEST POSITION.".format(t))
        self.best_bat(t)

        # STEP 8.2, compute the contigue number of the global optimum in the final solution and print the individual
        print("\tG-{} --> STEP-8.2 :: SELECTING THE GLOBAL BEST POSITION.".format(t))
        print("\nSOLUTIONS::\n")
        self.x_best.contigs_number(self.scores)
        print(self.x_best)
        print("------------")
