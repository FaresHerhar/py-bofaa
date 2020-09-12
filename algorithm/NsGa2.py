from typing import List, Set
from random import sample, randint, random

from models.Solution import Solution


class NsGa2:

    @staticmethod
    def select_cross_solutions(population: List[Solution], cross_over_propability: int) -> List[int]:
        """This function use the binary tournament selection method to select the layouts
        for crossover and mutation i.e to generate the parent population PP.
        It uses a tousize of 2, i.e selecting two solution randomly, comparing then picking,
        untill the number of the solutions it equals to the poolsize, where the
        poolsize is (cross_over_propability * the size of the population.)

        ...

        Parameters
        ----------
        cross_over_probability: float
            The probability of operating a crossover.
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

        # To always make sure we have faircross pairs.
        if (pool_size % 2) != 0:
            pool_size -= 1

        while selection_counter != pool_size:
            # Randomly select two solution from the population, and
            # since we're using indexes, its easier to use integers.
            first_selection = randint(0, population_size - 1)
            second_selection = sample(
                [i for i in range(population_size) if i != first_selection], 1)[0]

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

    @staticmethod
    def crossover(population: List[Solution], selection: List[int], hash_values: Set[int], generation_counter: int) -> List[Solution]:
        """This function if for operating the cross over operation on the selection pool solutions
        in order to for new child solution from two parents.
        The functions uses the double point crossover, while checking the validity and the existance
        of the solution before.

        ...

        Parameters
        ----------
        population: list
            A list of solutions.
        selectoin: list
            A list of int(solutions indexes) for the selection pool.
        hash_values: set
            A set of int, that contains the hash values of the already exists solutions.
        generation_counter: int
            An int, that represents the generation number that the solution will be created at.


        Returns
        -------
        list
            A list of solutions, that represents the new created solotions from the crossover process.
        """
        # carry the cross over childs.
        cross_childs = list()
        # The number of the fragments, whch equals to the lenght of the solution.
        g_len = len(population[0].genome)

        # Parcour the solutions by pair, step equels to 2.
        for p in range(0, len(selection) - 1, 2):
            # Copyt the solution picked genomes, to avoid mutability damage.
            p_1 = population[selection[p]].genome.copy()
            p_2 = population[selection[p + 1]].genome.copy()

            # Generate to random points, and make sure they are not equal.
            point_1 = randint(0, g_len - 1)
            point_2 = sample([i for i in range(g_len) if i != point_1], 1)[0]

            # Make sure the first point is less than the second one.
            if point_1 > point_2:
                point_1, point_2 = point_2, point_1

            # Create children.
            c_1 = p_1[:point_1]+p_2[point_1:point_2]+p_1[point_2:]
            c_2 = p_2[:point_1]+p_1[point_1:point_2]+p_2[point_2:]

            for sol in [c_1, c_2]:
                # Check if it is a valid solution.
                if len(set(sol)) == g_len:
                    # We can't calculate the hash value of mutable objects.
                    hash_val = hash(tuple(sol))

                    # Check if the solution already exists.
                    if hash_val not in hash_values:
                        hash_values.add(hash_val)
                        cross_childs.append(
                            Solution(sol, generation=generation_counter))

        # we should return the hash values to update it, in the main function
        return cross_childs

    @staticmethod
    def mutation(population: List[Solution], selection: List[int], hash_values: Set[int], mutation_probability: float, generation_counter: int) -> List[Solution]:
        """
        This function if for operating the mutation operation on the selection pool solutions
        in order to for new child solution from mutating one parent.
        The functions uses the swaping mutation type, we randomly pick two point, and swap them;
        while checking the validity and the existance of the solution before.

        ...

        Parameters
        ----------
        population: list
            A list of solutions.
        selectoin: list
            A list of int(solutions indexes) for the selection pool.
        hash_values: set
            A set of int, that contains the hash values of the already exists solutions.
        mutation_probability: float
            A float, to determine the mutation rate.
        generation_counter: int
            An int, that represents the generation number that the solution will created at.


        Returns
        -------
        list
            A list of solutions, that represents the new created solotions from the crossover process.
        """
        # carry the cross over childs.
        mutation_childs = list()
        # The number of the fragments, whch equals to the lenght of the solution.
        g_len = len(population[0].genome)

        for p in selection:
            if random() < mutation_probability:
                # Copy the solution picked genomes, to avoid mutability damage.
                sol = population[p].genome.copy()

                # Generate to random points, and make sure they are not equal.
                point_1 = randint(0, g_len - 1)
                point_2 = sample(
                    [i for i in range(g_len) if i != point_1], 1)[0]

                # Do swap mutation.
                sol[point_2], sol[point_1] = sol[point_1], sol[point_2]
                hash_val = hash(tuple(sol))

                # Check if the solution already exists.
                if hash_val not in hash_values:
                    hash_values.add(hash_val)
                    mutation_childs.append(
                        Solution(sol, generation=generation_counter))

        return mutation_childs
