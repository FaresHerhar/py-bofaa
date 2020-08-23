from typing import List


class Solution:
    """This is a Solution class that represent a Solution for the
    dna sequencing problem.

    ...

    Attributes
    ----------
    genome: list
        The list of the Fragments that represents a solution.
    genome_size: int
        The lenght of the genome(the solution array).
    oaf: int
        The value of the first objective function, overlaping adjacent
        fragments.
    odf: int
        The value of the second objective function, overlaping distinct
        fragments.
    generation: int
        The number(index) of the generation, that can tell at what
        generation the solution was created.
    rank: int
        The rank if the solution, i.e the fonts that the solution belongs to,
        and it changes from a generation to another.
    crowding_distance: float, optional
        The crowding distance of the solution that changes, from a generation
        to a generation, based on the font that it belongs to.
    contigs: int
        The number of contigs.

    Methods
    -------
    __init__(oaf=-1, odf=-1, generation=-1): None
        The constructor.
    __str__: str
        The print formating method.
    """

    def __init__(self, genome, generation=-1):
        """The constructor.

        ...

        Parameters
        ----------
        genome: list
            The list of the Fragments that represents a solution.
        genration: int, optional
            The number(index) of the generation, that can tell at what
            generation the solution was created. Durring the run time,
            once the Solution is computed, the generation will be assigned.

        Rturns
        ------
        None
        """

        self.genome = genome
        self.genome_size = len(self.genome)
        self.generation = generation
        self.oaf = -1
        self.odf = -1
        self.rank = -1
        self.crowding_distance = -1
        self.contigs = -1

    def __str__(self):
        """This method returns the formating print format, to print out
        a solution, while all the details all printed.
        """

        out = "* Genome:: {}\n* Genome size:: {}\n* OAF::{}\n* ODF:: {}\n* Rank:: {}\n* Crowding distance:: {}\n* Contigs number:: {}\n* Generation:: {}"
        return out.format(self.genome, self.genome_size, self.oaf, self.odf, self.rank, self.crowding_distance, self.contigs, self.generation)

    def oaf_objective(self, scores: List[List[int]]) -> None:
        """It is  the first objective function, Overlaping Adjacent Fragments.

        ...

        Parameters
        ----------
        scores: list
            A list of list(matrix) of int, that contains the overlaping scores.

        Returns
        -------
        None
        """

        self.oaf = 0
        # Can't explain, take a look at the research paper(/papers)
        for i in range(self.genome_size - 1):
            self.oaf += scores[self.genome[i]][self.genome[i + 1]] * 2

    def odf_objective(self, scores: List[List[int]]) -> None:
        """It is  the second objective function, Overlaping Distant Fragments.

        ...

        Parameters
        ----------
        scores: list
            A list of list(matrix) of int, that contains the overlaping scores.

        Returns
        -------
            None
        """

        self.odf = 0
        # Can't explain, take a look at the research paper(/papers)
        for i in range(self.genome_size - 2):
            p = i
            for j in range(i + 2, self.genome_size):
                self.odf += ((j - p) *
                             scores[self.genome[i]][self.genome[j]]) * 2

    def contigs_number(self, scores: List[List[int]]) -> None:
        """It is the function that calculates the number of contigs in a solution.
        For a genome of N fragments, initially NC=1 and for i=1, 2, ...,N-1,
        if score[i, i+1]=0 then NC=NC+1.

        ...

        Parameters
        ----------
        scores: list
            A list of list(matrix) of int, that contains the overlaping scores.

        Returns
        -------
        None
        """
        self.contigs = 1

        # If a score between to fragments, is less than a score condition calulated
        # we increment the number of the contigs
        for index in range(0, self.genome_size - 1):
            if scores[self.genome[index]][self.genome[index + 1]] == 0:
                self.contigs += 1
