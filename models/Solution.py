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

    Methods
    -------
    __init__(oaf=-1, odf=-1, generation=-1): None
        The constructor.
    __str__: str
        The print formating method.
    """

    def __init__(self, genome, oaf=-1, odf=-1, generation=-1):
        """The constructor

        ...

        Parameters
        ----------
        genome: list
            The list of the Fragments that represents a solution.
        oaf: int, optional
            The value of the first objective function, overlaping adjacent
            fragments. Durring the run time this variable will be
            assigned by the first objective function.
        odf: int, optional
            The value of the second objective function, overlaping distinct
            fragments.  Durring the run time this variable will be
            assigned by the second objective function.
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
        self.oaf = oaf
        self.odf = odf
        self.generation = generation

    def __str__(self):
        """This method returns the formating print format, to print out
        a solution, while all the details all printed.
        """

        out = "* Genome:: {}\n* Genome size:: {}\n* OAF::{}\n* ODF:: {}\n* Generation:: {}"
        return out.format(self.genome, self.genome_size, self.oaf, self.odf, self.generation)

    def __lt__(self, other):
        """This function is for comparing an object with another.
        This comparing is for ascendant ordering.
        """
        return self.oaf < other.oaf