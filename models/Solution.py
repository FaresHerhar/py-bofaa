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
    rank:
        The rank if the solution, i.e the fonts that the solution belongs to,
        and it changes from a generation to another.
    crowding_distance: float, optional
        The crowding distance of the solution that changes, from a generation
        to a generation, based on the font that it belongs to.

    Methods
    -------
    __init__(oaf=-1, odf=-1, generation=-1): None
        The constructor.
    __str__: str
        The print formating method.
    """

    def __init__(self, genome, oaf=-1, odf=-1, generation=-1, rank=-1, crowding_distance=-1):
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
        rank: int, optional
            The rank if the solution, i.e the fonts that the solution belongs to.
        crowding_distance: float, optional
            The crowding distance of the solution that changes, from a generation
            to a generation, based on the font that it belongs to.

        Rturns
        ------
        None
        """

        self.genome = genome
        self.genome_size = len(self.genome)
        self.oaf = oaf
        self.odf = odf
        self.generation = generation
        self.rank = rank
        self.crowding_distance = crowding_distance

    def __str__(self):
        """This method returns the formating print format, to print out
        a solution, while all the details all printed.
        """

        out = "* Genome:: {}\n* Genome size:: {}\n* OAF::{}\n* ODF:: {}\n* Rank:: {}\n* Crowding distance:: {}\n* Generation:: {}"
        return out.format(self.genome, self.genome_size, self.oaf, self.odf, self.rank, self.crowding_distance, self.generation)
