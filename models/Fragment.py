class Fragment:
    """This is a Fragment class that represent a Fragment for the
    dna sequencing problem.
    I could noave not used this class, but for the records, and keeping
    track of every data, and to make things make more sense, I defined it.

    ...

    Attributes
    ----------
    sequence: str
        This is the DNA sequence of the fragment.
    lenght: int
        The lenght of the read DNA sequence.
    index: int
        When reading the file, to extract fragments I keep the index
        of each fragments in the file. Not necessary, neither
        for the execution, nor the presentation,
        but for debugging(the solution's validity) needs. 
    orientation: bool
        Durring the run time, we need to know the orientation of the fragment
        from left to right if True, False else.

    Methods
    -------
    __init__(orientation=True): None
        The constructor.
    __str__: str
        The print formating method.
    inverse:
        For inversing the DNA sequence string.

    """

    def __init__(self, sequence, lenght, index):
        """The constructor

        ...

        Parameters
        ----------
        sequence: str
            This is the DNA sequence of the fragment.
        lenght: int
            The lenght of the read DNA sequence.
        index: int
            When reading the file, to extract fragments I keep the index
            of each fragments in the file. Not necessary, neither
            for the execution, nor the presentation,
            but for debugging(the solution's validity) needs. 


        Rturns
        ------
        None
        """

        self.sequence = sequence
        self.lenght = lenght
        self.index = index
        self.orientation = True

    def inverse(self):
        """The inverse method inverses the DNA sequence string."""
        
        self.sequence = self.sequence[::-1]
        self.orientation = not self.orientation

    def __str__(self):
        """This method returns the formating print format, to print out
        a fragment, while all the details all printed.
        """

        out = "* Sequence:: {}\n* The Lenght:: {}\n* The Index:: {}\n"
        if self.orientation:
            out += "* From Left to Right."
        else:
            out += "* From Right to Left."

        return out.format(self.sequence, self.lenght, self.index)
