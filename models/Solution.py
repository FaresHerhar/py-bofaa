class Solution:
    def __init__(self):
        self.genome = list()
        self.genome_size = -1
        self.oaf = -1
        self.odf = -1

    def __init__(self, genome, oaf=-1, odf=-1):
        self.genome = genome
        self.genome_size = len(self.genome)
        self.oaf = oaf
        self.odf = odf

    def __str__(self):
        out = "* Genome:: {}\n* Genome size:: {}\n* OAF::{}\n* ODF:: {}"
        return out.format(self.genome, self.genome_size, self.oaf, self.odf)