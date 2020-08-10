class Fragment:
    def __init__(self):
        self.sequence = ""
        self.lenght = -1
        self.index = -1
        self.orientation = True

    def __init__(self, sequence, lenght, index, orientation=True):
        self.sequence = sequence
        self.lenght = lenght
        self.index = index
        self.orientation = True

    def inverse(self):
        self.sequence = self.sequence[::-1]
        self.orientation = not self.orientation

    def __str__(self):
        out = "* Sequence:: {}\n* The Lenght:: {}\n* The Index:: {}\n"
        if self.orientation:
            out += "* From Left to Right."
        else:
            out += "* From Right to Left."

        return out.format(self.sequence, self.lenght, self.index)
