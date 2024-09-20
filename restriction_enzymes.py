from Bio.Seq import Seq


class RestrictionEnzyme(Seq):
    def __init__(self,
                 enzyme_name:str,
                 enzyme_seq:str,
                 head_after_digestion:str,
                 tail_after_digestion:str
                 ):
        super().__init__(enzyme_seq)
        self.enz_name = enzyme_name
        self.enz_seq = enzyme_seq
        self.head_add = - len(tail_after_digestion)
        self.tail_add = len(head_after_digestion)
        self.length = len(enzyme_seq)

    def complement(self, inplace=False):
        return str(Seq(self.enz_seq).complement())

