from Bio.Seq import Seq
from restriction_enzymes import RestrictionEnzyme
from utils import find_occurrences_between_patterns

dna_sequence = "AAGCTTATGCGAATTCCGGAAGCTTTGAATTCCGAATTCGAATTCGAATTCGAATTCGAATTC"


class Fragment:
    def __init__(self, head: int, tail: int, re1: str = None, re2: str = None):
        self.head = head  # co-ords
        self.tail = tail  # co-ords
        self.re1 = re1
        self.re2 = re2


def basic_re_cut(
        full_sequence: str,
        fragment: Fragment,
        enzyme_list: list[RestrictionEnzyme],
        which_cut: int
):
    fragments_list = []
    for re in enzyme_list:
        c_enzyme = re.complement()
        occurrences = find_occurrences_between_patterns(full_sequence, c_enzyme)
        for coord_start, coord_end in occurrences:
            fragments_list.append(
                Fragment(coord_start, coord_end)
            )


class ParentSeq(Seq):
    def __init__(self, str_seq: str, h_roi: int, t_roi: int):
        super().__init__(str_seq)
        self.full_sequence = str_seq
        self.re1_cut = False
        self.re2_cut = False

        self.re1_fragments = []
        self.re2_fragments = []

        self.h_roi = h_roi
        self.t_roi = t_roi

    def re1_cut(self, enzyme_list):
        fragment = Fragment(0, len(self.full_sequence))
        return basic_re_cut(
            full_sequence=self.full_sequence,
            fragment=fragment,
            enzyme_list=enzyme_list,
            which_cut=1
        )

    def re2_cut(self):
        fragments = basic_re_cut()
        pass
