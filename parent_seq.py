from Bio.Seq import Seq
from restriction_enzymes import RestrictionEnzyme
from utils import find_occurrences_between_patterns, find_left_and_right_most_occurrences

dna_sequence = "AAGCTTATGCGAATTCCGGAAGCTTTGAATTCCGAATTCGAATTCGAATTCGAATTCGAATTC"


class Fragment:
    def __init__(self,
                 head: int,
                 tail: int,
                 re1: RestrictionEnzyme = None,
                 re2: RestrictionEnzyme = None,
                 re2_head: bool = None
                 ):
        self.head = head  # co-ords
        self.tail = tail  # co-ords
        self.re1 = re1
        self.re2 = re2

        # Needed in primers, to determine if the head was cut by RE1 or by RE2
        self.re2_head = re2_head

    def size(self):
        return abs(self.tail - self.head)

def basic_re_cut(
        full_sequence: str,
        fragment: Fragment,
        enzyme_list: list[RestrictionEnzyme],
        which_cut: int
):
    partial_sequence = full_sequence[fragment.head:fragment.tail + 1]
    fragments_list = []
    for rest_enz in enzyme_list:
        complement_enzyme = rest_enz.complement()
        if which_cut == 1:
            occurrences = find_occurrences_between_patterns(partial_sequence, complement_enzyme)
            for coord_start, coord_end in occurrences:
                fragments_list.append(
                    Fragment(coord_start+rest_enz.head_add, coord_end+rest_enz.tail_add, re1=rest_enz, re2=None)
                )
        elif which_cut == 2:
            occurrences = find_left_and_right_most_occurrences(partial_sequence, complement_enzyme)
            if occurrences:
                coord_start, coord_end = occurrences[0]
                new_head = fragment.head + coord_start
                new_tail = fragment.head + coord_end + rest_enz.tail_add
                fragments_list.append(
                    Fragment(new_head, new_tail, re1=fragment.re1, re2=rest_enz, re2_head=False)
                )

                coord_start, coord_end = occurrences[1]
                new_head = fragment.head + coord_start + rest_enz.head_add
                new_tail = fragment.head + coord_end
                fragments_list.append(
                    Fragment(new_head, new_tail, re1=fragment.re1, re2=rest_enz, re2_head=True)
                )

    return fragments_list


class ParentSeq(Seq):
    def __init__(self, str_seq: str, h_roi: int, t_roi: int):
        super().__init__(str_seq)
        self.full_sequence = str_seq

        self.re1_fragments = []
        self.re2_fragments = []

        self.re1_options = []
        self.re2_options = []

        self.h_roi = h_roi
        self.t_roi = t_roi

    def re1_cut(self, enzyme_list):
        fragment = Fragment(0, len(self.full_sequence) - 1)
        self.re1_fragments_ = basic_re_cut(
            full_sequence=self.full_sequence,
            fragment=fragment,
            enzyme_list=enzyme_list,
            which_cut=1
        )
        for frag in self.re1_fragments_:
            if self.h_roi <= frag.head <self.t_roi or self.h_roi < frag.tail <= self.t_roi:
                self.re1_fragments.append(frag)
        self.re1_options = [frag.re1 for frag in self.re1_fragments]

    def re2_cut(self, enzyme_list):
        for fragment in self.re1_fragments:
            self.re2_fragments_ += basic_re_cut(
                full_sequence=self.full_sequence,
                fragment=fragment,
                enzyme_list=enzyme_list,
                which_cut=2
            )
        for frag in self.re2_fragments_:
                if self.h_roi <= frag.head <self.t_roi or self.h_roi < frag.tail <= self.t_roi:
                    if 200 <= fragment.size() <= 1500:
                        self.re2_fragments.append(frag)
        self.re2_options = [frag.re2 for frag in self.re2_fragments]
