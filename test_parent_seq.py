import unittest
from parent_seq import ParentSeq
from restriction_enzymes import RestrictionEnzyme

# tests (1, 2, 4, 6) should be returned, (3, 5, 7) shouldn't be.
ex_dna1 = "GCGA_test1_TTCG_test2_GCGA_test3_GCGA"
"""     ROI=   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               4                           32
"""
ex_dna2 = "GCGA_test4_TTCG_test5_TTCG_test6_GCGA_test7_GCGA"
"""     ROI=   ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
               4                               38
"""

ex_re1 = RestrictionEnzyme(
    enzyme_name="example_re1",
    enzyme_seq="CGCT",
    head_after_digestion="C",
    tail_after_digestion="CT"
)
ex_re2 = RestrictionEnzyme(
    enzyme_name="example_re1",
    enzyme_seq="AAGC",
    head_after_digestion="AAG",
    tail_after_digestion=""
)


class TestParentSeq(unittest.TestCase):
    def test_valid_input(self):
        # Test with valid input
        ps1 = ParentSeq(
            str_seq=ex_dna1,
            h_roi=4,
            t_roi=32
        )
        ps2 = ParentSeq(
            str_seq=ex_dna1,
            h_roi=4,
            t_roi=38
        )
        enzymes_list1 = [ex_re1]
        enzymes_list2 = [ex_re2]
        ps1.re1_cut(enzymes_list1)
        ps1.re2_cut(enzymes_list2)
        ps2.re1_cut(enzymes_list1)
        ps2.re2_cut(enzymes_list2)

        print("###########################ex_dna1###########################")
        for fragment in ps1.re2_fragments:
            print(ex_dna1[fragment.head:fragment.tail +1])
        print("###########################ex_dna2###########################")
        for fragment in ps2.re2_fragments:
            print(ex_dna2[fragment.head:fragment.tail +1])


# This runs the test when you run the file
if __name__ == '__main__':
    unittest.main()
