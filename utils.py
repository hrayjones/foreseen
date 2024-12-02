import re
import pytest
from fixtures import short_dna_example, short_pattern_example


def find_strings_between_patterns(string: str, pattern: str):
    if len(re.findall(pattern, string)) <= 1:
        return []  # only strings flanked by the pattern should be returned
    else:
        pattern_coordinates = []
        for match in re.finditer(pattern, string):
            pattern_coordinates.append(
                (match.start(), match.end())
            )
        
        occurrences_coordinates = []
        for idx in range(len(pattern_coordinates)):
            if idx < len(pattern_coordinates) - 1:
                occurrences_coordinates.append(
                    (pattern_coordinates[idx][1], pattern_coordinates[idx+1][0]-1)
                )

        return occurrences_coordinates

def find_left_and_right_most_strings(string: str, pattern: str):  # strings before first pattern occurence and strings after last pattern occurence
    occurrences_coordinates = []  # test for number of occurences?
    if pattern in string:
        occurrences_coordinates.append((0, string.find(pattern)-1))
        occurrences_coordinates.append((string.rfind(pattern) + len(pattern), len(string)-1))

    return occurrences_coordinates


# Tests
def test_find_strings_between_patterns(short_dna_example, short_pattern_example):
    occ1, occ2 = find_strings_between_patterns(short_dna_example, short_pattern_example)
    first_occurence = short_dna_example[occ1[0]:occ1[1]+1]
    second_occurence = short_dna_example[occ2[0]:occ2[1]+1]
    print(f"First Occurence = {first_occurence}")
    print(f"Second Occurence = {second_occurence}")
    assert (first_occurence == "ACGTG" and \
        second_occurence == "AGAGAGAGAG")
    
def test_find_strings_between_patterns_one_match_only(short_dna_example, short_pattern_example):
    modified_dna = short_pattern_example.join(short_dna_example.split(short_pattern_example)[0:2])
    print(f"Modified DNA, only one match = {modified_dna}")
    assert len(find_strings_between_patterns(modified_dna, short_pattern_example)) == 0

def test_find_strings_between_patterns_no_matches(short_dna_example, short_pattern_example):
    modified_dna = short_pattern_example.join(short_dna_example.split(short_pattern_example)[0:1]) 
    print(f"Modified DNA, no matches = {modified_dna}")
    assert len(find_strings_between_patterns(modified_dna, short_pattern_example)) == 0

def test_find_left_and_right_most_strings(short_dna_example, short_pattern_example):
    left_coord, right_coord = find_left_and_right_most_strings(short_dna_example, short_pattern_example)
    assert short_dna_example[left_coord[0]:left_coord[1]+ 1] == "A" and \
        short_dna_example[right_coord[0]:right_coord[1]+ 1] == "AC"