import re

ex_string = "AGTCACGTGGTCAGAGAGAGAGGTCA"
ex_pattern = "GTC"


def find_occurrences_between_patterns(string: str, pattern: str):
    newest_endpoint = 999999999999999
    occurrences_coordinates = []
    for match in re.finditer(pattern, string):
        newest_startpoint = match.start()
        if newest_startpoint > newest_endpoint:
            occurrences_coordinates.append((newest_endpoint, newest_startpoint - 1))
        newest_endpoint = match.end()

    return occurrences_coordinates


def find_left_and_right_most_occurrences(string: str, pattern: str):
    occurrences_coordinates = []
    occurrences_coordinates.append((0, string.find(pattern) - 1))
    occurrences_coordinates.append((string.rfind(pattern) + len(pattern), len(string) - 1))

    return occurrences_coordinates
