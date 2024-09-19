import re


ex_string = "AGTCACGTGGTCAGAGAGAGAGGTCA"
ex_pattern = "GTC"


def find_occurrences_between_patterns(string:str, pattern:str):
    newest_endpoint = 999999999999999
    occurences_coordinates = []
    for match in re.finditer(pattern, string):
        newest_startpoint = match.start()
        if newest_startpoint > newest_endpoint:
            occurences_coordinates.append((newest_endpoint, newest_startpoint-1))
        newest_endpoint = match.end()

    return occurences_coordinates


