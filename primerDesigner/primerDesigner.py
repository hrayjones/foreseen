"""primerDesigner.py
   For expanded sequence and the subsequent digested fragments will output a set of pairs of primers
   Run as "primerDesigner.py [fragment_data] [expanded_sequence]"
"""
import os
import sys
import argparse

# Functions ---------------------------------

def convertFragmentData(data, sequence):

    # NOTE PARSING WILL CHANGE DEPENDING ON OUTPUT OF PREVIOUS STEP
    # Split data into different parts
    start_RE = data.split()[0]
    if start_RE = 1:
        start_state = "reading"
        end_state = "non_reading"
    else:
        start_state = "non_reading"
        end_state = "reading"
    start_point = data.split()[1]
    end_point = data.split()[2]
    # Convert to sequence and assign which end is reading or non reading for each fragment
    frag_sequence = sequence[[start_point]:[end_point]]
    f = {
            "frag_sequence":frag_sequence,
            "start_point":start_point,
            "end_point":end_point,
            "start_state":start_state,
            "end_state":end_state
        }
    return f

def defineSearchSpace:

def primerDesigner:

def primerMatcher:

def primerSelector(data, sequence):

    primerSet = []

    for fragment in data: 
        # Convert Data from step 3 to sequence
        fragmentData = convertFragmentData(data, sequence)

        # Define search space
        forward_search_space = defineSearchSpace(
            fragmentData[frag_sequence], 
            fragmentData[end_point], 
            fragmentData[end_state]
            )
        reverse_search_space = defineSearchSpace(
            fragmentData[frag_sequence], 
            fragmentData[start_point], 
            fragmentData[start_state]
            )
        
        # Select Forward and Reverse Primer
        forwardPrimers = primerDesign(
            fragment[frag_sequence], 
            forward_search_space, 
            forward
        )
        reversePrimers = primerDesign(
            fragment[frag_sequence], 
            reverse_search_space, 
            reverse
        )

        # Pair up forward and reverse primers
        matchedPrimers = primerMatcher(forwardPrimers, reversePrimers)

        primerSet = # Append matched primers for set ... unsure how best to setup data structure for this

    return primerSet
    
# Code Body ---------------------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('data')
    parser.add_argument('sequence')
    args = parser.parse_args()
    primers = primerSelector(args.data, args.sequence)