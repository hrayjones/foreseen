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


def defineSearchSpace(sequence, read_length, start_state="reading"):
    """
    Defines the Reverse Primer and Forward Primer Search Space coordinates and assigns labels 
    based on the provided read length within the sequence. The start position is always set to 0.
    
    Parameters:
    sequence (str): The genetic sequence (e.g., DNA sequence) as a string.
    read_length (int): The read length, which must be between 75 and 150 bp.
    start_state (str): Indicates whether the sequence starts with RE1 or RE2. 
                       Accepts "RE1" or "RE2".
    
    Returns:
    dict: A dictionary containing the coordinates and sequence slices for Reverse Primer and Forward Primer Search Space,
          and primer labels.

    ### TODO
    ### Add the case where we want to have our genetic variant between the primer and the cut site. This can be either for RE1 or RE2
    
    """
    
    # Validate read length to ensure it's between 75 and 150 bp (recommended sequencing lengths for 4C-seq)
    if not (75 <= read_length <= 150):
        raise ValueError("Read length must be between 75 and 150 bp.")
    
    # Set the start position to 0
    start = 0
    end = len(sequence)
    
    # Define Reverse Primer Search Space End and Forward Primer Search Space Start based on the start_state
    if start_state == "reading":
        Reverse_primer_search_space_end = start + (read_length - 40)
        forward_primer_search_space_start = end - 150
        reverse_primer_label = "Reading Primer"
        forward_primer_label = "Non-reading Primer"
    elif start_state == "non_reading":
        Reverse_primer_search_space_end = start + 150
        forward_primer_search_space_start = end - (read_length - 40)
        reverse_primer_label = "Non-reading Primer"
        forward_primer_label = "Reading Primer"
    else:
        raise ValueError('Invalid start_state. Must be "reading" or "non_reading".')
    
    # Ensure the Reverse Primer Search Space does not exceed the end position
    if Reverse_primer_search_space_end > end:
        Reverse_primer_search_space_end = end

   # Define coordinates for Reverse and Forward Primer Search Spaces
    reverse_primer_coordinates = (0, Reverse_primer_search_space_end)
    forward_primer_coordinates = (forward_primer_search_space_start, end)
    
    # Extract the Reverse Primer and Forward Primer Search Space
    reverse_primer_search_seq = sequence[0:Reverse_primer_search_space_end]
    forward_primer_search_seq = sequence[forward_primer_search_space_start:end]
    
    # Return details
    result = {
        "reverse_primer_coordinates": reverse_primer_coordinates,
        "forward_primer_coordinates": forward_primer_coordinates,
        "reverse_primer_label": reverse_primer_label,
        "forward_primer_label": forward_primer_label,
        "reverse_primer_search_seq": reverse_primer_search_seq,
        "forward_primer_search_seq": forward_primer_search_seq
    }
   return result

# Example usage
sequence = "ATGCGTACCGGTTAGCTAGGCTAGCTAGCTAGGCTA"  # Example sequence
read_length = 100    # Input read length (between 75 and 150)
start_state = "reading"  # Set to "reading" or "non_reading", set by convertFragmentData

def primerDesigner:

def primerMatcher:

def primerSelector(data, sequence):

    primerSet = []

    for fragment in data: 
        # Convert Data from step 3 to sequence
        fragmentData = convertFragmentData(data, sequence)

        # Define search space
        search_spaces = defineSearchSpace(
            fragmentData[frag_sequence], 
            userInput[read_length], ### this needs to come from user input
            fragmentData[start_state]
            )
        
        # Select Forward and Reverse Primer
        forwardPrimers = primerDesign(
            search_spaces["forward_primer_search_seq"], 
            forward
        )
        reversePrimers = primerDesign(
            search_spaces["reverse_primer_search_seq"], 
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
