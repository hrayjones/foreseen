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


def defineSearchSpace(sequence, read_length, start_state="reading", mode="basic", SNP_location=None):

   # NOTE 1: ONLY BASIC MODE WILL WORK UNTIL SNP LOCATION WRT FRAGMENT ENDS IS BROUGHT FROM UPSTREAM PIPELINE 
   # NOTE 2: WE NEED TO BRING IN USER INPUT OF SEQUENCING READ LENGTH HERE ("read_length" argument)
   
    """
    Defines the Reverse Primer and Forward Primer Search Space coordinates and assigns labels 
    based on the provided read length within the sequence. The start position is always set to 0.
    
    Parameters:
    sequence (str): The genetic sequence (e.g., DNA sequence) as a string.
    read_length (int): The read length, which must be between 75 and 150 bp.
    start_state (str): Indicates whether the start of the sequence is "reading" or "non_reading".
                       Accepts "reading" or "non_reading".
    mode (str): Mode in which the function runs. Accepts "basic", "allele-specific", and "allele-aware".
                Default is "basic".
    SNP_location (int): The position of the SNP within the sequence. Required if mode is "allele-aware".
    
    Returns:
    dict: A dictionary containing the coordinates and sequence slices for Reverse Primer and Forward Primer Search Space,
          and primer labels.
    """
    # Validate read length to ensure it's between 75 and 150 bp
    if not (75 <= read_length <= 150):
        raise ValueError("Read length must be between 75 and 150 bp.")
    
    # Set the start position to 0
    start = 0
    end = len(sequence)
    
    # Validate start state
    if start_state not in ["reading", "non_reading"]:
        raise ValueError('Invalid start_state. Must be "reading" or "non_reading".')

    # Validate mode
    if mode not in ["basic", "allele-specific", "allele-aware"]:
        raise ValueError('Invalid mode. Must be "basic", "allele-specific", or "allele-aware".')

    # Validate SNP_location if mode is "allele-aware"
    if mode == "allele-aware":
        if SNP_location is None:
            raise ValueError('SNP_location is required when mode is "allele-aware".')
        # Check if SNP is within 100 bp of the start or the end of the sequence
        if SNP_location > 100 and SNP_location < (end - 100):
            raise ValueError("SNP must be within 100 bp of the start or end of the sequence.")

    # Define default primer search space end positions based on the start_state
    if start_state == "reading":
        reverse_primer_search_space_end = start + (read_length - 40)
        forward_primer_search_space_start = end - 150
        reverse_primer_label = "Reading Primer"
        forward_primer_label = "Non-reading Primer"
    else:  # If the start_state is "non_reading"
        reverse_primer_search_space_end = start + 150
        forward_primer_search_space_start = end - (read_length - 40)
        reverse_primer_label = "Non-reading Primer"
        forward_primer_label = "Reading Primer"
    
    # Adjust primer search spaces if mode is "allele-aware"
    if mode == "allele-aware":
        if SNP_location <= 100:
            # SNP is within 100bp of the start; adjust the reverse primer search space to start just after SNP
            reverse_primer_search_space_start = SNP_location + 1
            reverse_primer_search_space_end = min(reverse_primer_search_space_end, end)
        elif SNP_location >= (end - 100):
            # SNP is within 100bp of the end; adjust the forward primer search space to end just before SNP
            forward_primer_search_space_start = max(forward_primer_search_space_start, 0)
            forward_primer_search_space_end = SNP_location
            forward_primer_search_coordinates = (forward_primer_search_space_start, forward_primer_search_space_end)
        else:
            # SNP not near boundaries; keep default settings
            reverse_primer_search_space_start = 0
            forward_primer_search_coordinates = (forward_primer_search_space_start, end)

    # Define coordinates for Reverse and Forward Primer Search Spaces
    reverse_primer_search_coordinates = (reverse_primer_search_space_start, reverse_primer_search_space_end)
    forward_primer_search_coordinates = (forward_primer_search_space_start, end)
    
    # Extract the Reverse Primer and Forward Primer Search Space
    reverse_primer_search_seq = sequence[reverse_primer_search_space_start:reverse_primer_search_space_end]
    forward_primer_search_seq = sequence[forward_primer_search_space_start:end]
    
    # Return details
    result = {
        "reverse_primer_search_coordinates": reverse_primer_search_coordinates,
        "forward_primer_search_coordinates": forward_primer_search_coordinates,
        "reverse_primer_label": reverse_primer_label,
        "forward_primer_label": forward_primer_label,
        "reverse_primer_search_seq": reverse_primer_search_seq,
        "forward_primer_search_seq": forward_primer_search_seq
    }
    
    return result

# Example usage
genetic_sequence = "ATGCGTACCGGTTAGCTAGGCTAGCTAGCTAGGCTA"  # Example sequence
read_length = 100    # Input read length (between 75 and 150)
start_state = "reading"  # Set to "reading" or "non_reading" based on the sequence start
mode = "allele-aware"  # Choose from "basic", "allele-specific", or "allele-aware"
SNP_location = 80    # Position of the SNP within the sequence (required for "allele-aware" mode)



def primerDesigner:

def primerMatcher(forward_primers, reverse_primers, forward_template_length):
    """
    Matches forward and reverse primers from two different runs of Primer3, scoring them based on product size,
    melting temperature similarity, and GC content similarity.

    Args:
        forward_primers (dict): A dictionary with forward primer IDs as keys and a dictionary of properties 
                                (sequence, penalty, melting temperature, GC content, start, end) as values.
        reverse_primers (dict): A dictionary with reverse primer IDs as keys and a dictionary of properties 
                                (sequence, penalty, melting temperature, GC content, start, end) as values.
        forward_template_length (int): The length of the template sequence for the forward primer.

    Returns:
        list: A ranked list of paired primers, scored and sorted.
    """
    def calculate_product_size(forward_end, reverse_start, forward_template_length):
        # Calculate the product size based on the positions of the primers on the template
        return (forward_template_length - forward_end) + reverse_start
    
    def score_melting_temp(Tm_f, Tm_r):
        # Score the similarity of melting temperatures between forward and reverse primers (closer is better)
        return abs(Tm_f - Tm_r)
    
    def score_gc_content(gc_f, gc_r):
        # Score the similarity of GC content between forward and reverse primers (closer is better)
        return abs(gc_f - gc_r)
    
    def score_product_size(product_size):
        # Score product size, smaller product sizes score better (penalty increases with size)
        return product_size
    
    paired_primers = []

    # Iterate through all possible pairs of forward and reverse primers
    for f_primer, f_props in forward_primers.items():
        for r_primer, r_props in reverse_primers.items():
            # Calculate the product size
            product_size = calculate_product_size(
                f_props['end'], 
                r_props['start'], 
                forward_template_length
            )
            
            # Calculate scores based on melting temperature, GC content, and product size
            tm_score = score_melting_temp(f_props['melting_temp'], r_props['melting_temp'])
            gc_score = score_gc_content(f_props['gc_content'], r_props['gc_content'])
            product_size_score = score_product_size(product_size)
            
            # Combine scores (the lower the score, the better)
            total_score = tm_score + gc_score + product_size_score
            
            # Append the primer pair and its scores
            paired_primers.append({
                'forward_primer_id': f_primer,
                'forward_primer_sequence': f_props['sequence'],
                'reverse_primer_id': r_primer,
                'reverse_primer_sequence': r_props['sequence'],
                'product_size': product_size,
                'tm_score': tm_score,
                'gc_score': gc_score,
                'product_size_score': product_size_score,
                'total_score': total_score
            })

    # Sort the primer pairs by total score (ascending)
    ranked_primers = sorted(paired_primers, key=lambda x: x['total_score'])
    
    return ranked_primers

# Example input
forward_primers = {
    'F1': {
        'sequence': 'AGCTAGCTAGCTAGCTAGCT', 
        'penalty': 1.0, 
        'melting_temp': 60.5, 
        'gc_content': 50.0, 
        'start': 10, 
        'end': 30
    },
    'F2': {
        'sequence': 'CGTACGTACGTACGTACGTA', 
        'penalty': 0.8, 
        'melting_temp': 62.0, 
        'gc_content': 52.0, 
        'start': 20, 
        'end': 40
    },
}

reverse_primers = {
    'R1': {
        'sequence': 'TGCATGCATGCATGCATGCA', 
        'penalty': 1.2, 
        'melting_temp': 61.0, 
        'gc_content': 51.0, 
        'start': 100, 
        'end': 120
    },
    'R2': {
        'sequence': 'ATCGATCGATCGATCGATCG', 
        'penalty': 0.9, 
        'melting_temp': 63.5, 
        'gc_content': 53.0, 
        'start': 110, 
        'end': 130
    },
}

forward_template_length = 150

# Running the function
ranked_pairs = match_primers(forward_primers, reverse_primers, forward_template_length)

# Output the ranked primer pairs
for pair in ranked_pairs:
    print(pair)

   

def primerSelector(data, sequence):

    primerSet = []

    for fragment in data: 
        # Convert Data from step 3 to sequence
        fragmentData = convertFragmentData(data, sequence)

        # Define search space
        search_spaces = defineSearchSpace(
            fragmentData[frag_sequence], 
            fragmentData[start_state], 
            fragmentData[SNP_location], # this needs to come from earlier in pipeline: location of SNP in relation to fragment ends
            userInput[read_length], # needs to come from user
            userInput[mode], # needs to come from user
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
