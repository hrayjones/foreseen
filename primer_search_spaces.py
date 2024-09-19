def define_search_space(sequence, read_length, starts_with_re1=True):
    """
    Defines the Reverse Primer and Forward Primer Search Space coordinates and assigns labels 
    based on the provided read length within the sequence. The start position is always set to 0.
    
    Parameters:
    sequence (str): The genetic sequence (e.g., DNA sequence) as a string.
    read_length (int): The read length, which must be between 75 and 150 bp.
    starts_with_re1 (bool): A boolean indicating if the sequence starts with RE1. 
                            If False, the sequence starts with RE2.
    
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
    
    # Define Reverse Primer Search Space End and Forward Primer Search Space Start based on the start condition
    if starts_with_re1:
        Reverse_primer_search_space_end = start + (read_length - 40)
        forward_primer_search_space_start = end - 150
        reverse_primer_label = "Reading Primer"
        forward_primer_label = "Non-reading Primer"
    else:  # If the sequence starts with RE2
        Reverse_primer_search_space_end = start + 150
        forward_primer_search_space_start = end - (read_length - 40)
        reverse_primer_label = "Non-reading Primer"
        forward_primer_label = "Reading Primer"
    
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
genetic_sequence = "ATGCGTACCGGTTAGCTAGGCTAGCTAGCTAGGCTA"  # Example sequence
read_length = 100    # Input read length (between 75 and 150)
starts_with_re1 = True  # Set to True if sequence starts with RE1, False if it starts with RE2

# Querying the Reverse Primer and Forward Primer Search Space Coordinates and Labels
search_results = define_search_space(genetic_sequence, read_length, starts_with_re1)
print(search_results)
