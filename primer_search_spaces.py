def define_search_space(sequence, start, end, read_length, starts_with_re1=True):
    """
    Defines the Reverse Primer and Forward Primer Search Space coordinates and assigns labels 
    based on the provided start and end positions within the sequence.
    
    Parameters:
    sequence (str): The genetic sequence (e.g., DNA sequence) as a string.
    start (int): The starting position of the search space defined by RE1 or RE2 (0-indexed).
    end (int): The ending position of the search space defined by RE1 or RE2 (0-indexed).
    read_length (int): The read length, which must be between 75 and 150 bp.
    starts_with_re1 (bool): A boolean indicating if the sequence starts with RE1. 
                            If False, the sequence starts with RE2.
    
    Returns:
    dict: A dictionary containing the coordinates for Reverse Primer and Forward Primer Search Space,
          and primer labels.
    """
    # Validate read length to ensure it's between 75 and 150 bp
    if not (75 <= read_length <= 150):
        raise ValueError("Read length must be between 75 and 150 bp.")
    
    # Validate the start and end positions
    if start < 0 or end > len(sequence) or start >= end:
        raise ValueError("Invalid start or end positions in the sequence.")
    
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
    
    # Ensure the Reverse Primer Search Space does not exceed the specified end position
    if Reverse_primer_search_space_end > end:
        Reverse_primer_search_space_end = end
    
    # Extract the Reverse Primer Search Space
    reverse_primer_space = sequence[start:Reverse_primer_search_space_end]
    
    # Define coordinates for Reverse and Forward Primer Search Spaces
    reverse_primer_coordinates = (0, Reverse_primer_search_space_end)
    forward_primer_coordinates = (forward_primer_search_space_start, len(sequence))
    
    # Return details
    result = {
        "reverse_primer_coordinates": reverse_primer_coordinates,
        "forward_primer_coordinates": forward_primer_coordinates,
        "reverse_primer_label": reverse_primer_label,
        "forward_primer_label": forward_primer_label,
        "reverse_primer_space": reverse_primer_space
    }
    
    return result

# Example usage
genetic_sequence = "ATGCGTACCGGTTAGCTAGGCTAGCTAGCTAGGCTA"  # Example sequence
start_position = 5   # Known start position of RE1 or RE2
end_position = 35    # Known end position of RE1 or RE2
read_length = 100    # Input read length (between 75 and 150)
starts_with_re1 = True  # Set to True if sequence starts with RE1, False if it starts with RE2

# Querying the Reverse Primer and Forward Primer Search Space Coordinates and Labels
search_results = define_search_space(genetic_sequence, start_position, end_position, read_length, starts_with_re1)
print(search_results)
