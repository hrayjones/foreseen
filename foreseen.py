# Example DNA sequence (for testing purposes, should be replaced with real sequences)
dna_sequence = "AAGCTTATGCGAATTCCGGAAGCTTTGAATTCCGAATTCGAATTCGAATTCGAATTCGAATTC"

# Define RE1 and RE2 panels (four-cutters and suitable restriction enzymes)
re1_panel = {
    "HindIII": "AAGCTT",
    "EcoRI": "GAATTC",
    "MspI": "CCGG",
    "HaeIII": "GGCC"
}

re2_panel = {
    "BamHI": "GGATCC",
    "XhoI": "CTCGAG",
    "NcoI": "CCATGG"
}

# Define the target region (example)
roi_start = 15
roi_end = 20

# Expansion size
expansion_size = 2000

# 1. Expand the target region +/- 2kb
expanded_start = max(0, roi_start - expansion_size)
expanded_end = roi_end + expansion_size

# Function to find cut positions of a restriction enzyme
def find_cut_positions(sequence, enzyme_site, start, end):
    cut_positions = []
    position = start
    region_sequence = sequence[start:end]
    while position < end:
        pos = region_sequence.find(enzyme_site, position - start)
        if pos == -1:
            break
        cut_positions.append(start + pos)
        position = start + pos + 1
    return cut_positions

# 2. Digest with all RE1 options
re1_options = {}
for re1_name, re1_site in re1_panel.items():
    cut_positions = find_cut_positions(dna_sequence, re1_site, expanded_start, expanded_end)
    if cut_positions:
        re1_options[re1_name] = cut_positions

if not re1_options:
    print("No RE1 enzymes cut in the expanded region. Stopping.")
    exit()

# 3. Digest remaining RE1 options with RE2 panel
valid_re1_re2_combinations = {}

for re1_name, re1_cuts in re1_options.items():
    for re2_name, re2_site in re2_panel.items():
        re2_cuts = find_cut_positions(dna_sequence, re2_site, expanded_start, expanded_end)
        if not re2_cuts:
            continue
        
        # 4. Check if RE2 cuts result in fragment sizes between 200 bp and 1500 bp
        for re1_cut in re1_cuts:
            for re2_cut in re2_cuts:
                fragment_size = abs(re2_cut - re1_cut)
                if 200 <= fragment_size <= 1500:
                    if re1_name not in valid_re1_re2_combinations:
                        valid_re1_re2_combinations[re1_name] = []
                    valid_re1_re2_combinations[re1_name].append((re2_name, re1_cut, re2_cut, fragment_size))

if not valid_re1_re2_combinations:
    print("No valid RE1-RE2 combinations found. Stopping.")
else:
    print("Valid RE1-RE2 combinations:")
    for re1_name, re2_data in valid_re1_re2_combinations.items():
        for re2_name, re1_cut, re2_cut, fragment_size in re2_data:
            print(f"RE1: {re1_name}, RE1 cut: {re1_cut}, RE2: {re2_name}, RE2 cut: {re2_cut}, Fragment size: {fragment_size} bp")
