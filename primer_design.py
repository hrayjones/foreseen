from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Restriction import AllEnzymes
from Bio import SeqIO

# from Bio.Seq import Seq
# from Bio.Restriction import RestrictionBatch, Analysis
# from Bio.Restriction import AllEnzymes

def design_primer_with_restriction_present(site_name, target_sequence):
    # Create a Biopython Seq object from the target sequence
    dna_seq = Seq(target_sequence)
    
    # Define the restriction enzyme(s)
    restriction_batch = RestrictionBatch([site_name])

    # Perform restriction analysis
    analysis = Analysis(restriction_batch, dna_seq)
    
    # Get the cut sites of the enzyme on the target sequence
    cut_sites = analysis.full()

    # Check if the enzyme cuts the sequence at any site
    enzyme = restriction_batch.get(site_name)
    recognition_site = enzyme.site
    
    # If the restriction enzyme cuts the sequence
    if enzyme in cut_sites:
        cut_positions = cut_sites[enzyme]
        print(f"Restriction enzyme {site_name} site is already present at positions: {cut_positions}")

        # Forward primer: start from restriction site and extend downstream
        cut_position = cut_positions[0]  # Using the first restriction site found
        primer_fwd = target_sequence[cut_position:cut_position + 20]  # 20 bp after restriction site
        
        # Reverse primer: start from the reverse complement and extend upstream
        primer_rev = str(dna_seq.reverse_complement()[cut_position:cut_position + 20])

    else:
        print(f"Restriction enzyme site {site_name} not found. Adding the site to the primers.")

        # Design the primers by adding restriction site manually if not present
        primer_fwd = recognition_site + target_sequence[:20]  # Add the restriction site to the 5' end of the forward primer
        primer_rev = recognition_site + str(dna_seq.reverse_complement()[:20])  # Add the restriction site to the 5' end of the reverse primer
    
    return primer_fwd, primer_rev

# Example usage
target_seq = "ATGCGTACGTAGCTGAATTCGATCGTAGCTAGTACG"
restriction_site = "EcoRI"  # Restriction enzyme

fwd_primer, rev_primer = design_primer_with_restriction_present(restriction_site, target_seq)

print("Forward Primer: ", fwd_primer)
print("Reverse Primer: ", rev_primer)
