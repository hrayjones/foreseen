from Bio.Seq import Seq
from Bio.Restriction import RestrictionBatch, Analysis
from Bio.Restriction import AllEnzymes
from Bio import SeqIO

def design_primer_with_restriction(site_name, target_sequence):
    # Create a Biopython Seq object from the target sequence
    dna_seq = Seq(target_sequence)
    
    # Define the restriction enzyme(s)
    restriction_batch = RestrictionBatch([site_name])

    # Perform restriction analysis
    analysis = Analysis(restriction_batch, dna_seq)
    
    # Check the presence of restriction enzyme sites in the sequence
    cut_sites = analysis.full()
    
    if cut_sites:
        print(f"Restriction enzyme site {site_name} is already present in the sequence. Use a different site.")
        return None
    
    # Get the restriction enzyme and its recognition site
    enzyme = restriction_batch.get(site_name)
    recognition_site = enzyme.site
    
    # Design the primers
    primer_fwd = recognition_site + target_sequence[:20]  # Add the restriction site to the 5' end of the forward primer
    primer_rev = recognition_site + str(dna_seq.reverse_complement()[:20])  # Add the restriction site to the 5' end of the reverse primer
    
    return primer_fwd, primer_rev

# Example usage
target_seq = "ATGCGTACGTAGCTAGCGTACGATCGTAGCTAGTACG"
restriction_site = "EcoRI"  # Restriction enzyme you want to use

fwd_primer, rev_primer = design_primer_with_restriction(restriction_site, target_seq)

if fwd_primer and rev_primer:
    print("Forward Primer: ", fwd_primer)
    print("Reverse Primer: ", rev_primer)
