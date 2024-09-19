"""expandSelection.py
   Searches the human_genome for a matching sequence of nucleotides and outputs an exteneded sequence
   Run as "expandSelection.py [inputSequence]"
"""
import os
import subprocess
import sys
import argparse
from Bio.Blast import NCBIXML
from Bio import SeqIO

#-----------------------------------

def run_blast(sequence, db, blast_type="blastn", output="blast_output.xml"):

    # Run Blast command
    blast = f"{blast_type} -query {sequence} -db {db} -out {output} -outfmt 5"
    subprocess.run(blast, shell=True, check=True)

    return output

def blast_parser(blast_xml):

    # Parse the outputed blast results
    with open(blast_xml) as result_handle:
        blast_records = NCBIXML.parse(result_handle)
        for blast_record in blast_records:
            if blast_record.alignments:
                top_alignment = blast_record.alignments[0]
                top_hsp = top_alignment.hsps[0]
                hit = {
                    "chromosome":top_alignment.hit_def.split()[0],
                    "start": top_hsp.sbjct_start,
                    "end": top_hsp.sbjct_end
                }
                return hit
    return None

def extend_sequence(hit, genome_fna):

    # Load genome sequence
    genome_seq = SeqIO.to_dict(SeqIO.parse(genome_fna, 'fasta'))

    # Get chromosome sequence
    chromosome = genome_seq[hit["chromosome"]]

    # Count sequence size
    size = hit["end"] - hit["start"]

    # Calculate extended coordinates
    if size > 1000:
        extended_start = max(0, hit["start"] - 1000)
        extended_end = min(len(chromosome), hit["end"] + 1000)
    else:
        extended_start = max(0, hit["start"] - 2000)
        extended_end = min(len(chromosome), hit["end"] + 2000)

    # Extract extended sequence
    extended_sequence = chromosome.seq[extended_start:extended_end]

    return extended_sequence

def search_and_expand_selection(input_sequence):

    # Define databases
    db="human_genome_db"
    genome_fna="human_genome.fna"
    output_file="results.txt"

    # Run BLAST search
    blast_search = run_blast(input_sequence, db)

    # Extract position details
    hit = blast_parser(blast_search)

    # Expand selection
    if hit:
        expanded_sequence = extend_sequence(hit, genome_fna)

        # Save the extended sequence to file
        with open(output_file, "w") as file:
            file.write(f">Expanded_sequence\n{expanded_sequence}\n")
        return expanded_sequence
    else:
        return "No hit found"

# Code Body -----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file')
    args = parser.parse_args()
    expanded_selection = search_and_expand_selection(args.file)
