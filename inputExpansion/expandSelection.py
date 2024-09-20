"""expandSelection.py
   Searches the human_genome for a matching sequence of nucleotides and outputs an exteneded sequence
   Run as "expandSelection.py [inputSequence] or [input_position] {optional arguments --snp_pos}"
   Input sequence can either be in a .fasta or .fna format
   Input position must be in format of "ChrX [start_pos] [end_pos]"

   Example standard input using a fasta or fna file
   python3 expandSelection.py file.fna

   Example of snp_pos aware using positional input argument
   python3 expandSelection.py "Chr7 50000 51000" --snp_pos 50005

   To run local the following must be in the running directory:
   human_genome_db - see GITHUB
   human_genome.fna
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
        extension_size =  hit["start"] - extended_start
        extended_end = min(len(chromosome), hit["end"] + 1000)
    else:
        extended_start = max(0, hit["start"] - 2000)
        extension_size =  hit["start"] - extended_start
        extended_end = min(len(chromosome), hit["end"] + 2000)

    # Extract extended sequence
    extended_sequence = chromosome.seq[extended_start:extended_end]

    return extended_sequence, extension_size

def search_and_expand_selection(input,snp_pos):

    # Define databases
    db="human_genome_db"
    genome_fna="human_genome.fna"
    output_file="results.txt"

    if input.endswith(".fna") or input.endswith(".fasta"):
        input_sequence = input
        # Run BLAST search
        blast_search = run_blast(input_sequence, db)
        # Extract position details
        hit = blast_parser(blast_search)
    else:
        input_string = input[3:]
        chromosome_no = input_string.split()[0]
        with open(genome_fna, 'r') as f:
            for line in f.readlines():
                if f"chromosome {chromosome_no}," in line:
                    chromosome_info = line
        hit = {
            "chromosome": chromosome_info.split()[0][1:],
            "start": int(input_string.split()[1]),
            "end": int(input_string.split()[2])
        }

    # Expand selection
    if hit:
        expanded_sequence, extension_size = extend_sequence(hit, genome_fna)

        # If snp_pos included
        if snp_pos != None:
            # Calc expanded position
            snp_pos_exp = snp_pos + extension_size
            # Save extended sequence and snp_pos
            with open(output_file, "w") as file:
                file.write(f">Expanded_sequence\n{expanded_sequence}\n>snp_pos_original\n{snp_pos}\n>snp_pos_expanded\n{snp_pos_exp}")
        else: 
            # Save the extended sequence to file
            with open(output_file, "w") as file:
                file.write(f">Expanded_sequence\n{expanded_sequence}\n")
        
        return expanded_sequence
    else:
        return "No hit found"

# Code Body -----------------------
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('input')
    parser.add_argument('--snp_pos', type=int,
                    help='An optional integer argument')
    args = parser.parse_args()
    expanded_selection = search_and_expand_selection(args.input, args.snp_pos)
