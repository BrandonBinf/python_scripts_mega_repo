#!/usr/bin/env python3
import os, argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

# Define the standard genetic code start codon and stop codons
START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}


''''
This function reads in a fasta file and returns a dictionary with chromosome IDs as keys and RNA sequences as values.
It converts the DNA sequences to RNA sequences by transcribing them.
It also extracts the chromosome ID from the description of each record.
dict looks like this:
{ "chr1": rna_seq1, "chr2": rna_seq2, ... }
'''
def read_sequence_from_file(filename):
	fasta_data = {}
	with open(filename, "r") as file: 
		for record in SeqIO.parse(file, "fasta"):
			annotation = record.description
			parts = annotation.split(",")
			chromosome = parts[1].strip()
			#print(chromosome)
			dna_seq = record.seq.upper()
			rna_seq = str(dna_seq.transcribe())
			fasta_data[chromosome] = rna_seq
		return fasta_data
      

'''
dict looks like this:
{    "chr1": {
        "description1": [m6a_site1, m6a_site2, ...],
        "description2": [m6a_site1, m6a_site2, ...], ...},
    "chr2": {
        "description1": [m6a_site1, m6a_site2, ...],
        "description2": [m6a_site1, m6a_site2, ...], ... },
    ...
'''
def read_bed(bed_file):
    CHROM_COL_INDEX = 0
    M6A_SITE_COL_INDEX = 1
    DESCRIPTION_COL_INDEX = 3
    m6a_dict = {} 
    bed_df = pd.read_csv(
        bed_file,
        sep='\t',
        skiprows=1,
        header=None )
    for row_tuple in bed_df.itertuples(index=False):
        chrom = row_tuple[CHROM_COL_INDEX]
        description = row_tuple[DESCRIPTION_COL_INDEX]
        m6a_site = int(row_tuple[M6A_SITE_COL_INDEX])
        if chrom not in m6a_dict:
            m6a_dict[chrom] = {}
        if description not in m6a_dict[chrom]:
            m6a_dict[chrom][description] = []
        m6a_dict[chrom][description].append(m6a_site)
    return m6a_dict
        
def read_mut_bed(mut_bed):
    """Read mutation coordinates from a BED file into a per-chromosome dictionary.
        dict looks like this:
        {
            "chr1": [(start1, end1), (start2, end2), ...],
            "chr2": [(start1, end1), (start2, end2), ...],
            ...
        }
    """
    mut_dict = {}
    bed_df = pd.read_csv(mut_bed, sep='\t', skiprows=1, header=None)
    for row_tuple in bed_df.itertuples(index=False):
        chrom = row_tuple[0]
        mut_start = int(row_tuple[1])
        mut_end = int(row_tuple[2])
        if chrom not in mut_dict:
            mut_dict[chrom] = []  
        mut_dict[chrom].append((mut_start, mut_end))
    total_mutations = sum(len(coords) for coords in mut_dict.values())
    print(f"Found {total_mutations} mutations in the file {mut_bed}")
    return mut_dict

'''This function checks if m6A sites overlap with mutations.
dict looks like this:
[
    {
        "chromosome": "chr1",
        "m6a_site": m6a_site_coord,
        "description": description,
        "mut_overlap": True/False,
        "mut_start": mut_start if overlap else None,
        "mut_end": mut_end if overlap else None
    },...
]'''

def check_mutation_for_m6a(m6a_sites, mut_coords):
    overlapping_mutations = []
    for m6a_chrom, data in m6a_sites.items():
        for description, list_of_m6a_coords in data.items():
            for m6a_site_coord in list_of_m6a_coords:
                # Assume no overlap initially for the current m6A site
                has_overlap = False
                # Check if this m6A site overlaps any mutation on the same chromosome
                if m6a_chrom in mut_coords:
                    for mut_start, mut_end in mut_coords[m6a_chrom]:
                        if mut_start <= m6a_site_coord < mut_end:
                            has_overlap = True
                            break
                # Append result only once per m6A site
                overlapping_mutations.append({
                    "chromosome": m6a_chrom,
                    "m6a_site": m6a_site_coord,
                    "description": description,
                    "mut_overlap": has_overlap,
                    "mut_start": mut_start if has_overlap else "NA",
                    "mut_end": mut_end if has_overlap else "NA"
                })

    print(f"Found {len(overlapping_mutations)} m6A sites processed.")
    return overlapping_mutations


'''This function finds open reading frames (ORFs) in a given RNA sequence.
It searches for ORFs in all three reading frames of the sequence.
It returns a dictionary where the keys are tuples of start and end coordinates of the ORFs. 
dict looks like this:
{(start1, end1): orf_seq1, (start2, end2): orf_seq2, ...}
where start and end are the coordinates of the ORF in the sequence, and orf_seq is the corresponding ORF sequence.
'''
def find_orfs_seqs(seq): #Look for ORFs in all three reading frames
	seq_len = len(seq)
	orfs = {}
	for frame in range(3):
		start = None
		for i in range(frame, seq_len-2, 3):
			codon = seq[i:i+3]
			if len(codon) < 3:
				break  
			if codon == START_CODON and start is None:
				start = i  # Mark the start of a new ORF
			if codon in STOP_CODONS and start is not None:
				orf_end = i +3
				orf_seq = seq[start:orf_end] 
				orfs[start, orf_end] = orf_seq
				start = None  # Reset to find the next ORF
	return orfs


'''This function finds all m6A motifs in the ORFs.
It searches for specific motifs in the sequences of the ORFs.
It returns a list of tuples containing the start and end coordinates of each motif found.

dict looks like this:
[(motif_start1, motif_end1), (motif_start2, motif_end2), ...]
where motif_start and motif_end are the coordinates of the motif in the sequence.

'''

def find_m6_motifs(orfs_data):
	motif_coords = []
	all_motifs = ["AAC", "GAC"]
	for coords, sequence in orfs_data.items():
		orf_start = coords[0]
		orf_end = coords[1]   
		for motif in all_motifs:
				index = 0 #search index
				while True:
					rel_motif_start = sequence.find(motif, index)
					if rel_motif_start == -1:
						break
					abs_motif_start = orf_start + rel_motif_start
					abs_motif_end = abs_motif_start + len(motif)
					motif_coords.append((abs_motif_start, abs_motif_end))
					index = rel_motif_start + 1 
	return motif_coords	
     
'''This function calculates the distance from m6A sites to the stop codon in each ORF.
It takes the ORFs data, m6A sites, and motif coordinates as input. 
It returns a list of dictionaries containing the results, including the m6A site coordinates, ORF start and end coordinates, distance to stop codon, length of ORF, stop codon sequence, and motif information if applicable.'''
def calc_dist2m6(orfs_data, mut_data, motif_coords):
    results = []

    for coords, sequence in orfs_data.items():
        stop_codon = sequence[-3:] if len(sequence) >= 3 else "N/A"
        orf_start, orf_end = coords
        orf_len = len(sequence)
        motifs_in_current_orf = []

        # Find all motifs within this ORF
        for motif_start, motif_end in motif_coords:
            if orf_start <= motif_start and motif_end <= orf_end:
                rel_motif_start = motif_start - orf_start
                rel_motif_end = motif_end - orf_start
                if 0 <= rel_motif_start < rel_motif_end <= len(sequence):
                    motif_seq = sequence[rel_motif_start:rel_motif_end]
                else:
                    motif_seq = "N/A_out_of_bounds"

                motifs_in_current_orf.append({
                    "sequence": motif_seq,
                    "start": motif_start,
                    "end": motif_end
                })

        # Check each m6A site against every motif in the ORF (allow overlaps)
        for m6a_entry in mut_data:
            chromosome = m6a_entry["chromosome"]
            m6a_site_coord = int(m6a_entry["m6a_site"])
            description = m6a_entry["description"]
            mut_overlap = m6a_entry["mut_overlap"]
            mut_start = m6a_entry["mut_start"]
            mut_end = m6a_entry["mut_end"]

            if orf_start <= m6a_site_coord <= orf_end:
                d2stop = orf_end - m6a_site_coord
                for motif_info in motifs_in_current_orf:
                    if motif_info["start"] <= m6a_site_coord < motif_info["end"]:
                        results.append({
                            "chromosome": chromosome,
                            "description": description,
                            "m6a site": m6a_site_coord,
                            "ORF Start": orf_start,
                            "ORF Stop": orf_end,
                            "distance to stop": d2stop,
                            "length of ORF": orf_len,
                            "stop codon": stop_codon,
                            "motif sequence": motif_info["sequence"],
                            "motif start": motif_info["start"],
                            "motif end": motif_info["end"],
                            "mut_overlap": mut_overlap,
                            "mut_start": mut_start,
                            "mut_end": mut_end
                        })
    return results



def main():
    # Set up argument parser
    parser = argparse.ArgumentParser(
        description="Find all stop codons coords in a genome fasta and process m6a sites."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        default="GCA_000002765.3_GCA_000002765_genomic.fna",
        help="Input FASTA file containing the DNA sequence (default: GCA_000002765.3_GCA_000002765_genomic.fna)"
    )
    parser.add_argument(
        "-b", "--bed_file",
        type=str,
        default="m6A_Ring_standardized.bed",
        help="Input BED file containing the m6a sites (default: m6A_Ring_standardized.bed)"
    )
    parser.add_argument(
        "-m", "--mut_file",
        type=str,
        default="MUT_uniq_standardized.bed",
        help="Name and location of UTR gff file (default: MUT_uniq_standardized.bed)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="scan4m6a_out.tsv",
        help="Name and location of desired output file (default: scan4m6a_out.tsv)"
    )

    args = parser.parse_args()

    # Get the file paths from parsed arguments
    input_fasta_path = args.input
    bed_file_path = args.bed_file
    output_file_path = args.output
    mut_file_path = args.mut_file


    # --- Input File Existence Checks ---
    if not os.path.exists(input_fasta_path):
        print(f"Error: The input FASTA file '{input_fasta_path}' does not exist.")
        exit(1) 

    if not os.path.exists(bed_file_path):
        print(f"Error: The BED file '{bed_file_path}' does not exist.")
        exit(1) 

    if not os.path.exists(mut_file_path):
        print(f"Error: The BED file '{mut_file_path}' does not exist.")
        exit(1) 

    output_directory = os.path.dirname(output_file_path)
    if output_directory and not os.path.exists(output_directory):
        print(f"Creating output directory: '{output_directory}'")
        try:
            os.makedirs(output_directory)
        except OSError as e:
            print(f"Error creating directory '{output_directory}': {e}")
            exit(1) 


    # Read the DNA sequence from the specified input file
    all_sequences = read_sequence_from_file(input_fasta_path)
    #print(f"Found {len(all_sequences)} chromosomes in the input file {input_fasta_path}")
    #print(all_sequences.keys())  # Print chromosome IDs for debugging

    m6a_sites = read_bed(bed_file_path)
    mut_coords = read_mut_bed(mut_file_path)
    overlapping_mutations = check_mutation_for_m6a(m6a_sites, mut_coords)
    print(Counter([entry["mut_overlap"] for entry in overlapping_mutations]))

    all_results = []

    # Process each chromosome one by one
    for chrom_id, rna_sequence in all_sequences.items():
        print(f"Processing {chrom_id}")
        orfs = find_orfs_seqs(rna_sequence)
        m6_motifs = find_m6_motifs(orfs)
    # Filter mut_data to just this chromosome for better performance
        chrom_specific_mut_data = [entry for entry in overlapping_mutations if entry["chromosome"] == chrom_id]
        results = calc_dist2m6(orfs, chrom_specific_mut_data, m6_motifs)
        if not results:
            print(f"No m6A hits found for {chrom_id}")
        else:
            all_results.extend(results)
    
    print(f"Total results collected: {len(all_results)}")
    df = pd.DataFrame(all_results)
    # Reorder columns so "chromosome" is first
    columns = ["chromosome"] + [col for col in df.columns if col != "chromosome"]
    df = df[columns]
    fin_df = df.drop_duplicates()
    fin_df.to_csv(output_file_path, sep='\t', index=False)

# Main entry point
if __name__ == "__main__":
    main()
