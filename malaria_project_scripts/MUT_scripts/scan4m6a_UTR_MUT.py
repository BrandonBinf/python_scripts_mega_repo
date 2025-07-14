#!/usr/bin/env python3

#fix imports
import os, argparse
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from collections import Counter

# Define the standard genetic code start codon and stop codons
START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}

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

def read_UTR_file(UTR_file):
    "Read UTR coordinates from a GFF file."
    utr_dict = {}
    bed_df = pd.read_csv(UTR_file, sep='\t', skiprows=1, header=None)
    for row_tuple in bed_df.itertuples(index=False):
        chrom = row_tuple[0]
        utr_start = int(row_tuple[3])
        utr_end = int(row_tuple[4])
        if chrom not in utr_dict:
            utr_dict[chrom] = []  
        utr_dict[chrom].append((utr_start, utr_end))
    return utr_dict

def read_bed(bed_file):
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

def pull_UTR_sequences(rna_sequence, utr_dict):
    UTR_sequences = {}
    for chromosome, coords in utr_dict.items():
        for start, end in coords:
            UTR_seq = rna_sequence[int(start):int(end)]
            if (start, end) not in UTR_sequences:
                UTR_sequences[(start, end)] = []
            UTR_sequences[(start, end)].append(UTR_seq)
    return UTR_sequences

def find_m6a_sites_UTR(UTR_sequences, m6a_dict):
    m6a_sites_list = []
    for coords, utr_seqs in UTR_sequences.items():
        utr_start, utr_end = coords
        for utr_seq in utr_seqs:
            for m6a_chrom, data in m6a_dict.items():
                for description, m6a_coords in data.items():
                    for m6a_site_coord in m6a_coords:
                        if utr_start <= m6a_site_coord < utr_end:
                            m6a_sites_list.append({
                                "chromosome": m6a_chrom,
                                "utr_start": utr_start,
                                "utr_end": utr_end,
                                "description": description,
                                "m6a_site": m6a_site_coord,
                                "utr_sequence": utr_seq
                            })
    return m6a_sites_list

def check_mutation_for_m6a(m6a_sites_by_UTR, mut_dict):
    overlapping_mutations = []
    for data in m6a_sites_by_UTR:
        utr_start = data["utr_start"]
        utr_end = data["utr_end"]
        m6a_chrom = data["chromosome"]
        m6a_site_coord = data["m6a_site"]
        description = data["description"]
        utr_seq = data["utr_sequence"]
        has_overlap = False
        if m6a_chrom in mut_dict:
            for mut_start, mut_end in mut_dict[m6a_chrom]:
                if mut_start <= m6a_site_coord < mut_end:
                    has_overlap = True
                    mut_start = mut_start
                    mut_end = mut_end
                    break
            overlapping_mutations.append({
                "chromosome": m6a_chrom, "utr_start": utr_start,"utr_end": utr_end, "m6a_site": m6a_site_coord, "description": description,
                "utr_sequence": utr_seq,  "mut_overlap": has_overlap, "mut_start": mut_start if has_overlap else "NA", "mut_end": mut_end if has_overlap else "NA"
            })
    return overlapping_mutations

def find_m6a_motifs_in_utr(overlapping_mutations):
    final_results = {}
    all_motifs = ["AAACA", "AAACC", "AAACU", "AGACA", "AGACC", "AGACU",
                  "GAACA", "GAACC", "GAACU", "GGACA", "GGACC", "GGACU"]  # RRACH
    for entry in overlapping_mutations:
        sequence = entry["utr_sequence"]
        utr_start = entry["utr_start"]
        utr_end = entry["utr_end"]
        coords = (utr_start, utr_end)
        description = entry["description"]
        m6a_site = entry["m6a_site"]
        mut_overlap = entry["mut_overlap"]
        mut_start = entry["mut_start"]
        mut_end = entry["mut_end"]
        # Search for motifs in the UTR sequence
        for motif in all_motifs:
            index = 0  # search index
            while True:
                rel_motif_start = sequence.find(motif, index)
                if rel_motif_start == -1:
                    break
                abs_motif_start = utr_start + rel_motif_start
                abs_motif_end = abs_motif_start + len(motif)
                # Only keep motifs that contain the m6a site
                if abs_motif_start <= m6a_site < abs_motif_end:
                    if coords not in final_results:
                        final_results[coords] = []
                    final_results[coords].append({
                        "chromosome": entry["chromosome"],
                        "utr_start": utr_start,
                        "utr_end": utr_end,
                        "description": description,
                        "m6a_site": m6a_site,
                        "motif": motif,
                        "motif_start": abs_motif_start,
                        "motif_end": abs_motif_end,
                        "mut_overlap": mut_overlap,
                        "mut_start": mut_start if mut_overlap else "NA",
                        "mut_end": mut_end if mut_overlap else "NA"
                    })
                index = rel_motif_start + 1
    return final_results
     
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
        help="Name and location of MUT bed file (default: MUT_uniq_standardized.bed)"
    )
    parser.add_argument(
        "-u", "--utr_file",
        type=str,
        default="Llinas_3-UTR_standard.gff",
        help="Name and location of UTR gff file (default: Llinas_3-UTR_standard.gff)"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="scan4m6a_out_UTR.tsv",
        help="Name and location of desired output file (default: scan4m6a_out_UTR.tsv)"
    )

    args = parser.parse_args()

    # Get the file paths from parsed arguments
    input_fasta_path = args.input
    bed_file_path = args.bed_file
    output_file_path = args.output
    mut_file_path = args.mut_file
    UTR_file_path = args.utr_file


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
    m6a_sites = read_bed(bed_file_path)
    UTR_coords = read_UTR_file(UTR_file_path)
    mut_dict = read_mut_bed(mut_file_path)
    all_results = []
    # Process each chromosome one by one
    for chrom_id, rna_sequence in all_sequences.items():
        print(f"Processing chromosome: {chrom_id}")
        if chrom_id not in UTR_coords:
            continue
        # 1. Pull UTR sequences for this chromosome
        chrom_utr_coords = {chrom_id: UTR_coords[chrom_id]}
        UTR_sequences = pull_UTR_sequences(rna_sequence, chrom_utr_coords)
        # 2. Find m6a sites that fall inside UTRs
        m6a_data = find_m6a_sites_UTR(UTR_sequences, {chrom_id: m6a_sites.get(chrom_id, {})})
        # 3. Check if these m6a sites overlap with mutations
        overlaps = check_mutation_for_m6a(m6a_data, mut_dict)
        # 4. Identify motifs that contain the m6a site
        motif_results = find_m6a_motifs_in_utr(overlaps)
        for entries in motif_results.values():
            all_results.extend(entries)
    print(f"Total results collected: {len(all_results)}")
    df = pd.DataFrame(all_results)
    # Reorder columns so "chromosome" is first
    columns = ["chromosome"] + [col for col in df.columns if col != "chromosome"]
    df = df[columns]
    fin_df = df.drop_duplicates()
    fin_df.to_csv(output_file_path, sep='\t', index=False)

if __name__ == "__main__":
    main()