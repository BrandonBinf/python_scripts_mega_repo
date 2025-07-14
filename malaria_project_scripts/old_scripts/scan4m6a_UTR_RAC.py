#!/usr/bin/env python3

import os, argparse
import pandas as pd
from Bio import SeqIO

START_CODON = "AUG"
STOP_CODONS = {"UAA", "UAG", "UGA"}

def read_sequence_from_file(filename):
    fasta_data = {}
    with open(filename, "r") as file:
        for record in SeqIO.parse(file, "fasta"):
            annotation = record.description
            parts = annotation.split(",")
            chromosome = parts[1].strip()
            dna_seq = record.seq.upper()
            rna_seq = str(dna_seq.transcribe())
            fasta_data[chromosome] = rna_seq
    return fasta_data

def read_UTR_file(UTR_file):
    utr_dict = {}
    bed_df = pd.read_csv(UTR_file, sep='\t', skiprows=1, header=None)
    for row in bed_df.itertuples(index=False):
        chrom = row[0]
        utr_start = int(row[3])
        utr_end = int(row[4])
        if chrom not in utr_dict:
            utr_dict[chrom] = []
        utr_dict[chrom].append((utr_start, utr_end))
    return utr_dict

def pull_UTR_sequences(rna_sequences, UTR_coords):
    UTR_sequences = {}
    for chrom, coords in UTR_coords.items():
        if chrom not in rna_sequences:
            continue
        for start, end in coords:
            UTR_seq = rna_sequences[chrom][start:end]
            if chrom not in UTR_sequences:
                UTR_sequences[chrom] = []
            UTR_sequences[chrom].append((start, end, str(UTR_seq)))
    return UTR_sequences

def read_and_format_bed(bed_file):
    CHROM_COL_INDEX = 0
    M6A_SITE_COL_INDEX = 1
    DESCRIPTION_COL_INDEX = 3
    formatted_m6a_data = []
    bed_df = pd.read_csv(bed_file, sep='\t', skiprows=1, header=None)
    for row in bed_df.itertuples(index=False):
        chrom = row[CHROM_COL_INDEX]
        m6a_site = int(row[M6A_SITE_COL_INDEX])
        description = row[DESCRIPTION_COL_INDEX]
        formatted_m6a_data.append({
            "chromosome": chrom,
            "m6a_site": m6a_site,
            "description": description,
        })
    print(f"Found {len(formatted_m6a_data)} m6A sites processed.")
    return formatted_m6a_data

def find_m6a_sites_UTR(UTR_sequences, m6a_data):
    m6a_sites_by_utr = {}
    for chrom, utr_list in UTR_sequences.items():
        for start, end, utr_seq in utr_list:
            for entry in m6a_data:
                if entry["chromosome"] != chrom:
                    continue  # Only match m6A to UTRs on the same chromosome
                site = entry["m6a_site"]
                if start <= site < end:
                    coords = (chrom, start, end)
                    if coords not in m6a_sites_by_utr:
                        m6a_sites_by_utr[coords] = []
                    m6a_sites_by_utr[coords].append({
                        "chromosome": chrom,
                        "utr_start": start,
                        "utr_end": end,
                        "description": entry["description"],
                        "m6a_site": site,
                        "utr_sequence": utr_seq
                    })
    return m6a_sites_by_utr

def find_m6a_motifs_in_utr(m6a_data):
    final_results = {}
    motifs = ["AAC", "GAC"]
    for coords, entries in m6a_data.items():
        for entry in entries:
            chrom = entry["chromosome"]
            seq = entry["utr_sequence"]
            start = entry["utr_start"]
            end = entry["utr_end"]
            m6a_site = entry["m6a_site"]
            desc = entry["description"]

            for motif in motifs:
                index = 0
                while True:
                    rel_start = seq.find(motif, index)
                    if rel_start == -1:
                        break
                    abs_start = start + rel_start
                    abs_end = abs_start + len(motif)
                    if abs_start <= m6a_site < abs_end:
                        if coords not in final_results:
                            final_results[coords] = []
                        final_results[coords].append({
                            "chromosome": chrom,
                            "utr_start": start,
                            "utr_end": end,
                            "description": desc,
                            "m6a_site": m6a_site,
                            "motif": motif,
                            "motif_start": abs_start,
                            "motif_end": abs_end
                        })
                    index = rel_start + 1
    return final_results

def main():
    parser = argparse.ArgumentParser(
        description="Find all stop codon coords in a genome fasta and process m6a sites."
    )
    parser.add_argument(
        "-i", "--input",
        type=str,
        default="/home/bbendick/malaria_work/GCA_000002765.3_GCA_000002765_genomic.fna",
        help="Input FASTA file containing the DNA sequence"
    )
    parser.add_argument(
        "-b", "--bed_file",
        type=str,
        default="/home/bbendick/malaria_work/m6A_Ring_standardized.bed",
        help="Input BED file containing the m6a sites"
    )
    parser.add_argument(
        "-u", "--utr_file",
        type=str,
        default="/home/bbendick/malaria_work/Llinas_3-UTR_standard.gff",
        help="Name and location of UTR GFF file"
    )
    parser.add_argument(
        "-o", "--output",
        type=str,
        default="/home/bbendick/malaria_work/scan4m6a_UTR.tsv",
        help="Name and location of desired output file"
    )

    args = parser.parse_args()

    # File checks
    for path in [args.input, args.bed_file, args.utr_file]:
        if not os.path.exists(path):
            print(f"Error: File not found - {path}")
            return

    output_directory = os.path.dirname(args.output)
    if output_directory and not os.path.exists(output_directory):
        try:
            os.makedirs(output_directory)
        except OSError as e:
            print(f"Error creating directory '{output_directory}': {e}")
            return

    # Load all input data
    all_sequences = read_sequence_from_file(args.input)
    utr_coords = read_UTR_file(args.utr_file)
    m6a_sites = read_and_format_bed(args.bed_file)

    # Container for all matching results
    all_results = []

    # Process each chromosome separately
    for chrom_id, rna_sequence in all_sequences.items():
        print(f"Processing chromosome: {chrom_id}")
        if chrom_id not in utr_coords:
            continue

        # Filter UTRs and m6A to this chromosome only
        chrom_utr_coords = {chrom_id: utr_coords[chrom_id]}
        chrom_rna_seq = {chrom_id: rna_sequence}

        # Correct function calls:
        utr_sequences = pull_UTR_sequences(chrom_rna_seq, chrom_utr_coords)
        m6a_by_utr = find_m6a_sites_UTR(utr_sequences, m6a_sites)
        motif_results = find_m6a_motifs_in_utr(m6a_by_utr)

        for entries in motif_results.values():
            all_results.extend(entries)

    # Save results
    print(f"Total results collected: {len(all_results)}")
    if not all_results:
        print("No motifs found overlapping m6A sites in UTRs.")
        return

    df = pd.DataFrame(all_results)
    df = df[["chromosome"] + [col for col in df.columns if col != "chromosome"]]
    df.drop_duplicates().to_csv(args.output, sep='\t', index=False)

if __name__ == "__main__":
    main()
