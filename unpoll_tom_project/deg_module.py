#!/usr/bin/python3
import os
import pandas as pd

"""
    Extracts relevant information (gene ID, LFC, padj, baseMean, and direction) from a DataFrame.

    Args:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        dict: A dictionary where keys are gene IDs and values are dictionaries
              containing 'LFC', 'padj', 'baseMean', and 'Direction'.
"""
# Set the first column (gene ID) as index
def gather_data(df):
    gene_dict = {}
    genes = df.index.tolist() #get geneid
    LFC = df["log2FoldChange"].tolist()
    LFCSE = df["lfcSE"]
    padj = df["padj"].tolist()
    pval = df["pvalue"].tolist()
    basemean = df["baseMean"].tolist()
    direction = df["Direction"].tolist()
    for i, gene in enumerate(genes):
        gene_dict[gene] = {"LFC": LFC[i], "lfcSE": LFCSE[i], "padj": padj[i], "pvalue": pval[i], "baseMean": basemean[i], "Direction": direction[i]}
    return gene_dict

"""
    Finds genes in a TSV file that match a specified pattern across given columns.

    Args:
        file_path (str): Path to the TSV file.
        columns (list of str): List of column names to check for the pattern.
        pattern (list of int): The pattern to match (e.g., [0, 1, 1, 1]).

    Returns:
        list: A list of gene IDs that match the pattern.
              Returns an empty list if no genes match or if there's an error.
"""
def find_genes_matching_pattern(file_path, columns, pattern):
    try:
        # Read the TSV file into a pandas DataFrame
        df = pd.read_csv(file_path, sep='\t', index_col=0)

        # Check if the specified columns exist in the DataFrame
        if not all(col in df.columns for col in columns):
            print(f"Error: One or more of the required columns ({columns}) were not found in the file.")
            return []

        # Create a boolean mask for the pattern match
        pattern_match = True  # Initialize to True for the first column
        for i, col in enumerate(columns):
            pattern_match &= (df[col] == pattern[i])

        # Get the gene IDs that match the pattern
        matching_genes = df[pattern_match].index.tolist()
        return matching_genes

    except FileNotFoundError:
        print(f"Error: The file was not found at the specified path: {file_path}")
        return []
    except Exception as e:
        print(f"An error occurred: {e}")
        return []


def combine_deg_data(dataframes, variety_names, common_genes):
    """
    Combines pandas DataFrames of differentially expressed genes into a single DataFrame,
    containing only the overlapping genes across the specified varieties,
    considering only unique gene entries.

    Args:
        dataframes (list of pd.DataFrame): A list of DataFrames for the varieties to combine.
        variety_names (list of str): A list of names for the varieties.

    Returns:
        pd.DataFrame: A new DataFrame containing the overlapping genes with columns
                        formatted as 'variety_baseMean', 'variety_LFC', 'variety_padj',
                        and 'variety_direction' for each variety. Returns an empty DataFrame
                        if no common genes are found.

    Raises:
        TypeError: If any input is not a pandas DataFrame.
        ValueError: If the number of variety names does not match the number of input DataFrames.
    """
    # Ensure input validity
    if not all(isinstance(df, pd.DataFrame) for df in dataframes):
        raise TypeError("All inputs must be pandas DataFrames.")
    if len(variety_names) != len(dataframes):
        raise ValueError("The number of variety names must match the number of input DataFrames.")

    # Extract and combine data
    data_dicts = [gather_data(df) for df in dataframes]
    data_dict_combined = {}

    for gene in common_genes:
        row_data = {}
        for i, variety_name in enumerate(variety_names):
            data_dict = data_dicts[i]
            if gene in data_dict:
                row_data[f'{variety_name}_baseMean'] = data_dict[gene]['baseMean']
                row_data[f'{variety_name}_LFC'] = data_dict[gene]['LFC']
                row_data[f"{variety_name}_lfcSE"] = data_dict[gene]["lfcSE"]
                row_data[f'{variety_name}_padj'] = data_dict[gene]['padj']
                row_data[f'{variety_name}_pvalue'] = data_dict[gene]['pvalue']
                row_data[f'{variety_name}_direction'] = data_dict[gene]['Direction']
        data_dict_combined[gene] = row_data

    # Create final merged DataFrame
    combined_df = pd.DataFrame.from_dict(data_dict_combined, orient='index')
    return combined_df