#!/usr/bin/python3
import pandas as pd

def find_genes_matching_pattern(file_path, columns, pattern):
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

if __name__ == "__main__":
    # Specify the path to your TSV file
    file_path = 'input4complexUpset_3hr.tsv'  # Replace 'your_file.tsv' with the actual path

    # Specify the columns and pattern
    required_columns = ['Heinz 3hr 25Cvs37C', 'Mal 3hr 25Cvs37C', 'Nag 3hr 25Cvs37C', 'Tam 3hr 25Cvs37C']
    pattern_to_match = [0, 1, 1, 1]

    # Find the genes matching the pattern
    matching_genes = find_genes_matching_pattern(file_path, required_columns, pattern_to_match)

    # Print the results
    if matching_genes:
        print(f"Genes matching the pattern {pattern_to_match}:")
        print(matching_genes)
        print(f"Number of genes matching the pattern: {len(matching_genes)}")
    else:
        print(f"No genes found matching the pattern {pattern_to_match}.")