#!/usr/bin/python3
import deg_module as dm
import pandas as pd  



def main():

    # Specify the path to your TSV file
    upset_file = './DEGs/8hr/DEGs_heinz37C_8hr/input4complexUpset_8hr_vs_vs_Heinz_8hr_37C.tsv'
    # Specify the columns and pattern
    required_columns = ['Heinz 8hr 25C vs Heinz 37C', 'Mal 8hr 25C vs Heinz 37C', 'Nag 8hr 25C vs Heinz 37C', 'Tam 8hr 25C vs Heinz 37C']
    pattern_to_match = [0, 1, 1, 1]
    # Find the genes matching the pattern
    matching_genes = dm.find_genes_matching_pattern(upset_file, required_columns, pattern_to_match)
    print(len(matching_genes))

    # initialize dfs
    Heinz_df = pd.read_csv("./DEGs/8hr/DEGs_heinz37C_8hr/Heinz_DEG_8hr_25C_vs_Heinz_8hr_37C.tsv", sep="\t", index_col = 0)
    Mal_df = pd.read_csv("./DEGs/8hr/DEGs_heinz37C_8hr/Mal_DEG_8hr_25C_vs_Heinz_8hr_37C.tsv", sep="\t", index_col = 0)
    Nag_df = pd.read_csv("./DEGs/8hr/DEGs_heinz37C_8hr/Nag_DEG_8hr_25C_vs_Heinz_8hr_37C.tsv", sep="\t", index_col = 0)
    Tam_df = pd.read_csv("./DEGs/8hr/DEGs_heinz37C_8hr/Tam_DEG_8hr_25C_vs_Heinz_8hr_37C.tsv", sep="\t", index_col = 0)

    """
    # Run for all varieties
    df_list = [Heinz_df, Mal_df, Nag_df, Tam_df]
    varieties = ["Heinz", "Mal", "Nag", "Tam"]
    results = dm.combine_deg_data(df_list, varieties, matching_genes)
    results.index = results.index.str.replace("gene:", "")
    result_table = results.to_csv("./results/tables_8hr/genes_vs_H37_8hr/HMNT_genes_H37C_8hr.tsv", sep='\t')
    print(results)
    
    """

    #run for only some
    df_list = [Mal_df, Nag_df, Tam_df]
    varieties = ["Mal", "Nag", "Tam"]
    results = dm.combine_deg_data(df_list, varieties, matching_genes)
    results.index = results.index.str.replace("gene:", "")
    result_table = results.to_csv("./results/tables_8hr/genes_vs_H37_8hr/MNT_genes_H37_8hr.tsv", sep='\t')
    print(results)


    """

    #run for just one variety
    df_list = [Tam_df]
    varieties = ["Tam"]
    results = dm.combine_deg_data(df_list, varieties, matching_genes)
    results.index = results.index.str.replace("gene:", "")
    result_table = results.to_csv("./results/tables_8hr/genes_vs_H37_8hr/T_genes_H37_8hr.tsv", sep='\t')
    print(results)
    """

# Main entry point
if __name__ == "__main__":
    main()