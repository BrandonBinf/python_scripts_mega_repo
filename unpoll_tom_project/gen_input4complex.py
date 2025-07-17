"This code is used to generate inputs for complex upset plots in R"

import os
import pandas as pd
work_dir = r"C:\Code_projects\Work_stuff\R_work\Kelsie_NEXTFLOW\Unpollinated_deseq\Complex_Upset"
os.chdir(work_dir)

def gather_data(list):
    needed_info_list = []
    for df in list:
        gene_dict = {}
        genes = df.iloc[:, 0].index.tolist() #get geneid 
        LFC = df["log2FoldChange"].tolist() #get LFC
        padj = df["padj"].tolist() #get padj
        for i, gene in enumerate(genes):
            gene_dict[gene] = {"LFC": LFC[i], "padj": padj[i]} #creates a nested dict, gene id is key, values are two other dicts with LFC and padj
        needed_info_list.append(gene_dict) #make massive list of all needed pieces
    return needed_info_list

def set_up_table(needed_info):
    data = {"Heinz 3hr 25Cvs37C": [], "Mal 3hr 25Cvs37C": [], "Nag 3hr 25Cvs37C": [], "Tam 3hr 25Cvs37C": []}
    treatment_names = list(data.keys()) # Get the treatment names in order
    #print(treatment_names)
    
    for i, dict_item in enumerate(needed_info): #enumerate to get the order
        treatment = treatment_names[i] #get the treatment name based on the index.
        gene_results = {} #dictionary for the gene results of each treatment.
        for gene, values in dict_item.items():
            LFC_val = values["LFC"]
            padj_val = values["padj"]
            if (LFC_val > 1 or LFC_val < -1) and padj_val < 0.05: #find DE genes
                gene_results[gene] = 1
            else:
                gene_results[gene] = 0
        data[treatment].append(gene_results) #append the dictionary to the list corresponding to the treatment.
    return data
                
def make_binary_table(data):
    all_genes = set()
    for treatment_data in data.values():
        if treatment_data: #check if the list is not empty
            all_genes.update(treatment_data[0].keys())

    all_genes = sorted(list(all_genes)) #sort the genes for consistent order

    df = pd.DataFrame(index=all_genes)

    for treatment, treatment_data in data.items():
        if treatment_data: #check if the list is not empty
            gene_results = treatment_data[0]
            df[treatment] = df.index.map(lambda gene: gene_results.get(gene, 0)) #fill with 0 if gene is absent.
        else:
            df[treatment] = 0 #fill with zeros if there is no gene data for the treatment.

    df.to_csv('input4complexUpset_3hr.tsv', sep='\t', index=True)
    return df
    


def main():
    #initialize dfs
    Heinz_df = pd.read_csv(".//DEGs_3hr//Heinz_DEG_3hr_25C_vs_37C.tsv", sep="\t")
    Mal_df = pd.read_csv(".//DEGs_3hr//Mal_DEG_3hr_25C_vs_37C.tsv", sep="\t")
    Nag_df = pd.read_csv(".//DEGs_3hr//Nag_DEG_3hr_25C_vs_37C.tsv", sep="\t")
    Tam_df = pd.read_csv(".//DEGs_3hr//Tam_DEG_3hr_25C_vs_37C.tsv", sep="\t")

#organize dfs into a list
    df_list = [Heinz_df, Mal_df, Nag_df, Tam_df]
    id_LFC_padj_dict = gather_data(df_list)
    binary_data = set_up_table(id_LFC_padj_dict)
    fin_results = make_binary_table(binary_data)
    
# Main entry point
if __name__ == "__main__":
    main()
