#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import Counter

parser = argparse.ArgumentParser(description="Script to select all the orthogroups containing two or more broken genes")
parser.add_argument("--broken_genes", "-b", required=True, metavar="broken_genes", help="Paralogous genes sequentially annotated on the same chr and strand. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--orthogroups", "-o", required=True, metavar="orthogroups", help="Complete gene orthogroups. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species identifier as reported in broken_genes and orthogroups inputs")
parser.add_argument("--filtered_orthogroups", "-f", required=True, metavar="filtered_orthogroups", help="Path to output file. Filtered orthogroups containing two or more broken genes. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--broken_pairs", "-p", required=True, metavar="broken_pairs", help="Path to output file. Pairs of broken genes.")

##### Read arguments #####
args = parser.parse_args()
my_broken_genes = args.broken_genes
my_orthogroups = args.orthogroups
my_species = args.species
my_filtered_orthogroups = args.filtered_orthogroups
my_broken_pairs = args.broken_pairs

##### Main #####
broken_genes_df = pd.read_table(str(my_broken_genes), sep="\t", header=None, names=["OG_ID", "Species", "chr:strand", "gene1", "gene1_comb", "gene2", "gene2_comb"])
parsed_orthogroups = pd.read_table(str(my_orthogroups), sep="\t", header=None, names=["OG_ID", "Species", "GeneID", "GeneName"])
#Create dictionary with key=GeneID, value=OG_ID
geneID_OG_ID_dict = pd.Series(parsed_orthogroups.OG_ID.values, index=parsed_orthogroups.GeneID).to_dict()
#Create dictionary with key=OG_ID, value=number of broken pieces
OG_gene_counts_dict = dict(Counter(list(broken_genes_df["OG_ID"])))
#select all the OG_IDs of orthogroups containing 2 or more broken gene parts
OGs_to_save = [element for element in list(OG_gene_counts_dict.keys()) if OG_gene_counts_dict[element] >= 1]
#Filter only for orthogroups with given number of broken genes
final_df = broken_genes_df.loc[broken_genes_df.OG_ID.isin(OGs_to_save)]
final_df = final_df[["OG_ID", "gene1", "gene2"]]
final_df.to_csv(str(my_broken_pairs), sep="\t", na_rep="NA", index=False, header=False) #save to file
#print the genes one by one
all_genes = list(final_df["gene1"])+list(final_df["gene2"])
single_genes_df = pd.DataFrame({"GeneID" : all_genes})
single_genes_df["Species"] = my_species
single_genes_df["OG_ID"] = list(single_genes_df["GeneID"].map(geneID_OG_ID_dict))
single_genes_df = single_genes_df.loc[:,["OG_ID", "Species", "GeneID"]]
#save to file
single_genes_df.to_csv(str(my_filtered_orthogroups), sep="\t", header=False, index=False, na_rep="NA")
