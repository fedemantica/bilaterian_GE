#!/usr/bin/env python3

import argparse
import pandas as pd
import itertools

parser = argparse.ArgumentParser(description="Script to generate all pairs between broken genes and relative human orthologs.")
parser.add_argument("--broken_genes", "-b", required=True, metavar="broken_genes", help="Paralogous genes sequentially annotated on the same chr and strand. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Complete gene orthogroups. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species identifier as reported in broken_genes and orthogroups inputs")
parser.add_argument("--human_id", "-hs", required=True, metavar="human_id", help="Human species identifier as reported in  orthogroups inputs")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All possible pairs of broken genes and relative human orthologs.")

##### Read arguments #####
args = parser.parse_args()
my_broken_genes = args.broken_genes
my_orthogroups = args.orthogroups
my_species = args.species
my_human_id = args.human_id
my_output = args.output

##### Main #####
#Read inputs
broken_genes_df = pd.read_table(str(my_broken_genes), sep="\t", header=None, names=["OG_ID", "Species", "GeneID"])
orthogroups_df = pd.read_table(str(my_orthogroups), sep="\t", header=None, names=["OG_ID", "Species", "GeneID", "GeneName"])

#Select OG_IDs of interesting orthogroups (i.e. containing broken genes)
OGs_to_save = list(set(list(broken_genes_df["OG_ID"])))
#Subset general orthogroups for the human genes in the interesting orthogroups
filtered_orthogroups_df = orthogroups_df.loc[(orthogroups_df.OG_ID.isin(OGs_to_save)) & (orthogroups_df.Species==str(my_human_id))]
combined_df = pd.concat([broken_genes_df, filtered_orthogroups_df])
#Group dataframe by OG_ID
my_grouped_df = combined_df.groupby("OG_ID")
final_df = pd.DataFrame(columns=["OG_ID", "broken_gene", "human_gene"])
#Cycle on each orthogroup
for name, group in my_grouped_df:
  human_genes = list(group.loc[group.Species==str(my_human_id)]["GeneID"]) #list of the human genes in the orthogroup
  broken_genes = list(group.loc[group.Species==my_species]["GeneID"])
  #Generate all possible combinations between the members of the two lists
  all_combs = [(x,y) for x in broken_genes for y in human_genes]
  for element in all_combs:
    element_df = pd.DataFrame({"OG_ID" : [name], "broken_gene" : [element[0]], "human_gene" : [element[1]]})
    final_df = pd.concat([final_df, element_df])
#Save to file
final_df.to_csv(str(my_output), sep="\t", index=False, header=False, na_rep="NA")
