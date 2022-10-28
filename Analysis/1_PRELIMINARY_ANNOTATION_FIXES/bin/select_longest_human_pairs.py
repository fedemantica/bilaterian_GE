#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to select the longest human orthologs for each broken gene.")
parser.add_argument("--human_pairs", "-hp", required=True, metavar="human_pairs", help="File containing all possible pairs of broken genes and relative orthologs. Format is col1=OG_ID, col2=GeneID_broken, col3=GeneID_human")
parser.add_argument("--prot_length", "-p", required=True, metavar="prot_length", help="File containing length information for human reference proteins. Format is col1=GeneID, col2=Length")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All possible pairs of broken genes and relative human orthologs.")

##### Read arguments #####
args = parser.parse_args()
my_human_pairs = args.human_pairs
my_prot_length = args.prot_length
my_output = args.output

##### Main #####
#Read inputs
gene_pairs_df = pd.read_table(str(my_human_pairs), sep="\t", index_col=False, header=None, names=["OG_ID", "GeneID", "Hsa_geneID"])
prot_length_df = pd.read_table(str(my_prot_length), sep="\t", index_col=False, header=None, names=["GeneID", "Prot_length"])
#Create dictionary {geneID : prot_length}
prot_length_dict = pd.Series(prot_length_df.Prot_length.values, index=prot_length_df.GeneID)
grouped_df = gene_pairs_df.groupby("OG_ID")
#Generate final df
final_df = pd.DataFrame(columns=list(gene_pairs_df.columns))
#Select the longest protein for each group (OG_ID)
for name, group in grouped_df:
  hsa_genes = list(set(list(group["Hsa_geneID"])))
  length_list = [prot_length_dict[element] for element in hsa_genes]
  my_max = max(length_list)
  my_prot_length_dict = {key:prot_length_dict[key] for key in hsa_genes}
  longest_gene = [key for key in list(my_prot_length_dict.keys()) if my_prot_length_dict[key]==my_max][0]

  group_df = group.loc[group.Hsa_geneID==longest_gene]
  #Concatenate to final df
  final_df = pd.concat([final_df, group_df])
#Save final_df to file
final_df.to_csv(str(my_output), sep="\t", header=False, index=False, na_rep="NA")
