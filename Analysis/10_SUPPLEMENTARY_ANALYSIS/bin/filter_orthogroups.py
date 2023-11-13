#!/usr/bin/env python3

import argparse
import pandas as pd
from collections import Counter
import re

parser = argparse.ArgumentParser(description="Script to filter the gene orthogroups for a maximum number of genes. Bigger orthogroups are considered only if the proportions between the represented species are balanced")
parser.add_argument("--input", "-i", required=True, metavar="input_file", help="Orthogroup file, with format col1=OrthogoupID, col2=Species, col3=GeneID")
parser.add_argument("--max_genes", "-g", required=True, metavar="max_genes", help="Maximum number of genes allowed in the orthogroup")
parser.add_argument("--max_proportion", "-p", required=True, metavar="max_proportion", help="Maximum number of genes allowed in the orthogroup")
parser.add_argument("--output", "-o", required=True, metavar="output_file", help="Two columns file contaning with col1=ProteinID, col2=ProteinLength")

#Read arguments
args = parser.parse_args()
input_file = args.input
max_genes = args.max_genes
max_proportion = args.max_proportion
output_file = args.output

####################################
######### MAIN #####################
####################################

OG_df = pd.read_table(input_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#Get list of orthogroups
all_OG_list = list(OG_df["OG_ID"])
#Count number of genes in orthogroups
OG_gene_count_dict = dict(Counter(all_OG_list))
#Isolate all OGs with less than max_genes genes
saved_OGs = [key for key in list(OG_gene_count_dict.keys()) if OG_gene_count_dict[key] <= int(max_genes)]
saved_OGs_df = OG_df.loc[OG_df["OG_ID"].isin(saved_OGs)]
#Isolate all OGs with more than max_genes genes
unsaved_OGs = [key for key in list(OG_gene_count_dict.keys()) if OG_gene_count_dict[key] > int(max_genes)]
#Save the OGs with more than max_genes genes if the max proportion among the species is <= max_proportion
unsaved_OGs_df = OG_df.loc[OG_df["OG_ID"].isin(unsaved_OGs)]
unsaved_OGs_grouped_df = unsaved_OGs_df.groupby("OG_ID")
for OG_ID, group in unsaved_OGs_grouped_df:
  tot_genes = group.shape[0]
  species_count_dict = dict(Counter(list(group["Species"])))
  all_proportions_list = [count/tot_genes for count in list(species_count_dict.values())]
  group_max_proportion = max(all_proportions_list)
  if group_max_proportion <= float(max_proportion) and tot_genes <= 120: #in any case it has to have less than 120. Otherwise it is just crazy.
    saved_OGs_df = pd.concat([saved_OGs_df, group])
  else:
    print(OG_ID)

#Save to output file
saved_OGs_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
