#!/usr/bin/env python3

import argparse
import pandas as pd
import collections

parser = argparse.ArgumentParser(description="Script to separate orthogroups containing chimeric genes from orthogroups containing potetial fused genes.")
parser.add_argument("--input", "-i", required=True, metavar="input", help="File containing gene orthologs. Format is: col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--output_chimeric", "-oc", required=True, metavar="output_chimeric", help="Path to output file containing chimeric orthogroups.")
parser.add_argument("--output_fused", "-of", required=True, metavar="output_fused", help="Path to output file containing potentially fused orthogroups.")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output_chimeric = args.output_chimeric
my_output_fused = args.output_fused

##### Main #####
#Read input
OG_df = pd.read_table(str(my_input), sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#Isolate genes shared between two or more OGs
chimeric_genes = [element for element, my_count in collections.Counter(list(OG_df["GeneID"])).items() if my_count >= 2]
#Generate a list of OG_IDs containing the chimeric genes
chimeric_OGs = list(set(list(OG_df[OG_df["GeneID"].isin(chimeric_genes)]["OG_ID"])))
#Filter the original df with those OG_IDs
chimeric_OGs_df = OG_df[OG_df["OG_ID"].isin(chimeric_OGs)]
#Add chimeric gene label
chimeric_OGs_df["Chimeric_label"] = ["chimeric" if element in chimeric_genes else "not_chimeric" for element in list(chimeric_OGs_df["GeneID"])]
#Isolate OGs with chimeric genes coming from more than one species
chimeric_OGs_list = list(chimeric_OGs_df[chimeric_OGs_df["Chimeric_label"]=="chimeric"][["OG_ID", "Species"]].drop_duplicates()["OG_ID"])
#In the same OG, there are chimeric genes from at least two species
pot_fused_OG_list = [element for element, my_count in collections.Counter(chimeric_OGs_list).items() if my_count >=2]
#The genes from the different species have to be shared between the same OGs
pot_fused_chimeric_genes_df = chimeric_OGs_df[(chimeric_OGs_df["Chimeric_label"]=="chimeric") & (chimeric_OGs_df["OG_ID"].isin(pot_fused_OG_list))]
#Isolate the OG sharing the exact same chimeric genes from different species
grouped_df = pot_fused_chimeric_genes_df.groupby("GeneID")
chimeric_gene_OG_dict = {}
for gene, group in grouped_df:
  chimeric_gene_OG_dict[gene] = tuple(sorted(list(group["OG_ID"])))
duplicated_OGs_full = [list(element) for element, my_count in collections.Counter(list(chimeric_gene_OG_dict.values())).items() if my_count >=2 and len(element)>=2]
#Get a dicionary with the chimeric genes contained by each potentially fused OG        
grouped_OG_df = chimeric_OGs_df[chimeric_OGs_df["Chimeric_label"]=="chimeric"].groupby("OG_ID")
OG_chimeric_gene_dict = {}
for OG, group in grouped_OG_df:
  OG_chimeric_gene_dict[OG] = tuple(sorted(list(group["GeneID"])))
duplicated_OGs_double = [element for element in duplicated_OGs_full if len(element)==2 and OG_chimeric_gene_dict[element[0]]==OG_chimeric_gene_dict[element[1]]]
#Make sure that the chimeric genes shared between more orthogroups are not the single chimeric genes in others
grouped_df = chimeric_OGs_df[chimeric_OGs_df["Chimeric_label"]=="chimeric"].groupby("GeneID") #Build a dictionary where ALL chimeric OG (not only the potentially fused) are considered
chimeric_gene_OG_dict_all = {}
for gene, group in grouped_df:
  chimeric_gene_OG_dict_all[gene] = tuple(sorted(list(group["OG_ID"])))

duplicated_OGs_flatten = []
for OG_pair in duplicated_OGs_double:
  OG_gene_list = OG_chimeric_gene_dict[OG_pair[0]]
  if len([element for element in OG_gene_list if list(chimeric_gene_OG_dict_all[element])==list(OG_pair)]) == len(OG_gene_list):
    duplicated_OGs_flatten = duplicated_OGs_flatten+OG_pair

final_chimeric_OGs_df = chimeric_OGs_df[~(chimeric_OGs_df["OG_ID"].isin(duplicated_OGs_flatten))]
final_pot_fused_chimeric_OGs_df = chimeric_OGs_df[chimeric_OGs_df["OG_ID"].isin(duplicated_OGs_flatten)]
#Save to output file
final_chimeric_OGs_df.to_csv(str(my_output_chimeric), sep="\t", index=False, header=True, na_rep="NA")
final_pot_fused_chimeric_OGs_df.to_csv(str(my_output_fused), sep="\t", index=False, header=True, na_rep="NA")
