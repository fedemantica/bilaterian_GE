#!/usr/bin/env python3

import argparse
import math
import pandas as pd

parser = argparse.ArgumentParser(description="Script to filter out genes from orthogroups based on min_expr")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Corrected and filtered gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--expr_table_dir", "-ed", required=True, metavar="expr_table_dir", help="Path to directory containing the tissue average expression tables (in TPMs) by species")
parser.add_argument("--expr_file_suffix", "-es", required=True, metavar="expr_file_suffix", help="Suffix for expression files to upload")
parser.add_argument("--clade_species", "-s", required=True, metavar="clade_species", help="Comma separated list of all the species represented in the orthogroups. They will be used to read in the expression tables")
parser.add_argument("--min_expr", "-m", required=True, metavar="min_expr", help="Minimum expression (in TPMs) in at least one tissue")
parser.add_argument("--output_orthogroups", "-oo", required=True, metavar="output_orthogroups", help="Path to output file with orthogroups")
parser.add_argument("--output_lowly_expressed", "-ol", required=True, metavar="output_lowly_expressed", help="Path to output file with saved lowly expressed genes")
parser.add_argument("--output_removed_genes", "-or", required=True, metavar="output_removed_genes", help="Path to output file with removed genes")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
expr_table_dir = args.expr_table_dir
expr_file_suffix = args.expr_file_suffix
clade_species = args.clade_species.split(",")
min_expr = args.min_expr
output_orthogroups_file = args.output_orthogroups
output_lowly_expressed_file = args.output_lowly_expressed
output_removed_genes_file = args.output_removed_genes

##################################
###### READ INPUTS ###############
##################################
#Read orthogroups table
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])

#Read in expression tables and create a unique dataframe (fill with NAs in case of missing values)
all_species_expr_df = pd.DataFrame()
for species in clade_species:
  species_expr_file = expr_table_dir + "/" + species + expr_file_suffix
  species_expr_df = pd.read_table(species_expr_file, sep="\t", header=0, index_col=0)
  all_species_expr_df = pd.concat([all_species_expr_df, species_expr_df])

all_species_expr_df["max_expr"] = all_species_expr_df.apply(lambda x: max(x), axis=1)

##################################
###### MAIN ######################
##################################
#Create dictionary with key=geneID, value=max_expr across tissues
all_species_geneID_max_expr_dict = pd.Series(all_species_expr_df.max_expr.values, index=all_species_expr_df.index).to_dict()

#Filter out genes absent from GTF somehow
not_annot_genes_list = [gene for gene in list(orthogroups_df["GeneID"]) if gene not in list(all_species_geneID_max_expr_dict.keys())]
orthogroups_df = orthogroups_df.loc[~orthogroups_df["GeneID"].isin(not_annot_genes_list)]
removed_genes_list = removed_genes_list + [gene+"_annot_problem" for gene in not_annot_genes_list]

#Cycle by orthogroups
orthogroups_grouped_df = orthogroups_df.groupby("OG_ID")
final_orthogroups_df = pd.DataFrame()
lowly_expr_genes_list = [] #Initialize list of lowly expressed genes
removed_genes_list = [] #Initialize list of removed genes

for OG_ID, group in orthogroups_grouped_df:
  duplicated_species = list(set(list([species for species in list(group["Species"]) if list(group["Species"]).count(species) >= 2])))
  #for species in duplicated_species:
  for species in [element for element in duplicated_species if element in clade_species]:
    genes_to_remove = []
    geneID_max_expr_dict = {}
    species_genes = list(group.loc[group["Species"]==species]["GeneID"])
    for gene in species_genes:
      if all_species_geneID_max_expr_dict[gene] < float(min_expr):
        genes_to_remove.append(gene)
        geneID_max_expr_dict[gene] = all_species_geneID_max_expr_dict[gene]
    if len(genes_to_remove) == len(species_genes): #If none of the species genes has a maximum expression higher than the specified cutoff.
      gene_to_keep = [gene for gene in genes_to_remove if geneID_max_expr_dict[gene] == max(list(geneID_max_expr_dict.values()))][0] #keep the gene with the highest maximum expr
      lowly_expr_genes_list.append(gene_to_keep) #save it in a list of problematic genes
      genes_to_remove = [element for element in genes_to_remove if element != lowly_expr_genes_list] #Update genes to remove list
    #remove the genes from the group
    group = group.loc[~group["GeneID"].isin(genes_to_remove)]  #Remove all the genes which do not pass the cutoff
    if len(genes_to_remove) > 0:
       removed_genes_list = removed_genes_list + genes_to_remove 
  final_orthogroups_df = pd.concat([final_orthogroups_df, group])

#Write to file
final_orthogroups_df.to_csv(output_orthogroups_file, sep="\t", index=False, header=False, na_rep="NA")
pd.Series(lowly_expr_genes_list).to_csv(output_lowly_expressed_file, sep="\t", index=False, header=False, na_rep="NA")
pd.Series(removed_genes_list).to_csv(output_removed_genes_file, sep="\t", index=False, header=False, na_rep="NA")
