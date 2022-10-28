#!/usr/bin/env python3

import argparse
import pandas as pd
import math
from scipy.stats import spearmanr
import itertools
import numpy as np

parser = argparse.ArgumentParser(description="Script to compute expression correlation between all pairs of genes from different species in the same gene orthogroup")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--expr_table_dir", "-e", required=True, metavar="expr_table_dir", help="Path to directory containing the tissue median/average expression tables by species")
parser.add_argument("--expr_file_suffix", "-ef", required=True, metavar="expr_file_suffix", help="Suffix for expression files to be uploaded")
parser.add_argument("--clade_species", "-c", required=True, metavar="clade_species", help="Comma separated list of all the species represented in the orthogroups. They will be used to read in the expression tables")
parser.add_argument("--tissues", "-t", required=True, metavar="tissues", help="Comma separated list of all tissues")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
expr_table_dir = args.expr_table_dir
expr_file_suffix = args.expr_file_suffix
clade_species_string = args.clade_species
tissue_string = args.tissues
output_file = args.output

##################################
###### READ INPUTS ###############
##################################
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
clade_species_list = clade_species_string.split(",")
all_tissues = tissue_string.split(",")
print(all_tissues)
#Generate geneID-species dict
geneID_species_dict = pd.Series(orthogroups_df.Species.values, index=orthogroups_df.GeneID).to_dict()

#Read the expression tables in. Join them in a unique table, adding columns of NAs when one tissue is missing.
joint_expr_table = pd.DataFrame()
for species in clade_species_list:
  species_expr_table = pd.read_table(expr_table_dir+"/" + species + expr_file_suffix, sep="\t", header=0, index_col=0)
  if any(tissue not in list(species_expr_table.columns.values) for tissue in all_tissues): #check if there are any missing tissues from the table
    missing_tissues = [tissue for tissue in all_tissues if tissue not in list(species_expr_table.columns.values)]
    for tissue in missing_tissues:
      species_expr_table[tissue] = np.nan
  #reorder entries
  species_expr_table = species_expr_table[all_tissues]
  #join to final table
  joint_expr_table = pd.concat([joint_expr_table, species_expr_table])

##################################
###### MAIN ######################
##################################
final_df = pd.DataFrame()

#group dataframe and cycle on groups
grouped_df = orthogroups_df.groupby("OG_ID")
for OG_ID, group in grouped_df:
  #Generate all possible pair of genes between genes belonging
  all_genes = list(group["GeneID"])
  all_gene_pairs = list(itertools.combinations(all_genes, 2))
  #Filter out pairs where genes come from the same species
  filt_gene_pairs = [gene for gene in all_gene_pairs if geneID_species_dict[gene[0]] != geneID_species_dict[gene[1]]]
  filt_gene_pairs = [[gene[0], gene[1]] for gene in filt_gene_pairs]

  for pair in filt_gene_pairs:
    geneID1 = pair[0]
    geneID2 = pair[1]
    species1 = geneID_species_dict[geneID1]
    species2 = geneID_species_dict[geneID2]
    #check if the gene expression has been quantified
    if geneID1 not in list(joint_expr_table.index.values) or geneID2 not in list(joint_expr_table.index.values):
      expr_cor = np.nan
    else:
      geneID1_expr = list(joint_expr_table.loc[geneID1,:])
      geneID2_expr = list(joint_expr_table.loc[geneID2,:])
      if any([math.isnan(element) for element in geneID1_expr+geneID2_expr]):
        #identify the nan tissue and remove it
        problematic_indexes = [index for index in list(range(0, len(geneID1_expr))) if math.isnan(geneID1_expr[index])] + [index for index in list(range(0, len(geneID2_expr))) if math.isnan(geneID2_expr[index])]
        geneID1_expr = [geneID1_expr[position] for position in list(range(0, len(geneID1_expr))) if position not in problematic_indexes]
        geneID2_expr = [geneID2_expr[position] for position in list(range(0, len(geneID2_expr))) if position not in problematic_indexes]     
      #Compute correlation
      expr_cor = spearmanr(geneID1_expr, geneID2_expr)[0]
      #Join to final_df
    final_df = pd.concat([final_df, pd.DataFrame({"OG_ID" : [OG_ID, OG_ID], "Species1" : [species1, species2], "Species2" : [species2, species1], "GeneID1" : [geneID1, geneID2], "GeneID2" : [geneID2, geneID1], "Pearons_cor" : [expr_cor, expr_cor]})])

#Save to output file
final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
