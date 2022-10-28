#!/usr/bin/env python3

import argparse
import pandas as pd
import re
from collections import Counter
import numpy as np
import math 

parser = argparse.ArgumentParser(description="Script to select the best hits for each species among all the genes included in a gene orthogroup")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Gene orthogroups (in pricinple corrected and reclustered), with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--seq_sim", "-s", required=True, metavar="seq_sim", help="Files containing pairwise sequence similarities between genes in the same gene orthogroup")
parser.add_argument("--expr_cor", "-e", required=True, metavar="expr_cor", help="Files containing the pairwise pearson expression correlations between genes in the same gene orthogroup")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output directory")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
seq_sim_file = args.seq_sim
expr_cor = args.expr_cor
output_file = args.output

##################################
###### READ INPUTS ###############
##################################
#Orthogroups
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#Generate geneID-species dict
geneID_species_dict = pd.Series(orthogroups_df.Species.values, index=orthogroups_df.GeneID).to_dict()
#seq sim dataframe
seq_sim_df = pd.read_table(seq_sim_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species1", "Species2", "GeneID1", "GeneID2", "Seq_sim"])
#expr cor dataframe
expr_cor_df = pd.read_table(expr_cor, sep="\t", index_col=False, header=None, names=["OG_ID", "Species1", "Species2", "GeneID1", "GeneID2", "Expr_cor"])

##################################
###### MAIN ######################
##################################

group_df = orthogroups_df.groupby("OG_ID")
final_df = pd.DataFrame()
for OG_ID, group in group_df:
  print(OG_ID)
  general_OG_ID = re.sub("\\..*", "", OG_ID)
  final_group = group #Create variable for the end of the script
  group_genes = list(group["GeneID"])
  #Check if there are any species with more than 2 genes
  all_species = list(group["Species"])
  species_count_dict = dict(Counter(all_species))
  duplicated_species = [species for species in list(species_count_dict.keys()) if species_count_dict[species]>=2]

  if len(duplicated_species) > 0: #If there are any duplicated species
    for species in duplicated_species:
      species_genes = list(group.loc[group["Species"]==species]["GeneID"])
    
      #This is in the case we do not have the average seq similarities and average expression correlations at the level of reclustered orthogroups
      #There are only pairwise correspondences between genes in the same orthogroup
      subset_seq_sim_df = seq_sim_df.loc[(seq_sim_df["OG_ID"]==general_OG_ID) & (seq_sim_df["GeneID1"].isin(species_genes)) & (seq_sim_df["GeneID2"].isin(group_genes))] #get all pairwise sequence similarities for all species genes 
      subset_seq_sim_grouped_df = subset_seq_sim_df.groupby("GeneID1")
      subset_seq_sim_dict = {}
      for gene, sub_group in subset_seq_sim_grouped_df:
        avg_seq_sim = np.nanmean(sub_group["Seq_sim"])
        subset_seq_sim_dict[gene] = avg_seq_sim
      #Get the maximum seq similarities among the genes from the same species present in the gene orthogroup
      all_seq_sims = list(subset_seq_sim_dict.values())

      ######## Very ugly temporal piece of code, for one problematic gene orthogroups ########
      #if general_OG_ID=="GF_005484":
      #  best_hit_gene = species_genes[0]
      #  continue
      ########################################################################################

      max_seq_sim = max(all_seq_sims)
      #If there are more genes, get them all
      max_seq_sim_genes = [gene for gene in species_genes if subset_seq_sim_dict[gene] == max_seq_sim]
      #If there is only one genes with maximum sequence similarity, and that sequence similarity is at least 0.2 higher than the second highest sequence similarity, select the relative gene as best hit.
      all_seq_sims.remove(max_seq_sim)
      counter = 1
      best_hit_candidates = max_seq_sim_genes
      while counter == 1:
        second_highest_seq_sim = max(all_seq_sims)
        if max_seq_sim - second_highest_seq_sim >= 0.2:
          counter = 2
        else:
          best_hit_candidates = best_hit_candidates + [gene for gene in species_genes if subset_seq_sim_dict[gene] == second_highest_seq_sim]
          all_seq_sims.remove(second_highest_seq_sim)
          if len(all_seq_sims) == 0:
            counter = 2
      #if there is only one best hit candidate, select that as best hit:
      if len(best_hit_candidates) == 1:
        best_hit_gene = best_hit_candidates[0]
      else:
        #here we need to check the average expresssion correlation.
        subset_expr_cor_df = expr_cor_df.loc[(expr_cor_df["OG_ID"]==general_OG_ID) & (expr_cor_df["GeneID1"].isin(best_hit_candidates)) & (expr_cor_df["GeneID2"].isin(group_genes))]
        subset_expr_cor_grouped_df = subset_expr_cor_df.groupby("GeneID1")
        subset_expr_cor_dict = {}
        for gene, sub_group in subset_expr_cor_grouped_df:
          avg_expr_cor = np.nanmean(sub_group["Expr_cor"])
          subset_expr_cor_dict[gene] = avg_expr_cor
        all_expr_cors = list(subset_expr_cor_dict.values())
        max_expr_cor = max(all_expr_cors)
        print(max_expr_cor)
        if math.isnan(max_expr_cor):
          best_hit_gene = list(subset_expr_cor_dict.keys())[0] #This is temporal. Take a random gene if the expression correlation is nan.
        else:
          best_hit_gene = [gene for gene in list(subset_expr_cor_dict.keys()) if subset_expr_cor_dict[gene] == max_expr_cor][0] #in case there is exactly the same expression correlation, just select a random one.
    #join to the final dataframe
      final_group = final_group.loc[final_group["Species"] != species] #remove all genes from the examined species.
      final_group = pd.concat([final_group, group.loc[group["GeneID"]==best_hit_gene]]) #re-add the entry for the best-hit gene.
  
  final_df = pd.concat([final_df, final_group])

final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
