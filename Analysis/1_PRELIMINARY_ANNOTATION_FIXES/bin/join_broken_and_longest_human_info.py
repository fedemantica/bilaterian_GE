#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to join the info from the reciprocal comparison of broken genes pairs/groups and the comparison with their respective longest human ortholog.")
parser.add_argument("--broken_info", "-b", required=True, metavar="broken_info", help="File containing info from the reciprocal comparisons of broken genes pairs/groups.")
parser.add_argument("--human_info", "-hi", required=True, metavar="human_info", help="File containing info from the comparisons of each broken pair/group with the respective longest human ortholog.")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All combined information")

##### Read arguments #####
args = parser.parse_args()
my_broken_info = args.broken_info
my_human_info = args.human_info
my_output = args.output

##### Main #####
#Read inputs
broken_genes_df = pd.read_table(str(my_broken_info), sep="\t", index_col=False, header=0)
longest_human_df = pd.read_table(str(my_human_info), sep="\t", index_col=False, header=0)
#Create a dictionary with key=geneID1;geneID2, value=sim_score1;sim_score2
broken_genes_df["Gene_pairs"] = broken_genes_df["GeneID1"]+";"+broken_genes_df["GeneID2"]
broken_genes_df["Sim_scores"] = pd.Series([str(element) for element in list(broken_genes_df["GeneID1_sim_score"])])+";"+pd.Series([str(element) for element in list(broken_genes_df["GeneID2_sim_score"])])
gene_pair_sim_score_dict = pd.Series(broken_genes_df.Sim_scores.values, index=broken_genes_df.Gene_pairs).to_dict()
#Add columns with gene pair and sim scores to human longest df
longest_human_df["Gene_pairs"] = longest_human_df["GeneID1"]+";"+longest_human_df["GeneID2"]
longest_human_df["Sim_scores"] = longest_human_df["Gene_pairs"].map(gene_pair_sim_score_dict)
#separate sim scores
longest_human_df["GeneID1_sim_score"] = [element.split(";")[0] for element in list(longest_human_df["Sim_scores"])]
longest_human_df["GeneID2_sim_score"] = [element.split(";")[1] for element in list(longest_human_df["Sim_scores"])]
#Select columns for final dataframe
final_df = longest_human_df[["OG_ID", "Gene_pairs", "Hsa_geneID", "GeneID1_sim_score", "GeneID2_sim_score", "Match_overlap_proportion"]]
#Save df to file
final_df.to_csv(str(my_output), sep="\t", index=False, header=True, na_rep="NA")
