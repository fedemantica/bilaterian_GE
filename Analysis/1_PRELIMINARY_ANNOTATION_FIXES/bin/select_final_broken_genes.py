#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to select the final broken genes per species.")
parser.add_argument("--input", "-i", required=True, metavar="input", help="File containing all combined info of broken gene comparisons (reciprocal and vs the respective longest human ortholog).")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All combined information")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Main #####
more_broken_df = pd.read_table(str(my_input), sep="\t", index_col=False, header=None, names=["OG_ID", "Gene_pairs", "Hsa_gene", "GeneID1_sim_score", "GeneID2_sim_score", "Match_overlap_proportion"])
#Join genes on the "more" side
more_broken_df["GeneID1"] = [element.split(";")[0] for element in list(more_broken_df["Gene_pairs"])]
more_broken_df["GeneID2"] = [element.split(";")[1] for element in list(more_broken_df["Gene_pairs"])]
more_broken_grouped_df = more_broken_df.groupby("OG_ID")
final_df = pd.DataFrame()
for name, group in more_broken_grouped_df:
  if group.shape[0] >= 2:
    group = group[["GeneID1", "GeneID2"]]
    geneID1_list = list(group["GeneID1"])
    geneID2_list = list(group["GeneID2"])
    common_genes = [gene for gene in geneID1_list if gene in geneID2_list]
    if len(common_genes) >= 1:
      for gene in common_genes:
        #Join the two lines
        subset_group_all = pd.concat([group[group["GeneID1"]==gene].reset_index(drop=True), group[group["GeneID2"]==gene].reset_index(drop=True)], axis=1)
        all_genes = list(set(subset_group_all.values.tolist()[0]))
        final_group = pd.Series([all_genes[0]+";"+all_genes[1]+";"+all_genes[2]]).rename("Gene_pairs")
    else:
      group["Gene_pairs"] = group["GeneID1"]+";"+group["GeneID2"]
      final_group = group["Gene_pairs"]
  else:
      group["Gene_pairs"] = group["GeneID1"]+";"+group["GeneID2"]
      final_group = group["Gene_pairs"]
  final_df = pd.concat([final_df, final_group])

#Save to final file
final_df.to_csv(str(my_output), sep="\t", index=False, header=False, na_rep="NA")
