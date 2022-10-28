#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np
import re

parser = argparse.ArgumentParser(description="Script to compute average expression correlation between each gene all the genes (from other species) in the same gene orthogroups")
parser.add_argument("--input", "-i", required=True, metavar="input", help="File with pairwise expression correlations between all genes within an orthogroup. The format is col1=OG_ID, col2=Species1, col3=Species2, col4=GeneID1, col5=GeneID2, col6=expr_cor")
parser.add_argument("--measure", "-m", required=True, metavar="measure", help="average or median. It indicates the summarizing measure to use")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
input_file = args.input
my_measure = args.measure
output_file = args.output

##################################
###### MAIN ######################
##################################
orthogroups_df = pd.read_table(input_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species1", "Species2", "GeneID1", "GeneID2", "Expr_cor"])
#Generate column with OG_ID;GeneID1 and group df by it
orthogroups_df["OG_ID;GeneID1"] = [str(element[0])+";"+str(element[1]) for element in zip(list(orthogroups_df["OG_ID"]), list(orthogroups_df["GeneID1"]))]
grouped_df = orthogroups_df.groupby("OG_ID;GeneID1")
final_df = pd.DataFrame()
for ID, group in grouped_df:
  OG_ID = list(group["OG_ID"])[0]
  GeneID = list(group["GeneID1"])[0]
  if my_measure == "average":
    average_expr_cor = np.mean(group["Expr_cor"])
  elif my_measure == "median":
    average_expr_cor = np.median(group["Expr_cor"])
  final_df = pd.concat([final_df, pd.DataFrame({"OG_ID" : [OG_ID], "GeneID" : [GeneID], "Avg_expr_cor" : [average_expr_cor]})])

#Save to file
final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
