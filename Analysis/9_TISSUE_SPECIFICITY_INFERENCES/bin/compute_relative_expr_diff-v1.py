#!/usr/bin/env python3

import argparse
import pandas as pd
import math
import numpy as np
from collections import defaultdict

parser = argparse.ArgumentParser(description="Script to get the relative expression difference between the expression of the top two tissues and the third one")
parser.add_argument("--input", "-i", required=True, nargs="+",  metavar="input", help="List of files (one per species) with scaled gene expression across tissues")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##############################
##### Read arguments #########
##############################
args = parser.parse_args()
my_input_list = args.input
my_output = args.output

############################
##### Define function ######
############################
def get_sorted_tissues_list(x):
 highest_to_lowest_tissues = list(reversed([value for _, value in sorted(zip(x, tissues))]))
 #Add last common tissue
 highest_to_lowest_tissues = highest_to_lowest_tissues + ["None"] 
 return(highest_to_lowest_tissues)

###########################
##### Main ################
###########################
#Initialize final dataframe
super_final_df = pd.DataFrame()

#Cycle on all files containing the tissue scaled expression
for my_input in my_input_list:
  my_df = pd.read_table(my_input, sep="\t", header=0, index_col=0)
  tissues = list(my_df.columns.values)
  
  #Get dictionary with key=geneID, value=ordered_tissues
  geneID_ordered_tissues_dict = my_df.apply(get_sorted_tissues_list, axis=1).to_dict()
  #Get nested dictionary: dict["geneID"]["Tissue"] = Expr. This is a very elegant command
  #Add expression to None
  my_df["None"] = my_df.min(axis=1)
  geneID_tissue_expr_dict = my_df.transpose().to_dict()
  
  #Pivot to long df
  my_df["GeneID"] = list(my_df.index.values)
  my_long_df = pd.melt(my_df, id_vars="GeneID", var_name="Tissue", value_name="Rel_Expr")
  my_long_df = my_long_df.loc[my_long_df["Tissue"]!="None"]
  #Add a column with the difference in expr from that tissue to the following
  my_long_df["Expr_diff"]= [geneID_tissue_expr_dict[element[0]][element[1]] - geneID_tissue_expr_dict[element[0]][geneID_ordered_tissues_dict[element[0]][geneID_ordered_tissues_dict[element[0]].index(element[1])+1]] for element in list(zip(list(my_long_df["GeneID"]), list(my_long_df["Tissue"])))]
  #Subset to only the two top tissues for each genes
  my_long_grouped_df = my_long_df.groupby("GeneID")
  final_df = pd.DataFrame()
  for name, group in my_long_grouped_df:
    final_group = group.loc[group["Tissue"].isin(geneID_ordered_tissues_dict[name][0:2])]
    final_df = pd.concat([final_df, final_group])
  #Get the expr diff from two following tissues
  final_df["Second_expr_diff"] = [geneID_tissue_expr_dict[element[0]][element[1]] - geneID_tissue_expr_dict[element[0]][geneID_ordered_tissues_dict[element[0]][geneID_ordered_tissues_dict[element[0]].index(element[1])+2]] for element in list(zip(list(final_df["GeneID"]), list(final_df["Tissue"])))]
  
  #For the top tissue, select the second value. For the second top tissue, select the third value.
  grouped_final_df = final_df.groupby("GeneID")
  for name, group in grouped_final_df:
    group["Final_expr_diff"] = [element[2] if geneID_ordered_tissues_dict[name][0]==element[0] else element[1] for element in zip(list(group["Tissue"]), list(group["Expr_diff"]), list(group["Second_expr_diff"]))]
    super_final_df = pd.concat([super_final_df, group])

##### Remove the useless columns from the final dataframe
super_final_df = super_final_df[["GeneID", "Tissue", "Final_Expr_diff"]]

###########################
##### Save ################
###########################
super_final_df.to_csv(my_output, sep="\t", index=False, header=False, na_rep="NA")
