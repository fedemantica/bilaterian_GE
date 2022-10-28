#!/usr/bin/env python3

import argparse
import pandas as pd
import math
import numpy as np

parser = argparse.ArgumentParser(description="Script to compute the relative tissue specificity")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Logged expression table with genes on rows and tissues on columns")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Define function #####
def compute_relative_ts(x):
  valid_tissue_num = len([element for element in x if not math.isnan(element)])
  #valid_tissues = [element for element in x if not math.isnan(element)]
  valid_tissues = [0 if element <= 0 else element for element in x]
  #if len([element for element in x if math.isnan(element)]) == len(x):
  if valid_tissue_num == 0:
    final_values = pd.Series([np.nan for element in x])
  #elif max(valid_tissues) <= math.log2(2): #require the maximum to be at least 1 TPM (to which we added 1 when computing the log2)
  #  final_values = pd.Series([np.nan for element in x])
  else:
    total_expr = np.nansum(valid_tissues)
    final_values = pd.Series([element/total_expr for element in valid_tissues])
  return(final_values)

##### Main #####
my_df = pd.read_table(my_input, sep="\t", header=0, index_col=0)
colnames = list(my_df.columns.values)
new_colnames_dict = {element : colnames[element] for element in list(range(0,len(colnames)))}
#header: Tissues
final_df = my_df.apply(compute_relative_ts, axis=1)
final_df = final_df.rename(columns=new_colnames_dict)
final_df.to_csv(my_output, sep="\t", header=True, index=True, na_rep="NA")
