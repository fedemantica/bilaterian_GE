#!/usr/bin/env python3

import argparse
import pandas as pd
import math
import numpy as np

parser = argparse.ArgumentParser(description="Script to get the tissue associated with the tissue specificity")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Scaled expression across TS genes (Tau cutoff)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Define function #####
def get_tissue_with_ts(x):
  tissue_expr_dict = {x[index] : tissues[index] for index in list(range(0,len(x)))}
  max_expr = max(x) #get maximum expression
  second_max_expr = max([element for element in x if element != max_expr]) #get second maximum expression
  max_expr_diff = max_expr - second_max_expr #compute the difference
  if np.isnan(max_expr_diff) == False:
    if max_expr_diff > 0.1:
      final_tissue = tissue_expr_dict[max_expr]
    elif max_expr_diff <= 0.1:
      ordered_tissues = [tissue_expr_dict[max_expr], tissue_expr_dict[second_max_expr]]
      ordered_tissues.sort()
      final_tissue = ordered_tissues[0]+";"+ordered_tissues[1]
  else:
    final_tissue="NA"
  return(final_tissue)

##### Main #####
my_df = pd.read_table(my_input, sep="\t", header=0, index_col=0)
tissues = list(my_df.columns.values)
#header: Tissues
final_df = my_df.apply(get_tissue_with_ts, axis=1)
final_df.to_csv(my_output, sep="\t", header=False, index=True, na_rep="NA")
