#!/usr/bin/env python3

import argparse
import pandas as pd
import math
import numpy as np

parser = argparse.ArgumentParser(description="Script to get the tissue associated with the tissue specificity")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Scaled expression across TS genes (Tau cutoff)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### In this case: if the difference between the first and the second >= 0.10 and FC >= 1.7: select the first tissue only
##### If not: check the difference between the second and the third. If this difference is >= 0.15, select the first two tissues
##### If not: select no tissue at all

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Define function #####
def get_tissue_with_ts(x):
  tissue_expr_dict = {x[index] : tissues[index] for index in list(range(0,len(x)))}
  #max_expr = max(x) #get maximum expression
  x = [element for element in x if ~np.isnan(element)]
  if (len(x)==0):
    final_tissue="NA"
  else:
    max_expr = sorted(x, reverse=True)[0]
    second_max_expr = sorted(x, reverse=True)[1]
    third_max_expr = sorted(x, reverse=True)[2]
    #Compute the differences
    max_expr_diff = max_expr - second_max_expr
    second_max_expr_diff = second_max_expr - third_max_expr
    #Check the cutoffs
    if second_max_expr == 0: #If the expression of the second highest tissue is zero
      final_tissue = tissue_expr_dict[max_expr]
      final_tissue_paralog = final_tissue
    else:
      if max_expr_diff >= 0.10 and max_expr/second_max_expr >= 1.7:
        final_tissue = tissue_expr_dict[max_expr]
        final_tissue_paralog = final_tissue
      else:
        if second_max_expr_diff >= 0.15:
          if max_expr != second_max_expr:
            ordered_tissues = [tissue_expr_dict[max_expr], tissue_expr_dict[second_max_expr]]
            ordered_tissues.sort()
            final_tissue = ordered_tissues[0]+";"+ordered_tissues[1]
            final_tissue_paralog = final_tissue
          else: #If the two top tissues have exactly the same values
            ordered_tissues = [tissues[index] for index in list(range(0,len(x))) if x[index]==max_expr]
            ordered_tissues.sort()
            final_tissue = ordered_tissues[0]+";"+ordered_tissues[1]
            final_tissue_paralog = final_tissue
        else:
          #select anyways the top two tissue
          ordered_tissues = [tissue_expr_dict[max_expr], tissue_expr_dict[second_max_expr]]
          ordered_tissues.sort()
          final_tissue = "None"
          final_tissue_paralog = ordered_tissues[0]+";"+ordered_tissues[1]  
  #else:
  #  final_tissue="NA"
  return([final_tissue, final_tissue_paralog])

##### Main #####
my_df = pd.read_table(my_input, sep="\t", header=0, index_col=0)
tissues = list(my_df.columns.values)
#header: Tissues
final_df = my_df.apply(get_tissue_with_ts, axis=1)
final_df = pd.DataFrame(final_df, index=final_df.index.values)
final_df = pd.DataFrame(final_df[0].to_list(), index=final_df.index.values)
final_df.to_csv(my_output, sep="\t", header=False, index=True, na_rep="NA")
