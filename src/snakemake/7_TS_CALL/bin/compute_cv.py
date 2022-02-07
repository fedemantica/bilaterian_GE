#!/usr/bin/env python3

import argparse
import pandas as pd
import math
from statistics import stdev
from statistics import mean

parser = argparse.ArgumentParser(description="Script to compute Tau (a tissue specificity measure ranging 0-1)")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Logged expression table with genes on rows and tissues on columns")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Define function #####
def compute_cv(x):
  valid_tissue_num = len([element for element in x if not math.isnan(element)])
  valid_tissues = [element for element in x if not math.isnan(element)]
  valid_tissues = [element if element >= 0 else 0 for element in valid_tissues]
  #if len([element for element in x if math.isnan(element)]) == len(x):
  if valid_tissue_num == 0:
    my_cv = "NA"
  elif max(valid_tissues) <= 1:
    my_cv = "NA"
  else:
    my_cv = stdev(valid_tissues)/mean(valid_tissues)
  return(my_cv)

##### Main #####
my_df = pd.read_table(my_input, sep="\t", header=0, index_col=0)
#header: Tissues
final_df = my_df.apply(compute_cv, axis=1)
final_df.to_csv(my_output, sep="\t", header=False, index=True, na_rep="NA")