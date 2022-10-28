#!/usr/bin/env python3

import argparse
import pandas as pd
import math
import numpy as np

parser = argparse.ArgumentParser(description="Script to compute Tau (a tissue specificity measure ranging 0-1)")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Logged expression table with genes on rows and tissues on columns")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Define function #####
def compute_tau(x):
  valid_tissue_num = len([element for element in x if not math.isnan(element)])
  valid_tissues = [element for element in x if not math.isnan(element)]
  #valid_tissues = [element if element >= 0 else 0 for element in valid_tissues] #This line was needed only for the SVA
  if valid_tissue_num == 0:
    my_tau = np.nan
  elif max(valid_tissues) <= math.log2(5): #require the maximum to be at least 5 TPM
    my_tau = np.nan
  else:
    my_tau = sum([1-element/max(valid_tissues) for element in valid_tissues])/(valid_tissue_num-1)
  return(my_tau)

##### Main #####
my_df = pd.read_table(my_input, sep="\t", header=0, index_col=0)
#header: Tissues
final_df = my_df.apply(compute_tau, axis=1)
#Drop lines with NA
#final_df.dropna(inplace=True)
final_df.to_csv(my_output, sep="\t", header=False, index=True, na_rep="NA")
