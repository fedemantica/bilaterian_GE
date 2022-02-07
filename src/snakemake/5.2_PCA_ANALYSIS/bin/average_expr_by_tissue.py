#!/usr/bin/env python3

import argparse
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description="Script to compute the average or median expression across samples from the same tissue")
parser.add_argument("--expr", "-e", required=True, metavar="expression_table", help="Expression tables with genes on rows and metasamples on columns")
parser.add_argument("--metadata", "-md", required=True, metavar="metadata", help="Metadata table for the species of interest")
parser.add_argument("--measure", "-ms", required=True, metavar="metadata", help="average or median. It indicates the summarizing measure to use")
parser.add_argument("--output", "-o", required=True, metavar="output_file", help="Two columns file contaning with col1=ProteinID, col2=ProteinLength")

#Read arguments
args = parser.parse_args()
expr_file = args.expr
metadata_file = args.metadata
my_measure = args.measure
output_file = args.output

#######################################
######### READ INPUTS #################
#######################################
expr_df = pd.read_table(expr_file, sep="\t", index_col=0, header=0)
metadata_df = pd.read_table(metadata_file, sep="\t", index_col=False, header=0)


#######################################
######### MAIN ########################
#######################################
all_tissues = list(set(list(metadata_df["Tissue"]))) #generate tissue list
#cycle on all tissues
final_df = pd.DataFrame()
for tissue in all_tissues:
  tissue_sample_list = list(set(list(metadata_df.loc[metadata_df["Tissue"]==tissue]["Metasample"])))
  tissue_df = expr_df.loc[:,tissue_sample_list]
  if my_measure == "average":
    tissue_measure_series = tissue_df.apply(np.nanmean, axis=1)
  if my_measure == "median":
    tissue_measure_series = tissue_df.apply(np.nanmedian, axis=1)
  tissue_measure_df = pd.DataFrame(tissue_measure_series).rename(columns={0 : tissue})
  final_df = pd.concat([final_df, tissue_measure_df], axis=1)

#Save to output
final_df.to_csv(output_file, sep="\t", header=True, index=True, na_rep="NA")
