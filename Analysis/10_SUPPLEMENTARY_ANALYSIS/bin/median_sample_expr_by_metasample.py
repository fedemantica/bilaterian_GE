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
all_metasamples = list(set(list(metadata_df["Metasample"]))) #generate metasample list
#cycle on all metasamples
final_df = pd.DataFrame()
for metasample in all_metasamples:
  metasample_sample_list = list(set(list(metadata_df.loc[metadata_df["Metasample"]==metasample]["Sample"])))
  metasample_df = expr_df.loc[:,metasample_sample_list]
  if my_measure == "average":
    metasample_measure_series = metasample_df.apply(np.nanmean, axis=1)
  if my_measure == "median":
    metasample_measure_series = metasample_df.apply(np.nanmedian, axis=1)
  metasample_measure_df = pd.DataFrame(metasample_measure_series).rename(columns={0 : metasample})
  final_df = pd.concat([final_df, metasample_measure_df], axis=1)

#Save to output
final_df.to_csv(output_file, sep="\t", header=True, index=True, na_rep="NA")
