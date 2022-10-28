#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import csv

parser = argparse.ArgumentParser(description="Script to correct the human reference GTF for those cases in which the selected isoform was not the one with the longest CDS.")
parser.add_argument("--original_gtf", "-og", required=True, metavar="original_gtf", help="Original reference gtf to be corrected (one isoform per gene)")
parser.add_argument("--longestCDS_gtf", "-lg", required=True, metavar="longestCDS_gtf", help="Reference gtf (one isoform per gene) where the isoform with the longest CDS is selected")
parser.add_argument("--entries_to_correct", "-ec", required=True, metavar="entries_to_correct", help="")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file: corrected reference gtf")

###### Read arguments
args = parser.parse_args()
original_gtf_file = args.original_gtf
longestCDS_gtf_file = args.longestCDS_gtf
entries_to_correct_file = args.entries_to_correct
output_file = args.output

#####################################
######### READ INPUTS ###############
#####################################
original_gtf_df = pd.read_table(original_gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
longestCDS_gtf_df = pd.read_table(longestCDS_gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
entries_to_correct_df = pd.read_table(entries_to_correct_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"]) #NB: this only contains human entries

######################################
########## CORRECT GTF ###############
######################################
#Add geneID column
original_gtf_df["GeneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(original_gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part]
longestCDS_gtf_df["GeneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(longestCDS_gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part]
#Isolate genes to be modified
entries_to_correct_list = list(entries_to_correct_df["GeneID"])
#Remove genes to be modified from the original gtf
subsetted_original_gtf_df = original_gtf_df.loc[~original_gtf_df["GeneID"].isin(entries_to_correct_list)]
#Isolate the entries to be modified from the longestCDS gtf and add them to the original gtf
subsetted_longestCDS_gtf_df = longestCDS_gtf_df.loc[longestCDS_gtf_df["GeneID"].isin(entries_to_correct_list)]
final_gtf_df = pd.concat([subsetted_original_gtf_df, subsetted_longestCDS_gtf_df])
final_gtf_df = final_gtf_df.drop(columns=["GeneID"])

######################################
######### SAVE TO OUTPUT #############
######################################
final_gtf_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA", quoting=csv.QUOTE_NONE)
