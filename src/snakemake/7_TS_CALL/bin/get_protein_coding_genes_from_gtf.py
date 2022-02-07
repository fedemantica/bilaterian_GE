#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to get protein_coding genes from a gtf file")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Gtf file")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Main #######
#Read input
gtf_df = pd.read_table(my_input, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
gtf_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part] #add geneID as a separate field
CDS_gtf_df = gtf_df.loc[gtf_df["type"]=="CDS"]
#Isolate protein coding genes
protein_coding_genes_df = pd.Series(list(set(list(CDS_gtf_df["geneID"]))))
#Save to file
protein_coding_genes_df.to_csv(my_output, sep="\t", header=False, index=False, na_rep="NA")
