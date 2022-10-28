#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to retrieve start-stop genomic coordinates by gene")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Reference GTF (i.e. just one representative isoform for each gene)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_output = args.output

##### Main #####
#Read gtf in and add geneID column
my_gtf_df = pd.read_table(my_input, sep="\t", names=["chr", "source", "feature", "start", "stop", "score", "strand", "phase", "attributes"])
intermediate_gene_list = [re.sub(";.*", "", element) for element in list(my_gtf_df["attributes"])]
my_gtf_df["geneID"] = pd.Series([re.sub(" ", "", re.sub("gene_id ", "", re.sub('"', '', element))) for element in intermediate_gene_list])
#Filter only for the exons entry and subset columns
my_filtered_gtf = my_gtf_df.loc[my_gtf_df.feature=="exon"]
my_filtered_gtf = my_gtf_df.loc[:,["geneID","chr","start","stop"]]
#group by geneID and select lower start and higher stop
#I am just looking at DNA, so here is actually the same independently from the strand [start is always lower than stop]
final_df = pd.DataFrame(columns=list(my_filtered_gtf.columns.values))
my_grouped_df = my_filtered_gtf.groupby("geneID")
for name, group in my_grouped_df:
  my_chr = list(group["chr"])[0]
  my_start = min(list(group["start"]))
  my_stop = max(list(group["stop"]))
  group_df = pd.DataFrame({"geneID" : [name], "chr" : [my_chr], "start" : [my_start], "stop" : [my_stop]})
  final_df = pd.concat([final_df, group_df])
#save to output file
final_df.to_csv(my_output, sep="\t", index=False, header=False)
