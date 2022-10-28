#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to isolate paralogous genes that are annotated on the same chr and on the same strand")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Paralogous genes from the same species. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--chr_strand_info", "-c", required=True, metavar="chr_strand_info", help="Chr and strand info by gene")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##### Read arguments #####
args = parser.parse_args()
my_input = args.input
my_chr_strand_info = args.chr_strand_info
my_output = args.output

##### Main #####
#Read inputs
orthogroups_df = pd.read_table(str(my_input), sep="\t", header=None, names=["OG_ID", "Species", "GeneID", "GeneName"])
chr_strand_info_df = pd.read_table(str(my_chr_strand_info), sep="\t", header=None, names=["GeneID", "Chr", "Strand"])
#create dictionaries for chr and strand
my_chr_dict = pd.Series(chr_strand_info_df.Chr.values, index=chr_strand_info_df.GeneID).to_dict()
my_strand_dict = pd.Series(chr_strand_info_df.Strand.values, index=chr_strand_info_df.GeneID).to_dict()
#group by orthogroups ID
final_df = pd.DataFrame(columns=list(orthogroups_df.columns.values))
grouped_orthogroups_df = orthogroups_df.groupby("OG_ID")
for name, group in grouped_orthogroups_df:
  genes = list(group["GeneID"])
  genes_to_save = []
  strands = [my_strand_dict[element] for element in genes]
  chrs = [my_chr_dict[element] for element in genes]
  #Maybe not the most logical way. But for each gene, I check if there is at least another gene annotated on the same chr/strand. If yes, I save the gene to output.
  for gene in genes:
    if len([element for element in strands if element==my_strand_dict[gene]]) >= 2:
      if len([element for element in chrs if element==my_chr_dict[gene]]) >= 2:
        genes_to_save.append(gene)
  if len(genes_to_save) >= 2:
    new_group = group.loc[group.GeneID.isin(genes_to_save)]
    final_df = pd.concat([final_df, new_group])

#Save to file
final_df.to_csv(str(my_output), sep="\t", header=False, index=False)
