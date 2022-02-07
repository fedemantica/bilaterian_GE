#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to get gmt (Gene Matrix Transposed) format starting from the output of a GO transfer (filtered for the species of interest)")
parser.add_argument("--trasferred_GO", "-GO", required=True, metavar="retina_exons", help="Output of a GO transfer, with col1=Gene_ID, col2=GO_ID, col3=description. No Header.")
parser.add_argument("--output", "-o", required=True, metavar="retina_psi", help="Output file of get_tisAS_v6.pl: dPSIs_Retina*")

###### Read arguments ############
args = parser.parse_args()
trasferred_GO_file = args.trasferred_GO
output_file = args.output

transferred_GO_df = pd.read_table(trasferred_GO_file, sep="\t", index_col=False, header=None, names=["Gene_ID", "GO_ID", "description"])
GO_ID_description_df = transferred_GO_df[["GO_ID", "description"]].drop_duplicates()
GO_ID_description_dict = pd.Series(GO_ID_description_df.description.values, index=GO_ID_description_df.GO_ID).to_dict()
GO_ID_Gene_ID_df = transferred_GO_df[["GO_ID", "Gene_ID"]]
GO_ID_Gene_ID_grouped_df = GO_ID_Gene_ID_df.groupby("GO_ID")

for GO, group in GO_ID_Gene_ID_grouped_df:
  description = GO_ID_description_dict[GO]
  GO_df = pd.concat([pd.Series([GO, description]), pd.Series(group.transpose().iloc[1,:])]).to_frame().T
  GO_df.to_csv(output_file, sep="\t", index=False, header=False, mode="a", na_rep="NA")
