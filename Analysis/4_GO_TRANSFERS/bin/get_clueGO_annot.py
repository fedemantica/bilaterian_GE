#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to get clueGO annotation in a format useful to run the GO transfer and containing Ensembl geneIDs")
parser.add_argument("--input_annot", "-i", required=True, metavar="input_annot", help="File with clueGO annotation")
parser.add_argument("--key", "-k", required=True, metavar="key", help="File with correspondences between Entrez and Ensembl IDs, also provided by clueGO")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Output file with Ensembl IDs ready to be used for the GO transfers")
parser.add_argument("--output_levels", "-ol", required=True, metavar="output_levels", help="Output file with Ensembl IDs ready to be used for the GO transfers; GO level is included")

###### Read arguments ############
args = parser.parse_args()
input_annot_file = args.input_annot
key_file = args.key
output_file = args.output
output_levels_file = args.output_levels

###### Read input files ##########
annot_df = pd.read_table(str(input_annot_file), sep="\t", index_col=False, header=0, names=["GOterm", "Levels", "GOdef", "Genes"])
key_df = pd.read_table(str(key_file), sep="\t", index_col=False, header=0) #first two cols: UniqueID#EntrezGeneID, EnsemblGeneID

###### main #####################
#get dictionary EntrezID:EnsemblID.
#NB: EntrezID are strings.
entrez_ensembl_dict = pd.Series(key_df.EnsemblGeneID.values, index=key_df["UniqueID#EntrezGeneID"]).to_dict()

#define functions
def get_gene_list(my_str):
  final_list = [re.sub(".*:", "", element) for element in list(my_str.split("|"))]
  return(final_list)

def get_level_list(my_str):
  final_list = [element for element in list(my_str.split(","))]
  return(final_list)

#Just explode the dataframe
annot_df["Genes"] = [get_gene_list(str(element)) for element in list(annot_df["Genes"])]  
exploded_df = annot_df.explode("Genes")
exploded_df["EnsemblID"] = exploded_df["Genes"].map(entrez_ensembl_dict)
final_df = exploded_df[["EnsemblID", "GOterm", "GOdef"]]
final_df = final_df.dropna(subset=["EnsemblID"]) #Remove rows with NAs in the EnsemblID.
final_df["EnsemblID"] = [re.sub("\\|.*", "", str(element)) for element in list(final_df["EnsemblID"])] #In case of EnsemblIDs divided by pipe, just take the first one. Maybe not ideal, but it can be a start.

#Just explode the dataframes but keeping the levels
exploded_df["Levels"] = [get_level_list(str(element)) for element in list(exploded_df["Levels"])]
exploded_levels_df = exploded_df.explode("Levels")
exploded_levels_df["EnsemblID"] = exploded_levels_df["Genes"].map(entrez_ensembl_dict)
levels_final_df = exploded_levels_df[["EnsemblID", "GOterm", "GOdef", "Levels"]]
levels_final_df = levels_final_df.dropna(subset=["EnsemblID"])
levels_final_df["EnsemblID"] = [re.sub("\\|.*", "", str(element)) for element in list(levels_final_df["EnsemblID"])]

#save to output file
final_df.to_csv(str(output_file), sep="\t", index=False, header=False, na_rep = "NA")
levels_final_df.to_csv(str(output_levels_file), sep="\t", index=False, header=False, na_rep = "NA")
