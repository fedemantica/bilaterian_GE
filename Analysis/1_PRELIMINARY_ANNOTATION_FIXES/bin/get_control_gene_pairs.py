#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to select pairs of controls genes to identify cutoffs able to distinguish between true and false potentially broken genes.")
parser.add_argument("--broken_genes", "-b", required=True, metavar="broken_genes", help="Paralogous genes sequentially annotated on the same chr and strand. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Complete gene orthogroups. Format is col1=OG_ID, col2=Species, col3=GeneID")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species identifier as reported in broken_genes and orthogroups inputs")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All pairs of control genes.")

##### Read arguments #####
args = parser.parse_args()
my_broken_genes = args.broken_genes
my_orthogroups = args.orthogroups
my_species = args.species
my_output = args.output

##### Main #####
#Read inputs
orthogroups_df = pd.read_table(str(my_orthogroups), sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
PBG_df = pd.read_table(str(my_broken_genes), sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
PBG_orthogroups = list(set(list(PBG_df["OG_ID"]))) #select the PBG orthogroups
#Subset orthogroups_df only for species of interest
species_orthogroups_df = orthogroups_df.loc[orthogroups_df.Species==my_species]
#Select the orthogroups where there are exactly 2 genes per species which are NOT potentially broken genes.
all_species_orthogroups = list(species_orthogroups_df["OG_ID"])
orthogroups_with_pairs = list(set([element for element in all_species_orthogroups if all_species_orthogroups.count(element) == 2]))
no_PBG_orthogroups = [element for element in orthogroups_with_pairs if element not in PBG_orthogroups]
final_species_orthogroups = species_orthogroups_df.loc[species_orthogroups_df.OG_ID.isin(no_PBG_orthogroups)]
#Put the two genes one next to the other
final_df = pd.DataFrame(columns=["OG_ID", "GeneID_1", "GeneID_2"])
my_grouped_df = final_species_orthogroups.groupby("OG_ID")
for name, group in my_grouped_df:
  my_geneID_1 = list(group["GeneID"])[0]
  my_geneID_2 = list(group["GeneID"])[1]
  group_df = pd.DataFrame({"OG_ID" : [name], "GeneID_1" : [my_geneID_1], "GeneID_2" : [my_geneID_2]})
  final_df = pd.concat([final_df, group_df])
#Save to file
final_df.to_csv(str(my_output), sep="\t", header=False, index=False, na_rep="NA")
