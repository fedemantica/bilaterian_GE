#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to combine the sequence similarity and overlap information for potentially broken genes.")
parser.add_argument("--broken_pairs", "-b", required=True, metavar="broken_pairs", help="Pairs of broken genes for the species of interest. Format is col1=OG_ID, col2=GeneID1, col2=GeneID2")
parser.add_argument("--sim_scores", "-s", required=True, metavar="sim_scores", help="File with similarity scores information. Format is col1=GeneID1, col2=GeneID2, col=sim_score")
parser.add_argument("--overlap", "-ov", required=True, metavar="overlap", help="File with overlap information. Format is col1=GeneID1, col2=GeneID2, col=sim_score")
parser.add_argument("--species", "-sp", required=True, metavar="species", help="Species identifier as reported in broken_genes and orthogroups inputs")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All pairs of control genes.")

##### Read arguments #####
args = parser.parse_args()
my_broken_pairs = args.broken_pairs
my_sim_scores = args.sim_scores
my_overlap = args.overlap
my_species = args.species
my_output = args.output

##### Main #####
#Read input
broken_input_df = pd.read_table(str(my_broken_pairs), sep="\t", index_col=False, header=None, names=["OG_ID", "GeneID1", "GeneID2", "GeneName"])
sim_scores_df = pd.read_table(str(my_sim_scores), sep="\t", header=None, names=["GeneID1", "GeneID2", "Sim_score"])
overlap_df = pd.read_table(str(my_overlap), sep="\t", header=None, names=["GeneID1", "GeneID2", "Overlap"])
#Create dictionary with each of the genes first.
sim_scores_df["Sim_score"] = pd.Series([round(element, 4) for element in list(sim_scores_df["Sim_score"])]) #round to 4 significant digits
sim_scores_dict = pd.Series(sim_scores_df.Sim_score.values, index=sim_scores_df.GeneID1).to_dict()
overlap_df["Overlap"] = pd.Series([round(element, 4) for element in list(overlap_df["Overlap"])]) #round to 4 significant digits
overlap_dict = pd.Series(overlap_df.Overlap.values, index=overlap_df.GeneID1).to_dict()
#Translate GeneID1 and GeneID2 with the respective sim_score/overlap from the same dict
broken_input_df["GeneID1_sim_score"] = broken_input_df["GeneID1"].map(sim_scores_dict)
broken_input_df["GeneID2_sim_score"] = broken_input_df["GeneID2"].map(sim_scores_dict)
broken_input_df["GeneID1_overlap"] = broken_input_df["GeneID1"].map(overlap_dict)
broken_input_df["GeneID2_overlap"] = broken_input_df["GeneID2"].map(overlap_dict)
#Add species
broken_input_df["Species"] = pd.Series([my_species for element in list(range(broken_input_df.shape[0]))])
#Rearrange columns
broken_input_df = broken_input_df[["OG_ID", "Species", "GeneID1", "GeneID2", "GeneID1_sim_score", "GeneID2_sim_score", "GeneID1_overlap", "GeneID2_overlap"]]
#Save to file
broken_input_df.to_csv(str(my_output), sep="\t", index=False, header=True, na_rep="NA")
