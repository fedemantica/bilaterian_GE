#!/usr/bin/env python3

import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Script to combine the overlap info obtained by matching each broken gene pair to their longest human ortholog.")
parser.add_argument("--broken_pairs", "-b", required=True, metavar="broken_pairs", help="Pairs of broken genes for the species of interest. Format is col1=OG_ID, col2=GeneID1, col2=GeneID2")
parser.add_argument("--human_pairs", "-hp", required=True, metavar="human_pairs", help="Pairs of broken genes associated with the relative longest human ortholog. Format is col1=GeneID1, col2=GeneID2, col3=Human_geneID")
parser.add_argument("--interval_overlap", "-iv", required=True, metavar="interval_overlap", help="File with overlap information with the longest human ortholog. Format is col1=GeneID1, col2=Human_geneID, col3=Overlap_start, col4=Overlap_stop, col5=Overlap_coords. Coords refer to the human protein sequence")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file. All pairs of control genes.")

##### Read arguments #####
args = parser.parse_args()
my_broken_pairs = args.broken_pairs
my_human_pairs = args.human_pairs
my_interval_overlap = args.interval_overlap
my_output = args.output

##### Main #####
#Read inputs
broken_pairs_df = pd.read_table(str(my_broken_pairs), sep="\t", index_col=False, header=None, names=["OG_ID", "GeneID1", "GeneID2"])
longest_human_df = pd.read_table(str(my_human_pairs), sep="\t", index_col=False, header=None, names=["OG_ID", "GeneID1", "GeneID2"])
intervals_df = pd.read_table(str(my_interval_overlap), sep="\t", index_col=False, header=None, names=["GeneID", "Hsa_geneID", "Start", "Stop", "Interval", "Protein_length"])

#Filter the broken_genes_df
orthogroups_to_save = list(longest_human_df["OG_ID"])
broken_pairs_df = broken_pairs_df.loc[broken_pairs_df.OG_ID.isin(orthogroups_to_save)]

#Create dictionaries to add information to the final df
geneID_overlap_dict = pd.Series(intervals_df.Interval.values, index=intervals_df.GeneID).to_dict()
geneID_start_dict = pd.Series(intervals_df.Start.values, index=intervals_df.GeneID).to_dict()
geneID_stop_dict = pd.Series(intervals_df.Stop.values, index=intervals_df.GeneID).to_dict()
#NB: this is the geneID of the other species but the length of the human protein as value
protein_len_dict = pd.Series(intervals_df.Protein_length.values, index=intervals_df.GeneID).to_dict()
geneID_hsa_gene = pd.Series(intervals_df.Hsa_geneID.values, index=intervals_df.GeneID).to_dict()

broken_pairs_df["Hsa_geneID"] = broken_pairs_df["GeneID1"].map(geneID_hsa_gene)
broken_pairs_df["Hsa_prot_len"] = broken_pairs_df["GeneID1"].map(protein_len_dict)
#These are the ones I will print. I also need them separated because of stupid python
broken_pairs_df["GeneID1_match"] = broken_pairs_df["GeneID1"].map(geneID_overlap_dict)
broken_pairs_df["GeneID2_match"] = broken_pairs_df["GeneID2"].map(geneID_overlap_dict)
broken_pairs_df["GeneID1_match_start"] = broken_pairs_df["GeneID1"].map(geneID_start_dict)
broken_pairs_df["GeneID1_match_stop"] = broken_pairs_df["GeneID1"].map(geneID_stop_dict)
broken_pairs_df["GeneID2_match_start"] = broken_pairs_df["GeneID2"].map(geneID_start_dict)
broken_pairs_df["GeneID2_match_stop"] = broken_pairs_df["GeneID2"].map(geneID_stop_dict)
#I need to add this line because there are some chimeric genes (how shitty is that??) that end up in more orthogroups. But each gene ID can be repeated only once as a key.
broken_pairs_df = broken_pairs_df.dropna()
# Remove cases where the overlap in one of the genes is 0 to begin with. These will not be included in the final dataframe (Jan 2025)
broken_pairs_df = broken_pairs_df.loc[~(((broken_pairs_df["GeneID1_match_start"]==0) & (broken_pairs_df["GeneID1_match_stop"]==0)) | ((broken_pairs_df["GeneID2_match_start"]==0) & (broken_pairs_df["GeneID2_match_stop"]==0)))]

my_geneID1_match_list = [range(int(start), int(stop)) for start, stop in zip(list(broken_pairs_df["GeneID1_match_start"]), list(broken_pairs_df["GeneID1_match_stop"]))]
my_geneID2_match_list = [range(int(start), int(stop)) for start, stop in zip(list(broken_pairs_df["GeneID2_match_start"]), list(broken_pairs_df["GeneID2_match_stop"]))]
#Compute the range overlap.
my_match_overlap = [range(max(my_geneID1_match_list[x][0], my_geneID2_match_list[x][0]), min(my_geneID1_match_list[x][-1], my_geneID2_match_list[x][-1])+1) for x in list(range(broken_pairs_df.shape[0]))]
my_match_overlap_len = [len(element) for element in my_match_overlap]

#Add to final_df
broken_pairs_df["Match_overlap"] = my_match_overlap_len
######## Compute the match_overlap_proportion based on the shortest match (changed on 26/03/2021)
broken_pairs_df["GeneID1_match_length"] = broken_pairs_df["GeneID1_match_stop"] - broken_pairs_df["GeneID1_match_start"]
broken_pairs_df["GeneID2_match_length"] = broken_pairs_df["GeneID2_match_stop"] - broken_pairs_df["GeneID2_match_start"]
broken_pairs_df["Higher_match"] = [element[0] if float(element[0]) >= float(element[1]) else element[1] for element in list(zip(list(broken_pairs_df["GeneID1_match_length"]), list(broken_pairs_df["GeneID2_match_length"])))]
broken_pairs_df["Match_overlap_proportion"] = broken_pairs_df["Match_overlap"]/broken_pairs_df["Higher_match"]

##############
final_df = broken_pairs_df.loc[:,["OG_ID", "GeneID1", "GeneID2", "Hsa_geneID", "Hsa_prot_len", "GeneID1_match", "GeneID2_match", "Match_overlap", "Match_overlap_proportion"]]
#Write final df to file
final_df.to_csv(str(my_output), sep="\t", header=True, index=False, na_rep="NA")
