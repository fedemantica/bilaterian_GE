#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to merge the gene orthogroups when known human onhologs have not been assigned to the same orthogroup.")
parser.add_argument("--orthogroups", "-og", required=True, metavar="orthogroups", help="Parsed and corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--human_orthologs", "-ho", required=True, metavar="human_orthologs", help="Ensembl one2one correspondences between human and mouse or human and cow (in case the mouse ortholog is missing)")
parser.add_argument("--onhologs", "-on", required=True, metavar="onhologs", help="File containing known human onholog groups, with col1=Onholog_groupID, col2=GeneID, col3=GeneName")
parser.add_argument("--output_human", "-oh", required=True, metavar="output_human", help="Path to output file: extra human entries")
parser.add_argument("--output_orthogroups", "-oo", required=True, metavar="output_orthogroups", help="Path to output file: enriched and corrected orthogroups")
parser.add_argument("--output_merged", "-om", required=True, metavar="output_merged", help="Path to output file: correspondences between old and new orthogroup IDs")
parser.add_argument("--output_mixed", "-omix", required=True, metavar="output_mixed", help="Path to output file: problematic entries to fix manually (mixed up onhologs)")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
human_orthologs_file = args.human_orthologs
onhologs_file = args.onhologs
output_human_file = args.output_human
output_orthogroups_file = args.output_orthogroups
output_merged_file = args.output_merged
output_mixed_file = args.output_mixed

#####################################
######### READ INPUTS ###############
#####################################

#Read inputs
##Debugging
#orthogroups_df = pd.read_table("corrected_orthogroups.txt", sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#human_orthologs_df = pd.read_table("~/projects/bilaterian_GE/data/DB/ensembl_orthologs/Hs2-Mm2_Bt2_one2one_orthologs_NOredundant", sep="\t", index_col=False, header=None, names=["Human_GeneID", "Ortholog_GeneID", "Ortholog_Status"])
#onhologs_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/DB/onholog_info/Ohnologs-Hs2-V6-formatted.tab", sep="\t", index_col=False, header=None, names=["Onholog_ID", "GeneID", "GeneName"])
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
human_orthologs_df = pd.read_table(human_orthologs_file, sep="\t", index_col=False, header=None, names=["Human_GeneID", "Ortholog_GeneID", "Ortholog_Status"])
onhologs_df = pd.read_table(onhologs_file, sep="\t", index_col=False, header=None, names=["Onholog_ID", "GeneID", "GeneName"])

######################################
##### ENRICH HUMAN ORHTOLOGY #########
######################################

#Filter human_orthologs_df for human entries not included in the orthogroups
conserved_human_geneIDs = list(orthogroups_df.loc[orthogroups_df["Species"]=="Hs2"]["GeneID"])
excluded_human_orthologs_df = human_orthologs_df.loc[~human_orthologs_df["Human_GeneID"].isin(conserved_human_geneIDs)]
#Generate a dictionary with key=mouse/cow_geneID, value=orthogroup_ID
mouse_cow_orthogroups_df = orthogroups_df.loc[orthogroups_df["Species"].isin(["Mm2", "Bt2"])]
mouse_cow_geneID_OGID_dict = pd.Series(mouse_cow_orthogroups_df.OG_ID.values, index=mouse_cow_orthogroups_df.GeneID).to_dict()
#Associate the missing human ortholog with a orthogroupID, if possible, and generate a dataframe with the extra human entries
human_to_ortholog_geneID_dict = pd.Series(excluded_human_orthologs_df.Ortholog_GeneID.values, index=excluded_human_orthologs_df.Human_GeneID).to_dict()
extra_human_orthogroups_df = pd.DataFrame()
unconserved_df = pd.DataFrame() #This is just for debugging
for human_gene in list(human_to_ortholog_geneID_dict.keys()):
  ortholog = human_to_ortholog_geneID_dict[human_gene]
  if ortholog in list(mouse_cow_geneID_OGID_dict.keys()):
    OG_ID = mouse_cow_geneID_OGID_dict[ortholog]
    extra_human_orthogroups_df = pd.concat([extra_human_orthogroups_df, pd.DataFrame({"OG_ID" : [OG_ID], "Species" : ["Hs2"], "GeneID" : [human_gene]})])
  else:
    unconserved_df = pd.concat([unconserved_df, pd.DataFrame({"Human_ID" : [human_gene], "Ortholog" : [ortholog]})]) #This is just for debugging

#Merge the extra human entries dataframe with the original one and reorder
merged_orthogroups_df = pd.concat([orthogroups_df, extra_human_orthogroups_df]).sort_values(by=["OG_ID", "Species"])

######################################
##### MERGE ONHOLOGS ORTHOGROUPS #####
######################################

#Generate dictionary with key=geneID, value=OGID
geneID_OGID_dict = pd.Series(merged_orthogroups_df.OG_ID.values, index=merged_orthogroups_df.GeneID).to_dict()
onhologs_df["OG_ID"] = onhologs_df["GeneID"].map(geneID_OGID_dict)

#Isolate all onholog groups with different OGID
presubsetted_onholog_df = onhologs_df[["Onholog_ID", "OG_ID"]].dropna().drop_duplicates(keep="first", inplace=False)
subsetted_onholog_df = presubsetted_onholog_df[presubsetted_onholog_df.duplicated(subset=["Onholog_ID"], keep=False)]
#Generate a list with the OnhologIDs for all those cases where the Onhologs end up in different groups.
duplicated_onhologs = list(set(list(subsetted_onholog_df["Onholog_ID"])))
duplicated_onhologs.sort() #sort the list

#Generate a common OGID starting from the highest one still available
highest_OGID = max(list(merged_orthogroups_df["OG_ID"])) #This works even with the alphanumeric string
OGID_length = len(re.sub("GF_", "", highest_OGID))
highest_OGID_num = int(re.sub("GF_", "", highest_OGID))
new_OGIDs = [highest_OGID_num+i for i in list(range(1, len(duplicated_onhologs)+1))] #Generate a alphanumeric ID for each groups of orthogroups to be merged.
new_OGIDs = ["GF_" + "0"*(OGID_length-len(str(element)))+str(element) for element in new_OGIDs]

#Generate dictionary with key=original_OGID, value=common_OGID; all other cases: key=original_OGID, value=original_OGID (fake dict entries)
onhologID_OGID_dict = pd.Series(new_OGIDs, index=duplicated_onhologs).to_dict()
subsetted_onholog_df["new_OGID"] = subsetted_onholog_df["Onholog_ID"].map(onhologID_OGID_dict)
mixed_up_onhologs_df = subsetted_onholog_df[subsetted_onholog_df.duplicated(subset=["OG_ID"], keep=False)].sort_values(by=["OG_ID"]) #Get the duplicated OG_ID entries, corresponding to cases where different onholog groups are mixed
final_onhologs_df = subsetted_onholog_df.drop_duplicates(subset=["OG_ID"], keep=False).sort_values(by=["new_OGID"]) #Remove the entries with duplicated OG_IDs 
OGID_newOG_ID_dict = pd.Series(final_onhologs_df.new_OGID.values, index=final_onhologs_df.OG_ID)

#Translate the original_geneID with the new ones
merged_orthogroups_df["new_OGID"] = [OG_ID if OG_ID not in list(final_onhologs_df["OG_ID"]) else OGID_newOG_ID_dict[OG_ID] for OG_ID in list(merged_orthogroups_df["OG_ID"])]
final_df = merged_orthogroups_df[["new_OGID", "Species", "GeneID"]]

######################################
######### SAVE TO OUTPUT #############
######################################
extra_human_orthogroups_df.to_csv(output_human_file, sep="\t", index=False, header=False, na_rep="NA") #Save to file separately the extra human entries (I will use them to modify the reference GTF)
final_df.to_csv(output_orthogroups_file, sep="\t", index=False, header=False, na_rep="NA") #Save the enriched and corrected orthogroups to file
final_onhologs_df.to_csv(output_merged_file, sep="\t", index=False, header=True, na_rep="NA") #Save the correspondences between old and new geneIDs to file
mixed_up_onhologs_df.to_csv(output_mixed_file, sep="\t", index=False, header=True, na_rep="NA") #Save to output those cases where different onhologs have been mixed up (to fix manually)
