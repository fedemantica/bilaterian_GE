#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to correct the parsed orthogroups by replacing the broken and chimeric original geneIDs with the fixed ones")
parser.add_argument("--orthogroups", "-og", required=True, metavar="species", help="Parsed gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--geneIDs", "-g", required=True, metavar="geneIDs", help="Output of generate_fixed_geneIDs.py, col1=original_ID, col2=new_ID, col3=category (broken || chimeric)")
parser.add_argument("--unresolved_chimeric", "-uc", required=True, metavar="unresolved_chimeric", help="File with all the chimeric genes that could not be splitted into new genes. The last field report the gene orthogroup with the longest aligned fragment, corresponding to the one where the gene should be saved")
parser.add_argument("--resolved_chimeric", "-rc", required=True, metavar="resolved_chimeric", help="File with all the chimeric genes that could be splitted. The last field is in the format OG1;new_OD1|OG2;new_ID2, where the new_IDs correspond to the geneIDs assigned to the novel splitted genes matching each orthogroup")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
geneIDs_file = args.geneIDs #This need to be collapsed for all the species
unresolved_chimeric_file = args.unresolved_chimeric #This need to be collapsed for all the species
resolved_chimeric_file = args.resolved_chimeric #This need to be collapsed for all the species
output_file = args.output

##################################
###### READ INPUTS ###############
##################################
#Header: ['OG_ID', 'species', 'geneID']
orthogroup_df = pd.read_table(str(orthogroups_file), sep="\t", index_col=False, header=None, names=["OG_ID", "species", "geneID"])
#Header: ['geneID', 'new_IDs', 'category']
geneIDs_df = pd.read_table(str(geneIDs_file), sep="\t", index_col=False, header=0)
#Header: ['species', 'chimeric_geneID', 'orthogroup_ids', 'first-last_aligned_ex', 'first-last_aligned_aa', 'ex:start|stop_coord', 'chimeric_class', 'ex_overlap']
unresolved_chimeric_df = pd.read_table(str(unresolved_chimeric_file), sep="\t", index_col=False, header=0)
#Header: ['species', 'chimeric_geneID', 'orthogroup_ids', 'first-last_aligned_ex', 'first-last_aligned_aa', 'ex:start|stop_coord', 'chimeric_class']
resolved_chimeric_df = pd.read_table(str(resolved_chimeric_file), sep="\t", index_col=False, header=0)

### debugging
#orthogroup_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version1/parsed_orthogroups.txt", sep="\t", index_col=False, header=None, names=["OG_ID", "species", "geneID"])
#geneIDs_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version1/corrected_gtfs/all_species-new_geneIDs.txt", sep="\t", index_col=False, header=0)
#unresolved_chimeric_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version1/corrected_gtfs/all_species-unresolved_chimeric_genes.tab", sep="\t", index_col=False, header=0)
#resolved_chimeric_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version1/corrected_gtfs/all_species-resolved_chimeric_genes.tab", sep="\t", index_col=False, header=0)

##################################
###### PRE-MANIPULATION ##########
##################################

#### Modify the resolved chimeric df to add the new chimeric geneIDs
#create a dictionary with key = chimericID, value=newIDs
chimeric_geneIDs_df = geneIDs_df.loc[geneIDs_df["category"]=="chimeric"]
reverse_geneIDs_dict = pd.Series(chimeric_geneIDs_df.new_IDs.values, index=chimeric_geneIDs_df.geneID).to_dict()
resolved_chimeric_df["new_IDs"] = resolved_chimeric_df["chimeric_geneID"].map(reverse_geneIDs_dict)

resolved_chimeric_dict_df = resolved_chimeric_df[["chimeric_geneID", "orthogroup_ids", "new_IDs"]]
resolved_chimeric_dict_df["OG1"] = [element.split(";")[0] for element in list(resolved_chimeric_dict_df["orthogroup_ids"])]
resolved_chimeric_dict_df["OG2"] = [element.split(";")[1] for element in list(resolved_chimeric_dict_df["orthogroup_ids"])]
resolved_chimeric_dict_df["chimeric_ID1"] = [element.split(";")[0] for element in list(resolved_chimeric_dict_df["new_IDs"])]
resolved_chimeric_dict_df["chimeric_ID2"] = [element.split(";")[1] for element in list(resolved_chimeric_dict_df["new_IDs"])]

#Create dictionary to translate the chimeric_geneID;OG into the new_ID (a chimeric ID)
OG1_dict = {element[0]+";"+element[1] : element[2] for element in list(zip(resolved_chimeric_dict_df["chimeric_geneID"], resolved_chimeric_dict_df["OG1"], resolved_chimeric_dict_df["chimeric_ID1"]))}
OG2_dict = {element[0]+";"+element[1] : element[2] for element in list(zip(resolved_chimeric_dict_df["chimeric_geneID"], resolved_chimeric_dict_df["OG2"], resolved_chimeric_dict_df["chimeric_ID2"]))}
resolved_chimeric_dict = {}
resolved_chimeric_dict.update(OG1_dict)
resolved_chimeric_dict.update(OG2_dict)

#### Modify the unresolved chimeric to determine the OG with the longest overlap
#this code works, but it is not super easy to read...
unresolved_chimeric_df["first_fragment"] = [float(element.split(";")[0].split("-")[1])-float(element.split(";")[0].split("-")[0]) for element in list(unresolved_chimeric_df["first-last_aligned_aa"])]
unresolved_chimeric_df["second_fragment"] = [float(element.split(";")[1].split("-")[1])-float(element.split(";")[1].split("-")[0]) for element in list(unresolved_chimeric_df["first-last_aligned_aa"])]
unresolved_chimeric_df["longest_fragment"] = ["OG1" if list(unresolved_chimeric_df["first_fragment"])[index] > list(unresolved_chimeric_df["second_fragment"])[index] else "OG2" for index in list(unresolved_chimeric_df.index.values)]
unresolved_chimeric_df["selected_orthogroup"] = [list(unresolved_chimeric_df["orthogroup_ids"])[index].split(";")[0] if list(unresolved_chimeric_df["longest_fragment"])[index]=="OG1" else list(unresolved_chimeric_df["orthogroup_ids"])[index].split(";")[1] for index in list(unresolved_chimeric_df.index.values)]
#create a list with chimericID;OG_to_save
chimeric_gene_selected_OG_list = [element[0]+";"+element[1] for element in zip(list(unresolved_chimeric_df["chimeric_geneID"]), list(unresolved_chimeric_df["selected_orthogroup"]))]

##################################
###### MAIN ######################
##################################

#####################
##### Broken ########
#####################
#subset orthogroup to only broken genes
broken_genes_raw = list(geneIDs_df.loc[geneIDs_df["category"]=="broken"]["geneID"])
broken_genes = [part for element in broken_genes_raw for part in element.split(";")]
broken_orthogroups_df = orthogroup_df.loc[orthogroup_df["geneID"].isin(broken_genes)]

#generate dictionary to translate the broken geneIDs to the corrected geneID
broken_geneIDs_df = geneIDs_df.loc[geneIDs_df["category"]=="broken"]
geneIDs_dict = pd.Series(broken_geneIDs_df.new_IDs.values, index=broken_geneIDs_df.geneID).to_dict()
#broken_genes_keys = [key for key in list(geneIDs_dict.keys()) if ";" in key]
for key in list(geneIDs_dict.keys()):
  for num in list(range(0, len(key.split(";")))):
    geneIDs_dict[key.split(";")[num]] = geneIDs_dict[key]

#translate with new geneID and remove duplicated entries
broken_orthogroups_df["geneID"] = broken_orthogroups_df["geneID"].map(geneIDs_dict)
broken_orthogroups_df = broken_orthogroups_df.drop_duplicates() #remove duplicated entries

#####################
##### Chimeric ######
#####################
#subset only by unresolved chimeric genes
unresolved_chimeric_genes = list(unresolved_chimeric_df["chimeric_geneID"])
unresolved_chimeric_orthogroups_df = orthogroup_df.loc[orthogroup_df["geneID"].isin(unresolved_chimeric_genes)]

#unresolved_chimeric_orthogroups_df = unresolved_chimeric_orthogroups_df.loc[unresolved_chimeric_orthogroups_df["OG_ID"]]
#create an entry with geneID;OGID
unresolved_chimeric_orthogroups_df["geneID;OG_ID"] = [element[0]+";"+element[1] for element in zip(list(unresolved_chimeric_orthogroups_df["geneID"]), list(unresolved_chimeric_orthogroups_df["OG_ID"]))]
#keep only the entry where the OG_ID corresponds to the longest fragment
selected_unresolved_chimeric_orthogroups_df = unresolved_chimeric_orthogroups_df.loc[unresolved_chimeric_orthogroups_df["geneID;OG_ID"].isin(chimeric_gene_selected_OG_list)]

#subset only by resolved chimeric genes
resolved_chimeric_genes = list(resolved_chimeric_df["chimeric_geneID"])
#create OG:chimeric entry -> translate with the relative new chimeric geneID
resolved_chimeric_orthogroups_df = orthogroup_df.loc[orthogroup_df["geneID"].isin(resolved_chimeric_genes)]
resolved_chimeric_orthogroups_df["geneID;OG_ID"] = [element[0]+";"+element[1] for element in list(zip(list(resolved_chimeric_orthogroups_df["geneID"]), list(resolved_chimeric_orthogroups_df["OG_ID"])))]
resolved_chimeric_orthogroups_df["geneID"] = resolved_chimeric_orthogroups_df["geneID;OG_ID"].map(resolved_chimeric_dict)
#08/09/21: remove all the entries whose geneID translation is NA (i.e. chimeric genes shared between more than 2 orthogroups, which can only be splitted only in two parts. The original entries for the other parts have to be remove)
resolved_chimeric_orthogroups_df = resolved_chimeric_orthogroups_df.dropna(subset=["geneID"])

#####################
#### Rejoin parts ###
#####################
#Isolate healthy dataframe
all_messy_genes = broken_genes + unresolved_chimeric_genes + resolved_chimeric_genes
healthy_df = orthogroup_df.loc[~(orthogroup_df["geneID"].isin(all_messy_genes))] #select all entries where the genes are neither broken nor chimeric
#join parts
final_orthogroup_df = pd.concat([healthy_df, broken_orthogroups_df, selected_unresolved_chimeric_orthogroups_df, resolved_chimeric_orthogroups_df])
final_orthogroup_df = final_orthogroup_df.drop(columns=["geneID;OG_ID"])
#reorder dataframe
final_orthogroup_df = final_orthogroup_df.sort_values(by=["OG_ID", "species"])

### Save to file ####
final_orthogroup_df.to_csv(output_file, sep="\t", header=False, index=False, na_rep="NA")
