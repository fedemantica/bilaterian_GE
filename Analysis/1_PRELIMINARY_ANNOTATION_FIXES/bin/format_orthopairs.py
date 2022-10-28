#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to recluster the gene orthogroups based on the orthopairs connections between a subset of species")
parser.add_argument("--orthogroups", "-og", required=True, metavar="parsed_orthogroups", help="Parsed gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--corrected_orthogroups", "-cg", required=True, metavar="corrected_orthogroups", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--orthopairs", "-op", required=True, metavar="orthopairs", help="Broccoli orthopairs but in the format col1=geneID1, col2=geneID2")
parser.add_argument("--geneIDs", "-g", required=True, metavar="geneIDs", help="File with correspondence between the original gene ID and the new geneID for the broken and chimeric genes")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
corrected_orthogroups_file = args.corrected_orthogroups
orthopairs_file = args.orthopairs
geneIDs_file = args.geneIDs
output_file = args.output

##################################
###### MAIN ######################
##################################
###Read inputs
#Orthopairs
orthopairs_df = pd.read_table(orthopairs_file, sep="\t", index_col=False, header=None, names=["geneID1", "geneID2"])
orthopairs_df["geneID1"] = [re.sub(".*\\|", "", element) for element in list(orthopairs_df["geneID1"])]
orthopairs_df["geneID2"] = [re.sub(".*\\|", "", element) for element in list(orthopairs_df["geneID2"])]
#Orthogroups
orthogroups_df = pd.read_table(orthogroups_file, sep="\t", index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
corrected_orthogroups_df = pd.read_table(corrected_orthogroups_file, index_col=False, header=None, names=["OG_ID", "Species", "GeneID"])
#Header is: ["geneID", "new_IDs", "category"]
geneIDs_df = pd.read_table(geneIDs_file, sep="\t", index_col=False, header=0)

##################################
###### BROKEN ####################
##################################
#Isolate original_broken_geneIDs
nested_broken_geneIDs = list(geneIDs_df.loc[geneIDs_df["category"]=="broken"]["geneID"])
nested_broken_gene_IDs = [element.split(";") for element in nested_broken_geneIDs]
broken_gene_IDs = [item for sublist in nested_broken_gene_IDs for item in sublist]
#Filter orthopairs_df: col1 || col2 in original broken geneIDs
first_broken_orthopairs_df = orthopairs_df.loc[orthopairs_df["geneID1"].isin(broken_gene_IDs)]
second_broken_orthopairs_df = orthopairs_df.loc[orthopairs_df["geneID2"].isin(broken_gene_IDs)]
#Move all broken genes to col2
first_broken_orthopairs_df["geneID3"] = first_broken_orthopairs_df["geneID1"]
first_broken_orthopairs_df["geneID1"] = first_broken_orthopairs_df["geneID2"]
first_broken_orthopairs_df["geneID2"] = first_broken_orthopairs_df["geneID3"]
first_broken_orthopairs_df = first_broken_orthopairs_df.drop(columns=["geneID3"])
broken_orthopairs_df = pd.concat([first_broken_orthopairs_df, second_broken_orthopairs_df])

#Create dictionary with key=original_broken_geneID, value=new_broken_geneID
broken_geneIDs_df = geneIDs_df.loc[geneIDs_df["category"]=="broken"]
broken_geneIDs_dict = pd.Series(broken_geneIDs_df.new_IDs.values, index=broken_geneIDs_df.geneID).to_dict()
for key in list(broken_geneIDs_dict.keys()):
  for num in list(range(0, len(key.split(";")))):
    broken_geneIDs_dict[key.split(";")[num]] = broken_geneIDs_dict[key]

#Translate column2 with new_broken_geneID from the previously generated dictionary
broken_orthopairs_df["new_broken_ID"] = broken_orthopairs_df["geneID2"].map(broken_geneIDs_dict)
broken_orthopairs_df = broken_orthopairs_df[["geneID1", "new_broken_ID"]]
broken_orthopairs_df = broken_orthopairs_df.rename(columns={"new_broken_ID" : "geneID2"})

##################################
###### CHIMERIC ##################
##################################
#Isolate original chimeric geneIDs
chimeric_gene_IDs = list(geneIDs_df.loc[geneIDs_df["category"]=="chimeric"]["geneID"])
#Filter orthopairs_df: col1 || col2 in original_chimeric_geneIDs
first_chimeric_orthopairs_df = orthopairs_df.loc[orthopairs_df["geneID1"].isin(chimeric_gene_IDs)]
second_chimeric_orthopairs_df = orthopairs_df.loc[orthopairs_df["geneID2"].isin(chimeric_gene_IDs)]
#Move all chimeric geneIDs to col2
first_chimeric_orthopairs_df["geneID3"] = first_chimeric_orthopairs_df["geneID1"]
first_chimeric_orthopairs_df["geneID1"] = first_chimeric_orthopairs_df["geneID2"]
first_chimeric_orthopairs_df["geneID2"] = first_chimeric_orthopairs_df["geneID3"]
first_chimeric_orthopairs_df = first_chimeric_orthopairs_df.drop(columns=["geneID3"])
chimeric_orthopairs_df = pd.concat([first_chimeric_orthopairs_df, second_chimeric_orthopairs_df])

#Starting from orthogroups_df, create a dictionary with key=geneID, value=OG_ID
geneID_OGID_dict = pd.Series(orthogroups_df.OG_ID.values, index=orthogroups_df.GeneID).to_dict()
#Add OG_ID to each row based on the value present in col1 (pair of the chimeric gene which is not itself chimeric)
chimeric_orthopairs_df["OG_ID"] = chimeric_orthopairs_df["geneID1"].map(geneID_OGID_dict)
#Add a complete ID composed of OG_ID;geneID2 (This is used to map the new chimeric gene).
chimeric_orthopairs_df["CompleteID"] = [element[0]+";"+element[1] for element in zip(list(chimeric_orthopairs_df["OG_ID"]), list(chimeric_orthopairs_df["geneID2"]))]
#Generate chimeric geneID datafame
chimeric_geneIDs_df = geneIDs_df.loc[geneIDs_df["category"]=="chimeric"]
chimeric_geneIDs_df["new_IDs"] = [element.split(";") for element in list(chimeric_geneIDs_df["new_IDs"])]
chimeric_geneIDs_df = chimeric_geneIDs_df.explode("new_IDs") #Split the columns with the new geneIDs
#Add column with the gene orthogroup to the chimeric geneID dataframe
nested_new_chimeric_geneIDs = list(geneIDs_df.loc[geneIDs_df["category"]=="chimeric"]["new_IDs"])
nested_new_chimeric_geneIDs = [element.split(";") for element in nested_new_chimeric_geneIDs]
new_chimeric_geneIDs = [item for sublist in nested_new_chimeric_geneIDs for item in sublist]
chimeric_corrected_orthogroups_df = corrected_orthogroups_df.loc[corrected_orthogroups_df["GeneID"].isin(new_chimeric_geneIDs)] #I need to start from the corrected orthogroups
new_chimeric_OGID_dict = pd.Series(chimeric_corrected_orthogroups_df.OG_ID.values, index=chimeric_corrected_orthogroups_df.GeneID)
chimeric_geneIDs_df["OG_ID"] = chimeric_geneIDs_df["new_IDs"].map(new_chimeric_OGID_dict)
#remove the entries with NA, corresponding to unresolved chimeric genes 
chimeric_geneIDs_df = chimeric_geneIDs_df.dropna()

#Starting from the geneID df, create a dictionary with key=OG_ID;geneID, value=old chimeric geneID 
chimeric_geneIDs_df["CompleteID"] = [element[0]+";"+element[1] for element in zip(list(chimeric_geneIDs_df["OG_ID"]), list(chimeric_geneIDs_df["geneID"]))]
complete_ID_new_chimeric_dict = pd.Series(chimeric_geneIDs_df.new_IDs.values, index=chimeric_geneIDs_df.CompleteID).to_dict()

#Add new chimeric geneID to the orthopairs_df
chimeric_orthopairs_df["new_chimeric"] = chimeric_orthopairs_df["CompleteID"].map(complete_ID_new_chimeric_dict)
#Filter out columns with NA (these again correspond to unresolved chimeric genes)
chimeric_orthopairs_df = chimeric_orthopairs_df.dropna()
#Filter out uninteresting columns
chimeric_orthopairs_df = chimeric_orthopairs_df[["geneID1", "new_chimeric"]]
chimeric_orthopairs_df = chimeric_orthopairs_df.rename(columns = {"new_chimeric" : "geneID2"})

##################################
###### HEALTHY GENES #############
##################################
#Select all the genes in the orthopairs df where neither geneID1 nor geneID2 are broken/chimeric
all_problematic_genes = chimeric_gene_IDs + broken_gene_IDs
first_healthy_orthopairs_df = orthopairs_df.loc[~(orthopairs_df["geneID1"].isin(all_problematic_genes))]
second_healthy_orthopairs_df = orthopairs_df.loc[~(orthopairs_df["geneID2"].isin(all_problematic_genes))]
healthy_orthopairs_df = pd.concat([first_healthy_orthopairs_df, second_healthy_orthopairs_df])

##################################
###### JOIN PARTS ################
##################################
all_new_orthopairs_df = pd.concat([broken_orthopairs_df, chimeric_orthopairs_df, healthy_orthopairs_df])

##################################
###### ADD SPECIES ###############
##################################
#Get a dictionary with key=GeneID and value=Species (starting from the corrected gene orthogroups)
geneID_species_dict = pd.Series(corrected_orthogroups_df.Species.values, index=corrected_orthogroups_df.GeneID).to_dict()
#Add the species to the relative genes
all_new_orthopairs_df["species1"] = all_new_orthopairs_df["geneID1"].map(geneID_species_dict)
all_new_orthopairs_df["species2"] = all_new_orthopairs_df["geneID2"].map(geneID_species_dict)
#Repeat the df, so that all the genes from the same specses actually appear in the first column
all_new_orthopairs_df = all_new_orthopairs_df[["species1", "geneID1", "species2", "geneID2"]]
all_new_orthopairs_alt_df = all_new_orthopairs_df[["species2", "geneID2", "species1", "geneID1"]]
all_new_orthopairs_alt_df.columns = ["species1", "geneID1", "species2", "geneID2"] #not the safest operation ever, but I want the genes from all species to be repeated
#Concatenate in a unique dataframe
final_df = pd.concat([all_new_orthopairs_df, all_new_orthopairs_alt_df])
final_df = final_df.dropna() #remove entries relative to original broken genes, which do not have species translation.

##################################
###### SAVE OUTPUT ###############
##################################
#remove duplicates
final_df = final_df.drop_duplicates()
#remove entries with species1 == species2
final_df = final_df.loc[final_df["species1"] != final_df["species2"]]
#save to output file
final_df.to_csv(output_file, sep="\t", header=False, index=False, na_rep="NA")
