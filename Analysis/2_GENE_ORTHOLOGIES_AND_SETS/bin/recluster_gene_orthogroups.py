#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to recluster the gene orthogroups based on the orthopairs connections between a subset of species")
parser.add_argument("--orthogroups", "-og", required=True, metavar="species", help="Corrected gene orthogroups, with col1=orthogroupID, col2=species, col3=geneID")
parser.add_argument("--orthopairs", "-op", required=True, metavar="geneIDs", help="Broccoli orthopairs but in the format col1=species1, col2=geneID1, col=species2, col4=geneID2")
parser.add_argument("--clade_species", "-s", required=True, metavar="species", help="Comma separated string containing the selected species. Ex: 'Hs2,Mm2,Bt2'")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
orthogroups_file = args.orthogroups
orthopairs_file = args.orthopairs
clade_species_string = args.clade_species
output_file = args.output

##################################
###### MAIN ######################
##################################
#Read inputs
orthogroups_df = pd.read_table(str(orthogroups_file), sep="\t", header=None, names=["OG_ID", "species", "geneID"])
orthopairs_df = pd.read_table(str(orthopairs_file), sep="\t", header=None, names=["species1", "geneID1", "species2", "geneID2"])
clade_species_list = clade_species_string.split(",")

#Integrate orthogroups information
geneID_OG_ID_dict = pd.Series(orthogroups_df.OG_ID.values, index=orthogroups_df.geneID).to_dict() #Generate dictionary with key=geneID, value=OG_ID
orthopairs_df["OG_ID1"] = orthopairs_df["geneID1"].map(geneID_OG_ID_dict) #Add OG_ID to orthopairs df based on geneID1
orthopairs_df = orthopairs_df[["OG_ID1", "species1", "geneID1", "species2", "geneID2"]] #Reorder orthopairs df
orthopairs_df["OG_ID2"] = orthopairs_df["geneID2"].map(geneID_OG_ID_dict) #Add OG_ID to orthopairs df based on geneID2

#Filter for the orthopairs where both genes end up in the same orthogroups.
#This is not always the case, because chimeric genes have been assigned to only one orthogroup.
#Some of the original chimeric orthopairs will therefore be invalid.
orthopairs_df = orthopairs_df.loc[orthopairs_df["OG_ID1"] == orthopairs_df["OG_ID2"]]
orthopairs_df = orthopairs_df.drop(columns="OG_ID2")
orthopairs_df = orthopairs_df.rename(columns={"OG_ID1" : "OG_ID"})
#Filter for orthopairs involving genes between the selected subset of species
orthopairs_filt_df = orthopairs_df.loc[(orthopairs_df.species1.isin(clade_species_list)) & (orthopairs_df.species2.isin(clade_species_list))]
                
#Initialize final dataframe with the same column as the orthogroups df
final_df = pd.DataFrame(columns=list(orthogroups_df.columns.values))
#Group the orthopairs df
grouped_df = orthopairs_filt_df.groupby("OG_ID")
#Cycle on the orthopairs df
for name, group in grouped_df:
  n=1 #Initialize subcluster number
  original_group = group

  while group.shape[0]>=1:
    index="0"+str(n)
    first_gene = group.iloc[0,:].loc["geneID1"] #select a first random gene
    second_genes = list(original_group.loc[original_group.geneID1==first_gene]["geneID2"]) #Select all the direct orthologs of the first random gene
    third_genes = list(original_group.loc[original_group.geneID1.isin(second_genes)]["geneID2"]) #Select all the direct orthologs of the direct orthologs of the first random gene
    all_genes = list(set([first_gene] + second_genes + third_genes)) #Join all geneIDs in a single list
    #Repeat cycle until no new genes are retrieved
    cycle_value = 1
    while cycle_value ==1:
      new_genes_1 = list(original_group.loc[original_group.geneID1.isin(all_genes)]["geneID2"]) 
      new_genes_2 = list(original_group.loc[original_group.geneID2.isin(all_genes)]["geneID1"])
      ##
      new_genes_3 = new_genes_1 + new_genes_2
      new_genes_3 = list(set(new_genes_3))
      #new_genes_3 = list(set(new_genes_1, new_genes_2))
      new_genes = [gene for gene in new_genes_3 if gene not in all_genes]
      if len(new_genes) == 0: #If there are no more new connected genes
        cycle_value = 0 #Exit the loop
      else:
        all_genes = all_genes + new_genes

    #Filter group only for selected genes
    final_group = orthogroups_df.loc[orthogroups_df.geneID.isin(all_genes)] 
    final_group["OG_ID"] = [element+"."+index for element in list(final_group["OG_ID"])]
    ##### This is an extra filter. Make sure to save to output only the orthogroups with a minimum number of species
    final_group_species = list(set(list(final_group["species"])))
    if len(final_group_species) >= 6: #Here I am requiring at least 6 species (meaning at least 6 Vertebrates or 6 Insects)
      final_df = pd.concat([final_df, final_group]) #Concatenate to final df
    #####
    group = group.loc[~group.geneID1.isin(all_genes)] #Filter out the selected genes from the group dataframe
    n = n+1 #generate another sub-cluster

#Save to output file
final_df.to_csv(output_file, sep="\t", header=False, index=False, na_rep="NA")
