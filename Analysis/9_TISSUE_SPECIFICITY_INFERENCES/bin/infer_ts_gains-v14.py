#!/usr/bin/env python3

import argparse
import json
import pandas as pd
import math
import re
from collections import Counter
import numpy as np

parser = argparse.ArgumentParser(description="Script to infer tissue-specificity gains along the phylogeny based on Tau")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Complete bilaterian conserved orthogroups (at least one gene per species) where I select the highest Tau")
parser.add_argument("--up_tau_cutoff", "-u", required=True, metavar="up tau cutoff", help="Stricter cutoff to infer TS gain")
parser.add_argument("--low_tau_cutoff", "-l", required=True, metavar="low tau cutoff", help="Milder cutoff to infer TS gain")
parser.add_argument("--deut_species", "-d", required=True, metavar="deuterostome species", help="Species on the deuterostome side of the tree")
parser.add_argument("--prot_species", "-p", required=True, metavar="protostome species", help="Species on the protostome side of the tree")
parser.add_argument("--ancestral_clades", "-a", required=True, metavar="protostome species", help="Nodes on both sides of the tree where an inferred gain will be mapped back to Bilateria")
parser.add_argument("--node_name_dict", "-n", required=True, metavar="node name dict", help="Json-formatted dictionary with key=earliest diverging species and value=node name")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")
parser.add_argument("--query_tissue", "-q", required=True, metavar="query tissue", help="Tissue for which the tissue-specificity gains are evaluated")

##########################
##### Read arguments #####
##########################
args = parser.parse_args()
my_input = args.input
up_tau_cutoff = float(args.up_tau_cutoff)
low_tau_cutoff = float(args.low_tau_cutoff)
species_to_classify_deut = args.deut_species.split(",")
species_to_classify_prot = args.prot_species.split(",")
ancestral_clades = args.ancestral_clades.split(",")
node_name_dict = json.loads(args.node_name_dict)
query_tissue = str(args.query_tissue)
my_output = args.output

###########################
##### Define function #####
###########################

#Transform to external arguments
def get_node_TS_OG(species1, species2, outer_species, OG_to_classify):
  #NB: species1 is a list
  #NB: species2 is a string
  node_OGs = []
  for OG_ID in OG_to_classify:
    #if the gene in the tested species is tissue-specific, above the tau cutoff and above the expr cutoff
    if (species_tau_dict[species2][OG_ID] >= up_tau_cutoff and query_tissue in species_tissue_dict[species2][OG_ID] and species_expr_dict[species2][OG_ID]=="YES_EXPR"):
      #inner_species_TS = [species for species in species1 if species_tau_dict[species][OG_ID] >= up_tau_cutoff]
      #outer_species_TS = [species for species in outer_species if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
      #Changed on 13/07: also check the tissue association
      inner_species_TS = [species for species in species1 if species_tau_dict[species][OG_ID] >= low_tau_cutoff and query_tissue in species_tissue_dict[species][OG_ID]]
      #Consider outer species only if they pass the expr cutoff
      outer_species_TS = [species for species in outer_species if species_tau_dict[species][OG_ID] >= low_tau_cutoff and query_tissue in species_tissue_dict[species][OG_ID] and species_expr_dict[species][OG_ID]=="YES_EXPR"]
      conserved_species = [species for species in species1 if ~np.isnan(species_tau_dict[species][OG_ID])] #count only the species that actually have a gene
      #if len(inner_species_TS) >= math.floor(len(species1)/2): #round down the 50%
      if len(inner_species_TS) >= math.floor(len(conserved_species)/2):
        if len(outer_species_TS) == 0:
          node_OGs = node_OGs + [OG_ID]
          OG_to_classify = [element for element in OG_to_classify if element != OG_ID] #remove the OG_IDs from the ones to classify
  return([node_OGs, OG_to_classify])   

def get_node_TS_OG_unclassified(species_to_check, OG_to_classify, node_name_dict, result_dict):
  final_OG_to_classify = OG_to_classify
  for OG_ID in OG_to_classify:
    #final_OG_to_classify = OG_to_classify
    #TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= up_tau_cutoff]
    #original_semi_TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
    TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= up_tau_cutoff and query_tissue in species_tissue_dict[species][OG_ID] and species_expr_dict[species][OG_ID]=="YES_EXPR"]
    original_semi_TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff and query_tissue in species_tissue_dict[species][OG_ID]]
    #### Exclude all species "higher" than the highest TS_species if they have NO_EXPR
    max_TS_species_indexes = max([idx for idx in range(len(species_to_check)) if species_to_check[idx] in TS_species])
    original_semi_TS_species = [species for species in original_semi_TS_species if species in species_to_check[0:max_TS_species_indexes+1] or species_expr_dict[species][OG_ID]=="YES_EXPR"]
    #######
    semi_TS_species = original_semi_TS_species
    #semi_TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
    if len(TS_species) > 0: #if there is at least one gene with tau >= up_cutoff and passing the expression cutoff
      if len(original_semi_TS_species) == 1:
        result_dict[OG_ID] = semi_TS_species[0]
      else: #if there are at least 2 genes with tau >= low_cutoff
        my_round = 0
        while len(semi_TS_species) >= 2: #while there are more than 2 genes with tau >= low_cutoff
          first_species = semi_TS_species[0]
          second_species = semi_TS_species[1]
          first_species_index = [idx for idx in range(len(species_to_check)) if species_to_check[idx] == first_species][0] #isolate index of species1 phylogenetic position
          second_species_index = [idx for idx in range(len(species_to_check)) if species_to_check[idx] == second_species][0] #isolate index of species2 phylogenetic position
          not_conserved_intermediate = len([species for species in species_to_check[first_species_index:second_species_index] if np.isnan(species_tau_dict[species][OG_ID])]) #Count the not-conserved species in between 
          if len(semi_TS_species) > 2: #if three are still more than 2 species:
            if second_species_index - first_species_index - not_conserved_intermediate <= 3: #if species1 and species2 are separated by 3 or less positions (quite phylogenetically close)
              semi_TS_species = semi_TS_species[1:] #remove first species and move on to the next iteration
              continue
            else: #if species1 and species2 are separated by more than 3 phylogenetic positions
              my_round = my_round+1 #update round to use in label
              if first_species == original_semi_TS_species[0]: #consider species1 as a tissue-specificity gain
                result_dict[OG_ID+"."+str(my_round)] = first_species
              else:
                result_dict[OG_ID+"."+str(my_round)] = node_name_dict[first_species]
              semi_TS_species = semi_TS_species[1:] #update the semi-ts species
              continue
          elif len(semi_TS_species) == 2: #if there are exactly 2 species left
            my_round = my_round+1
            if second_species_index - first_species_index - not_conserved_intermediate <= 3: #if species1 and species2 are separated by less than 3 phylogenetic positions, consider it as a unique tissue-specificity gain
              result_dict[OG_ID+"."+str(my_round)] = node_name_dict[second_species]
            else: #if they are more distant, consider them as two separated tissue-specificity gains
              if first_species == original_semi_TS_species[0]:
                result_dict[OG_ID+"."+str(my_round)] = first_species
              else:
                result_dict[OG_ID+"."+str(my_round)] = node_name_dict[first_species]
              my_round = my_round+1 #update round to use in label
              result_dict[OG_ID+"."+str(my_round)] = second_species 
            semi_TS_species = []
    final_OG_to_classify = [element for element in final_OG_to_classify if element != OG_ID]
  return([result_dict, final_OG_to_classify])
          
##########################
##### Main ###############
##########################

##### Read inputs ########
#orthogroups_df = pd.read_table(my_input, sep="\t", header=0, index_col=False)
orthogroups_df = pd.read_table(my_input, sep="\t", header=None, index_col=False, names=["OG_ID", "Species", "GeneID", "Tau", "Tissue", "Expr_cutoff"])
#I think for my purposes it's alright to replace NAs with zeros. Those genes I can't say if are tissue-specific or not
orthogroups_df = orthogroups_df.fillna(0)
#Header is: ["OG_ID", "Species", "Tau", "GeneID"]
#orthogroups_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/ts_gains_losses/All_version/Bilateria/Bilateria_conserved_orthogroups-COMPLETE-best_tau.txt", sep="\t", header=0, index_col=False)
grouped_df = orthogroups_df.groupby("OG_ID")

##### Initialize dictionaries ######
TS_num_deut_dict = {}
pot_TS_num_deut_dict = {} #potentially TS
TS_num_prot_dict = {}
pot_TS_num_prot_dict = {} #potentially TS
species_tau_dict = {}
species_tissue_dict = {}
species_expr_dict = {}
#Dictionaries with key=OG_ID, value=tissue_specific gain inference
OG_TS_classification_dict_deut = {}
OG_TS_classification_dict_prot = {}

#Initialize dictionary for all species
for species in list(set(list(orthogroups_df["Species"]))):
  species_tau_dict[species] = {}
  species_tissue_dict[species] = {}
  species_expr_dict[species] = {}

for OG_ID, group in grouped_df:
  #Deuterostoma branch
  deut_group = group.loc[group["Species"].isin(species_to_classify_deut)]
  TS_num = deut_group.loc[(deut_group["Tau"] >= up_tau_cutoff) & (deut_group["Tissue"].str.contains(query_tissue)) & (deut_group["Expr_cutoff"]=="YES_EXPR")].shape[0]
  pot_TS_num = deut_group.loc[(deut_group["Tau"] >= low_tau_cutoff) & (deut_group["Tau"] < up_tau_cutoff) & (deut_group["Tissue"].str.contains(query_tissue))].shape[0]
  #TS_num = len([element for element in list(deut_group["Tau"]) if element >= up_tau_cutoff]) #Count number of species with Tau >= up_cutoff
  #pot_TS_num = len([element for element in list(deut_group["Tau"]) if element >= low_tau_cutoff and element < up_tau_cutoff]) #Count number of species with Tau >= low_cutoff and < up_cutoff
  TS_num_deut_dict[OG_ID] = TS_num
  pot_TS_num_deut_dict[OG_ID] = pot_TS_num
  #Protostoma branch
  prot_group = group.loc[group["Species"].isin(species_to_classify_prot)]
  TS_num = prot_group.loc[(prot_group["Tau"] >= up_tau_cutoff) & (prot_group["Tissue"].str.contains(query_tissue)) & (prot_group["Expr_cutoff"]=="YES_EXPR")].shape[0]
  pot_TS_num = prot_group.loc[(prot_group["Tau"] >= low_tau_cutoff) & (prot_group["Tau"] < up_tau_cutoff) & (prot_group["Tissue"].str.contains(query_tissue))].shape[0]
  #TS_num = len([element for element in list(prot_group["Tau"]) if element >= up_tau_cutoff]) #Count number of species with Tau >= up_cutoff
  #pot_TS_num = len([element for element in list(prot_group["Tau"]) if element >= low_tau_cutoff and element < up_tau_cutoff]) #Count number of species with Tau >= low_cutoff and < up_cutoff
  TS_num_prot_dict[OG_ID] = TS_num
  pot_TS_num_prot_dict[OG_ID] = pot_TS_num
  #### Work on Tau dictionary: key=[species][orthogroup_ID], value=Tau
  for species in list(group["Species"]):
    species_tau_dict[species][OG_ID] = list(group.loc[group["Species"]==species]["Tau"])[0]
  #Add NA as an entry for the species where the gene is not conserved
  all_species = species_to_classify_deut + species_to_classify_prot
  not_conserved_species = [species for species in all_species if species not in list(group["Species"])]
  for species in not_conserved_species:
    species_tau_dict[species][OG_ID] = np.nan
  #### Work on Tissue dictionary: key=[species][orthogroup_ID], value=associated_tissue
  for species in list(group["Species"]):
    species_tissue_dict[species][OG_ID] = list(group.loc[group["Species"]==species]["Tissue"])[0]
    species_expr_dict[species][OG_ID] = list(group.loc[group["Species"]==species]["Expr_cutoff"])[0]
  for species in not_conserved_species:
    species_tissue_dict[species][OG_ID] = "None"
    species_expr_dict[species][OG_ID] = "NO_EXPR"

#Just select all the orthogroups
#OG_to_classify = list(TS_num_dict.keys())
OG_to_classify = list(set(list(orthogroups_df["OG_ID"])))


##### Work on OGs with one TS or no TS  #######
#select OGs where there is only one gene with Tau >= up_cutoff
#and there is no other gene with Tau >= low_cutoff
single_TS_OGs_deut = [OG_ID for OG_ID in OG_to_classify if TS_num_deut_dict[OG_ID]==1 and pot_TS_num_deut_dict[OG_ID]==0]
single_TS_OGs_prot = [OG_ID for OG_ID in OG_to_classify if TS_num_prot_dict[OG_ID]==1 and pot_TS_num_prot_dict[OG_ID]==0]
#Deuterostoma branch
for OG_ID in single_TS_OGs_deut:
  TS_species = list(orthogroups_df.loc[(orthogroups_df["OG_ID"]==OG_ID) & (orthogroups_df["Species"].isin(species_to_classify_deut)) & (orthogroups_df["Tau"] >= up_tau_cutoff) & (orthogroups_df["Tissue"].str.contains(query_tissue)) & (orthogroups_df["Expr_cutoff"]=="YES_EXPR")]["Species"])[0]
  OG_TS_classification_dict_deut[OG_ID] = TS_species

#Protostoma branch
for OG_ID in single_TS_OGs_prot:
  TS_species = list(orthogroups_df.loc[(orthogroups_df["OG_ID"]==OG_ID) & (orthogroups_df["Species"].isin(species_to_classify_prot)) & (orthogroups_df["Tau"] >= up_tau_cutoff) & (orthogroups_df["Tissue"].str.contains(query_tissue)) & (orthogroups_df["Expr_cutoff"]=="YES_EXPR")]["Species"])[0]
  OG_TS_classification_dict_prot[OG_ID] = TS_species


###### This will not happen in this script by definition
#Select OGs where none of the genes has Tau >= 0.75
#no_TS_OGs = [OG_ID for OG_ID in OG_to_classify if TS_num_deut_dict[OG_ID]==0] #there are no TS genes at all.
#Initialize the OGs to classify for the deuterostoma branch.
#single_TS_OG_deut = [OG_ID for OG_ID in OG_to_classify if TS_num_deut_dict[OG_ID] == 1 and pot_TS_num_deut_dict[OG_ID]==0]
#OG_to_classify_deut = [OG_ID for OG_ID in OG_to_classify if OG_ID not in no_TS_OGs and OG_ID not in single_TS_OGs_deut]
deuterostome_OGs = list(orthogroups_df.loc[(orthogroups_df["Species"].isin(species_to_classify_deut)) & (orthogroups_df["Tau"] >= up_tau_cutoff) & (orthogroups_df["Tissue"].str.contains(query_tissue)) & (orthogroups_df["Expr_cutoff"]=="YES_EXPR")]["OG_ID"])
OG_to_classify_deut = [OG_ID for OG_ID in OG_to_classify if OG_ID not in single_TS_OGs_deut]
OG_to_classify_deut = [OG_ID for OG_ID in OG_to_classify_deut if OG_ID in deuterostome_OGs]

print("Number of orthogroups on the deuterostome with more than one species with tissue-specificity")
print(len(OG_to_classify_deut))


####################################
####### FIRST BRANCH ###############
####################################
##### Work on the first pair #######
species1 = species_to_classify_deut[0]
species2 = species_to_classify_deut[1]
species_to_check = species_to_classify_deut

first_node_OGs = []
for OG_ID in OG_to_classify_deut:
  #Filter based on the species that are present in that orthogroup
  #species_to_check = [species for species in species_to_classify_deut if species in list(orthogroups_df.loc[orthogroups_df["OG_ID"]==OG_ID]["Species"])] #this should maintain the order
  #species1 = species_to_check[0]
  #species2 = species_to_check[1]
  #if at least one species has Tau >= up_cutoff and both speciese have tau >= low_cutoff:
  if (query_tissue in species_tissue_dict[species1][OG_ID] and query_tissue in species_tissue_dict[species2][OG_ID]): #Check that both species are associated with the query tissue
    if (species_tau_dict[species1][OG_ID] >= up_tau_cutoff and species_expr_dict[species1][OG_ID]=="YES_EXPR" and species_tau_dict[species2][OG_ID] >= low_tau_cutoff) or (species_tau_dict[species2][OG_ID] >= up_tau_cutoff and species_expr_dict[species2][OG_ID]=="YES_EXPR" and species_tau_dict[species1][OG_ID] >= low_tau_cutoff): #Check if the taus if these two species pass the cutoffs
      species_to_check = [species for species in species_to_check if species != species1 and species != species2]
      if all(species_tau_dict[species][OG_ID] < low_tau_cutoff for species in species_to_check): #if none of the other species has tau >= low_cutoff
        first_node_OGs = first_node_OGs + [OG_ID]
      else: #if some of the other species has tau >= low_cutoff, but all in other tissues:
        high_tau_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
        tissue_high_tau_species = [species for species in high_tau_species if query_tissue in species_tissue_dict[species][OG_ID]]
        if len(tissue_high_tau_species)==0:
          first_node_OGs = first_node_OGs + [OG_ID]

OG_to_classify_deut = [element for element in OG_to_classify_deut if element not in first_node_OGs]
for OG_ID in first_node_OGs:
  OG_TS_classification_dict_deut[OG_ID] = node_name_dict[species1+";"+species2]       

##### Work on the rest of the clade #######
species1 = [species1] #Transform to list
for species in species_to_classify_deut[2:]:
  node_name = node_name_dict[species] #isolate node name from relative dictionary
  species1 = species1 + [species2] #this comes from the previous iteration
  species2 = species 
  outer_species = [element for element in species_to_check if element != species2 and element not in species1] #select all species outside the boundary
  node_TS_OG_res = get_node_TS_OG(species1, species2, outer_species, OG_to_classify_deut)
  node_TS_OGs = node_TS_OG_res[0]
  OG_to_classify_deut = node_TS_OG_res[1]
  for OG_ID in node_TS_OGs:
    OG_TS_classification_dict_deut[OG_ID] = node_name

####################################
####### SECOND BRANCH ##############
####################################

#Reinitialize the orthogroups to classify.
OG_to_classify = list(set(list(orthogroups_df["OG_ID"])))
#select OGs where none of the genes has Tau >= 0.75
#no_TS_OGs = [OG_ID for OG_ID in OG_to_classify if TS_num_prot_dict[OG_ID]==0]
#Initialize the OGs to classify for the protostoma branch.
#single_TS_OG_prot = [OG_ID for OG_ID in OG_to_classify if TS_num_prot_dict[OG_ID]==1 and pot_TS_num_prot_dict[OG_ID]==0]
#OG_to_classify_prot = [OG_ID for OG_ID in OG_to_classify if OG_ID not in no_TS_OGs and OG_ID not in single_TS_OGs_prot]
#Remove the OGs that are conserved in protostomes but do not have any TS gene
protostome_OGs = list(orthogroups_df.loc[(orthogroups_df["Species"].isin(species_to_classify_prot)) & (orthogroups_df["Tau"] >= up_tau_cutoff) & (orthogroups_df["Tissue"].str.contains(query_tissue)) & (orthogroups_df["Expr_cutoff"]=="YES_EXPR")]["OG_ID"])
OG_to_classify_prot = [OG_ID for OG_ID in OG_to_classify if OG_ID not in single_TS_OGs_prot]
OG_to_classify_prot = [OG_ID for OG_ID in OG_to_classify_prot if OG_ID in protostome_OGs]

##### Work on the first pair #######
species1 = species_to_classify_prot[0]
species2 = species_to_classify_prot[1]
species_to_check = species_to_classify_prot

first_node_OGs = []
for OG_ID in OG_to_classify_prot:
  if (query_tissue in species_tissue_dict[species1][OG_ID] and query_tissue in species_tissue_dict[species2][OG_ID]): #Check that both species are associated with the query tissue
    if (species_tau_dict[species1][OG_ID] >= up_tau_cutoff and species_expr_dict[species1][OG_ID]=="YES_EXPR" and species_tau_dict[species2][OG_ID] >= low_tau_cutoff) or (species_tau_dict[species2][OG_ID] >= up_tau_cutoff and species_expr_dict[species2][OG_ID]=="YES_EXPR" and species_tau_dict[species1][OG_ID] >= low_tau_cutoff):
      species_to_check = [species for species in species_to_check if species != species1 and species != species2]
      if all(species_tau_dict[species][OG_ID] < low_tau_cutoff for species in species_to_check):
        first_node_OGs = first_node_OGs + [OG_ID]
      else: #if some of the other species has tau >= low_cutoff, but all in other tissues:
        high_tau_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
        tissue_high_tau_species = [species for species in high_tau_species if query_tissue in species_tissue_dict[species][OG_ID]]
        if len(tissue_high_tau_species)==0:
          first_node_OGs = first_node_OGs + [OG_ID]

OG_to_classify_prot = [element for element in OG_to_classify_prot if element not in first_node_OGs]
for OG_ID in first_node_OGs:
  OG_TS_classification_dict_prot[OG_ID] = node_name_dict[species1+";"+species2]       

##### Work on the rest of the clade #######
species1 = [species1]
for species in species_to_classify_prot[2:]:
  node_name = node_name_dict[species] #isolate node name from relative dictionary
  species1 = species1 + [species2] #this comes from the previous iteration
  species2 = species 
  outer_species = [element for element in species_to_check if element != species2 and element not in species1] #select all species outside the boundary
  node_TS_OG_res = get_node_TS_OG(species1, species2, outer_species, OG_to_classify_prot)
  node_TS_OGs = node_TS_OG_res[0]
  OG_to_classify_prot = node_TS_OG_res[1]
  for OG_ID in node_TS_OGs:
    OG_TS_classification_dict_prot[OG_ID] = node_name

####################################
####### BILATERIA CONS  ############
####################################
#Among the classified OG_IDs, check how many are Bilaterian:
ancestral_deut = [OG_ID for OG_ID in list(OG_TS_classification_dict_deut.keys()) if OG_TS_classification_dict_deut[OG_ID] in ancestral_clades]
ancestral_prot = [OG_ID for OG_ID in list(OG_TS_classification_dict_prot.keys()) if OG_TS_classification_dict_prot[OG_ID] in ancestral_clades]
bilaterian_OG_IDs = [OG_ID for OG_ID in ancestral_deut if OG_ID in ancestral_prot]
#Change the classification for the identified orthogroups
for OG_ID in bilaterian_OG_IDs:
  OG_TS_classification_dict_deut[OG_ID] = "Bilateria"
  OG_TS_classification_dict_prot[OG_ID] = "Bilateria"

classified_OGs = list(set(list(OG_TS_classification_dict_prot.keys()) + list(OG_TS_classification_dict_deut.keys())))
#OG_to_classify_deut = [OG_ID for OG_ID in OG_to_classify_deut if OG_ID not in no_TS_OGs and OG_ID not in classified_OGs]

#prova = dict(Counter(list(OG_TS_classification_dict.values())))
deut_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_deut, orient="index").rename(columns={0:"TS_gain"})
deut_final_df["OG_ID"] = list(deut_final_df.index)
deut_final_df = deut_final_df.reset_index(drop=True)

prot_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_prot, orient="index").rename(columns={0:"TS_gain"})
prot_final_df["OG_ID"] = list(prot_final_df.index)
prot_final_df = prot_final_df.reset_index(drop=True)

#The drop_duplicates is necessary because of the Bilateria classification
strict_final_df = pd.concat([deut_final_df, prot_final_df]).drop_duplicates()
strict_final_df["Category"] = "STRICT"


####################################
####### UNCLASSIFIED OGs  ##########
####################################
#For all the others, I will assign the gain to the last common ancestor of all species with Tau >= low_cutoff.
#This will be done by evaluating the two branches separately.
OG_TS_classification_dict_deut_LASSO = {}
OG_TS_classification_dict_prot_LASSO = {}

#LASSO classification on the deuterostome node
deut_node_TS_OG_LASSO_res = get_node_TS_OG_unclassified(species_to_classify_deut, OG_to_classify_deut, node_name_dict, OG_TS_classification_dict_deut_LASSO)
OG_TS_classification_dict_deut_LASSO = deut_node_TS_OG_LASSO_res[0]
OG_to_classify_deut = deut_node_TS_OG_LASSO_res[1] 

#LASSO classification on the protostome node
prot_node_TS_OG_LASSO_res = get_node_TS_OG_unclassified(species_to_classify_prot, OG_to_classify_prot, node_name_dict, OG_TS_classification_dict_prot_LASSO)
OG_TS_classification_dict_prot_LASSO = prot_node_TS_OG_LASSO_res[0]
OG_to_classify_prot = prot_node_TS_OG_LASSO_res[1]

#Check how many are missing
missing_common = [element for element in OG_to_classify_prot if element in OG_to_classify_deut]
if len(missing_common) == 0: #each OG has been classified on at least one branch.
  print("All orthogroups were classified")
else:
  print(str(len(missing_common))+" orthogroups failed to be classified")
  print(missing_common)

######## BILATERIAN inference ############
ancestral_deut = [re.sub("\\..", "", OG_ID) for OG_ID in list(OG_TS_classification_dict_deut_LASSO.keys()) if OG_TS_classification_dict_deut_LASSO[OG_ID] in ancestral_clades]
ancestral_prot = [re.sub("\\..", "", OG_ID) for OG_ID in list(OG_TS_classification_dict_prot_LASSO.keys()) if OG_TS_classification_dict_prot_LASSO[OG_ID] in ancestral_clades]
bilaterian_OG_IDs = [OG_ID for OG_ID in ancestral_deut if OG_ID in ancestral_prot]
#change the classification for the identified orthogroups
for OG_ID in bilaterian_OG_IDs:
  OG_TS_classification_dict_deut_LASSO[OG_ID] = "Bilateria"
  OG_TS_classification_dict_prot_LASSO[OG_ID] = "Bilateria"

#Remove the keys from the dictionary when the value is in ancestral clades and the key is in bilaterian_OG_IDs
for OG_ID in list(OG_TS_classification_dict_deut_LASSO.keys()):
  mod_ID = re.sub("\\..", "", OG_ID)
  if mod_ID in bilaterian_OG_IDs and "." in OG_ID:
    del OG_TS_classification_dict_deut_LASSO[OG_ID]
for OG_ID in list(OG_TS_classification_dict_prot_LASSO.keys()):
  mod_ID = re.sub("\\..", "", OG_ID)
  if mod_ID in bilaterian_OG_IDs and "." in OG_ID:
    del OG_TS_classification_dict_prot_LASSO[OG_ID]

##### Get final dataframes
deut_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_deut_LASSO, orient="index").rename(columns={0:"TS_gain"})
deut_final_df["OG_ID"] = list(deut_final_df.index)
deut_final_df = deut_final_df.reset_index(drop=True)

prot_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_prot_LASSO, orient="index").rename(columns={0:"TS_gain"})
prot_final_df["OG_ID"] = list(prot_final_df.index)
prot_final_df = prot_final_df.reset_index(drop=True)


lasso_final_df = pd.concat([deut_final_df, prot_final_df]).drop_duplicates()
lasso_final_df["Category"] = "LASSO"
#Added 18/07
lasso_final_df["OG_ID"] = [re.sub("\\..","", element) for element in list(lasso_final_df["OG_ID"])]

#Added 10/08
#Remove species-specific inferred gains in the LASSO when the inferred gain does not have tau >= 0.75 and passes the expr_cutoff
grouped_lasso_final_df = lasso_final_df.groupby(["OG_ID", "TS_gain"])
for name, group in grouped_lasso_final_df: #name = ("OG_ID", "TS_gain")
  if name[1] in species_to_classify_deut or name[1] in species_to_classify_prot: #If the gain inference is species-specific.
    if species_expr_dict[name[1]][name[0]] == "NO_EXPR" or species_tau_dict[name[1]][name[0]] < up_tau_cutoff: #If the entry does not pass the expression cutoff and/or does not have tau >= 0.75
      lasso_final_df = lasso_final_df.loc[~((lasso_final_df["OG_ID"]==name[0]) & (lasso_final_df["TS_gain"]==name[1]))] #Remove entry from lasso final dataframe

##########################################
##########################################
##########################################
#Merge the inferences from the two classifications
OG_TS_classification_dict_deut_MERGED = {**OG_TS_classification_dict_deut, **OG_TS_classification_dict_deut_LASSO}
OG_TS_classification_dict_prot_MERGED = {**OG_TS_classification_dict_prot, **OG_TS_classification_dict_prot_LASSO}
#Get the ancestral bilaterian that comes from the two classifications
ancestral_deut_MERGED = [re.sub("\\..", "", OG_ID) for OG_ID in list(OG_TS_classification_dict_deut_MERGED.keys()) if OG_TS_classification_dict_deut_MERGED[OG_ID] in ancestral_clades]
ancestral_prot_MERGED = [re.sub("\\..", "", OG_ID) for OG_ID in list(OG_TS_classification_dict_prot_MERGED.keys()) if OG_TS_classification_dict_prot_MERGED[OG_ID] in ancestral_clades]
bilaterian_OG_IDs_MERGED = [OG_ID for OG_ID in ancestral_deut_MERGED if OG_ID in ancestral_prot_MERGED]

#Remove the merged entries from the strict and lasso dataframes
#For the STRICT, just remove the entries
strict_final_df = strict_final_df.loc[~strict_final_df["OG_ID"].isin(bilaterian_OG_IDs_MERGED)]
#For the LASSO, remove the entries when the value is in ancestral_clades
print(lasso_final_df.loc[lasso_final_df["OG_ID"]=="GF_010889"])
lasso_final_df = lasso_final_df.loc[~((lasso_final_df["OG_ID"].isin(bilaterian_OG_IDs_MERGED)) & (lasso_final_df["TS_gain"].isin(ancestral_clades)))]
print(lasso_final_df.loc[lasso_final_df["OG_ID"]=="GF_010889"])
#Generate the merged dataframe
merged_final_df = pd.DataFrame({"OG_ID" : bilaterian_OG_IDs_MERGED, "TS_gain" : ["Bilateria"]*len(bilaterian_OG_IDs_MERGED), "Category" : ["MERGED"]*len(bilaterian_OG_IDs_MERGED)})


##########################################
############# JOIN #######################
##########################################
final_df = pd.concat([strict_final_df, lasso_final_df, merged_final_df])
final_df["OG_ID"] = [re.sub("\\..","", element) for element in list(final_df["OG_ID"])]
final_df = final_df[["OG_ID", "TS_gain", "Category"]]

###### Consider the Testis exception: If a gene is classified both as Euteleostomi and Protostoma in Testis, infer the bilaterian ancestor.
if query_tissue == "Testis":
  grouped_final_df = final_df.groupby("OG_ID")
  for OG_ID, group in grouped_final_df:
    inferred_gains = list(group["TS_gain"])
    if "Euteleostomi" in inferred_gains and "Protostoma" in inferred_gains:
      #Remove original entry (I have to consider the cases where there might be gains in other nodes.)
      final_df = final_df.loc[~((final_df["OG_ID"]==OG_ID) & (final_df["TS_gain"].isin(["Euteleostomi", "Protostoma"])))]
      #Add new entry
      final_df = pd.concat([final_df, pd.DataFrame({"OG_ID" : [OG_ID], "TS_gain" : ["Bilateria"], "Category" : ["MERGED_Testis"]})])
  #reorder dataframe
  final_df = final_df.sort_values(by="OG_ID")


####################################
####### SAVE TO FILE  ##############
####################################
#Remove .round from the first column
final_df.to_csv(str(my_output), sep="\t", header=True, index=False, na_rep="NA")
