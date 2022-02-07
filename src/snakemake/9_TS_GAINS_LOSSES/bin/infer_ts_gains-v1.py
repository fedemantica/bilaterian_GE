#!/usr/bin/env python3

import argparse
import pandas as pd
import math
import re
from collections import Counter

parser = argparse.ArgumentParser(description="Script to infer tissue-specificity gains along the phylogeny based on Tau")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Complete bilaterian conserved orthogroups (at least one gene per species) where I select the highest Tau")
parser.add_argument("--up_tau_cutoff", "-u", required=True, metavar="up tau cutoff", help="Stricter cutoff to infer TS gain")
parser.add_argument("--low_tau_cutoff", "-l", required=True, metavar="low tau cutoff", help="Milder cutoff to infer TS gain")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

##########################
##### Read arguments #####
##########################
args = parser.parse_args()
my_input = args.input
up_tau_cutoff = float(args.up_tau_cutoff)
low_tau_cutoff = float(args.low_tau_cutoff)
my_output = args.output
#up_tau_cutoff = 0.75
#low_tau_cutoff = 0.60

###########################
##### Define function #####
###########################

species_to_classify_deut = ["Hs2", "Mm2", "Bt2", "Mdo", "Gga", "Xtr", "Dre", "Cmi", "Bla", "Sp2"]
species_to_classify_prot = ["Dme", "Eba", "Aae", "BmA", "Tca", "Ame", "Bge", "Cdi", "Sma", "Obi"]
node_name_dict = {"Hs2":"Hs2", "Hs2;Mm2":"Euarchontoglires", "Mm2":"Euarchontoglires", "Bt2":"Eutheria", "Mdo":"Mammalia", "Gga":"Amniota", "Xtr":"Tetrapoda", "Dre":"Euteleostomi", "Cmi":"Vertebrata", "Bla":"Chordata", "Sp2":"Deuterostoma", "Dme":"Dme","Dme;Eba":"Cyclorrapha", "Eba":"Cyclorrapha", "Aae":"Diptera", "BmA":"Panorpidae", "Tca":"Oligoneoptera", "Ame":"Holometabola", "Bge":"Neoptera", "Cdi":"Insecta", "Sma":"Arthropoda", "Obi":"Protostoma"}
#ancestral_clades = ["Chordata", "Deuterostoma", "Vertebrata", "Arthropoda", "Protostoma", "Insecta"]
ancestral_clades = ["Chordata", "Deuterostoma", "Arthropoda", "Protostoma"]

def get_node_TS_OG(species1, species2, outer_species, OG_to_classify):
  #NB: species1 is a list
  #NB: species2 is a string
  node_OGs = []
  for OG_ID in OG_to_classify:
    if (species_tau_dict[species2][OG_ID] >= up_tau_cutoff):
      inner_species_TS = [species for species in species1 if species_tau_dict[species][OG_ID] >= up_tau_cutoff]
      outer_species_TS = [species for species in outer_species if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
      if len(inner_species_TS) >= math.floor(len(species1)/2): #round down the 50%
        if len(outer_species_TS) == 0:
          node_OGs = node_OGs + [OG_ID]
          OG_to_classify = [element for element in OG_to_classify if element != OG_ID] #remove the OG_IDs from the ones to classify
  return([node_OGs, OG_to_classify])   

def get_node_TS_OG_unclassified(species_to_check, OG_to_classify, node_name_dict, result_dict):
  for OG_ID in OG_to_classify:
    TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= up_tau_cutoff]
    original_semi_TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff] 
    semi_TS_species = original_semi_TS_species
    #semi_TS_species = [species for species in species_to_check if species_tau_dict[species][OG_ID] >= low_tau_cutoff]
    if len(TS_species) > 0: #if there is at least one gene with 0.75
      if len(original_semi_TS_species) == 1:
        result_dict[OG_ID] = semi_TS_species[0]
      else: #if there is at least another gene with 0.60
      #if len(semi_TS_species) > 1: #and there is at least another gene with 0.60
        my_round = 0
        while len(semi_TS_species) >= 2:
          first_species = semi_TS_species[0]
          second_species = semi_TS_species[1]
          first_species_index = [idx for idx in range(len(species_to_check)) if species_to_check[idx] == first_species][0]
          second_species_index = [idx for idx in range(len(species_to_check)) if species_to_check[idx] == second_species][0]

          if len(semi_TS_species) > 2:
            if second_species_index - first_species_index <= 3:
              semi_TS_species = semi_TS_species[1:]
              continue
            else:
              my_round = my_round+1
              if first_species == original_semi_TS_species[0]:
                result_dict[OG_ID+"."+str(my_round)] = first_species
              else:
                result_dict[OG_ID+"."+str(my_round)] = node_name_dict[first_species]
              semi_TS_species = semi_TS_species[1:]
              continue
          elif len(semi_TS_species) == 2:
            my_round = my_round+1
            if second_species_index - first_species_index <= 3:
              result_dict[OG_ID+"."+str(my_round)] = node_name_dict[second_species]
            else:
              if first_species == original_semi_TS_species[0]:
                result_dict[OG_ID+"."+str(my_round)] = first_species
              else:
                result_dict[OG_ID+"."+str(my_round)] = node_name_dict[first_species]
                my_round = my_round+1
                result_dict[OG_ID+"."+str(my_round)] = second_species 
            semi_TS_species = []
      OG_to_classify = [element for element in OG_to_classify if element != OG_ID]
  return([result_dict, OG_to_classify])
          
##########################
##### Main ###############
##########################

##### Read inputs ########
orthogroups_df = pd.read_table(my_input, sep="\t", header=0, index_col=False)
#Header is: ["OG_ID", "Species", "Tau", "GeneID"]
#orthogroups_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/ts_gains_losses/All_version/Bilateria/Bilateria_conserved_orthogroups-COMPLETE-best_tau.txt", sep="\t", header=0, index_col=False)

##### Initialize dictionaries ######
grouped_df = orthogroups_df.groupby("OG_ID")
TS_num_dict = {}
not_TS_num_dict = {}
species_tau_dict = {}
#OG_TS_classification_dict = {}
OG_TS_classification_dict_deut = {}
OG_TS_classification_dict_prot = {}

for species in list(set(list(orthogroups_df["Species"]))):
  species_tau_dict[species] = {}

for OG_ID, group in grouped_df:
  TS_num = len([element for element in list(group["Tau"]) if element >= up_tau_cutoff]) #Count number of species with Tau >= 0.75
  not_TS_num = len([element for element in list(group["Tau"]) if element < low_tau_cutoff]) #Count number of species with Tau <= 0.60
  TS_num_dict[OG_ID] = TS_num
  not_TS_num_dict[OG_ID] = not_TS_num
  for species in list(group["Species"]):
    species_tau_dict[species][OG_ID] = list(group.loc[group["Species"]==species]["Tau"])[0]

OG_to_classify = list(TS_num_dict.keys())

##### Work on OGs with one TS or no TS  #######
#select OGs where there is only one gene with Tau >= 0.75
single_TS_OGs = [OG_ID for OG_ID in OG_to_classify if TS_num_dict[OG_ID] == 1]
for OG_ID in single_TS_OGs: #add entries to the classification dictionary
  TS_species = list(orthogroups_df.loc[(orthogroups_df["OG_ID"]==OG_ID) & (orthogroups_df["Tau"] >= up_tau_cutoff)]["Species"])[0]
  OG_TS_classification_dict_deut[OG_ID] = TS_species
  OG_TS_classification_dict_prot[OG_ID] = TS_species

#select OGs where none of the genes has Tau >= 0.75
no_TS_OGs = [OG_ID for OG_ID in OG_to_classify if TS_num_dict[OG_ID] == 0]
OG_to_classify = [OG_ID for OG_ID in OG_to_classify if OG_ID not in single_TS_OGs and OG_ID not in no_TS_OGs]

####################################
####### FIRST BRANCH ###############
####################################
##### Work on the first pair #######
species1 = species_to_classify_deut[0]
species2 = species_to_classify_deut[1]
#species_to_check = list(orthogroups_df.loc[orthogroups_df["OG_ID"]==OG_ID]["Species"])
species_to_check = species_to_classify_deut

first_node_OGs = []
for OG_ID in OG_to_classify:
  if (species_tau_dict[species1][OG_ID] >= up_tau_cutoff and species_tau_dict[species2][OG_ID] >= low_tau_cutoff) or (species_tau_dict[species2][OG_ID] >= up_tau_cutoff and species_tau_dict[species1][OG_ID] >= low_tau_cutoff):
    species_to_check = [species for species in species_to_check if species != species1 and species != species2]
    if all(species_tau_dict[species][OG_ID] < low_tau_cutoff for species in species_to_check):
      first_node_OGs = first_node_OGs + [OG_ID]
      OG_to_classify = [element for element in OG_to_classify if element != OG_ID] #remove the OG_IDs from the ones to classify

for OG_ID in first_node_OGs:
  OG_TS_classification_dict_deut[OG_ID] = node_name_dict[species1+";"+species2]       

##### Work on the rest of the clade #######
species1 = [species1]
for species in species_to_classify_deut[2:]:
  node_name = node_name_dict[species] #isolate node name from relative dictionary
  species1 = species1 + [species2] #this comes from the previous iteration
  species2 = species 
  outer_species = [element for element in species_to_check if element != species and element not in species1] #select all species outside the boundary
  node_TS_OG_res = get_node_TS_OG(species1, species2, outer_species, OG_to_classify)
  node_TS_OGs = node_TS_OG_res[0]
  OG_to_classify = node_TS_OG_res[1]
  for OG_ID in node_TS_OGs:
    OG_TS_classification_dict_deut[OG_ID] = node_name

####################################
####### SECOND BRANCH ##############
####################################

#Reinitialize the orthogroups to classify.
OG_to_classify = list(TS_num_dict.keys())
OG_to_classify = [OG_ID for OG_ID in OG_to_classify if OG_ID not in single_TS_OGs and OG_ID not in no_TS_OGs]

##### Work on the first pair #######
species1 = species_to_classify_prot[0]
species2 = species_to_classify_prot[1]
#species_to_check = list(orthogroups_df.loc[orthogroups_df["OG_ID"]==OG_ID]["Species"])
species_to_check = species_to_classify_prot

first_node_OGs = []
for OG_ID in OG_to_classify:
  if (species_tau_dict[species1][OG_ID] >= up_tau_cutoff and species_tau_dict[species2][OG_ID] >= low_tau_cutoff) or (species_tau_dict[species2][OG_ID] >= up_tau_cutoff and species_tau_dict[species1][OG_ID] >= low_tau_cutoff):
    species_to_check = [species for species in species_to_check if species != species1 and species != species2]
    if all(species_tau_dict[species][OG_ID] < low_tau_cutoff for species in species_to_check):
      first_node_OGs = first_node_OGs + [OG_ID]
      OG_to_classify = [element for element in OG_to_classify if element != OG_ID] #remove the OG_IDs from the ones to classify

for OG_ID in first_node_OGs:
  OG_TS_classification_dict_prot[OG_ID] = node_name_dict[species1+";"+species2]       

##### Work on the rest of the clade #######
species1 = [species1]
for species in species_to_classify_prot[2:]:
  node_name = node_name_dict[species] #isolate node name from relative dictionary
  species1 = species1 + [species2] #this comes from the previous iteration
  species2 = species 
  outer_species = [element for element in species_to_check if element != species and element not in species1] #select all species outside the boundary
  node_TS_OG_res = get_node_TS_OG(species1, species2, outer_species, OG_to_classify)
  node_TS_OGs = node_TS_OG_res[0]
  OG_to_classify = node_TS_OG_res[1]
  for OG_ID in node_TS_OGs:
    OG_TS_classification_dict_prot[OG_ID] = node_name

####################################
####### BILATERIA CONS  ############
####################################
#Of the remaining OG_IDs with TS, check how many are Bilaterian:
ancestral_deut = [OG_ID for OG_ID in list(OG_TS_classification_dict_deut.keys()) if OG_TS_classification_dict_deut[OG_ID] in ancestral_clades]
ancestral_prot = [OG_ID for OG_ID in list(OG_TS_classification_dict_prot.keys()) if OG_TS_classification_dict_prot[OG_ID] in ancestral_clades]
bilaterian_OG_IDs = [OG_ID for OG_ID in ancestral_deut if OG_ID in ancestral_prot]
#change the classification for the identified orthogroups
for OG_ID in bilaterian_OG_IDs:
  OG_TS_classification_dict_deut[OG_ID] = "Bilateria"
  OG_TS_classification_dict_prot[OG_ID] = "Bilateria"

classified_OGs = list(set(list(OG_TS_classification_dict_prot.keys()) + list(OG_TS_classification_dict_deut.keys())))
#Still 455 without proper classification
OG_to_classify = [OG_ID for OG_ID in OG_to_classify if OG_ID not in no_TS_OGs and OG_ID not in classified_OGs]

#prova = dict(Counter(list(OG_TS_classification_dict.values())))
deut_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_deut, orient="index").rename(columns={0:"TS_gain"})
deut_final_df["OG_ID"] = list(deut_final_df.index)
deut_final_df = deut_final_df.reset_index(drop=True)

prot_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_prot, orient="index").rename(columns={0:"TS_gain"})
prot_final_df["OG_ID"] = list(prot_final_df.index)
prot_final_df = prot_final_df.reset_index(drop=True)

strict_final_df = pd.concat([deut_final_df, prot_final_df]).drop_duplicates()
strict_final_df["Category"] = "STRICT"

#res_deut = pd.DataFrame.from_dict(OG_TS_classification_dict_deut, orient="index")
#res_subset_deut = res_deut.loc[res_deut[0].isin(ancestral_clades)]
#res_prot = pd.DataFrame.from_dict(OG_TS_classification_dict_prot, orient="index")
#res_subset_prot = res_prot.loc[res_prot[0].isin(ancestral_clades)]
#res_OG_IDS_list = list(res_subset_deut.index) + list(res_subset_prot.index)
#res_counts_dict = dict(Counter(res_OG_IDS_list))
#bilaterian_OG_IDs = [OG_ID for OG_ID in list(res_counts_dict.keys()) if res_counts_dict[OG_ID]==2]
#for OG_ID in bilaterian_OG_IDs:
#  OG_TS_classification_dict_deut[OG_ID] = "Bilateria"
#  OG_TS_classification_dict_prot[OG_ID] = "Bilateria"


####################################
####### UNCLASSIFIED OGs  ##########
####################################

#For all the others, I will assign the gain to the last common ancestor of all species with Tau>=0.65.
#This will be done by evaluating the two branches separately.
OG_TS_classification_dict_deut_LASSO = {}
OG_TS_classification_dict_prot_LASSO = {}

#LASSO classification on the deuterostome node
deut_node_TS_OG_LASSO_res = get_node_TS_OG_unclassified(species_to_classify_deut, OG_to_classify, node_name_dict, OG_TS_classification_dict_deut_LASSO)
OG_TS_classification_dict_deut_LASSO = deut_node_TS_OG_LASSO_res[0]
OG_to_classify_deut = deut_node_TS_OG_LASSO_res[1] 

#LASSO classification on the protostome node
prot_node_TS_OG_LASSO_res = get_node_TS_OG_unclassified(species_to_classify_prot, OG_to_classify, node_name_dict, OG_TS_classification_dict_prot_LASSO)
OG_TS_classification_dict_prot_LASSO = prot_node_TS_OG_LASSO_res[0]
OG_to_classify_prot = prot_node_TS_OG_LASSO_res[1]

#Check how many are missing
missing_common = [element for element in OG_to_classify_prot if element in OG_to_classify_deut]
if len(missing_common) == 0: #each OG has been classified on at least one branch.
  print("All orthogroups were classified")
else:
  print(str(len(missing_common))+" orthogroups failed to be classified")

######## BILATERIAN inference ############
ancestral_deut = [OG_ID for OG_ID in list(OG_TS_classification_dict_deut_LASSO.keys()) if OG_TS_classification_dict_deut_LASSO[OG_ID] in ancestral_clades]
ancestral_prot = [OG_ID for OG_ID in list(OG_TS_classification_dict_prot_LASSO.keys()) if OG_TS_classification_dict_prot_LASSO[OG_ID] in ancestral_clades]
bilaterian_OG_IDs = [OG_ID for OG_ID in ancestral_deut if OG_ID in ancestral_prot]
#change the classification for the identified orthogroups
for OG_ID in bilaterian_OG_IDs:
  OG_TS_classification_dict_deut_LASSO[OG_ID] = "Bilateria"
  OG_TS_classification_dict_prot_LASSO[OG_ID] = "Bilateria"

##### Get final dataframes
deut_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_deut_LASSO, orient="index").rename(columns={0:"TS_gain"})
deut_final_df["OG_ID"] = list(deut_final_df.index)
deut_final_df = deut_final_df.reset_index(drop=True)

prot_final_df = pd.DataFrame.from_dict(OG_TS_classification_dict_prot_LASSO, orient="index").rename(columns={0:"TS_gain"})
prot_final_df["OG_ID"] = list(prot_final_df.index)
prot_final_df = prot_final_df.reset_index(drop=True)

lasso_final_df = pd.concat([deut_final_df, prot_final_df]).drop_duplicates()
lasso_final_df["Category"] = "LASSO"

####################################
####### JOIN and SAVE  #############
####################################
final_df = pd.concat([strict_final_df, lasso_final_df])
#Remove .round from the first column
final_df["OG_ID"] = [re.sub("\\..","", element) for element in list(final_df["OG_ID"])]
final_df = final_df[["OG_ID", "TS_gain", "Category"]]
final_df.to_csv(str(my_output), sep="\t", header=True, index=False, na_rep="NA")
