#!/usr/bin/env python3

import argparse
import json
import pandas as pd

parser = argparse.ArgumentParser(description="Script to infer tissue-specificity gains along the phylogeny based on Tau")
parser.add_argument("--all_orthogroups_input", "-a", required=True, metavar="all orthogroups input", help="Complete bilaterian conserved orthogroups with all associated tissue and expression cutoff information")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Complete bilaterian conserved orthogroups (at least one gene per species) where I select the highest Tau")
parser.add_argument("--ts_gains_input", "-t", required=True, metavar="ts gains input", help="Inferred tissue-specificity gains for the query tissue")
parser.add_argument("--tissue_differences_input", "-tdi", required=True, nargs="+", metavar="tissue differences input", help="File containing for each species and tissue the difference in relative expr between that tissue and the one with immediately lower expr")
parser.add_argument("--low_tau_cutoff", "-l", required=True, metavar="low tau cutoff", help="Strong cutoff to infer TS losses")
parser.add_argument("--up_tau_cutoff", "-u", required=True, metavar="up tau cutoff", help="Strong cutoff to infer TS gain")
parser.add_argument("--deut_species", "-d", required=True, metavar="deuterostome species", help="Species on the deuterostome side of the tree")
parser.add_argument("--prot_species", "-p", required=True, metavar="protostome species", help="Species on the protostome side of the tree")
parser.add_argument("--node_species_dict", "-ns", required=True, metavar="node species dict", help="Json-formatted dictionary with key=node name and value=list of species belonging to that node")
parser.add_argument("--node_name_dict", "-nn", required=True, metavar="node name dict", help="Json-formatted dictionary with key=earliest diverging species and value=node name")
parser.add_argument("--deut_nodes", "-dn", required=True, metavar="deuterostome nodes", help="Nodes on the deuterostome side of the tree")
parser.add_argument("--prot_nodes", "-pn", required=True, metavar="protostome nodes", help="Nodes on the protostome side of the tree")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file with losses")
parser.add_argument("--output_gains", "-og", required=True, metavar="output gains", help="Path to output file with updated_gains")
parser.add_argument("--tissue_diff_cutoff", "-tdc", required=True, metavar="tissue diff cutoff", help="Cutoff on the maximum differences in expression proportion between the first/second tissue to call the loss")
parser.add_argument("--query_tissue", "-q", required=True, metavar="query tissue", help="Tissue for which the tissue-specificity gains are evaluated")

##########################
##### Read arguments #####
##########################
args = parser.parse_args()
all_orthogroups_input = args.all_orthogroups_input
my_input = args.input
ts_gains_input = args.ts_gains_input
tissue_differences_input = args.tissue_differences_input
low_tau_cutoff = float(args.low_tau_cutoff)
up_tau_cutoff = float(args.up_tau_cutoff)
species_to_classify_deut = args.deut_species.split(",")
species_to_classify_prot = args.prot_species.split(",")
deut_nodes = args.deut_nodes.split(",")
prot_nodes = args.prot_nodes.split(",")
node_species_dict = json.loads(args.node_species_dict)
node_name_dict = json.loads(args.node_name_dict)
query_tissue = str(args.query_tissue)
tissue_diff_cutoff = args.tissue_diff_cutoff
my_output = args.output
my_output_gains = args.output_gains

all_species = species_to_classify_deut + species_to_classify_prot

#Update deut and prot lists and dictionaries
deut_nodes = deut_nodes + ["Bilateria_deut"]
prot_nodes = prot_nodes + ["Bilateria_prot"]
node_species_dict["Bilateria_deut"] = node_species_dict["Deuterostoma"]
node_species_dict["Bilateria_prot"] = node_species_dict["Protostoma"]
         
##########################
##### Main ###############
##########################

##### Read inputs ########
#Upload all orthogroups (this contains the information about the tissue association independently from the tissue cutoff we apply)
all_orthogroups_df = pd.read_table(all_orthogroups_input, sep="\t", header=None, index_col=False, names=["OG_ID", "Species", "GeneID", "Tau", "Tissue", "Paralog_tissue", "Expr_cutoff"])
#Generate a dictionary where each gene is associated with the relative paralog tissue
geneID_paralog_tissue_dict = pd.Series(all_orthogroups_df.Paralog_tissue.values, index=all_orthogroups_df.GeneID).to_dict()

#Upload orthogroups
orthogroups_df = pd.read_table(my_input, sep="\t", header=None, index_col=False, names=["OG_ID", "Species", "GeneID", "Tau", "Tissue", "Expr_cutoff"])
#Add Paralog tissue info
orthogroups_df["Paralog_tissue"] = orthogroups_df["GeneID"].map(geneID_paralog_tissue_dict)
#I think for my purposes it's alright to replace NAs with zeros. Those genes I can't say if are tissue-specific or not
orthogroups_df = orthogroups_df.fillna(0)

###########################
###########################
#Upload table with the difference in relative expression between the first and the third tissue, and the second and the third tissue for all species
tissue_differences_df = pd.DataFrame()
for my_file in tissue_differences_input:
  species_tissue_differences_df = pd.read_table(my_file, sep="\t", header=None, index_col=False, names=["GeneID", "Tissue", "Relative_expr_diff"])
  tissue_differences_df = pd.concat([tissue_differences_df, species_tissue_differences_df])
#Subset to only query tissue
query_tissues_differences_df = tissue_differences_df.loc[tissue_differences_df["Tissue"]==query_tissue]
query_tissues_differences_df["Translation_key"] = [element[0]+";"+element[1] for element in zip(list(query_tissues_differences_df["GeneID"]), list(query_tissues_differences_df["Tissue"]))]

#Get a dictionary with key=GeneID and value=Relative_expr_diff
geneID_tissue_differences_dict = pd.Series(query_tissues_differences_df.Relative_expr_diff.values, index=query_tissues_differences_df.Translation_key).to_dict()
#Add the information to the orthogroup table; if the info is not present, NA will be added
orthogroups_df["Translation_key"] = [element+";"+query_tissue for element in list(orthogroups_df["GeneID"])]
orthogroups_df["Tissue_diff"] = orthogroups_df["Translation_key"].map(geneID_tissue_differences_dict)

###########################
###########################

grouped_df = orthogroups_df.groupby("OG_ID")
#Upload tissue-specific gains
#Header is: [OG_ID, TS_gain, Category]
ts_gains_df = pd.read_table(ts_gains_input, sep="\t", header=0, index_col=False)
ts_gains_final_df = ts_gains_df.copy()

########################################################
######### Duplicate the Bilateria entries ##############
########################################################
bilaterian_deut_ts_gains = ts_gains_df.loc[ts_gains_df["TS_gain"]=="Bilateria"]
bilaterian_deut_ts_gains["TS_gain"] = ["Bilateria_deut"]*bilaterian_deut_ts_gains.shape[0]
bilaterian_prot_ts_gains = ts_gains_df.loc[ts_gains_df["TS_gain"]=="Bilateria"]
bilaterian_prot_ts_gains["TS_gain"] = ["Bilateria_prot"]*bilaterian_prot_ts_gains.shape[0]

ts_gains_df = pd.concat([ts_gains_df.loc[ts_gains_df["TS_gain"]!="Bilateria"], bilaterian_deut_ts_gains, bilaterian_prot_ts_gains])
########################################################



#Initialize output dataframe
all_losses_df = pd.DataFrame()

#Create a list of tuples (OG_ID, node_gain, branch) for those cases with multiple ancestral gains on the same branch
multiple_ancestral_gains_list = []
for OG_ID, group in grouped_df:
  #Isolate the all "problematic species" for each orthogroup
  OG_ID_cutoff_excluded_species = list(group.loc[group["Tau"] <= low_tau_cutoff]["Species"]) #All species with Tau <= 0.45
  #All species with Tau > 0.45 but where the query tissue is not among the top 2 expressed tissues
  OG_ID_tissues_excluded_species = list(group.loc[(group["Tau"] > low_tau_cutoff) & (~group["Paralog_tissue"].str.contains(query_tissue))]["Species"])
  #All species with Tau > 0.45, associated with the query tissue but where the difference in relative expr between the first/second and the second/third tissue is <= the predefined cutoff.
  OG_ID_tissue_differences_excluded_species = list(group.loc[(group["Tau"] > low_tau_cutoff) & (group["Paralog_tissue"].str.contains(query_tissue)) & (group["Tissue_diff"] <= float(tissue_diff_cutoff))]["Species"])
  #OG_excluded_species = list(group.loc[((group["Tau"] <= low_tau_cutoff)) | ((group["Tau"] > low_tau_cutoff) & (~group["Tissue"].str.contains(query_tissue)))]["Species"])
  OG_excluded_species = OG_ID_cutoff_excluded_species + OG_ID_tissues_excluded_species + OG_ID_tissue_differences_excluded_species  

  #Isolate the node/species with gains:
  OG_ts_gains = list(ts_gains_df.loc[ts_gains_df["OG_ID"]==OG_ID]["TS_gain"])
  #For all gains that are not in single species
  node_gains = [gain for gain in OG_ts_gains if gain not in all_species]
  deut_gains_num = len([gain for gain in OG_ts_gains if gain in deut_nodes])
  prot_gains_num = len([gain for gain in OG_ts_gains if gain in prot_nodes])
  deut_species_gains = [gain for gain in OG_ts_gains if gain in species_to_classify_deut]
  prot_species_gains = [gain for gain in OG_ts_gains if gain in species_to_classify_prot]
  
  if len(node_gains) >= 1: 
    for node_gain in node_gains:
      node_type = list(ts_gains_df.loc[(ts_gains_df["OG_ID"]==OG_ID) & (ts_gains_df["TS_gain"]==node_gain)]["Category"])[0]
      inner_species = [species for species in node_species_dict[node_gain] if species in list(group["Species"])]
      if node_type == "STRICT": #if the gain type is STRICT:
        inferred_losses = [species for species in inner_species if species in OG_excluded_species]
      elif node_type == "LASSO" or node_type == "MERGED" or node_type == "MERGED_Testis": #if the gain type is LASSO:
        if (node_gain in deut_nodes and deut_gains_num == 1) or (node_gain in prot_nodes and prot_gains_num == 1):
          inferred_losses = [species for species in inner_species if species in OG_excluded_species]
        else:
          branch = "deut" if node_gain in deut_nodes else "prot"
          multiple_ancestral_gains_list = multiple_ancestral_gains_list + [(OG_ID, node_gain, branch)]
          continue

      ######################
      #Adjust the gains
      if (node_gain in deut_nodes and len(deut_species_gains)==1) or (node_gain in prot_nodes and len(prot_species_gains)==1):
        if node_gain in deut_nodes:
          lower_species_gains = [species for species in deut_species_gains if species in inner_species] #Isolate internal species gains
        else:
          lower_species_gains = [species for species in prot_species_gains if species in inner_species]
        if len(lower_species_gains) >= 1:
          ts_gains_final_df = ts_gains_final_df.loc[~((ts_gains_final_df["OG_ID"]==OG_ID) & (ts_gains_final_df["TS_gain"].isin(lower_species_gains)))]
      #This is because the species will be added again
      ######################

      if len(inferred_losses) >= 2: #You need at least two losses overall to go through this
        #Get the last node of the last ordered common species
        subset_inner_species = inner_species[inner_species.index(inferred_losses[0]):] #The inferred losses should be ordered. Get the inner species from the first lost (in the tree)
        consecutive_lost_species = [species for species in inferred_losses if inferred_losses.index(species) == subset_inner_species.index(species)]
        if len(consecutive_lost_species) >= 2: #If there are at least two consecutive species, call the ancestor
          temp_ancestral_loss =  node_name_dict[consecutive_lost_species[-1]]
          temp_inferred_losses = [temp_ancestral_loss]  + [species for species in inferred_losses if species not in node_species_dict[temp_ancestral_loss]]

          #Select the TS species lower than the first lost species
          TS_species = list(group.loc[(group["Tau"] >= up_tau_cutoff) & (group["Tissue"].str.contains(query_tissue)) & (group["Expr_cutoff"]=="YES_EXPR")]["Species"])
          potential_new_gains = [species for species in inner_species if inner_species.index(species) < inner_species.index(inferred_losses[0]) and species in TS_species]
          #Select the LCA of all these species
          if len(potential_new_gains) > 0:
            new_gain_num = 1
            if len(potential_new_gains) == 1:
              new_gain = potential_new_gains[0]
            else:
              new_gain = node_name_dict[potential_new_gains[-1]]
          else:
            new_gain = "None"
            new_gain_num = 0
  
          #If the number of new inferences is lower than the number of old inferences, save ancestral losses and ADD new gain
          if new_gain_num + len(temp_inferred_losses) < len(inferred_losses):
            #Difference between the inferred ancestral loss and the original gain
            ancestor_difference = inner_species.index(inner_species[-1]) - inner_species.index(consecutive_lost_species[-1])
            if OG_ID == "GF_003234":
              print(node_gain)
              print(temp_inferred_losses)
              print(inner_species[-1])
              print(consecutive_lost_species[-1])
              print(ancestor_difference)
            if ancestor_difference > 1:
              inferred_losses = temp_inferred_losses
              if new_gain != "None":
                ts_gains_final_df = pd.concat([ts_gains_final_df, pd.DataFrame({"OG_ID": [OG_ID], "TS_gain" : [new_gain], "Category" : [node_type]})])
              #print("For orthogroup " + OG_ID + " the gains will change:")
              #print("Original gain: " + node_gain)
              #print("New losses: " + str(inferred_losses))
              #print("New gains: " + new_gain)  

      #Join to final dataframe
      all_losses_df = pd.concat([all_losses_df, pd.DataFrame({"OG_ID" : [OG_ID]*len(inferred_losses), "TS_losses" : inferred_losses, "Category" : [node_type]*len(inferred_losses)})])

##### Take care of the cases with multiple ancestral gains
#print("These are the orthogroups with multiple tissue-specificity gains on the same branch")
#print(multiple_ancestral_gains_list)
multiple_ancestral_gains_OGs = list(set([element[0] for element in multiple_ancestral_gains_list]))
#print(multiple_ancestral_gains_OGs)

for OG_ID in multiple_ancestral_gains_OGs:
  TS_gains = [element[1] for element in multiple_ancestral_gains_list if element[0]==OG_ID]
  branch = [element[2] for element in multiple_ancestral_gains_list if element[0]==OG_ID][0]
  if branch == "deut":
    my_species = species_to_classify_deut
    group = orthogroups_df.loc[(orthogroups_df["OG_ID"]==OG_ID) & (orthogroups_df["Species"].isin(species_to_classify_deut))] #Subset the orthogroups to deut species
    highest_gain = [gain for gain in TS_gains if deut_nodes.index(gain)==max([deut_nodes.index(element) for element in TS_gains])][0] #Select the highest gain
    lowest_gain = [gain for gain in TS_gains if deut_nodes.index(gain)==min([deut_nodes.index(element) for element in TS_gains])][0] #Select the lowest gain
  elif branch == "prot":
    my_species = species_to_classify_prot
    group = orthogroups_df.loc[(orthogroups_df["OG_ID"]==OG_ID) & (orthogroups_df["Species"].isin(species_to_classify_prot))] #Subset the orthogroups to prot species
    highest_gain = [gain for gain in TS_gains if prot_nodes.index(gain)==max([prot_nodes.index(element) for element in TS_gains])][0] #Select the highest gain
    lowest_gain = [gain for gain in TS_gains if prot_nodes.index(gain)==min([prot_nodes.index(element) for element in TS_gains])][0] #Select the lowest
  else:
    print("Something is wrong")

  #Comments for this code are above
  OG_ID_cutoff_excluded_species = list(group.loc[group["Tau"] <= low_tau_cutoff]["Species"])
  OG_ID_tissues_excluded_species = list(group.loc[(group["Tau"] > low_tau_cutoff) & (~group["Paralog_tissue"].str.contains(query_tissue))]["Species"])
  OG_ID_tissue_differences_excluded_species = list(group.loc[(group["Tau"] > low_tau_cutoff) & (group["Paralog_tissue"].str.contains(query_tissue)) & (group["Tissue_diff"] <= float(tissue_diff_cutoff))]["Species"])
  OG_excluded_species = OG_ID_cutoff_excluded_species + OG_ID_tissues_excluded_species + OG_ID_tissue_differences_excluded_species

  #Isolate the inner species, which depend on the highest inferred gain
  inner_species = [species for species in node_species_dict[highest_gain] if species in list(group["Species"])] #This should be an ordered list
  inferred_losses = [species for species in inner_species if species in OG_excluded_species]
  int_branch_inferred_losses = [species for species in inferred_losses if my_species.index(species) > my_species.index(node_species_dict[lowest_gain][-1])] #Select only the ones that are higher then the highest lowest ancestor (I realize this comment is not super clear)

  if len(int_branch_inferred_losses) >= 2:
    #get the lowest inferred loss
    lowest_inferred_loss_index = min([my_species.index(species) for species in inferred_losses])
    int_inner_species = [species for species in inner_species if my_species.index(species) >= lowest_inferred_loss_index]
    #This is the same code as above
    uncommon_index_list = [index for index in list(range(0,len(int_branch_inferred_losses))) if int_branch_inferred_losses[index] != int_inner_species[index]]
    if len(uncommon_index_list) > 0:
      last_lost_species = int_branch_inferred_losses[uncommon_index_list[0]-1]
    else:
      last_lost_species = int_branch_inferred_losses[-1]
    #I am not 100% sure about this part yet
    ancestral_inferred_loss = node_name_dict[last_lost_species]
    inferred_losses = [ancestral_inferred_loss] + [species for species in int_branch_inferred_losses if species not in node_species_dict[ancestral_inferred_loss]] + [species for species in inferred_losses if species not in int_branch_inferred_losses]

  #Join to final dataframe
  all_losses_df = pd.concat([all_losses_df, pd.DataFrame({"OG_ID" : [OG_ID]*len(inferred_losses), "TS_losses" : inferred_losses, "Category" : [node_type]*len(inferred_losses)})])  

##########################
##### Save to file #######
##########################
#Header: OG_ID, loss, gain type
all_losses_df.to_csv(my_output, sep="\t", header=True, index=False, na_rep="NA")
#Order final gains dataframe and save to file
ts_gains_final_df.sort_values(by=["OG_ID"]).drop_duplicates().to_csv(my_output_gains, sep="\t", header=True, index=False, na_rep="NA")
