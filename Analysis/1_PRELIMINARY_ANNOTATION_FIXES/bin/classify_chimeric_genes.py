#!/usr/bin/env python3

import argparse
import pandas as pd
import re

parser = argparse.ArgumentParser(description="Script to categorize the chimeric genes based on the overlapping fragments with each orthogroup they belong to")
parser.add_argument("--input", "-i", required=True, metavar="input", help="Path to output of get_chimeric_aligned_region.py")
parser.add_argument("--species_dict", "-s", required=True, metavar="species_dict", help="File with col1=geneID, col2=species")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###############################
###### Read arguments #########
###############################
args = parser.parse_args()
my_input_file = args.input
my_species_dict_file = args.species_dict
my_output_file = args.output

##############################
###### Define functions ######
##############################

def select_longest_fragments(group):
  group = group.dropna() #select only elements which do not have nan. There are always at least 2 
  group["first_aligned_aa"] = group["first_aligned_aa"].astype(int) #change data type
  group["last_aligned_aa"] = group["last_aligned_aa"].astype(int)
  group["overlap_length"] = group["last_aligned_aa"]-group["first_aligned_aa"]
  #order by overlap_length and select only the first two entries
  two_longest_entries = group.sort_values(by=["overlap_length"], ascending=False).iloc[[0,1],]
  two_longest_entries = two_longest_entries.drop(columns=["overlap_length"])
  return(two_longest_entries)

def process_chimeric_entries(group):
    if all([element != "NA" for element in list(group["first_aligned_ex_num"])]): #There are valid overlapping fragments for all the exons
      #Change data type
      group["first_aligned_aa"] = group["first_aligned_aa"].astype(int)
      group["last_aligned_aa"] = group["last_aligned_aa"].astype(int)
      #Order by first aligned aa
      group = group.sort_values(by=["first_aligned_aa","last_aligned_aa"]).reset_index() #order by first aligned aminoacid
      OG1 = group.iloc[0,:]; OG2 = group.iloc[1,:]
      OG1_range = range(int(OG1["first_aligned_ex_num"]), int(OG1["last_aligned_ex_num"])+1) #I need to add 1 so that the exon is actually considered
      OG2_range = range(int(OG2["first_aligned_ex_num"]), int(OG2["last_aligned_ex_num"])+1)
      ex_intersection = set(OG1_range).intersection(OG2_range)
      #If the intersection is equal to one of the ranges -> Class=COMPLETE_OVERLAP
      if ex_intersection == set(OG1_range) or ex_intersection == set(OG2_range):
        chimeric_class = "COMPLETE_OVERLAP"
      #If the intersection is of one element (exon)
      else:
        if len(ex_intersection) == 1:
          AA1_range = range(int(OG1["first_aligned_aa"]), int(OG1["last_aligned_aa"])+1)
          AA2_range = range(int(OG2["first_aligned_aa"]), int(OG2["last_aligned_aa"])+1)
          aa_intersection = set(AA1_range).intersection(AA2_range)
          #If there is intersection between the AA cooords -> Class=SAME_EXON-OVERLAP
          if len(aa_intersection) >=1:
            chimeric_class = "SAME_EXON-OVERLAP"
            #Else (no intersection between the AA coords) -> Class=SAME_EXON-NO_OVERLAP
          else:
            chimeric_class = "SAME_EXON-NO_OVERLAP"
        #If the intersection is of more elements (exons) -> Class=DIFFERENT_EXONS-OVERLAP
        elif len(ex_intersection) > 1:
          chimeric_class = "DIFFERENT_EXONS-OVERLAP"
        elif len(ex_intersection) == 0: #Else (there is no intersection between the two ranges)
          if max([int(ex) for ex in list(group["first_aligned_ex_num"])]) - min([int(ex) for ex in list(group["last_aligned_ex_num"])]) == 1: #If diff between ranges==1 -> Class=DIFFERENT_EXONS-CONTINUOS
            chimeric_class = "DIFFERENT_EXONS-CONTINUOS"
          else:
            chimeric_class = "DIFFERENT_EXONS-NO_CONTINUOS"
    else: #There are NAs instead of overlapping fragment coordinates
      OG1 = group.iloc[0,:]; OG2 = group.iloc[1,:]
      chimeric_class = "INVALID_ALN"
  
    #generate entry for the final dataframe
    group_entry = "%s\t%s\t%s;%s\t%s-%s;%s-%s\t%s-%s;%s-%s\t%s|%s;%s|%s\t%s\n" % (species, gene, OG1["orthogroup_id"], OG2["orthogroup_id"], OG1["first_aligned_ex_num"], OG1["last_aligned_ex_num"], OG2["first_aligned_ex_num"], OG2["last_aligned_ex_num"], str(OG1["first_aligned_aa"]), str(OG1["last_aligned_aa"]), str(OG2["first_aligned_aa"]), str(OG2["last_aligned_aa"]), OG1["first_aligned_ex"], OG1["last_aligned_ex"], OG2["first_aligned_ex"], OG2["last_aligned_ex"], chimeric_class)
    return(group_entry)

##############################
####### Upload inputs ########
##############################
#Input header: chimeric_geneID, orthogroup_id, first_aligned_aa, last_aligned_aa, first_aligned_ex, last_aligned_ex
input_df = pd.read_table(str(my_input_file), sep="\t", index_col=False, header=0)
#Create columns containing just the exon number
input_df["first_aligned_ex_num"] = [re.sub(":.*", "", element) for element in list(input_df["first_aligned_ex"])]
input_df["last_aligned_ex_num"] = [re.sub(":.*", "", element) for element in list(input_df["last_aligned_ex"])]
#Add species
geneID_species_df = pd.read_table(my_species_dict_file, sep="\t", index_col=False, header=None, names=["geneID", "species"])
geneID_species_dict = pd.Series(geneID_species_df.species.values, index=geneID_species_df.geneID).to_dict()
input_df["species"] = input_df["chimeric_geneID"].map(geneID_species_dict)

##############################
####### Main  ################
##############################
#open output filehandle
my_output = open(my_output_file, "w")
#Group dataframe by chimeric geneID
grouped_df = input_df.groupby("chimeric_geneID")
#print header
my_output.write("species\tchimeric_geneID\torthogroup_ids\tfirst-last_aligned_ex\tfirst-last_aligned_aa\tex:start|stop_coord\tchimeric_class\n")
for gene, group in grouped_df: #For each group
  #If the chimeric gene is divided only between 2 orthogroups:
  species = list(group["species"])[0] #isolate the species to save at the end
  if group.shape[0] == 2:
    group_entry = process_chimeric_entries(group)
  elif group.shape[0] > 2:
    group = select_longest_fragments(group)
    group_entry = process_chimeric_entries(group)   
  my_output.write(group_entry)
  
my_output.close()
