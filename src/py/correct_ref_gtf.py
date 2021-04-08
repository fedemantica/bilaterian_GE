#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import numpy as np
import math
import csv

parser = argparse.ArgumentParser(description="Script to correct broken and chimeric genes in the reference gtf")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species of interest")
parser.add_argument("--gtf", "-g", required=True, metavar="species", help="Reference gtf (only the reference transcript for each protein)")
parser.add_argument("--broken", "-b", required=True, metavar="suffix", help="One column file with semicolon separated groups of broken geneID to be joint")
parser.add_argument("--chimeric", "-c", required=True, metavar="length", help="Output of classify_chimeric_genes.py, with classification of the chimeric genes")
parser.add_argument("--IDs", "-i", required=True, metavar="input", help="File with col1=original gene IDs (of broken or chimeric genes) and col2=new geneID after correction")
parser.add_argument("--params_file", "-p", required=True, metavar="input", help="Params file with geneID infos for each species (suffix, length)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file for corrected GTF")
parser.add_argument("--output_brochi", "-ob", required=True, metavar="output_brochi", help="Path where to save only the subsetted gtf with the corrected broken and chimeric genes")
parser.add_argument("--output_unresolved", "-or", required=True, metavar="output_unresolved", help="Path to output file with unresolved chimeric genes")

###### Read arguments
args = parser.parse_args()
species = args.species
gtf_file = args.gtf
broken_genes_file = args.broken
chimeric_genes_file = args.chimeric
geneIDs_file = args.IDs
my_params_file = args.params_file
output_file = args.output
output_brochi = args.output_brochi
output_unresolved = args.output_unresolved

##################################
###### DEFINE FUNCTIONS ##########
##################################

#function to separate the gtf atttribute field in a list of tuples, preserving the original order of the entries
def separate_attributes(attribute_field): #input is a list of attribute entries (9th field of the GTF)
  final_list = []
  for field in attribute_field:
    entry_list = [[element.split(" ")[0], element.split(" ")[1]] for element in field.split("; ")] #remove the last element
    final_list = final_list+[entry_list]
  return(final_list)

def add_exon_number(attribute_field):
  final_list = []
  for field in attribute_field:
    complete_entry_list = [part for element in field.split(";")[:-1] for part in element.split(" ")]
    if "exon_number" in complete_entry_list:
      exon_number = int([re.sub(".*[ ]", "", re.sub('"', "", element)) for element in field.split(";") if "exon_number" in element][0])
    else:
      exon_number = np.nan
    final_list = final_list+[exon_number]
  return(final_list)

def add_gene_name(attribute_field):
  final_list = []
  for field in attribute_field:
    complete_entry_list = [part for element in field.split(";")[:-1] for part in element.split(" ")]
    if "gene_name" in complete_entry_list:
      gene_name = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in field.split(";") if "gene_name" in element][0]
    else:
      gene_name = "NoName"
    final_list = final_list+[gene_name]
  return(final_list)

def def_boundary_exon(chimeric_df): #input is a dataframe with header: species, chimeric_geneID, orthogroup_ids, first-last_aligned_ex, first-last_aligned_aa, ex:start|stop_coord, chimeric_class
  final_df = pd.DataFrame()
  #NB: the boundary exon belongs to the first gene_ID (<= kind of boundary)
  #Format of useful columns: 1-11;11-13      14-879;866-1275 1:1-35|11:793-887;11:793-887|13:974-1290
  ########### DIFFERENT_EXONS-CONTINUOS:
  my_df = chimeric_df[chimeric_df["chimeric_class"]=="DIFFERENT_EXONS-CONTINUOS"]
  my_df["boundary_ex_left"] = [int(element.split(";")[0].split("-")[1]) for element in list(my_df["first-last_aligned_ex"])]
  my_df["boundary_ex_right"] = my_df["boundary_ex_left"]+1
  final_df = pd.concat([final_df, my_df[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]])
  ########### SAME_EXON-NO_OVERLAP and SAME_EXON_OVERLAP: assign the exon to the geneID for which the overlap is higher
  my_df = chimeric_df[chimeric_df["chimeric_class"].isin(["SAME_EXON-NO_OVERLAP","SAME_EXON-OVERLAP"])]
  my_df["ex_aa_start"] = [int(element.split("|")[1].split(";")[0].split("-")[0].split(":")[1]) for element in my_df["ex:start|stop_coord"]]
  my_df["ex_aa_stop"] = [int(element.split("|")[1].split(";")[0].split("-")[1]) for element in my_df["ex:start|stop_coord"]] 
  my_df["aa_stop1"] = [int(element.split(";")[0].split("-")[1]) for element in my_df["first-last_aligned_aa"]]
  my_df["aa_start2"] = [int(element.split(";")[1].split("-")[0]) for element in my_df["first-last_aligned_aa"]]
  my_df["overlap1"] = my_df["aa_stop1"]-my_df["ex_aa_start"]
  my_df["overlap2"] = my_df["ex_aa_stop"]-my_df["aa_start2"]
  my_df["longest_overlap"] = ["1" if element[0] >= element[1] else "2" for element in list(zip(list(my_df["overlap1"]), list(my_df["overlap2"])))]
  #If the overlap with the first gene is higher, set the boundary exon, otherwise subtract one.
  my_df["boundary_ex_left"] = [int(element[0].split(";")[0].split("-")[1]) if element[1]=="1" else int(element[0].split(";")[0].split("-")[1])-1 for element in zip(list(my_df["first-last_aligned_ex"]), list(my_df["longest_overlap"]))]
  my_df["boundary_ex_right"] = my_df["boundary_ex_left"]+1
  #Add to final dataframe
  final_df = pd.concat([final_df, my_df[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]])  
  ########### DIFFERENT_EXONS-NO_CONTINUOS: here split in correspondence of the last exon which has coverage
  my_df = chimeric_df[chimeric_df["chimeric_class"]=="DIFFERENT_EXONS-NO_CONTINUOS"]
  my_df["boundary_ex_left"] = [element.split(";")[0].split("-")[1] for element in list(my_df["first-last_aligned_ex"])]
  my_df["boundary_ex_right"] = [element.split(";")[1].split("-")[0] for element in list(my_df["first-last_aligned_ex"])]
  final_df = pd.concat([final_df, my_df[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]])
  ########### DIFFERENT_EXONS-OVERLAP: here in case the overlap is only 2 exons, we assign each overlapping to one of the genes
  my_df = chimeric_df[chimeric_df["chimeric_class"]=="DIFFERENT_EXONS-OVERLAP"]
  my_df["ex_stop1"] = [int(element.split(";")[0].split("-")[1]) for element in list(my_df["first-last_aligned_ex"])]
  my_df["ex_start2"] = [int(element.split(";")[1].split("-")[0]) for element in list(my_df["first-last_aligned_ex"])]
  my_df["ex_overlap"] = my_df["ex_stop1"]-my_df["ex_start2"]
  #select only the cases where the overlap is over one exon
  my_df_selected = my_df[my_df["ex_overlap"]==1]
  my_df_selected["boundary_ex_left"] = my_df_selected["ex_start2"]  #simply invert the boundaries
  my_df_selected["boundary_ex_right"] = my_df_selected["ex_stop1"]
  final_df = pd.concat([final_df, my_df_selected[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]]) 
  #separately save the DIFFERENT_EXONS-OVERLAP which can't be properly classified
  my_df_unselected = my_df[my_df["ex_overlap"]!=1]
  my_df_unselected = my_df_unselected.drop(columns=["ex_stop1", "ex_start2"]) #remove extra columns
  ######### COMPLETE_OVERLAP: this will not be corrected, I will act at the level of gene orthogroups.
  my_df = chimeric_df[chimeric_df["chimeric_class"]=="COMPLETE_OVERLAP"]
  my_df["ex_overlap"] = "complete"
  my_df_unselected = pd.concat([my_df_unselected, my_df])
  ######### INVALID_ALN
  my_df = chimeric_df[chimeric_df["chimeric_class"]=="INVALID_ALN"]
  my_df["ex_overlap"] = "invalid"
  my_df_unselected = pd.concat([my_df_unselected, my_df]) 
  return(final_df, my_df_unselected)

def modify_value_in_tuple(attribute_field, category, new_value):
  final_list = []
  for field in attribute_field:
    index = attribute_field.index(field)
    for part in field:
      if part[0] == category:
        sub_index = field.index(part)
        #field[sub_index][1] = new_value[index]
        field[sub_index][1] = f'"{new_value[index]}"' #this prints the string with surrounding quotes
    final_list = final_list+[field]
  return(final_list)

def add_entries_first_gene(first_gene_df, group, first_ex, last_ex, transcript_suffix):
  #gene entry and transcript entry: correct the right boundaries
  if "gene" in list(set(list(group["type"]))): #if there is a gene entry
    gene_entry = group[group["type"]=="gene"]
    gene_entry["new_geneID"] = [element.split(";")[0] for element in list(gene_entry["new_geneID"])] #select first geneID
    if str(list(gene_entry["strand"])[0]) == "+":
      gene_entry["stop"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["stop"])[0] #set stop from last exon
    elif str(list(gene_entry["strand"])[0]) == "-":
      gene_entry["start"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["start"])[0] #set start from first exon
    first_gene_df = pd.concat([first_gene_df, gene_entry])
  if "transcript" in list(set(list(group["type"]))): #if there is a transcript entry
    transcript_entry = group[group["type"]=="transcript"]
    transcript_entry["new_geneID"] = [element.split(";")[0] for element in list(transcript_entry["new_geneID"])] #select first geneID
    transcript_entry["new_transcriptID"] = transcript_entry["new_geneID"]+transcript_suffix
    #transcript_entry["new_transcriptID"] = [element.split(";")[0]+transcript_suffix for element in list(transcript_entry["new_geneID"])]
    transcript_entry["start"] = int(list(gene_entry["start"])[0])
    transcript_entry["stop"] = int(list(gene_entry["stop"])[0])
    first_gene_df = pd.concat([first_gene_df, transcript_entry])
    #if str(list(transcript_entry["strand"])[0]) == "+":
     # transcript_entry["stop"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["stop"])[0] #set stop from last exon
    #elif str(list(transcript_entry["strand"])[0]) == "-":
      #transcript_entry["start"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["start"])[0] #set stop from first exon
    #first_gene_df = pd.concat([first_gene_df, transcript_entry])          
  #The start codon has an exon number, so it is not necessary to add it.
  #if there is a start codon in the group, add it to the first_gene_df
  #if "start_codon" in list(set(list(group["type"]))):
   # start_codon_entry = group[group["type"]=="start_codon"]
    #first_gene_df = pd.concat([first_gene_df, start_codon_entry])
  return(first_gene_df)
 
def add_entries_second_gene(second_gene_df, group, first_ex, last_ex, transcript_suffix):
  #add gene entry
  if "gene" in list(set(list(group["type"]))): #if there is a gene entry
     gene_entry = group[group["type"]=="gene"]
     gene_entry["new_geneID"] = [element.split(";")[1] for element in list(gene_entry["new_geneID"])] #select the second geneID
     if str(list(gene_entry["strand"])[0]) == "+":
       gene_entry["start"] = list(second_gene_df[(second_gene_df["exon_number"]==first_ex) & (second_gene_df["type"]=="exon")]["start"])[0] #set start from first ex
       gene_entry["stop"] = list(second_gene_df[(second_gene_df["exon_number"]==last_ex) & (second_gene_df["type"]=="exon")]["stop"])[0] #set stop from last ex
     elif str(list(gene_entry["strand"])[0]) == "-":
        gene_entry["start"] = int(list(second_gene_df[(second_gene_df["exon_number"]==last_ex) & (second_gene_df["type"]=="exon")]["start"])[0]) #set start from last ex
        gene_entry["stop"] = int(list(second_gene_df[(second_gene_df["exon_number"]==first_ex) & (second_gene_df["type"]=="exon")]["stop"])[0]) #set stop from first ex
     second_gene_df = pd.concat([second_gene_df, gene_entry])
  #add transcript entry
  if "transcript" in list(set(list(group["type"]))):
    transcript_entry = group[group["type"]=="transcript"]
    transcript_entry["new_geneID"] = [element.split(";")[1] for element in list(transcript_entry["new_geneID"])] #select second geneID
    transcript_entry["new_transcriptID"] = transcript_entry["new_geneID"]+transcript_suffix
    #transcript_entry["new_transcriptID"] = [element.split(";")[1]+transcript_suffix for element in list(transcript_entry["new_geneID"])]
    transcript_entry["start"] = int(list(gene_entry["start"])[0])
    transcript_entry["stop"] = int(list(gene_entry["stop"])[0])
    second_gene_df = pd.concat([second_gene_df, transcript_entry])
  #Fix here the start of the second splitted gene in case the phase is not 0.
  return(second_gene_df)

def renumerate_exons(broken_exons_df, broken_parts, last_ex_tmp, first_broken_exons_df):
  final_df = first_broken_exons_df
  for part in broken_parts[1:]: #access all the broken genes after the first
    second_broken_exons_df = broken_exons_df[broken_exons_df["geneID"]==part] #subset by first gene
    ex_number_list = list(sorted(list(set(list(second_broken_exons_df["exon_number"])))))
    new_exon_number = list(range(last_ex_tmp+1, last_ex_tmp+1+len(ex_number_list)))
    ex_number_dict = {float(element) : new_exon_number[ex_number_list.index(element)] for element in ex_number_list}
    second_broken_exons_df["exon_number"] = second_broken_exons_df["exon_number"].map(ex_number_dict)
    last_ex_tmp = int(max(list(second_broken_exons_df["exon_number"])))
    final_df = pd.concat([final_df, second_broken_exons_df])
  return(final_df)

def add_entries_broken_genes(broken_exons_df, group, first_ex, last_ex):
  if "gene" in list(set(list(group["type"]))):
    gene_entry = group[group["type"]=="gene"].head(1) #select the first gene entry, but replace both start and stop just in case.
    if str(list(gene_entry["strand"])[0]) == "+":
      gene_entry["start"] = list(broken_exons_df[(broken_exons_df["exon_number"]==first_ex) & (broken_exons_df["type"]=="exon")]["start"])[0] #set start from first ex
      gene_entry["stop"] = list(broken_exons_df[(broken_exons_df["exon_number"]==last_ex) & (broken_exons_df["type"]=="exon")]["stop"])[0] #set stop from last ex
    elif str(list(gene_entry["strand"])[0]) == "-": #always select the min start and the max stop among all exon  entries
      gene_entry["start"] = min(list(broken_exons_df[broken_exons_df["type"]=="exon"]["start"]))
      gene_entry["stop"] = max(list(broken_exons_df[broken_exons_df["type"]=="exon"]["stop"]))
      #gene_entry["start"] = int(list(broken_exons_df[(broken_exons_df["exon_number"]==last_ex) & (broken_exons_df["type"]=="exon")]["start"])[0]) #set start from last ex
      #gene_entry["stop"] = int(list(broken_exons_df[(broken_exons_df["exon_number"]==first_ex) & (broken_exons_df["type"]=="exon")]["stop"])[0])      
    broken_exons_df = pd.concat([broken_exons_df, gene_entry])
  if "transcript" in list(set(list(group["type"]))):
    transcript_entry = group[group["type"]=="transcript"].head(1) #select the first transcript entry, but replace both start and stop just in case.
    transcript_entry["start"] = int(list(gene_entry["start"])[0])
    transcript_entry["stop"] = int(list(gene_entry["stop"])[0])
    broken_exons_df = pd.concat([broken_exons_df, transcript_entry])
  return(broken_exons_df) 

def rebuild_attribute_entry(mod_attribute_field):
  final_list = []
  for field in mod_attribute_field:
    new_entry = '; '.join([' '.join(element) for element in field])
    final_list = final_list + [new_entry]
  return(final_list)

##################################
###### READ INPUTS ###############
##################################

gtf_df = pd.read_table(gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
#gtf_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/DB/gtf/ref/BmA_annot-B.gtf", sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase","attribute"])
broken_genes_list = list(pd.read_table(broken_genes_file, sep="\t", index_col=False, header=None, names=["broken_genes"])["broken_genes"])
broken_gene_flatten_list = [part for element in broken_genes_list for part in element.split(";")]
#broken_genes_list = list(pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/broken_genes/BmA-selected_broken_genes.tab", sep="\t", index_col=False, header=None, names=["broken_genes"])["broken_genes"])
#Header: species, chimeric_geneID, orthogroup_ids, first-last_aligned_ex, first-last_aligned_aa, ex:start|stop_coord, chimeric_class
chimeric_genes_df = pd.read_table(chimeric_genes_file, sep="\t", index_col=False, header=0)
chimeric_genes_df = chimeric_genes_df[chimeric_genes_df["species"]==species] #subset for the species of interest
chimeric_genes_list = list(chimeric_genes_df["chimeric_geneID"])
#chimeric_genes_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/chimeric_proteins/classified_chimeric_genes.tab", sep="\t", index_col=False, header=0)
#Header: geneID, new_IDs, category
geneIDs_df = pd.read_table(geneIDs_file, sep="\t", index_col=False, header=0)
geneIDs_dict = pd.Series(geneIDs_df.new_IDs.values, index=geneIDs_df.geneID).to_dict()
reverse_geneID_dict = pd.Series(geneIDs_df.geneID.values, index=geneIDs_df.new_IDs).to_dict()
#correct broken genes entries
broken_genes_keys = [key for key in list(geneIDs_dict.keys()) if ";" in key]
for key in broken_genes_keys:
  for num in list(range(0, len(key.split(";")))):
    geneIDs_dict[key.split(";")[num]] = geneIDs_dict[key]
#geneIDs_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/corrected_gtfs/BmA_new_geneIDs.txt", sep="\t", index_col=False, header=0)

#Header: species, geneID_prefix, geneID_length, transcript_suffix, protein_suffix
params_df = pd.read_table(my_params_file, sep="\t", header=0, index_col=False)
transcript_suffix = str(list(params_df[params_df["species"]==species]["transcript_suffix"])[0])
protein_suffix = str(list(params_df[params_df["species"]==species]["protein_suffix"])[0])
#params_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/src/snakemake/1_PRELIMINARY_ANNOTATION_FIXES/geneID_params.txt", sep="\t", header=0, index_col=False)

print("%s: %d" % ("broken_genes", len(broken_gene_flatten_list))) #for debugging
print("%s: %d" % ("chimeric_genes", len(chimeric_genes_list))) #for debugging

##################################
####### HEALTHY GENES ############
##################################
#add geneID
gtf_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part]
brochi_genes_list = broken_gene_flatten_list + chimeric_genes_list
#filter out entries of brochi genes from the GTF
healthy_GTF_df = gtf_df[~(gtf_df["geneID"].isin(brochi_genes_list))] 
healthy_GTF_df = healthy_GTF_df.rename(columns={"attribute" : "attribute_mod"}) #rename field to match the broken and chimeric gtf dataframes
healthy_GTF_df = healthy_GTF_df[["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute_mod"]]

##################################
####### BROKEN GENES #############
##################################

#subset gtf to broken genes.
broken_GTF_df = gtf_df.loc[gtf_df["geneID"].isin(broken_gene_flatten_list)]
#add exon number column and transform the attribute field in a list of tuples
broken_GTF_df["exon_number"] = add_exon_number(list(broken_GTF_df["attribute"]))
broken_GTF_df["exon_number"] = broken_GTF_df["exon_number"].astype("Int64") #This is necessary to have integers with NaN
broken_GTF_df["gene_name"] = add_gene_name(list(broken_GTF_df["attribute"]))
broken_GTF_df["attribute_mod"] = separate_attributes(list(broken_GTF_df["attribute"])) 
#add new geneID column (use the dict)
broken_GTF_df["new_geneID"] = broken_GTF_df["geneID"].map(geneIDs_dict)

all_broken_gtf_df = pd.DataFrame() #initialize gtf_df for all broken genes
#groupby the new geneID column and cycle on the groups.
grouped_broken_GTF_df = broken_GTF_df.groupby("new_geneID")
for repaired_gene, group in grouped_broken_GTF_df:
  print(repaired_gene) #for debugging
  #group = broken_GTF_df[broken_GTF_df["new_geneID"]=="BGIBMGAB00011"]
  #derive variables
  broken_parts = reverse_geneID_dict[repaired_gene].split(";")
  #if list(group["strand"])[0] == "-": #revert the order of the broken parts if the strand is negative
    #broken_parts = broken_parts[::-1]
  new_gene_name = ';'.join([name for name in list(set(list(group["gene_name"]))) if name != "NoName"])
  #add transcriptID and proteinID columns to group.
  group["new_transcriptID"] = group["new_geneID"]+transcript_suffix
  group["new_proteinID"] = group["new_geneID"]+protein_suffix
  #broken_df_exons; select only entries with exons.
  broken_exons_df = group.dropna(subset=["exon_number"])

  #make sure exnos are ordered (by start and stop coords), and re-number them (the second gene will change).
  broken_exons_df = broken_exons_df.sort_values(by=["start", "stop"])
  first_broken_exons_df = broken_exons_df[broken_exons_df["geneID"]==broken_parts[0]] #subset by the first gene. I hope they are always ordered
  last_ex_tmp = int(max(list(first_broken_exons_df["exon_number"])))

  #renumber the exons (both exons and CDS)
  broken_exons_df = renumerate_exons(broken_exons_df, broken_parts, last_ex_tmp, first_broken_exons_df)
  first_ex = int(min(list(broken_exons_df["exon_number"]))) #this should be one in the majority of cases, but who knows
  last_ex = int(max(list(broken_exons_df["exon_number"])))
  #If there were start codons entries -> take only the start codon corresponding to new exon 1.
  #broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="start_codon") & (broken_exons_df["exon_number"] > 1))]
  #If there were stop codons entries -> take only the stop codon corresponding to the new last exon.
  if list(broken_exons_df["strand"])[0] == "+":
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="start_codon") & (broken_exons_df["exon_number"] > 1))] #only the start codon corresponding to new exon 1
    removed_stop_codons =  broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon") & (broken_exons_df["exon_number"] < last_ex)] #remove all the stop codons apart from the last
    #If there are removed stop codons, I need to fix the coordinates of the penultimate exon: not sure this should stay all through.
    if removed_stop_codons.shape[0] >= 1:
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="stop_codon") & (broken_exons_df["exon_number"] < last_ex))]
      last_stop_coords = list(removed_stop_codons["stop"])
      broken_exons_df["stop"] = [element if element not in last_stop_coords else element-3 for element in list(broken_exons_df["stop"])]
  elif list(broken_exons_df["strand"])[0] == "-":
    #only the start codon corresponding to the second gene
    removed_start_codons = broken_exons_df.loc[(broken_exons_df["type"]=="start_codon")].sort_values(by="start", ascending=False).iloc[1:]
    removed_stop_codons = broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon")].sort_values(by="start").iloc[1:] #remove all the stop codons apart from the first
    if removed_start_codons.shape[0] >= 1:
      removed_start_coord_start = list(removed_start_codons["start"])
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="start_codon") & (broken_exons_df["start"].isin(removed_start_coord_start)))]
    if removed_stop_codons.shape[0] >= 1:
      removed_stop_coord_start = list(removed_stop_codons["start"])
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="stop_codon") & (broken_exons_df["start"].isin(removed_stop_coord_start)))]
      last_start_coords = list(removed_stop_codons["start"])
      broken_exons_df["start"] = [element if element not in last_start_coords else element-3 for element in list(broken_exons_df["start"])]
  #if originally there were gene and transcript entries, take them and modify the start and stop.
  broken_exons_df = add_entries_broken_genes(broken_exons_df, group, first_ex, last_ex)

  #update the geneID, transcriptID and proteinID in loco
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "gene_id", list(broken_exons_df["new_geneID"]))
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "transcript_id", list(broken_exons_df["new_transcriptID"]))
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "protein_id", list(broken_exons_df["new_proteinID"]))
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "exon_number", list(broken_exons_df["exon_number"]))
  #modify source and gene name (new gene name=broken_name1;broken_name2...)
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "gene_source", ["brochi_pipe" for element in list(range(broken_exons_df.shape[0]))])
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "transcript_source", ["brochi_pipe" for element in list(range(broken_exons_df.shape[0]))])
  broken_exons_df["attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "gene_name", [new_gene_name for element in list(range(broken_exons_df.shape[0]))])
  #Rebuild the attribute field in the original format
  broken_exons_df["attribute_mod"] = rebuild_attribute_entry(list(broken_exons_df["attribute_mod"]))
  final_broken_exons_df = broken_exons_df[["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute_mod"]]
  all_broken_gtf_df = pd.concat([all_broken_gtf_df, final_broken_exons_df])

##################################
####### CHIMERIC GENES ###########
##################################

###Preprocessing: get exon boundaries for all categories
#chimeric_genes_to_correct_df = chimeric_genes_df[~(chimeric_genes_df["chimeric_class"]=="COMPLETE_OVERLAP")]
#print(chimeric_genes_to_correct_df)
#print(chimeric_genes_to_correct_df.shape)
#chimeric_genes_to_correct = list(chimeric_genes_to_correct_df["chimeric_geneID"])
chimeric_genes_to_correct = list(chimeric_genes_df["chimeric_geneID"])
#generate a boundary exons df
res = def_boundary_exon(chimeric_genes_df)
boundary_ex_df = res[0]
unresolved_chimeric_df = res[1]
#Save unresolved_df to file
unresolved_chimeric_df.to_csv(output_unresolved, sep="\t", index=False, header=True, na_rep="NA")

#create dictionary with correspondence between geneID - left boundary and geneID - right boundary
geneID_left_bound_dict = pd.Series(boundary_ex_df.boundary_ex_left.values, index=boundary_ex_df.chimeric_geneID).to_dict()
geneID_right_bound_dict = pd.Series(boundary_ex_df.boundary_ex_right.values, index=boundary_ex_df.chimeric_geneID).to_dict()

###Correct GTF
#subset GTF to chimeric genes to correct
chimeric_genes_to_correct = [element for element in chimeric_genes_to_correct if element not in list(unresolved_chimeric_df["chimeric_geneID"])] #remove "unersolved" genes (these are saved to file separately)
chimeric_genes_to_correct = [element for element in chimeric_genes_to_correct if element not in broken_gene_flatten_list] #remove those fucking cases where the gene is also broken

chimeric_GTF_df = gtf_df.loc[gtf_df["geneID"].isin(chimeric_genes_to_correct)]
#add columns with relevant information
#chimeric_GTF_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(chimeric_GTF_df["attribute"]) for part in element.split(";") if "gene_id" in part]
chimeric_GTF_df["exon_number"] = add_exon_number(list(chimeric_GTF_df["attribute"]))
chimeric_GTF_df["attribute_mod"] = separate_attributes(list(chimeric_GTF_df["attribute"])) #transform the attribute field in a list of tuples
#add new geneIDs as a separate column
chimeric_GTF_df["new_geneID"] = chimeric_GTF_df["geneID"].map(geneIDs_dict) #the new geneID for the chimeric gene is in the format "geneID1;geneID2"

all_chimeric_gtf_df = pd.DataFrame() #initialize gtf_df for all chimeric genes
#group by chimeric geneID and cycle on the groups
grouped_chimeric_GTF_df = chimeric_GTF_df.groupby("geneID")
for chimeric_gene, group in grouped_chimeric_GTF_df:
  print(chimeric_gene) #for debugging
  ### First gene: exon number <= boundary_ex_left
  #Test geneID: BGIBMGA001038
  #group = chimeric_GTF_df[chimeric_GTF_df["geneID"]=="BGIBMGA001038"]
  first_ex = min([int(element) for element in list(group["exon_number"]) if math.isnan(element) == False]) #I have to use min because sometimes the enumeration starts from 2
  last_ex = int(geneID_left_bound_dict[chimeric_gene]) #change data type
  print(last_ex); print(type(last_ex)) 
 
  group_exon_df = group.dropna(subset=["exon_number"])
  first_gene_df = group_exon_df[group_exon_df["exon_number"] <= last_ex]
  first_gene_df["exon_number"] = first_gene_df["exon_number"].astype("Int64") #this is necessary to have integers with NaN
  #generate the geneID, transcriptID and proteinID
  first_gene_df["new_geneID"] = [element.split(";")[0] for element in list(first_gene_df["new_geneID"])]
  first_gene_df["new_transcriptID"] = first_gene_df["new_geneID"]+transcript_suffix
  first_gene_df["new_proteinID"] = first_gene_df["new_geneID"]+protein_suffix
  #add entries for gene, transcript and start codon (if present in the original)
  first_gene_df = add_entries_first_gene(first_gene_df, group, first_ex, last_ex, transcript_suffix)
  #update the geneID, transcriptID and proteinID and exon number
  first_gene_df["attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "gene_id", list(first_gene_df["new_geneID"]))
  first_gene_df["attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "transcript_id", list(first_gene_df["new_transcriptID"]))
  first_gene_df["attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "protein_id", list(first_gene_df["new_proteinID"]))
  first_gene_df["attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "exon_number", list(first_gene_df["exon_number"]))
  #modify gene name (new gene name=chimeric_geneID_1)
  first_gene_df["attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "gene_name", [chimeric_gene+"_C1" for element in list(range(first_gene_df.shape[0]))])

  ### Second gene: exon number >= boundary_ex_right
  second_gene_df = group_exon_df[group_exon_df["exon_number"] >= int(geneID_right_bound_dict[chimeric_gene])]
  #renumber the exons (both exons and CDS)
  ex_number_list = list(sorted(list(set(list(second_gene_df["exon_number"])))))
  ex_number_dict = {float(element) : ex_number_list.index(element)+1 for element in ex_number_list}
  second_gene_df["exon_number"] = second_gene_df["exon_number"].map(ex_number_dict)
  second_gene_df["exon_number"] = second_gene_df["exon_number"].astype("Int64") #this is necessary to have integers with NaN

  #select first and last exon
  first_ex = 1 #this will always be one after re-numbering
  last_ex = max([int(element) for element in list(second_gene_df["exon_number"])])
  #generate the geneID, transcriptID and proteinID
  second_gene_df["new_geneID"] = [element.split(";")[1] for element in list(second_gene_df["new_geneID"])]
  second_gene_df["new_transcriptID"] = second_gene_df["new_geneID"]+transcript_suffix
  second_gene_df["new_proteinID"] = second_gene_df["new_geneID"]+protein_suffix
  #add entries for gene and transcript
  second_gene_df = add_entries_second_gene(second_gene_df, group, first_ex, last_ex, transcript_suffix)
  #update the geneID, transcriptID and proteinID
  second_gene_df["attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "gene_id", list(second_gene_df["new_geneID"]))
  second_gene_df["attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "transcript_id", list(second_gene_df["new_transcriptID"]))
  second_gene_df["attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "protein_id", list(second_gene_df["new_proteinID"]))
  second_gene_df["attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "exon_number", list(second_gene_df["exon_number"]))
  #modify gene name (new gene name=chimeric_geneID_1)
  second_gene_df["attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "gene_name", [chimeric_gene+"_C2" for element in list(range(second_gene_df.shape[0]))])

  #join dataframes and correct gene_source and transcript_source
  joint_chimeric_df = pd.concat([first_gene_df, second_gene_df])
  joint_chimeric_df["attribute_mod"] = modify_value_in_tuple(list(joint_chimeric_df["attribute_mod"]), "gene_source", ["brochi_pipe" for element in list(range(joint_chimeric_df.shape[0]))])
  joint_chimeric_df["attribute_mod"] = modify_value_in_tuple(list(joint_chimeric_df["attribute_mod"]), "transcript_source", ["brochi_pipe" for element in list(range(joint_chimeric_df.shape[0]))])

  joint_chimeric_df["attribute_mod"] = rebuild_attribute_entry(list(joint_chimeric_df["attribute_mod"]))
  final_joint_chimeric_df = joint_chimeric_df[["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute_mod"]] 
  all_chimeric_gtf_df = pd.concat([all_chimeric_gtf_df, final_joint_chimeric_df])

#############################################
########## JOIN ALL GTF PARTS ###############
#############################################
#save only the brochi_gtf to file
brochi_df = pd.concat([all_broken_gtf_df, all_chimeric_gtf_df])
brochi_df = brochi_df.sort_values(by=["chr","start", "stop", "type"]) #order gtf
brochi_df.to_csv(output_brochi, sep="\t", index=False, header=False, na_rep="NA", quoting=csv.QUOTE_NONE) #the csv.QUOTE_NONE avoids extra quotes aroung the attribute field

final_df = pd.concat([healthy_GTF_df, all_broken_gtf_df, all_chimeric_gtf_df]) #join all parts
#final_df = final_df.sort_values(by=["chr","start", "stop", "type"]) #order gtf (still need to figure out the exon-CDS order)
final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA", quoting=csv.QUOTE_NONE)  #save to file
