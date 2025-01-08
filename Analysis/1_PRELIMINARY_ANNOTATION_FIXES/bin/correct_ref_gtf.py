#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import numpy as np
import math
import csv
from itertools import groupby

parser = argparse.ArgumentParser(description="Script to correct broken and chimeric genes in the reference gtf")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species of interest")
parser.add_argument("--gtf", "-g", required=True, metavar="species", help="Reference gtf (only the reference transcript for each protein)")
parser.add_argument("--broken", "-b", required=True, metavar="suffix", help="One column file with semicolon separated groups of broken geneID to be joint")
parser.add_argument("--chimeric", "-c", required=True, metavar="length", help="Output of classify_chimeric_genes.py, with classification of the chimeric genes")
parser.add_argument("--IDs", "-i", required=True, metavar="input", help="File with col1=original gene IDs (of broken or chimeric genes) and col2=new geneID after correction")
parser.add_argument("--params_file", "-p", required=True, metavar="input", help="Params file with geneID infos for each species (suffix, length)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file for corrected GTF")
parser.add_argument("--output_brochi", "-ob", required=True, metavar="output_brochi", help="Path where to save only the subsetted gtf with the corrected broken and chimeric genes")
parser.add_argument("--output_unresolved", "-ou", required=True, metavar="output_unresolved", help="Path to output file with unresolved chimeric genes")
parser.add_argument("--output_resolved", "-or", required=True, metavar="output_resolved", help="Path to output file with 'resolved' chimeric genes")
parser.add_argument("--problematic_genes", "-pg", required=True, metavar="problematic_genes", help="Handmade generated file with genes that cannot be fixed")

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
output_resolved = args.output_resolved
problematic_genes_file = args.problematic_genes

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
  my_df = chimeric_df.loc[chimeric_df["chimeric_class"]=="DIFFERENT_EXONS-CONTINUOS"].copy()
  my_df.loc[:,"boundary_ex_left"] = [int(element.split(";")[0].split("-")[1]) for element in list(my_df["first-last_aligned_ex"])]
  my_df.loc[:,"boundary_ex_right"] = my_df["boundary_ex_left"]+1
  final_df = pd.concat([final_df, my_df[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]])
  ########### SAME_EXON-NO_OVERLAP and SAME_EXON_OVERLAP: assign the exon to the geneID for which the overlap is higher
  my_df = chimeric_df.loc[chimeric_df["chimeric_class"].isin(["SAME_EXON-NO_OVERLAP","SAME_EXON-OVERLAP"])].copy()
  my_df.loc[:,"ex_aa_start"] = [int(element.split("|")[1].split(";")[0].split("-")[0].split(":")[1]) for element in my_df["ex:start|stop_coord"]]
  my_df.loc[:,"ex_aa_stop"] = [int(element.split("|")[1].split(";")[0].split("-")[1]) for element in my_df["ex:start|stop_coord"]] 
  my_df.loc[:,"aa_stop1"] = [int(element.split(";")[0].split("-")[1]) for element in my_df["first-last_aligned_aa"]]
  my_df.loc[:,"aa_start2"] = [int(element.split(";")[1].split("-")[0]) for element in my_df["first-last_aligned_aa"]]
  my_df.loc[:,"overlap1"] = my_df["aa_stop1"]-my_df["ex_aa_start"]
  my_df.loc[:,"overlap2"] = my_df["ex_aa_stop"]-my_df["aa_start2"]
  my_df.loc[:,"longest_overlap"] = ["1" if element[0] >= element[1] else "2" for element in list(zip(list(my_df["overlap1"]), list(my_df["overlap2"])))]
  #If the overlap with the first gene is higher, set the boundary exon, otherwise subtract one.
  my_df.loc[:,"boundary_ex_left"] = [int(element[0].split(";")[0].split("-")[1]) if element[1]=="1" else int(element[0].split(";")[0].split("-")[1])-1 for element in zip(list(my_df["first-last_aligned_ex"]), list(my_df["longest_overlap"]))]
  my_df.loc[:,"boundary_ex_right"] = my_df["boundary_ex_left"]+1
  #Add to final dataframe
  final_df = pd.concat([final_df, my_df[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]])  
  ########### DIFFERENT_EXONS-NO_CONTINUOS: here split in correspondence of the last exon which has coverage
  my_df = chimeric_df.loc[chimeric_df["chimeric_class"]=="DIFFERENT_EXONS-NO_CONTINUOS"].copy()
  my_df.loc[:,"boundary_ex_left"] = [element.split(";")[0].split("-")[1] for element in list(my_df["first-last_aligned_ex"])]
  my_df.loc[:,"boundary_ex_right"] = [element.split(";")[1].split("-")[0] for element in list(my_df["first-last_aligned_ex"])]
  final_df = pd.concat([final_df, my_df[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]])
  ########### DIFFERENT_EXONS-OVERLAP: here in case the overlap is only 2 exons, we assign each overlapping to one of the genes
  my_df = chimeric_df.loc[chimeric_df["chimeric_class"]=="DIFFERENT_EXONS-OVERLAP"].copy()
  my_df.loc[:,"ex_stop1"] = [int(element.split(";")[0].split("-")[1]) for element in list(my_df["first-last_aligned_ex"])]
  my_df.loc[:,"ex_start2"] = [int(element.split(";")[1].split("-")[0]) for element in list(my_df["first-last_aligned_ex"])]
  my_df.loc[:,"ex_overlap"] = my_df["ex_stop1"]-my_df["ex_start2"]
  #select only the cases where the overlap is over one exon
  my_df_selected = my_df.loc[my_df["ex_overlap"]==1].copy()
  my_df_selected.loc[:,"boundary_ex_left"] = my_df_selected["ex_start2"]  #simply invert the boundaries
  my_df_selected.loc[:,"boundary_ex_right"] = my_df_selected["ex_stop1"]
  final_df = pd.concat([final_df, my_df_selected[["chimeric_geneID", "boundary_ex_left", "boundary_ex_right"]]]) 
  #separately save the DIFFERENT_EXONS-OVERLAP which can't be properly classified
  my_df_unselected = my_df.loc[my_df["ex_overlap"]!=1].copy()
  my_df_unselected = my_df_unselected.drop(columns=["ex_stop1", "ex_start2"]) #remove extra columns
  ######### COMPLETE_OVERLAP: this will not be corrected, I will act at the level of gene orthogroups.
  my_df = chimeric_df.loc[chimeric_df["chimeric_class"]=="COMPLETE_OVERLAP"].copy()
  my_df["ex_overlap"] = "complete"
  my_df_unselected = pd.concat([my_df_unselected, my_df])
  ######### INVALID_ALN
  my_df = chimeric_df.loc[chimeric_df["chimeric_class"]=="INVALID_ALN"].copy()
  my_df["ex_overlap"] = "invalid"
  my_df_unselected = pd.concat([my_df_unselected, my_df])
  ####### Generate a dataframe with all the selected entries to save to output
  selected_chimeric = list(final_df["chimeric_geneID"])
  my_df_selected = chimeric_df.loc[chimeric_df["chimeric_geneID"].isin(selected_chimeric)].copy() 
  return(final_df, my_df_unselected, my_df_selected)

def modify_value_in_tuple(attribute_field_raw, category, new_value):
  #group1["attribute_mod"] = group.attribute_mod.apply(list.copy)
  attribute_field = attribute_field_raw.copy()
  final_list = []
  index = 0
  for field_raw in attribute_field:
    field = field_raw.copy()
    for part in field:
      if part[0] == category:
        sub_index = field.index(part)
        field[sub_index][1] = f'"{new_value[index]}"' #this prints the string with surrounding quotes
    final_list = final_list+[field]
    index = index+1
  return(final_list)

def order_broken_genes(group, broken_parts):
  gene_start_coords_dict = {min(list(group.loc[(group["geneID"]==part) & (group["type"]=="CDS")]["start"])) : part for part in broken_parts} #consider only CDS entries. Gene entries are sometimes problematic.
  ordered_broken_genes = [gene_start_coords_dict[value] for value in sorted(list(gene_start_coords_dict.keys()))]
  #revert the order of the broken parts if the strand is negative
  if list(group["strand"])[0] == "-":
    ordered_broken_genes = ordered_broken_genes[::-1]
  return(ordered_broken_genes)

def add_entries_first_gene(first_gene_df, group, first_ex, last_ex, transcript_suffix):
  #gene entry and transcript entry: correct the right boundaries
  group1 = group.copy(deep=True)
  if "gene" in list(set(list(group1["type"]))): #if there is a gene entry
    gene_entry = group1.loc[group1["type"]=="gene"].copy(deep=True)
    gene_entry.loc[:,"new_geneID"] = [element.split(";")[0] for element in list(gene_entry["new_geneID"])] #select first geneID
    if str(list(gene_entry["strand"])[0]) == "+":
      gene_entry.loc[:,"stop"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["stop"])[0] #set stop from last exon
    elif str(list(gene_entry["strand"])[0]) == "-":
      gene_entry.loc[:,"start"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["start"])[0] #set start from first exon
    first_gene_df = pd.concat([first_gene_df, gene_entry])
  if "transcript" in list(set(list(group1["type"]))): #if there is a transcript entry
    transcript_entry = group1.loc[group1["type"]=="transcript"].copy(deep=True)
    transcript_entry.loc[:,"new_geneID"] = [element.split(";")[0] for element in list(transcript_entry["new_geneID"])] #select first geneID
    transcript_entry.loc[:,"new_transcriptID"] = transcript_entry["new_geneID"]+transcript_suffix
    transcript_entry.loc[:,"start"] = int(list(gene_entry["start"])[0])
    transcript_entry.loc[:,"stop"] = int(list(gene_entry["stop"])[0])
    first_gene_df = pd.concat([first_gene_df, transcript_entry])
  return(first_gene_df)

#"chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"
def add_entries_second_gene(second_gene_df, group, first_ex, last_ex, transcript_suffix):
  #all of this shit is necessary because pandas betrayed me.
  group2 = group.copy(deep=True)
  my_chr = list(second_gene_df["chr"])[0]
  my_db = list(second_gene_df["db"])[0]
  score = "."
  strand = list(second_gene_df["strand"])[0]
  phase = "."
  new_geneID = list(second_gene_df["new_geneID"])[0] #select the second geneID
  attribute_mod_gene = [['gene_id', new_geneID]]
  attribute_mod_transcript = [['gene_id', new_geneID], ['transcript_id', new_geneID+transcript_suffix]]
  if strand == "+":
    gene_entry_start = list(second_gene_df[(second_gene_df["exon_number"]==first_ex) & (second_gene_df["type"]=="exon")]["start"])[0] #set start from first ex
    gene_entry_stop = list(second_gene_df[(second_gene_df["exon_number"]==last_ex) & (second_gene_df["type"]=="exon")]["stop"])[0] #set stop from last ex
  elif strand == "-":
    gene_entry_start = int(list(second_gene_df[(second_gene_df["exon_number"]==last_ex) & (second_gene_df["type"]=="exon")]["start"])[0]) #set start from last ex
    gene_entry_stop = int(list(second_gene_df[(second_gene_df["exon_number"]==first_ex) & (second_gene_df["type"]=="exon")]["stop"])[0]) #set stop from first ex
  gene_entry = pd.DataFrame({"chr" : [my_chr], "db" : [my_db], "type" : ["gene"], "start" : [gene_entry_start], "stop" : [gene_entry_stop], "score" : [score], "strand" : [strand], "phase" : [phase], "attribute_mod" : [attribute_mod_gene], "new_geneID" : [new_geneID]}).rename(index={0 : 0.5})
  transcript_entry = pd.DataFrame({"chr" : [my_chr], "db" : [my_db], "type" : ["transcript"], "start" : [gene_entry_start], "stop" : [gene_entry_stop], "score" : [score], "strand" : [strand], "phase" : [phase], "attribute_mod" : [attribute_mod_transcript], "new_geneID" : [new_geneID], "new_transcriptID" : [new_geneID+transcript_suffix]}).rename(index={0 : 1.5})
  second_gene_df = pd.concat([second_gene_df, gene_entry, transcript_entry])
  return(second_gene_df)


def adjust_broken_phases(broken_exons_df_1, last_ex_first, broken_parts, phases_rest_dict):
  broken_exons_df = broken_exons_df_1.copy() #This is to avoid weird problems with the indexing
  strand = str(list(broken_exons_df["strand"])[0])
  first_gene = broken_parts[0]
  for second_gene in broken_parts[1:]: #access all the broken genes after the first
    first_ex_second = min(list(broken_exons_df[(broken_exons_df["geneID"]==second_gene) & (broken_exons_df["type"]=="CDS")]["exon_number"]))
    #### this is because I still have not renumbered
    last_CDS_first_gene = broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==first_gene)]
    first_CDS_second_gene = broken_exons_df.loc[(broken_exons_df["exon_number"]==first_ex_second) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==second_gene)]
    last_ex_first_phase = str(list(last_CDS_first_gene["phase"])[0]) #this works because we still have only the exons in the df
    first_ex_second_phase = str(list(first_CDS_second_gene["phase"])[0])
    phase_transition = last_ex_first_phase+"-"+first_ex_second_phase
    ##### corrected on 06/12/24
    first_ex_len = int(list(last_CDS_first_gene["stop"]-last_CDS_first_gene["start"])[0])+1
    #first_ex_len = int(list(last_CDS_first_gene["stop"]-last_CDS_first_gene["start"])[0])
    rest = first_ex_len%3 #Divide length of the last CDS exon of the first gene by 3
    if rest != phases_rest_dict[phase_transition]: #if the number of (3n + number) is different from the one required for that phase combination:
      if strand == "+":
        current_stop = list(broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==first_gene)]["stop"])[0]
        #current_start = list(broken_exons_df.loc[(broken_exons_df["exon_number"]==first_ex_second) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==second_gene)]["start"])[0]
        #cut the end of the last exon of the first gene so that it does not code for an interrupted nucleotide
        broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["geneID"]==first_gene), "stop"] = current_stop + rest - phases_rest_dict[phase_transition]        
        #the first exons of the second gene should always be in phase 0 (by definition)
        #current_stop = list(broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==first_gene)]["stop"])[0]
        #broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["geneID"]==first_gene), "stop"] = current_stop - rest + phases_rest_dict[phase_transition]
      elif strand == "-":
        current_start = list(broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==first_gene)]["start"])[0]
        broken_exons_df.loc[(broken_exons_df["exon_number"]==last_ex_first) & (broken_exons_df["geneID"]==first_gene), "start"] = current_start + rest - phases_rest_dict[phase_transition]
    #############
    #Remove all exon entries where the exon number is greater than last_ex_first(which is the last coding exon)
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["geneID"]==first_gene) & (broken_exons_df["exon_number"] > last_ex_first) & (broken_exons_df["type"]=="exon"))] #the exon filter should not be necessary, but just in case
    #############
    #update variables for next cycle.
    first_gene = second_gene
    last_ex_first = max(list(broken_exons_df.loc[(broken_exons_df["geneID"]==second_gene) & (broken_exons_df["type"]=="CDS")]["exon_number"])) 
  return(broken_exons_df)

def trim_3UTR(broken_exons_df, first_gene):
  broken_exons_df = broken_exons_df.copy()
  last_ex_first = max(list(broken_exons_df.loc[(broken_exons_df["type"]=="exon") & (broken_exons_df["geneID"]==first_gene)]["exon_number"]))
  last_CDS_first = max(list(broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==first_gene)]["exon_number"]))
  last_CDS_first_entry = broken_exons_df.loc[(broken_exons_df["exon_number"]==last_CDS_first) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==first_gene)]
  #select start and stop
  last_CDS_first_start = list(last_CDS_first_entry["start"])[0]
  last_CDS_first_stop = list(last_CDS_first_entry["stop"])[0]
  if last_CDS_first != last_ex_first:
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["geneID"]==first_gene) & (broken_exons_df["exon_number"]>=last_ex_first))] #remove all the 3'UTR exons
  #Make sure the coordinates are identical between last CDS and last exon
  broken_exons_df.loc[(broken_exons_df["geneID"]==first_gene) & (broken_exons_df["exon_number"]==last_CDS_first) & (broken_exons_df["type"]=="exon"), "start"] = last_CDS_first_start
  broken_exons_df.loc[(broken_exons_df["geneID"]==first_gene) & (broken_exons_df["exon_number"]==last_CDS_first) & (broken_exons_df["type"]=="exon"), "stop"] = last_CDS_first_stop
  return(broken_exons_df)

def trim_5UTR(broken_exons_df, last_gene):
  broken_exons_df = broken_exons_df.copy()
  first_ex_last = min(list(broken_exons_df.loc[(broken_exons_df["type"]=="exon") & (broken_exons_df["geneID"]==last_gene)]["exon_number"])) #first exon
  first_CDS_last = min(list(broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==last_gene)]["exon_number"])) #first CDS
  first_CDS_last_entry = broken_exons_df.loc[(broken_exons_df["exon_number"]==first_CDS_last) & (broken_exons_df["type"]=="CDS") & (broken_exons_df["geneID"]==last_gene)]
  first_CDS_last_start = list(first_CDS_last_entry["start"])[0]
  first_CDS_last_stop = list(first_CDS_last_entry["stop"])[0]
  if first_ex_last != first_CDS_last:
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["geneID"]==last_gene) & (broken_exons_df["exon_number"]<=first_ex_last))] #remove all 5'UTR exons
  #Make sure the coordinates are identical between first CDS and first exon
  broken_exons_df.loc[(broken_exons_df["geneID"]==last_gene) & (broken_exons_df["exon_number"]==first_CDS_last) & (broken_exons_df["type"]=="exon"), "start"] = first_CDS_last_start
  broken_exons_df.loc[(broken_exons_df["geneID"]==last_gene) & (broken_exons_df["exon_number"]==first_CDS_last) & (broken_exons_df["type"]=="exon"), "stop"] = first_CDS_last_stop
  return(broken_exons_df)

def trim_broken_UTRs(broken_exons_df, broken_parts):
  broken_exons_df = broken_exons_df.copy() #this is to avoid pandas warnings
  #correct 5'UTR of first broken gene
  first_gene = broken_parts[0]
  broken_exons_df = trim_3UTR(broken_exons_df, first_gene)
  #correct 3'UTR of last broken gene
  last_gene = broken_parts[-1] 
  broken_exons_df = trim_5UTR(broken_exons_df, last_gene)
  #correct both 5' and 3' UTR for intermediate genes
  if len(broken_parts) > 2:
    middle_genes = broken_parts[1:-1] #select only intermediate  pieces
    for gene in middle_genes:
      broken_exons_df = trim_3UTR(broken_exons_df, gene)
      broken_exons_df = trim_5UTR(broken_exons_df, gene)
  return(broken_exons_df) 

#### different version of fix_start_stop_codons: keep only the start_codon with the lowest start (+ strand) or with the highest stop ()
def fix_start_stop_codons(broken_exons_df, broken_parts):
  strand = list(broken_exons_df["strand"])[0]
  #positive strand
  if strand == "+":
    if "start_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]!=broken_parts[0]]["type"]): ### if there are start codons that are not in the first gene
    #if "start_codon" in list(broken_exons_df["type"]):
      #start codons
      if "start_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]==broken_parts[0]]["type"]):
        min_start_codon_start = min(list(broken_exons_df.loc[(broken_exons_df["type"]=="start_codon")]["start"])) #get min start codon entry
      else:
        min_start_codon_start = 0 #set arbitrary value
      removed_start_codons = broken_exons_df.loc[(broken_exons_df["start"] != min_start_codon_start) & (broken_exons_df["type"]=="start_codon")].copy() ## #isolate the removed start codons
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["start"] != min_start_codon_start) & (broken_exons_df["type"]=="start_codon"))] #remove the removed start codons
      #removed_start_codons = broken_exons_df.loc[(broken_exons_df["start"] != min_start_codon_start) & (broken_exons_df["type"]=="start_codon")] #isolate the removed start codons
      removed_start_codons_start = list(removed_start_codons["start"])
      broken_exons_df.loc[:,"start"] = [element if element not in removed_start_codons_start else element+3 for element in list(broken_exons_df["start"])]
    if "stop_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]!=broken_parts[-1]]["type"]): ### if there are stop codons that are not in the last gene
    #if "stop_codon" in list(broken_exons_df["type"]):
      #stop codons
      #### 06/12/24: I think I have to require that the stop codon is in the last gene
      #max_stop_codon_stop = max(list(broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon")]["stop"]))
      if "stop_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]==broken_parts[-1]]["type"]): #if the last gene has a stop codon, set that as maximum stop
        max_stop_codon_stop = max(list(broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon") & (broken_exons_df["geneID"]==broken_parts[-1])]["stop"]))
      else: ### remove all the stop codons
        max_stop_codon_stop = 0 #arbitrary number
      ###########
      removed_stop_codons = broken_exons_df.loc[((broken_exons_df["stop"] != max_stop_codon_stop) & (broken_exons_df["type"]=="stop_codon"))].copy() #isolate the removed stop codons
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["stop"] != max_stop_codon_stop) & (broken_exons_df["type"]=="stop_codon"))] #remove the removed stop codons
      removed_stop_codons_stop = list(removed_stop_codons["stop"])
      ##### 06/12/24
      broken_exons_df.loc[:,"stop"] = [element if element not in removed_stop_codons_stop else element-3 for element in list(broken_exons_df["stop"])]
  #negative strand
  if strand == "-":
    ##### also modified on 06/12/24
    #if "start_codon" in list(broken_exons_df["type"]):
    if "start_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]!=broken_parts[0]]["type"]):
      #start codons
      if "start_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]==broken_parts[0]]["type"]): #if there is an annotated start codon in the first gene
        max_start_codon_stop = max(list(broken_exons_df.loc[(broken_exons_df["type"]=="start_codon")]["stop"]))
      else:
        max_start_codon_stop = 0
      removed_start_codons = broken_exons_df.loc[(broken_exons_df["stop"] != max_start_codon_stop) & (broken_exons_df["type"]=="start_codon")]
      removed_start_codons_stop = list(removed_start_codons["stop"])
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["stop"] != max_start_codon_stop) & (broken_exons_df["type"]=="start_codon"))]
      #removed_start_codons = broken_exons_df.loc[(broken_exons_df["stop"] != max_start_codon_stop) & (broken_exons_df["type"]=="start_codon")]
      #removed_start_codons_stop = list(removed_start_codons["stop"])
      broken_exons_df.loc[:,"stop"] = [element if element not in removed_start_codons_stop else element-3 for element in list(broken_exons_df["stop"])]
    #if "stop_codon" in list(broken_exons_df["type"]):
    if "stop_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]!=broken_parts[-1]]["type"]):
      #stop codons
      if "stop_codon" in list(broken_exons_df.loc[broken_exons_df["geneID"]==broken_parts[-1]]["type"]): #if the last gene has a stop codon, set that as maximum stop
        min_stop_codon_start = min(list(broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon") & (broken_exons_df["geneID"]==broken_parts[-1])]["start"]))
      else:
        min_stop_codon_start = 0
      removed_stop_codons = broken_exons_df.loc[(broken_exons_df["start"] != min_stop_codon_start) & (broken_exons_df["type"]=="stop_codon")]
      removed_stop_codons_start = list(removed_stop_codons["start"])
      broken_exons_df = broken_exons_df.loc[~((broken_exons_df["start"] != min_stop_codon_start) & (broken_exons_df["type"]=="stop_codon"))]
      broken_exons_df.loc[:,"start"] = [element if element not in removed_stop_codons_start else element+3 for element in list(broken_exons_df["start"])]
  #change exon number of the stop_codon to the exon number of the last CDS exon
  last_CDS_exon = max(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["exon_number"])
  broken_exons_df.loc[broken_exons_df["type"]=="stop_codon", "exon_number"] = last_CDS_exon

  return(broken_exons_df)

def renumerate_exons(broken_exons_df, broken_parts, last_ex_first):
  final_df = broken_exons_df.loc[broken_exons_df["geneID"]==broken_parts[0]] #select first gene
  for part in broken_parts[1:]: #access all the broken genes after the first
    second_broken_exons_df = broken_exons_df.loc[broken_exons_df["geneID"]==part] #subset by gene
    second_broken_exons_df = second_broken_exons_df.copy()
    second_broken_exons_df.loc[:,"exon_number"] = second_broken_exons_df["exon_number"].astype("Int64")
    ex_number_list = list(sorted(list(set(list(second_broken_exons_df["exon_number"])))))
    new_exon_number = [int(x) for x in list(range(last_ex_first+1, last_ex_first+1+len(ex_number_list)))]
    ex_number_dict = {element : new_exon_number[ex_number_list.index(element)] for element in ex_number_list}
    second_broken_exons_df.loc[:,"exon_number"] = second_broken_exons_df["exon_number"].map(ex_number_dict)
    last_ex_first = int(max(list(second_broken_exons_df["exon_number"])))
    final_df = pd.concat([final_df, second_broken_exons_df])
  return(final_df)

#This function is an extra check that fixes the cases where there are overlapping exons between the broken genes
#(i.e. some Cterminal exons of geneID downstream of the Nterminal exons of gene2)
def renumerate_exons_overlapping_entries(broken_exons_df, strand):
  if strand == "+":
    broken_exons_df = broken_exons_df.sort_values(["start", "type"], ascending=[True,False])
  elif strand == "-":
    broken_exons_df = broken_exons_df.sort_values(["stop", "type"], ascending=[False,False])
  #remove duplicates for consecutive CDS and exon entries (physically overlapping exons from two different genes)
  broken_exons_df = broken_exons_df.loc[broken_exons_df["type"].shift(-1) != broken_exons_df["type"]]
  #order list of the exon numbers as they appear
  ex_number_list = [x[0] for x in groupby(list(broken_exons_df["exon_number"]))]
  new_exon_number = [int(x) for x in list(range(1, len(ex_number_list)+1))]
  ex_number_dict = {element : new_exon_number[ex_number_list.index(element)] for element in ex_number_list}   
  broken_exons_df.loc[:,"exon_number"] = broken_exons_df["exon_number"].map(ex_number_dict)
  return(broken_exons_df)

def add_entries_broken_genes(broken_exons_df, group, first_ex, last_ex):
  if "gene" in list(set(list(group["type"]))):
    gene_entry = group[group["type"]=="gene"].head(1) #select the first gene entry, but replace both start and stop just in case.
    if str(list(gene_entry["strand"])[0]) == "+":
      gene_entry["start"] = list(broken_exons_df[(broken_exons_df["exon_number"]==first_ex) & (broken_exons_df["type"]=="exon")]["start"])[0] #set start from first ex
      gene_entry["stop"] = list(broken_exons_df[(broken_exons_df["exon_number"]==last_ex) & (broken_exons_df["type"]=="exon")]["stop"])[0] #set stop from last ex
    elif str(list(gene_entry["strand"])[0]) == "-": #always select the min start and the max stop among all exon  entries
      gene_entry["start"] = min(list(broken_exons_df[broken_exons_df["type"]=="exon"]["start"]))
      gene_entry["stop"] = max(list(broken_exons_df[broken_exons_df["type"]=="exon"]["stop"]))
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

#adjust the phase of the first coding exon (if phase != 0)
#if strand == +: first_CDS_start = first_CDS_start + 1
#if strand == "-": first_CDS_stop = first_CDS_stop - 1
def adjust_chimeric_phases(second_gene_df):
  strand = list(second_gene_df["strand"])[0]
  first_CDS_exon = min(list(second_gene_df.loc[second_gene_df["type"]=="CDS"]["exon_number"]))
  first_CDS_exon_entry = second_gene_df.loc[(second_gene_df["type"]=="CDS") & (second_gene_df["exon_number"]==first_CDS_exon)]
  first_CDS_exon_phase = int(list(first_CDS_exon_entry["phase"])[0]) #only one value here in any case
  if first_CDS_exon_phase != 0:
    second_gene_df.loc[(second_gene_df["type"]=="CDS") & (second_gene_df["exon_number"]==first_CDS_exon), "phase"] = "0" #change the phase to 0
    if strand == "+":
      second_gene_df.loc[(second_gene_df["type"]=="CDS") & (second_gene_df["exon_number"]==first_CDS_exon), "start"] = first_CDS_exon_entry["start"]+first_CDS_exon_phase
    elif strand == "-":
      second_gene_df.loc[(second_gene_df["type"]=="CDS") & (second_gene_df["exon_number"]==first_CDS_exon), "stop"] = first_CDS_exon_entry["stop"]-first_CDS_exon_phase
  return(second_gene_df) 

### New function to be used only in case of overlapping genes
def fix_overlapping_broken_genes(broken_exons_df, broken_parts, strand, phases_current_dict):
  if strand == "+":
    ### sort dataframe based on coods
    broken_exons_df = broken_exons_df.sort_values(["start", "type"], ascending=[True,False])
    ### trim and remove intermediate UTR exons
    first_exon_CDS_start = min(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"]))
    last_exon_CDS_start = max(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"]))
    first_exon_CDS_stop = min(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"]))
    last_exon_CDS_stop = max(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"]))
    #for each CDS entry, find the corresponding exon entry and set the new boundaries (apart from first and last exon)
    CDS_entries_coords = [(element[0], element[1]) for element in zip(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"]), list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"]))]
    CDS_entries_starts = list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"])
    CDS_entries_stops = list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"])
    for CDS_entry in CDS_entries_coords:
      if CDS_entry[0] != first_exon_CDS_start or CDS_entry[0] != last_exon_CDS_start: #excluding first and last exons
        broken_exons_df.loc[((broken_exons_df["type"]=="exon") & ((broken_exons_df["start"]==CDS_entry[0]) | (broken_exons_df["stop"]==CDS_entry[1]))), "start"] = CDS_entry[0]
        broken_exons_df.loc[((broken_exons_df["type"]=="exon") & ((broken_exons_df["start"]==CDS_entry[0]) | (broken_exons_df["stop"]==CDS_entry[1]))), "stop"] = CDS_entry[1]

    ### Remove all exon entries that do not match a CDS and are not UTRs
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="exon") & ((broken_exons_df["start"] > first_exon_CDS_start) | (broken_exons_df["stop"] < first_exon_CDS_stop)) & ~(broken_exons_df["start"].isin(CDS_entries_starts)))].copy()

    ### Compare consecutive CDS if they come from different genes, and fix their phase
    for index in list(range(1,len(CDS_entries_coords))):
       previous_entry = broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index-1][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index-1][1]))].copy() 
       current_entry = broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index][1]))].copy()
       previous_gene = list(previous_entry["geneID"])[0]
       current_gene = list(current_entry["geneID"])[0] 
       if previous_gene != current_gene: ### check phases
         previous_phase = list(previous_entry["phase"])[0]
         current_phase = list(current_entry["phase"])[0]
         print(previous_entry)
         previous_len = list(previous_entry["stop"] - previous_entry["start"] +1)[0]
         current_len = list(current_entry["stop"]  - current_entry["start"] +1)[0]
         ### trim end of previous exon so that it matches the phase
         if previous_len%3 != phases_current_dict[previous_phase]:
           previous_stop = list(previous_entry["stop"])[0]
           if (previous_phase == 0 & previous_len%3==1) or (previous_phase == 1 & previous_len%3==0) or (previous_phase == 2 & previous_len%3==2):
             selected_len = 1
           elif (previous_phase == 0 & previous_len%3==2) or (previous_phase == 1 & previous_len%3==1) or (previous_phase == 2 & previous_len%3==0):
             selected_len = 1
           else:
             selected_len = 0
           ### correct length
           broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index-1][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index-1][1])), "stop"] = previous_stop - selected_len
         ### trim start of the current exon so that it matches the phase
         if current_len%3 != phases_current_dict[current_phase]:
           current_start = list(current_entry["start"])[0]
           if (previous_phase == 0 & current_len%3==1) or (current_phase == 1 & current_len%3==0) or (current_phase == 2 & current_len%3==2):
             selected_len = 1
           elif (current_phase == 0 & current_len%3==2) or (current_phase == 1 & current_len%3==1) or (current_phase == 2 & current_len%3==0):
             selected_len = 1
           else:
             selected_len = 0
           ### correct length AND phase
           broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index][1])), "start"] = current_start + selected_len
           broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index][1])), "phase"] = "0"

    #### remove intermediate start codon entries
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="start_codon") & (broken_exons_df["start"] != first_exon_CDS_start))]

    #### remove intermediate stop codon entries and correct relative exon/CDS entries
    removed_stop_codons = broken_exons_df.loc[((broken_exons_df["stop"] != last_exon_CDS_stop) & (broken_exons_df["type"]=="stop_codon"))].copy() #isolate the removed stop codons
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["stop"] != last_exon_CDS_stop) & (broken_exons_df["type"]=="stop_codon"))] #remove the removed stop codons
    removed_stop_codons_stop = list(removed_stop_codons["stop"])
    broken_exons_df.loc[:,"stop"] = [element if element not in removed_stop_codons_stop else element-3 for element in list(broken_exons_df["stop"])]

    #### re-enumerate exons
    broken_exons_df = broken_exons_df.copy() ## The dataframe should already be ordered
    broken_exons_df.loc[:,"exon_number"] = broken_exons_df["exon_number"].astype("Int64") # set data type to int
    #### Asssign exon numbers to exon entries. If the CDS/start_codon/stop_codon entries have the start/stop in common, assign the same exon number
    exon_entries_coords = [(element[0], element[1]) for element in zip(list(broken_exons_df.loc[broken_exons_df["type"]=="exon"]["start"]), list(broken_exons_df.loc[broken_exons_df["type"]=="exon"]["stop"]))]
    for index in list(range(0,len(exon_entries_coords))):
       broken_exons_df.loc[((broken_exons_df["type"]=="exon") & (broken_exons_df["start"]==exon_entries_coords[index][0]) & (broken_exons_df["stop"]==exon_entries_coords[index][1])), "exon_number"] = index+1
    # create start and stop dictionary with key=exon coords
    exon_entries_df = broken_exons_df.loc[broken_exons_df["type"]=="exon"].copy()
    exon_start_coords_dictionary = pd.Series(exon_entries_df.exon_number.values, index=exon_entries_df.start).to_dict()
    exon_stop_coords_dictionary = pd.Series(exon_entries_df.exon_number.values, index=exon_entries_df.stop).to_dict()
    broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & (broken_exons_df["start"].isin(exon_start_coords_dictionary)), "exon_number"] = broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & (broken_exons_df["start"].isin(exon_start_coords_dictionary))]["start"].map(exon_start_coords_dictionary)
    # for the CDS of the first exon, whose start might not match the exon start
    broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & ~(broken_exons_df["start"].isin(exon_start_coords_dictionary)), "exon_number"] = broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & ~(broken_exons_df["start"].isin(exon_start_coords_dictionary))]["stop"].map(exon_stop_coords_dictionary)
    # create start and stop dictionary with key=CDS coords
    CDS_entries_df = broken_exons_df.loc[broken_exons_df["type"]=="CDS"].copy()
    CDS_start_coords_dictionary = pd.Series(CDS_entries_df.exon_number.values, index=CDS_entries_df.start).to_dict()
    CDS_stop_coords_dictionary = pd.Series(CDS_entries_df.exon_number.values, index=CDS_entries_df.stop).to_dict()
    # add exon number to start and stop codons
    broken_exons_df.loc[(broken_exons_df["type"]=="start_codon"), "exon_number"] =  broken_exons_df.loc[broken_exons_df["type"]=="start_codon"]["start"].map(CDS_start_coords_dictionary) 
    broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon"), "exon_number"] =  broken_exons_df.loc[broken_exons_df["type"]=="stop_codon"]["stop"].map(CDS_stop_coords_dictionary)
    
   
  if strand == "-":
    ### sort dataframe based on coods
    broken_exons_df = broken_exons_df.sort_values(["stop", "type"], ascending=[False,False])
    ### trim and remove intermediate UTR exons
    first_exon_CDS_start = max(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"]))
    last_exon_CDS_start = min(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"]))
    first_exon_CDS_stop = max(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"]))
    last_exon_CDS_stop = min(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"]))
    #for each CDS entry, find the corresponding exon entry and set the new boundaries (apart from first and last exon)
    CDS_entries_coords = [(element[0], element[1]) for element in zip(list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"]), list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"]))]
    CDS_entries_starts = list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["start"])
    CDS_entries_stops = list(broken_exons_df.loc[broken_exons_df["type"]=="CDS"]["stop"])
    for CDS_entry in CDS_entries_coords:
      if CDS_entry[0] != first_exon_CDS_start or CDS_entry[0] != last_exon_CDS_start: #excluding first and last exons
        broken_exons_df.loc[((broken_exons_df["type"]=="exon") & ((broken_exons_df["start"]==CDS_entry[0]) | (broken_exons_df["stop"]==CDS_entry[1]))), "start"] = CDS_entry[0]
        broken_exons_df.loc[((broken_exons_df["type"]=="exon") & ((broken_exons_df["start"]==CDS_entry[0]) | (broken_exons_df["stop"]==CDS_entry[1]))), "stop"] = CDS_entry[1]

    ### Remove all exon entries that do not match a CDS and are not UTRs
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="exon") & ((broken_exons_df["start"] < first_exon_CDS_start) | (broken_exons_df["stop"] > first_exon_CDS_stop)) & ~(broken_exons_df["start"].isin(CDS_entries_starts)))]

    ### Compare consecutive CDS if they come from different genes, and fix their phase
    ### the first part is also identical between strands, in principle
    phases_current_dict = {"0":0,"1":2,"2":1}
    for index in list(range(1,len(CDS_entries_coords))):
       print(index)
       previous_entry = broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index-1][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index-1][1]))].copy() 
       current_entry = broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index][1]))].copy()
       previous_gene = list(previous_entry["geneID"])[0]
       current_gene = list(current_entry["geneID"])[0] 
       if previous_gene != current_gene: ### check phases
         previous_phase = list(previous_entry["phase"])[0]
         current_phase = list(current_entry["phase"])[0]
         previous_len = list((previous_entry["stop"]-previous_entry["start"])+1)[0]
         current_len = list((current_entry["stop"]-current_entry["start"])+1)[0]
         ### trim end of previous exon so that it matches the phase
         if previous_len%3 != phases_current_dict[previous_phase]:
           previous_start = list(previous_entry["start"])[0]
           if (previous_phase == 0 & previous_len%3==1) or (previous_phase == 1 & previous_len%3==0) or (previous_phase == 2 & previous_len%3==2):
             selected_len = 1
           elif (previous_phase == 0 & previous_len%3==2) or (previous_phase == 1 & previous_len%3==1) or (previous_phase == 2 & previous_len%3==0):
             selected_len = 1
           else:
             selected_len = 0
           ### correct length
           broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index-1][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index-1][1])), "start"] = previous_start + selected_len
         ### trim start of the current exon so that it matches the phase
         if current_len%3 != phases_current_dict[current_phase]:
           current_stop = list(current_entry["stop"])[0]
           if (previous_phase == 0 & current_len%3==1) or (current_phase == 1 & current_len%3==0) or (current_phase == 2 & current_len%3==2):
             selected_len = 1
           elif (current_phase == 0 & current_len%3==2) or (current_phase == 1 & current_len%3==1) or (current_phase == 2 & current_len%3==0):
             selected_len = 1
           else:
             selected_len = 0
           ### correct length AND phase
           broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index][1])), "stop"] = current_stop - selected_len
           broken_exons_df.loc[((broken_exons_df["type"]=="CDS") & (broken_exons_df["start"]==CDS_entries_coords[index][0]) & (broken_exons_df["stop"]==CDS_entries_coords[index][1])), "phase"] = "0"
         print(broken_exons_df)

    #### remove intermediate start codon entries
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["type"]=="start_codon") & (broken_exons_df["stop"] != first_exon_CDS_stop))]

    #### remove intermediate stop codon entries and correct relative exon/CDS entries
    removed_stop_codons = broken_exons_df.loc[((broken_exons_df["start"] != last_exon_CDS_start) & (broken_exons_df["type"]=="stop_codon"))].copy() #isolate the removed stop codons
    broken_exons_df = broken_exons_df.loc[~((broken_exons_df["start"] != last_exon_CDS_start) & (broken_exons_df["type"]=="stop_codon"))] #remove the removed stop codons
    removed_stop_codons_start = list(removed_stop_codons["start"])
    broken_exons_df.loc[:,"start"] = [element if element not in removed_stop_codons_start else element+3 for element in list(broken_exons_df["start"])]

    #### re-enumerate exons
    broken_exons_df = broken_exons_df.copy() ## The dataframe should already be ordered
    broken_exons_df.loc[:,"exon_number"] = broken_exons_df["exon_number"].astype("Int64") # set data type to int
    #### Asssign exon numbers to exon entries. If the CDS/start_codon/stop_codon entries have the start/stop in common, assign the same exon number
    exon_entries_coords = [(element[0], element[1]) for element in zip(list(broken_exons_df.loc[broken_exons_df["type"]=="exon"]["start"]), list(broken_exons_df.loc[broken_exons_df["type"]=="exon"]["stop"]))]
    for index in list(range(0,len(exon_entries_coords))):
       broken_exons_df.loc[((broken_exons_df["type"]=="exon") & (broken_exons_df["start"]==exon_entries_coords[index][0]) & (broken_exons_df["stop"]==exon_entries_coords[index][1])), "exon_number"] = index+1
    # create start and stop dictionary with key=exon coords
    exon_entries_df = broken_exons_df.loc[broken_exons_df["type"]=="exon"].copy()
    exon_start_coords_dictionary = pd.Series(exon_entries_df.exon_number.values, index=exon_entries_df.start).to_dict()
    exon_stop_coords_dictionary = pd.Series(exon_entries_df.exon_number.values, index=exon_entries_df.stop).to_dict()
    broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & (broken_exons_df["start"].isin(exon_start_coords_dictionary)), "exon_number"] = broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & (broken_exons_df["start"].isin(exon_start_coords_dictionary))]["start"].map(exon_start_coords_dictionary)
    # for the CDS of the last CDS exon whose start might not match the exon start
    broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & ~(broken_exons_df["start"].isin(exon_start_coords_dictionary)), "exon_number"] = broken_exons_df.loc[(broken_exons_df["type"]=="CDS") & ~(broken_exons_df["start"].isin(exon_start_coords_dictionary))]["stop"].map(exon_stop_coords_dictionary)
    # create start and stop dictionary with key=CDS coords
    CDS_entries_df = broken_exons_df.loc[broken_exons_df["type"]=="CDS"].copy()
    CDS_start_coords_dictionary = pd.Series(CDS_entries_df.exon_number.values, index=CDS_entries_df.start).to_dict()
    CDS_stop_coords_dictionary = pd.Series(CDS_entries_df.exon_number.values, index=CDS_entries_df.stop).to_dict()
    # add exon number to start and stop codons
    broken_exons_df.loc[(broken_exons_df["type"]=="start_codon"), "exon_number"] =  broken_exons_df.loc[broken_exons_df["type"]=="start_codon"]["stop"].map(CDS_stop_coords_dictionary) 
    broken_exons_df.loc[(broken_exons_df["type"]=="stop_codon"), "exon_number"] =  broken_exons_df.loc[broken_exons_df["type"]=="stop_codon"]["start"].map(CDS_start_coords_dictionary)
    print(broken_exons_df["exon_number"])
 
    ##### I still need to fix exon numbers 
  print(broken_exons_df) 
  return(broken_exons_df) 

##################################
###### READ INPUTS ###############
##################################

gtf_df = pd.read_table(gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
gtf_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part] #add geneID as a separate field
broken_genes_list = list(pd.read_table(broken_genes_file, sep="\t", index_col=False, header=None, names=["broken_genes"])["broken_genes"])
broken_gene_flatten_list = [part for element in broken_genes_list for part in element.split(";")]

#Header: species, chimeric_geneID, orthogroup_ids, first-last_aligned_ex, first-last_aligned_aa, ex:start|stop_coord, chimeric_class
chimeric_genes_df = pd.read_table(chimeric_genes_file, sep="\t", index_col=False, header=0)
chimeric_genes_df = chimeric_genes_df[chimeric_genes_df["species"]==species] #subset for the species of interest
chimeric_genes_list = list(chimeric_genes_df["chimeric_geneID"])

#Header: geneID, new_IDs, category
geneIDs_df = pd.read_table(geneIDs_file, sep="\t", index_col=False, header=0)
geneIDs_dict = pd.Series(geneIDs_df.new_IDs.values, index=geneIDs_df.geneID).to_dict()
reverse_geneID_dict = pd.Series(geneIDs_df.geneID.values, index=geneIDs_df.new_IDs).to_dict()
#correct broken genes entries
broken_genes_keys = [key for key in list(geneIDs_dict.keys()) if ";" in key]
for key in broken_genes_keys:
  for num in list(range(0, len(key.split(";")))):
    geneIDs_dict[key.split(";")[num]] = geneIDs_dict[key]

#Header: species, geneID_prefix, geneID_length, transcript_suffix, protein_suffix, transcript prefix, protein prefix
params_df = pd.read_table(my_params_file, sep="\t", header=0, index_col=False)
transcript_suffix = str(list(params_df[params_df["species"]==species]["transcript_suffix"])[0])
protein_suffix = str(list(params_df[params_df["species"]==species]["protein_suffix"])[0])
##added
transcript_prefix = str(list(params_df[params_df["species"]==species]["transcript_prefix"])[0])
protein_prefix = str(list(params_df[params_df["species"]==species]["protein_prefix"])[0])

print("%s: %d" % ("broken_genes", len(broken_gene_flatten_list))) #for debugging
print("%s: %d" % ("chimeric_genes", len(chimeric_genes_list))) #for debugging

#12/09/2021: read in problematic genes file. These genes can't be corrected and will have to be filtered out.
problematic_genes_list = list(pd.read_table(problematic_genes_file, sep="\t", index_col=False, header=None, names=["GeneIDs"])["GeneIDs"])
problematic_genes_flatten_list = [part for element in problematic_genes_list for part in element.split(";")]
broken_gene_flatten_list = [element for element in broken_gene_flatten_list if element not in problematic_genes_flatten_list]

##################################
####### BROKEN GENES #############
##################################

#subset gtf to broken genes.
broken_GTF_df = gtf_df.loc[gtf_df["geneID"].isin(broken_gene_flatten_list)]
#add exon number column and transform the attribute field in a list of tuples
broken_GTF_df.loc[:,"exon_number"] = add_exon_number(list(broken_GTF_df["attribute"]))
broken_GTF_df.loc[:,"exon_number"] = broken_GTF_df["exon_number"].astype("Int64") #This is necessary to have integers with NaN
broken_GTF_df.loc[:,"gene_name"] = add_gene_name(list(broken_GTF_df["attribute"]))
broken_GTF_df.loc[:,"attribute_mod"] = separate_attributes(list(broken_GTF_df["attribute"])) 
#add new geneID column (use the dict)
broken_GTF_df.loc[:,"new_geneID"] = broken_GTF_df["geneID"].map(geneIDs_dict)

#generate dictionary to adjust phases
#Add/subtract 2 to the stop/start of the first broken gene (depending on the strand) to pass to the previous phase (once 3n)
#Add/subtract 2 to the stop/start of the first broken gene (depending on the strand) to pass to the following phase (once 3n)
#leave unchanged when the phase is the same (once 3n)
phases_rest_dict = {"0-1":2, "1-2":2, "2-0":2,"1-0":1, "2-1":1, "0-2":1, "0-0":0, "1-1":0, "2-2":0}
phases_current_dict = {"0":0,"1":2,"2":1}

all_broken_gtf_df = pd.DataFrame() #initialize gtf_df for all broken genes
#groupby the new geneID column and cycle on the groups.
grouped_broken_GTF_df = broken_GTF_df.groupby("new_geneID")
for repaired_gene, group in grouped_broken_GTF_df:
  group = group.copy(deep=True)
  strand = list(group["strand"])[0]
  #derive variables
  broken_parts = reverse_geneID_dict[repaired_gene].split(";") #the genes are always ordered according to genomic coordinates
  #make sure that they are ordered based on the coding order (different in case of positive and negative strand)
  broken_parts = order_broken_genes(group, broken_parts)

  #add gene_name
  new_gene_name = ';'.join([name for name in list(set(list(group["gene_name"]))) if name != "NoName"])
  #add transcriptID and proteinID columns to group.
  group.loc[:,"new_transcriptID"] = [element+transcript_suffix for element in list(group["new_geneID"])]
  group.loc[:,"new_proteinID"] = [element+protein_suffix for element in list(group["new_geneID"])]
  #broken_df_exons; select only entries with exons.
  broken_exons_df = group.dropna(subset=["exon_number"])
  broken_exons_df = broken_exons_df.copy(deep=True)

  #make sure exons are ordered (by start and stop coords), and re-number them (the second gene will change).
  broken_exons_df = broken_exons_df.sort_values(by=["start", "stop"])
  last_ex_first_CDS = int(max(list(broken_exons_df[(broken_exons_df["geneID"]==broken_parts[0]) & (broken_exons_df["type"]=="CDS")]["exon_number"])))

###############
  #check if the last or the second gene overlap the first one. This is to fix weird overlapping cases.
  if strand == "+":
    last_ex_first_CDS_stop = int(max(list(broken_exons_df[(broken_exons_df["geneID"]==broken_parts[0]) & (broken_exons_df["type"]=="CDS")]["stop"]))) 
    first_ex_last_CDS_start = int(min(list(broken_exons_df[(broken_exons_df["geneID"]==broken_parts[1]) & (broken_exons_df["type"]=="CDS")]["start"]))) #added 07/12/24
    if first_ex_last_CDS_start <= last_ex_first_CDS_stop:
       overlapping_status = True
    else:
       overlapping_status = False
  if strand == "-":
    last_ex_first_CDS_stop = int(min(list(broken_exons_df[(broken_exons_df["geneID"]==broken_parts[0]) & (broken_exons_df["type"]=="CDS")]["stop"])))
    first_ex_last_CDS_start = int(max(list(broken_exons_df[(broken_exons_df["geneID"]==broken_parts[1]) & (broken_exons_df["type"]=="CDS")]["start"])))
    last_ex_first_CDS_start = int(min(list(broken_exons_df[(broken_exons_df["geneID"]==broken_parts[0]) & (broken_exons_df["type"]=="CDS")]["start"])))
    if first_ex_last_CDS_start >= last_ex_first_CDS_stop or first_ex_last_CDS_start >= last_ex_first_CDS_start:
      overlapping_status = True
    else:
      overlapping_status = False

###########
  if overlapping_status: #start of first CDS of second gene preceeds the last CDS of the first gene
    print(repaired_gene + " overlapping")
    broken_exons_df = fix_overlapping_broken_genes(broken_exons_df, broken_parts, strand, phases_current_dict)
##########
  else:
    #Adjust exon boundaries depending on phase combinations between last ex of first_broken and first_ex of second broken gene
    broken_exons_df = adjust_broken_phases(broken_exons_df, last_ex_first_CDS, broken_parts, phases_rest_dict)
    
    #Trim 3'UTR of the first broken gene and 5'UTR of the second broken gene (both for "intermediate") genes
    broken_exons_df = trim_broken_UTRs(broken_exons_df, broken_parts)
    #renumber the exons (both exons and CDS)
    broken_exons_df = renumerate_exons(broken_exons_df, broken_parts, last_ex_first_CDS)
    #Fix start and stop exons
    broken_exons_df = fix_start_stop_codons(broken_exons_df, broken_parts)
    #extra check necessary in case of overlapping exon entries between broken genes
    broken_exons_df = renumerate_exons_overlapping_entries(broken_exons_df, strand)
  
  first_ex = int(min(list(broken_exons_df["exon_number"]))) #this should be one in the majority of cases, but who knows
  last_ex = int(max(list(broken_exons_df["exon_number"])))
  #If originally there were gene and transcript entries, take them and modify the start and stop.
  broken_exons_df = add_entries_broken_genes(broken_exons_df, group, first_ex, last_ex)

  #update the geneID, transcriptID and proteinID in loco
  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "gene_id", list(broken_exons_df["new_geneID"]))
  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "transcript_id", list(broken_exons_df["new_transcriptID"]))
  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "protein_id", list(broken_exons_df["new_proteinID"]))
  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "exon_number", list(broken_exons_df["exon_number"]))

  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "gene_source", ["brochi_pipe" for element in list(range(broken_exons_df.shape[0]))])
  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "transcript_source", ["brochi_pipe" for element in list(range(broken_exons_df.shape[0]))])
  broken_exons_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(broken_exons_df["attribute_mod"]), "gene_name", [new_gene_name for element in list(range(broken_exons_df.shape[0]))])
  #Rebuild the attribute field in the original format
  broken_exons_df.loc[:,"attribute_mod"] = rebuild_attribute_entry(list(broken_exons_df["attribute_mod"]))
  final_broken_exons_df = broken_exons_df[["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute_mod"]]
  #order base on strand:
  custom_order = ["gene", "transcript", "exon", "CDS", "start_codon", "stop_codon"]
  final_broken_exons_df["type"] = pd.Categorical(final_broken_exons_df["type"], categories=custom_order, ordered=True)
  if strand == "+":
    final_broken_exons_df = final_broken_exons_df.sort_values(["start", "type"], ascending=[True,True])
  elif strand == "-":
    final_broken_exons_df = final_broken_exons_df.sort_values(["stop", "type"], ascending=[False,True])
  all_broken_gtf_df = pd.concat([all_broken_gtf_df, final_broken_exons_df])

##################################
####### CHIMERIC GENES ###########
##################################

###Preprocessing: get exon boundaries for all categories
chimeric_genes_to_correct = list(chimeric_genes_df["chimeric_geneID"])
#generate a boundary exons df
res = def_boundary_exon(chimeric_genes_df)
boundary_ex_df = res[0]
unresolved_chimeric_df = res[1]
resolved_chimeric_df = res[2]
#Save unresolved_df to file
unresolved_chimeric_df.to_csv(output_unresolved, sep="\t", index=False, header=True, na_rep="NA")
#Save resolved_df to file
resolved_chimeric_df.to_csv(output_resolved, sep="\t", index=False, header=True, na_rep="NA")


#create dictionary with correspondence between geneID - left boundary and geneID - right boundary
geneID_left_bound_dict = pd.Series(boundary_ex_df.boundary_ex_left.values, index=boundary_ex_df.chimeric_geneID).to_dict()
geneID_right_bound_dict = pd.Series(boundary_ex_df.boundary_ex_right.values, index=boundary_ex_df.chimeric_geneID).to_dict()

###Correct GTF
#subset GTF to chimeric genes to correct
chimeric_genes_to_correct = [element for element in chimeric_genes_to_correct if element not in list(unresolved_chimeric_df["chimeric_geneID"])] #remove "unresolved" genes (these are saved to file separately).
chimeric_genes_to_correct = [element for element in chimeric_genes_to_correct if element not in broken_gene_flatten_list] #remove those fucking cases where the gene is also broken.
#chimeric_genes_to_correct = [element for element in chimeric_genes_to_correct if element not in problematic_genes_flatten_list] #remove the problematic genes that I could not find a way to properly correct.

chimeric_GTF_df = gtf_df.loc[gtf_df["geneID"].isin(chimeric_genes_to_correct)]
#add columns with relevant information
#chimeric_GTF_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(chimeric_GTF_df["attribute"]) for part in element.split(";") if "gene_id" in part]
chimeric_GTF_df.loc[:,"exon_number"] = add_exon_number(list(chimeric_GTF_df["attribute"]))
chimeric_GTF_df["attribute_mod"] = separate_attributes(list(chimeric_GTF_df["attribute"])) #transform the attribute field in a list of tuples
#add new geneIDs as a separate column
chimeric_GTF_df.loc[:,"new_geneID"] = chimeric_GTF_df["geneID"].map(geneIDs_dict) #the new geneID for the chimeric gene is in the format "geneID1;geneID2"

all_chimeric_gtf_df = pd.DataFrame() #initialize gtf_df for all chimeric genes
#group by chimeric geneID and cycle on the groups
grouped_chimeric_GTF_df = chimeric_GTF_df.groupby("geneID")
for chimeric_gene, group in grouped_chimeric_GTF_df:
  strand = list(group["strand"])[0]
  print(chimeric_gene) #for debugging
  ### First gene: exon number <= boundary_ex_left
  #group = chimeric_GTF_df[chimeric_GTF_df["geneID"]=="BGIBMGA001038"]
  first_ex = min([int(element) for element in list(group["exon_number"]) if math.isnan(element) == False]) #I have to use min because sometimes the enumeration starts from 2
  last_ex = int(geneID_left_bound_dict[chimeric_gene]) #change data type
 
  group_exon_df = group.dropna(subset=["exon_number"])
  first_gene_df = group_exon_df[group_exon_df["exon_number"] <= last_ex]
  first_gene_df = first_gene_df.copy(deep=True)
  first_gene_df.loc[:,"exon_number"] = first_gene_df["exon_number"].astype("Int64") #this is necessary to have integers with NaN
  #generate the geneID, transcriptID and proteinID
  first_gene_df.loc[:,"new_geneID"] = [element.split(";")[0] for element in list(first_gene_df["new_geneID"])]
  first_gene_df.loc[:,"new_transcriptID"] = first_gene_df["new_geneID"]+transcript_suffix
  first_gene_df.loc[:,"new_proteinID"] = first_gene_df["new_geneID"]+protein_suffix
  #add entries for gene, transcript and start codon (if present in the original)
  first_gene_df = add_entries_first_gene(first_gene_df, group, first_ex, last_ex, transcript_suffix)
  #update the geneID, transcriptID and proteinID and exon number
  first_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "gene_id", list(first_gene_df["new_geneID"]))
  first_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "transcript_id", list(first_gene_df["new_transcriptID"]))
  first_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "protein_id", list(first_gene_df["new_proteinID"]))
  first_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "exon_number", list(first_gene_df["exon_number"]))
  #modify gene name (new gene name=chimeric_geneID_1)
  first_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "gene_name", [chimeric_gene+"_C1" for element in list(range(first_gene_df.shape[0]))])

  ### Second gene: exon number >= boundary_ex_right
  second_gene_df = group_exon_df[group_exon_df["exon_number"] >= int(geneID_right_bound_dict[chimeric_gene])]
  second_gene_df = second_gene_df.copy(deep=True)
  #renumber the exons (both exons and CDS)
  ex_number_list = list(sorted(list(set(list(second_gene_df["exon_number"])))))
  ex_number_dict = {float(element) : ex_number_list.index(element)+1 for element in ex_number_list}
  second_gene_df.loc[:,"exon_number"] = second_gene_df["exon_number"].map(ex_number_dict)
  second_gene_df.loc[:,"exon_number"] = second_gene_df["exon_number"].astype("Int64") #this is necessary to have integers with NaN

  #select first and last exon
  first_ex = 1 #this will always be one after re-numbering
  last_ex = max([int(element) for element in list(second_gene_df["exon_number"])])
  #generate the geneID, transcriptID and proteinID
  second_gene_df.loc[:,"new_geneID"] = [element.split(";")[1] for element in list(second_gene_df["new_geneID"])]
  second_gene_df.loc[:,"new_transcriptID"] = second_gene_df["new_geneID"]+transcript_suffix
  second_gene_df.loc[:,"new_proteinID"] = second_gene_df["new_geneID"]+protein_suffix
  #adjust the phase of the first coding exon (if phase != 0)
  second_gene_df = adjust_chimeric_phases(second_gene_df)
  #add entries for gene and transcript
  second_gene_df = add_entries_second_gene(second_gene_df, group, first_ex, last_ex, transcript_suffix)

  #update the geneID, transcriptID and proteinID
  second_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "gene_id", list(second_gene_df["new_geneID"]))
  second_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "transcript_id", list(second_gene_df["new_transcriptID"]))
  second_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "protein_id", list(second_gene_df["new_proteinID"]))
  second_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "exon_number", list(second_gene_df["exon_number"]))
  #modify gene name (new gene name=chimeric_geneID_1)
  second_gene_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(second_gene_df["attribute_mod"]), "gene_name", [chimeric_gene+"_C2" for element in list(range(second_gene_df.shape[0]))])

  #join dataframes and correct gene_source and transcript_source
  joint_chimeric_df = pd.concat([first_gene_df, second_gene_df])
  joint_chimeric_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(joint_chimeric_df["attribute_mod"]), "gene_source", ["brochi_pipe" for element in list(range(joint_chimeric_df.shape[0]))])
  joint_chimeric_df.loc[:,"attribute_mod"] = modify_value_in_tuple(list(joint_chimeric_df["attribute_mod"]), "transcript_source", ["brochi_pipe" for element in list(range(joint_chimeric_df.shape[0]))])
  #joint_chimeric_df["attribute_mod"] = modify_value_in_tuple(list(joint_chimeric_df["attribute_mod"]), "gene_id", list(joint_chimeric_df["new_geneID"]))

  joint_chimeric_df.loc[:,"attribute_mod"] = rebuild_attribute_entry(list(joint_chimeric_df["attribute_mod"]))
  final_joint_chimeric_df = joint_chimeric_df[["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute_mod"]]
  #order depending on strand
  custom_order = ["gene", "transcript", "exon", "CDS", "start_codon", "stop_codon"]
  final_joint_chimeric_df["type"] = pd.Categorical(final_joint_chimeric_df["type"], categories=custom_order, ordered=True)
  if strand == "+":
    final_joint_chimeric_df = final_joint_chimeric_df.sort_values(["start", "type"], ascending=[True,True])
  elif strand == "-":
    final_joint_chimeric_df = final_joint_chimeric_df.sort_values(["stop", "type"], ascending=[False,True])
  all_chimeric_gtf_df = pd.concat([all_chimeric_gtf_df, final_joint_chimeric_df])
  

##################################
####### HEALTHY GENES ############
##################################
#add geneID
brochi_genes_list = broken_gene_flatten_list + chimeric_genes_list
print(broken_gene_flatten_list)
print(chimeric_genes_list)
#Consider as "healthy" the chimeric genes that cannot be fixed and are not among the broken ones (another small exception)
brochi_genes_list = [gene for gene in brochi_genes_list if gene not in [element for element in list(unresolved_chimeric_df["chimeric_geneID"]) if element not in broken_gene_flatten_list]]
print(brochi_genes_list)
#filter out entries of brochi genes from the GTF
healthy_GTF_df = gtf_df[~(gtf_df["geneID"].isin(brochi_genes_list))] 
healthy_GTF_df = healthy_GTF_df.rename(columns={"attribute" : "attribute_mod"}) #rename field to match the broken and chimeric gtf dataframes
healthy_GTF_df = healthy_GTF_df[["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute_mod"]]

#############################################
########## JOIN ALL GTF PARTS ###############
#############################################
#save only the brochi_gtf to file
brochi_df = pd.concat([all_broken_gtf_df, all_chimeric_gtf_df])
#brochi_df = brochi_df.sort_values(by=["chr", "start", "type"], ascending=[True,True,False]) #order gtf
brochi_df.to_csv(output_brochi, sep="\t", index=False, header=False, na_rep="NA", quoting=csv.QUOTE_NONE) #the csv.QUOTE_NONE avoids extra quotes aroung the attribute field

final_df = pd.concat([healthy_GTF_df, all_broken_gtf_df, all_chimeric_gtf_df]) #join all parts
#final_df = final_df.sort_values(by=["chr", "start", "type"], ascending=[True,True,False]) #order gtf (still need to figure out the exon-CDS order)
final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA", quoting=csv.QUOTE_NONE)  #save to file
