#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import numpy as np
import math

parser = argparse.ArgumentParser(description="Script to correct broken and chimeric genes in the reference gtf")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species of interest")
parser.add_argument("--gtf", "-g", required=True, metavar="species", help="Reference gtf (only the reference transcript for each protein)")
parser.add_argument("--broken", "-b", required=True, metavar="suffix", help="One column file with semicolon separated groups of broken geneID to be joint")
parser.add_argument("--chimeric", "-c", required=True, metavar="length", help="Output of classify_chimeric_genes.py, with classification of the chimeric genes")
parser.add_argument("--IDs", "-i", required=True, metavar="input", help="File with col1=original gene IDs (of broken or chimeric genes) and col2=new geneID after correction")
parser.add_argument("--params_file", "-p", required=True, metavar="input", help="Params file with geneID infos for each species (suffix, length)")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file")

###### Read arguments
args = parser.parse_args()
species = args.species
gtf_file = args.gtf
broken_genes_file = args.broken
chimeric_genes_file = args.chimeric
geneIDs_file = args.IDs
output_file = args.output

###### Define functions
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
  return(final_df, my_df_unselected)

def modify_value_in_tuple(attribute_field, category, new_value):
  final_list = []
  for field in attribute_field:
    index = attribute_field.index(field)
    for part in field:
      if part[0] == category:
        sub_index = field.index(part)
        field[sub_index][1] = new_value[index]
    final_list = final_list+[field]
  return(final_list)

def add_entries_first_gene(first_gene_df, group, first_ex, last_ex, anscript_suffix):
  #gene entry and transcript entry: correct the right boundaries
  if "gene" in list(set(list(group["type"]))): #if there is a gene entry
    gene_entry = group[group["type"]=="gene"]
    gene_entry["new_geneID"] = [element.split(";")[0] for element in list(gene_entry["new_geneID"])] #select first geneID
    if str(list(gene_entry["strand"])[0]) == "+":
      gene_entry["stop"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["stop"])[0] #set stop from last exon
    elif str(list(gene_entry["strand"])[0]) == "-":
      gene_entry["start"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["start"])[0] #set start from first exon
    first_gene_df = pd.concat([first_gene_df, gene_entry])
  if "transcript" in list(set(list(group["type"]))): #if there is a gene entry
    transcript_entry = group[group["type"]=="transcript"]
    transcript_entry["new_geneID"] = [element.split(";")[0] for element in list(transcript_entry["new_geneID"])] #select first geneID
    transcript_entry["new_transcriptID"] = [element.split(";")[0]+transcript_suffix for element in list(transcript_entry["new_geneID"])]
    if str(list(transcript_entry["strand"])[0]) == "+":
      transcript_entry["stop"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["stop"])[0] #set stop from last exon
    elif str(list(transcript_entry["strand"])[0]) == "-":
      transcript_entry["start"] = list(first_gene_df[(first_gene_df["exon_number"]==last_ex) & (first_gene_df["type"]=="exon")]["start"])[0] #set stop from first exon
    first_gene_df = pd.concat([first_gene_df, transcript_entry])          
  #The start codon has an exon number, so it is not necessary to add it.
  #if there is a start codon in the group, add it to the first_gene_df
  #if "start_codon" in list(set(list(group["type"]))):
   # start_codon_entry = group[group["type"]=="start_codon"]
    #first_gene_df = pd.concat([first_gene_df, start_codon_entry])
  return(first_gene_df)
 

###### Read inputs
gtf_df = pd.read_table(gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase","attribute"])
#gtf_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/DB/gtf/ref/BmA_annot-B.gtf", sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase","attribute"])
broken_genes_list = list(pd.read_table(broken_genes_file, sep="\t", index_col=False, header=None, names=["broken_genes"])["broken_genes"])
#broken_genes_list = list(pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/broken_genes/BmA-selected_broken_genes.tab", sep="\t", index_col=False, header=None, names=["broken_genes"])["broken_genes"])
#Header: species, chimeric_geneID, orthogroup_ids, first-last_aligned_ex, first-last_aligned_aa, ex:start|stop_coord, chimeric_class
chimeric_genes_df = pd.read_table(chimeric_genes_file, sep="\t", index_col=False, header=0)
chimeric_genes_df = chimeric_genes_df[chimeric_genes_df["species"]==species] #subset for the species of interest
#chimeric_genes_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/chimeric_proteins/classified_chimeric_genes.tab", sep="\t", index_col=False, header=0)
#Header: geneID, new_IDs, category
geneIDs_df = pd.read_table(geneIDs_file, sep="\t", index_col=False, header=0)
#geneIDs_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/corrected_gtfs/BmA_new_geneIDs.txt", sep="\t", index_col=False, header=0)

#Header: species, geneID_prefix, geneID_length, transcript_suffix, protein_suffix
params_df = pd.read_table(my_params_file, sep="\t", header=0, index_col=False)
transcript_suffix = str(list(params_df[params_df["species"]==species]["transcript_suffix"])[0])
protein_suffix = str(list(params_df[params_df["species"]==species]["protein_suffix"])[0])
#params_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/src/snakemake/1_PRELIMINARY_ANNOTATION_FIXES/geneID_params.txt", sep="\t", header=0, index_col=False)

####### Chimeric genes
###Preprocessing: get exon boundaries for all categories
chimeric_genes_to_correct_df = chimeric_genes_df[chimeric_genes_df["chimeric_class"]!="COMPLETE_OVERLAP"]
#generate a boundary exons df
res = def_boundary_exon(chimeric_genes_to_correct_df)
boundary_ex_df = res[0]
unresolved_chimeric_df = res[1]
#create dictionary with correspondence between geneID - left boundary and geneID - right boundary
geneID_left_bound_dict = pd.Series(boundary_ex_df.boundary_ex_left.values, index=boundary_ex_df.chimeric_geneID).to_dict()
geneID_right_bound_dict = pd.Series(boundary_ex_df.boundary_ex_right.values, index=boundary_ex_df.chimeric_geneID).to_dict()

###Correct GTF
#subset GTF to chimeric genes to correct
chimeric_GTF_df = gtf_df.loc[gtf_df["geneID"].isin(chimeric_genes_to_correct)]
#add columns with relevant information
chimeric_GTF_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(chimeric_GTF_df["attribute"]) for part in element.split(";") if "gene_id" in part]
chimeric_GTF_df["exon_number"] = add_exon_number(list(chimeric_GTF_df["attribute"]))
chimeric_GTF_df["attribute_mod"] = separate_attributes(list(chimeric_GTF_df["attribute"])) #transform the attribute field in a list of tuples

#add new geneIDs as a separate column
geneIDs_dict = pd.Series(geneIDs_df.new_IDs.values, index=geneIDs_df.geneID).to_dict()
chimeric_GTF_df["new_geneID"] = chimeric_GTF_df["geneID"].map(geneIDs_dict) #the new geneID for the chimeric gene is in the format "geneID1;geneID2"


#group by chimeric geneID and cycle on the groups
grouped_chimeric_GTF_df = chimeric_GTF_df.groupby("geneID")
for chimeric_gene, group:
  ### First gene: exon number <= boundary_ex_left
  first_ex = min([int(element) for element in list(group["exon_number"]) if math.isnan(element) == False]) #I have to use min because sometimes the enumeration starts from 2
  last_ex = geneID_left_bound_dict[chimeric_gene]

  first_gene_df = group[group["exon_number"]<=last_ex]
  #generate the geneID, transcriptID and proteinID
  first_gene_df["new_geneID"] = [element.split(";")[0] for element in list(first_gene_df["new_geneID"])]
  first_gene_df["new_transcriptID"] = first_gene_df["new_geneID"]+transcript_suffix
  first_gene_df["new_proteinID"] = first_gene_df["new_geneID"]+protein_suffix
  #add entries for gene, transcript and start codon (if present in the original)
  first_gene_df = add_entries_first_gene(first_gene_df, group, first_ex, last_ex, transcript_suffix)
  #update the geneID, transcriptID and proteinID
  first_gene_df["attribute_mod"] = modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "gene_id", list(first_gene_df["new_geneID"]))
  first_gene_df["attribute_mod"] =  modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "transcript_id", list(first_gene_df["new_transcriptID"]))
  first_gene_df["attribute_mod"] =  modify_value_in_tuple(list(first_gene_df["attribute_mod"]), "protein_id", list(first_gene_df["new_proteinID"]))
 
  #Second gene: exon number >= boundary_ex_right
  second_gene_df = group[group["exon_number"]>=geneID_right_bound_dict[chimeric_gene]]
  #renumber the exons (both exons and CDS)
  #update geneID, transcriptID and proteinID
  #add new gene entry and transcript entry
