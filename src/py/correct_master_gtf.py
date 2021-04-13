#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import csv

parser = argparse.ArgumentParser(description="Script to correct broken and chimeric genes in the master gtf")
parser.add_argument("--species", "-s", required=True, metavar="species", help="Species of interest")
parser.add_argument("--gtf", "-g", required=True, metavar="gtf", help="Master gtf to be corrected")
parser.add_argument("--brochi_gtf", "-bg", required=True, metavar="brochi_gtf", help="Corrected reference gtf (subsetted to brochi genes)")
parser.add_argument("--IDs", "-i", required=True, metavar="input", help="File with col1=original gene IDs (of broken or chimeric genes) and col2=new geneID after correction")
parser.add_argument("--output", "-o", required=True, metavar="output", help="Path to output file for corrected master GTF")

###### Read arguments
args = parser.parse_args()
species = args.species
gtf_file = args.gtf
brochi_gtf_file = args.brochi_gtf
geneIDs_file = args.IDs
output_file = args.output

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

def modify_value_in_tuple(attribute_field, category, new_value):
  final_list = []
  index = 0
  for field in attribute_field:
    for part in field:
      if part[0] == category:
        sub_index = field.index(part)
        field[sub_index][1] = f'"{new_value[index]}"' #this prints the string with surrounding quotes
    final_list = final_list+[field]
    index = index+1
  return(final_list)

def rebuild_attribute_entry(mod_attribute_field):
  final_list = []
  for field in mod_attribute_field:
    new_entry = '; '.join([' '.join(element) for element in field])
    final_list = final_list + [new_entry]
  return(final_list)

##################################
###### READ INPUTS ###############
##################################
#master gtf
gtf_df = pd.read_table(gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
gtf_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part] #add geneID as a separate field
#corrected brochi gtf
brochi_gtf_df = pd.read_table(brochi_gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
brochi_gtf_df["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(brochi_gtf_df["attribute"]) for part in element.split(";") if "gene_id" in part]
####old-new geneID dictionary: Header: geneID, new_IDs, category
geneIDs_df = pd.read_table(geneIDs_file, sep="\t", index_col=False, header=0)
geneIDs_dict = pd.Series(geneIDs_df.new_IDs.values, index=geneIDs_df.geneID).to_dict()
#correct broken genes entries
broken_genes_keys = [key for key in list(geneIDs_dict.keys()) if ";" in key]
for key in broken_genes_keys:
  for num in list(range(0, len(key.split(";")))):
    geneIDs_dict[key.split(";")[num]] = geneIDs_dict[key]
####new-old dictionary
reverse_geneID_dict = pd.Series(geneIDs_df.geneID.values, index=geneIDs_df.new_IDs).to_dict()
chimeric_genes_keys = [key for key in list(reverse_geneID_dict .keys()) if ";" in key]
for key in chimeric_genes_keys:
  for num in list(range(0, len(key.split(";")))):
    reverse_geneID_dict[key.split(";")[num]] = reverse_geneID_dict[key]

#gtf_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/corrected_gtfs/Bm0_annot.gtf", sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase","attribute"])
#brochi_gtf_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/corrected_gtfs/BmA_annot-B-brochi_only.gtf", sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase", "attribute"])
#geneIDs_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/broccoli/BmA_version/corrected_gtfs/BmA_new_geneIDs.txt", sep="\t", index_col=False, header=0)
#gtf_df["attribute"] = separate_attributes(list(gtf_df["attribute"]))
#brochi_gtf_df["attribute"] = separate_attributes(list(brochi_gtf_df["attribute"]))

##################################
####### BROKEN GENES #############
##################################
#For the broken genes: I translate the previous geneID and I add the new reference transcript as an extra transcript
broken_genes = list(geneIDs_df[geneIDs_df["category"]=="broken"]["new_IDs"])
corrected_broken_genes = [element for element in broken_genes if element in list(brochi_gtf_df["geneID"])] #new gene ID of the repaired genes
broken_genes_to_repair = list(geneIDs_df[geneIDs_df["new_IDs"].isin(corrected_broken_genes)]["geneID"]) #genes to select in the master gtf (broken genes corresponding to repair genes)
broken_genes_to_repair = [part for element in broken_genes_to_repair for part in element.split(";")]
gtf_df = gtf_df.loc[~(gtf_df["geneID"].isin(broken_genes_to_repair))]

#filter gtf only on those genes. This is to save RAM.
#broken_genes_gtf = gtf_df.loc[gtf_df["geneID"].isin(broken_genes_to_repair)]
#remove gene entries
#broken_genes_gtf = broken_genes_gtf.loc[broken_genes_gtf["type"]!="gene"]
#broken_genes_gtf["attribute"] = separate_attributes(list(broken_genes_gtf["attribute"]))
#remove the filtered genes from the whole gtf.
#gtf_df = gtf_df.loc[~(gtf_df["geneID"].isin(broken_genes_to_repair))]

#add new geneID for broken genes
#broken_genes_gtf["geneID"] = broken_genes_gtf["geneID"].map(geneIDs_dict)
#replace geneID in the attribute field
#broken_genes_gtf["attribute"] = modify_value_in_tuple(list(broken_genes_gtf["attribute"]), "gene_id", list(broken_genes_gtf["geneID"])) #Replace gene ID of existing broken transcripts
#rebuild the attribute field
#broken_genes_gtf["attribute"] = rebuild_attribute_entry(list(broken_genes_gtf["attribute"]))

##################################
####### CHIMERIC GENES ###########
##################################
#For the chimeric genes: I remove the corrected chimeric genes from the master gtf and I only keep the new ones.
chimeric_genes = [part for element in list(geneIDs_df[geneIDs_df["category"]=="chimeric"]["new_IDs"]) for part in element.split(";")] #select all chimeric genes
corrected_chimeric_genes = [element for element in chimeric_genes if element in list(brochi_gtf_df["geneID"])] #filter only by the chimeric genes corrected in the brochi
chimeric_genes_to_repair = list(set([reverse_geneID_dict[element] for element in corrected_chimeric_genes])) #get the original geneID of the brochi genes
#remove the entries of the corrected repaired genes
gtf_df = gtf_df.loc[~(gtf_df["geneID"].isin(chimeric_genes_to_repair))]

#############################################
########## JOIN ALL GTF PARTS ###############
#############################################
#combine master and brochi gtf
#final_df = pd.concat([gtf_df, broken_genes_gtf, brochi_gtf_df]) #join all parts
print(gtf_df.shape)
print(brochi_gtf_df.shape)
final_df = pd.concat([gtf_df, brochi_gtf_df]) #join all parts
#remove  geneID column
final_df = final_df.drop(columns=["geneID"])
#final_df = final_df.sort_values(by=["chr", "start", "type"], ascending=[True,True,False]) #order gtf (still need to figure out the exon-CDS order)
final_df.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA", quoting=csv.QUOTE_NONE)  #save to file
