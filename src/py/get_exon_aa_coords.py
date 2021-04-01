#!/usr/bin/env python3

import argparse
from Bio import SeqIO
import pandas as pd
import glob
import re
import math

parser = argparse.ArgumentParser(description="Derive the portion of each chimeric gene actually aligned in each orthtogroup it belongs to")
parser.add_argument("--input_dir", "-i", required=True, metavar="input_dir", help="Directory containing all the multiple alignments between the chimeric genes and the respective outgroups")
parser.add_argument("--alignment", "-a", required=True, metavar="alignment", help="List of alignment files prefix in the format chimeric_geneID-orthogroupID")
parser.add_argument("--exon_positions_dir", "-e", required=True, metavar="exon_positions_dir", help="Directory containing all species files with the aminoacidic positions corresponding to each exon")
parser.add_argument("--overlap_stringency", "-s", required=True, metavar="overlap_stringency", help="Minumum percentage of aligned genes in OG required to define the overlap")
parser.add_argument("--output", "-o", required=True, metavar="output_dir", help="path to output file")

#Read arguments
args = parser.parse_args()
my_input_dir = args.input_dir
my_exon_positions_dir = args.exon_positions_dir
my_overlap_stringency = float(args.overlap_stringency)
my_output_file = args.output

gtf = pd.read_table("Hs2_annot-B.gtf", sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase","attribute"])
ref_prot_df = pd.read_table("/users/mirimia/fmantica/projects/bilaterian_GE/data/DB/ref_prot_info/Hs2_gene_ref_prot_dict", sep="\t", index_col=False, header=None, names=["geneID", "ref_prot_ID"])
#Filter for CDS exons only
gtf_CDS = gtf[(gtf["type"]=="CDS") & (gtf["attribute"].str.contains("protein_id"))]

#Add columns with relevant information
gtf_CDS["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", element.split(";")[0])) for element in list(gtf_CDS["attribute"])]
gtf_CDS["proteinID"] = [part for element in list(gtf_CDS["attribute"]) for part in element.split(";") if "protein_id" in part] #This ugly thing is necessary from when the proteinID is repeated in the GTF
gtf_CDS["proteinID"] = [re.sub(".*[ ]", "", re.sub('"', "", element)) for element in list(gtf_CDS["proteinID"])]
gtf_CDS["exon_number"]  = [int(re.sub(".*[ ]", "", re.sub('"', "", part))) for element in list(gtf_CDS["attribute"]) for part in element.split(";") if "exon_number" in part]
gtf_CDS["phase"] = pd.to_numeric(gtf_CDS["phase"])

#Filter only for ref proteins
gtf_CDS_ref = gtf_CDS[gtf_CDS["proteinID"].isin(list(ref_prot_df["ref_prot_ID"]))]

#groupby prot_ID
grouped_df = gtf_CDS_ref.groupby("proteinID")
final_df = pd.DataFrame()
for name, group in grouped_df:
  #order by exon number
  group = group.sort_values(["exon_number"])
  group["size"] = ((group["stop"]-group["start"]+1-group["phase"])/3).apply(lambda x: math.ceil(x))
  group["last_aa"] = group["size"].cumsum()
  group["first_aa"] = [1]+[element+1 for element in list(group["last_aa"])[0:len(list(group["last_aa"]))-1]]
  group["aa_coords"] = [str(element[0])+"-"+str(element[1]) for element in list(zip(list(group["first_aa"]), list(group["last_aa"])))]
  #group["aa_coords"] = [str(i) for i in list(group["first_aa"])]+"-"+[str(i) for i in list(group["last_aa"])]
  final_df = pd.concat([final_df, group])

#save to final output
final_df["ID"] = final_df["proteinID"]+"|"+final_df["geneID"]
final_df["start"] = final_df["start"].astype(str) #just change the data type
final_df["stop"] = final_df["stop"].astype(str) #just change the data type
final_df["coords"] = final_df["chr"]+":"+final_df["start"]+"-"+final_df["stop"]
final_df_ordered = final_df[["ID", "exon_number", "aa_coords", "coords", "strand"]]

#write to file
final_df_ordered.to_csv(my_output_file, sep="\t", index=False, header=False, na_rep="NA")


#test_df = gtf_CDS[gtf_CDS["proteinID"]=="ENSP00000349216"]
test_df = test_df.sort_values(["exon_number"])
#math.ceil is a function which returns the smallest integer which is greater than or equal to the input value
#In this way, when the phase is either 1 or 2, it goes to the first aminoacid of the following exon.
#test_df["size"] = ((test_df["stop"]-test_df["start"]+1)/3).apply(lambda x: math.ceil(x))
#test_df["size"] = ((test_df["stop"]-test_df["start"]+1-test_df["phase"])/3).apply(lambda x: math.ceil(x))
#test_df["last_aa"]=test_df["size"].cumsum()
#test_df["first_aa"] = [1]+[element+1 for element in list(test_df["last_aa"])[0:len(list(test_df["last_aa"]))-1]]
#order by exon number
#compute the size in terms of AA: stop-start+1
#compute the cumulative sum of the sizes (end of the aa range)
#select 1+end[0:penultimo] of sum_cum and set to start of the range


test = [re.sub(".*[ ]", "", re.sub('"', "", part[0]) for element in list(gtf["attribute"]) for part in element.split(";") if "protein_id" in part]

[part[0] for element in prova for part in element.split(";") if "e" in part]
