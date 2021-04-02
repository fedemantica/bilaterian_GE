#!/usr/bin/env python3

import argparse
import pandas as pd
import re
import math

parser = argparse.ArgumentParser(description="Script to compute the aminoacidid portion coded by each exon")
parser.add_argument("--gtf", "-g", required=True, metavar="gtf", help="gtf file")
parser.add_argument("--ref_prot", "-r", required=True, metavar="ref_prot", help="ref protein file with col1=geneID, col2=ref_prot_ID")
parser.add_argument("--output", "-o", required=True, metavar="output_dir", help="path to output file")

#Read arguments
args = parser.parse_args()
gtf_file = args.gtf
ref_prot_file = args.ref_prot
output_file = args.output

gtf = pd.read_table(gtf_file, sep="\t", index_col=False, header=None, names=["chr", "db", "type", "start", "stop", "score", "strand", "phase","attribute"])
ref_prot_df = pd.read_table(ref_prot_file, sep="\t", index_col=False, header=None, names=["geneID", "ref_prot_ID"])
#Filter for CDS exons only
gtf_CDS = gtf[(gtf["type"]=="CDS") & (gtf["attribute"].str.contains("protein_id"))]
#Change data type
gtf_CDS["start"] = gtf_CDS["start"].astype(int) #Change the data type
gtf_CDS["stop"] = gtf_CDS["stop"].astype(int) #Change the data type
gtf_CDS["phase"] = gtf_CDS["phase"].astype(int)

#Add columns with relevant information
gtf_CDS["geneID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(gtf_CDS["attribute"]) for part in element.split(";") if "gene_id" in part]
gtf_CDS["exon_number"]  = [int(re.sub(".*[ ]", "", re.sub('"', "", part))) for element in list(gtf_CDS["attribute"]) for part in element.split(";") if "exon_number" in part]
#This is necessary for when the proteinID field is repeated at the end of the GTF.
if len([element for element in  list(gtf_CDS["attribute"])[0].split(";") if "protein_id" in element]) == 2:
  gtf_CDS["proteinID"] = [re.sub(".*[ ]", "", re.sub('"', "", element.split(";")[-2])) for element in list(gtf_CDS["attribute"])] 
else:
  gtf_CDS["proteinID"] = [re.sub(".*[ ]", "", re.sub('"', "", part)) for element in list(gtf_CDS["attribute"]) for part in element.split(";") if "protein_id" in part]


#Filter only for ref proteins
gtf_CDS_ref = gtf_CDS[gtf_CDS["proteinID"].isin(list(ref_prot_df["ref_prot_ID"]))]

#groupby prot_ID
grouped_df = gtf_CDS_ref.groupby("proteinID")
final_df = pd.DataFrame()
for name, group in grouped_df:
  group = group.sort_values(["exon_number"]) #order by exon number
  group["size"] = ((group["stop"]-group["start"]+1-group["phase"])/3).apply(lambda x: math.ceil(x))
  group["last_aa"] = group["size"].cumsum()
  group["first_aa"] = [1]+[element+1 for element in list(group["last_aa"])[0:len(list(group["last_aa"]))-1]]
  group["aa_coords"] = [str(element[0])+"-"+str(element[1]) for element in list(zip(list(group["first_aa"]), list(group["last_aa"])))]
  final_df = pd.concat([final_df, group])

#save to final output
final_df["ID"] = final_df["proteinID"]+"|"+final_df["geneID"]
final_df["chr"] = final_df["chr"].astype(str)
final_df["start"] = final_df["start"].astype(str) #Change the data type
final_df["stop"] = final_df["stop"].astype(str) #Change the data type
final_df["coords"] = final_df["chr"]+":"+final_df["start"]+"-"+final_df["stop"]
final_df_ordered = final_df[["ID", "exon_number", "aa_coords", "coords", "strand"]]

#write to file
final_df_ordered.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
