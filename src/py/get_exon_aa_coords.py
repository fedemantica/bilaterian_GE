#!/usr/bin/env python3

import argparse
import pandas as pd
import glob
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
  group = group.sort_values(["exon_number"]) #order by exon number
  group["size"] = ((group["stop"]-group["start"]+1-group["phase"])/3).apply(lambda x: math.ceil(x))
  group["last_aa"] = group["size"].cumsum()
  group["first_aa"] = [1]+[element+1 for element in list(group["last_aa"])[0:len(list(group["last_aa"]))-1]]
  group["aa_coords"] = [str(element[0])+"-"+str(element[1]) for element in list(zip(list(group["first_aa"]), list(group["last_aa"])))]
  final_df = pd.concat([final_df, group])

#save to final output
final_df["ID"] = final_df["proteinID"]+"|"+final_df["geneID"]
final_df["start"] = final_df["start"].astype(str) #Change the data type
final_df["stop"] = final_df["stop"].astype(str) #Change the data type
final_df["coords"] = final_df["chr"]+":"+final_df["start"]+"-"+final_df["stop"]
final_df_ordered = final_df[["ID", "exon_number", "aa_coords", "coords", "strand"]]

#write to file
final_df_ordered.to_csv(output_file, sep="\t", index=False, header=False, na_rep="NA")
